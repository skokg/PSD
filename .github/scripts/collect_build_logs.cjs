/* Collect completed build/install/test job logs, including failed attempts. */
const fs = require('node:fs');
const path = require('node:path');
const crypto = require('node:crypto');

module.exports = async function collect({ github, context, core }) {
  const directory = 'release-bundle';
  const logsDirectory = path.join(directory, 'logs');
  fs.mkdirSync(logsDirectory, { recursive: true });
  const logPath = path.join(directory, 'workflow-logs.log');
  const runUrl = `${context.serverUrl}/${context.repo.owner}/${context.repo.repo}/actions/runs/${context.runId}`;
  fs.writeFileSync(logPath,
    `PSD build, installation and test logs\nRun: ${runUrl}\nCommit: ${context.sha}\n` +
    `Collection attempt: ${process.env.GITHUB_RUN_ATTEMPT}\n` +
    'Includes all available attempts of upstream jobs, including failures.\n');

  const targets = process.env.SELECTED_TARGET === 'all'
    ? ['linux-x86_64', 'windows-x86_64', 'macos-intel', 'macos-apple-silicon']
    : [process.env.SELECTED_TARGET];
  const versions = process.env.SELECTED_PYTHON === 'all'
    ? ['3.10', '3.11', '3.12', '3.13', '3.14', '3.15']
    : [process.env.SELECTED_PYTHON];
  const expected = new Set(targets.flatMap(target =>
    versions.map(version => `Python ${version} / ${target}`)));
  expected.add('Build universal wheel and source archive');
  expected.add('Rebuild and test source archive');
  expected.add('Minimum Python and runtime dependencies');
  const problems = [];
  let jobs = [];
  try {
    const available = await github.paginate(
      'GET /repos/{owner}/{repo}/actions/runs/{run_id}/jobs',
      { ...context.repo, run_id: context.runId, filter: 'all', per_page: 100 });
    jobs = available.filter(job => expected.has(job.name));
    jobs.sort((a, b) => a.name.localeCompare(b.name) || a.id - b.id);
  } catch (error) {
    problems.push(`Could not list jobs: ${error.message}`);
  }
  fs.writeFileSync(path.join(directory, 'JOBS.json'), JSON.stringify(jobs, null, 2) + '\n');
  const latest = new Map();
  for (const job of jobs) {
    const previous = latest.get(job.name);
    if (!previous || job.id > previous.id) latest.set(job.name, job);
    const heading = '\n' + '='.repeat(80) + '\n' + job.name +
      `\nJob ID: ${job.id}\nStarted: ${job.started_at}\n` +
      `Result: ${job.status} / ${job.conclusion}\nURL: ${job.html_url}\n` + '='.repeat(80) + '\n';
    fs.appendFileSync(logPath, heading);
    // Skipped jobs have no runner log; their status is recorded explicitly.
    if (job.conclusion === 'skipped') {
      fs.appendFileSync(logPath, 'Job skipped; no execution log exists.\n');
      continue;
    }
    let text;
    let lastError;
    // GitHub may need a short delay to expose logs after a job completes.
    for (let attempt = 0; attempt < 4; attempt++) {
      try {
        const response = await github.request(
          'GET /repos/{owner}/{repo}/actions/jobs/{job_id}/logs',
          { ...context.repo, job_id: job.id });
        text = typeof response.data === 'string'
          ? response.data : Buffer.from(response.data).toString('utf8');
        break;
      } catch (error) {
        lastError = error;
        const retryable = [404, 410, 429].includes(error.status) || error.status >= 500;
        if (attempt === 3 || !retryable) break;
        await new Promise(resolve => setTimeout(resolve, 1000 * (2 ** attempt)));
      }
    }
    if (text === undefined) {
      const message = `Could not download log for ${job.name} (job ${job.id}): ${lastError.message}`;
      problems.push(message);
      fs.appendFileSync(logPath, 'LOG DOWNLOAD FAILED: ' + message + '\n');
    } else {
      const filename = `${job.id}-${job.name.replace(/[^a-zA-Z0-9._-]+/g, '_')}.log`;
      fs.writeFileSync(path.join(logsDirectory, filename), heading + text + '\n');
      fs.appendFileSync(logPath, text + '\n');
    }
  }
  for (const name of expected) {
    const job = latest.get(name);
    if (!job || job.status !== 'completed' || job.conclusion !== 'success') {
      problems.push(`${name}: ${job ? `${job.status} / ${job.conclusion}` : 'missing job'}`);
    }
  }
  const packages = fs.readdirSync(directory).filter(name => name.endsWith('.whl') || name.endsWith('.tar.gz')).sort();
  if (packages.filter(name => name.endsWith('.whl')).length !== 1) problems.push('Expected one universal wheel.');
  if (packages.filter(name => name.endsWith('.tar.gz')).length !== 1) problems.push('Expected one source archive.');
  const report = `Bundle status: ${problems.length ? 'INCOMPLETE / CHECKS FAILED' : 'COMPLETE / ALL CHECKS PASSED'}\n` +
    `Run: ${runUrl}\nCommit: ${context.sha}\n` +
    `Selected target: ${process.env.SELECTED_TARGET}; Python: ${process.env.SELECTED_PYTHON}\n` +
    `Expected jobs: ${expected.size}; job attempts found: ${jobs.length}\n` +
    `Packages: ${packages.join(', ') || '(none)'}\n` +
    (problems.length ? problems.join('\n') + '\nDo not treat these packages as fully tested.\n' : '') +
    'The collector job itself is still running; its findings are recorded in this report.\n';
  fs.writeFileSync(path.join(directory, 'BUILD_STATUS.txt'), report);
  fs.appendFileSync(logPath, '\n' + report);

  // Include checksums for packages, status, metadata, and combined/individual logs.
  const files = [...packages, 'workflow-logs.log', 'BUILD_STATUS.txt', 'JOBS.json',
    ...fs.readdirSync(logsDirectory).sort().map(name => `logs/${name}`)];
  const checksums = files.map(name => crypto.createHash('sha256')
    .update(fs.readFileSync(path.join(directory, name))).digest('hex') + '  ' + name);
  fs.writeFileSync(path.join(directory, 'SHA256SUMS.txt'), checksums.join('\n') + '\n');
  await core.summary.addHeading('Distributions and detailed logs')
    .addRaw('Download **psd-python-build-results** for packages, **workflow-logs.log**, ' +
      'individual job logs, **JOBS.json**, **BUILD_STATUS.txt**, and **SHA256SUMS.txt**.\n\n')
    .addCodeBlock(report, 'text').write();
  if (problems.length) core.setFailed('Checks failed or logs are incomplete; see BUILD_STATUS.txt in the results ZIP.');
};
