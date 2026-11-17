"""Run the README usage example against the installed library."""
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest


class ExampleTests(unittest.TestCase):
    def test_synthetic_example(self):
        readme = Path(__file__).resolve().parents[1] / 'README.md'
        usage = readme.read_text(encoding='utf-8').split('## Usage\n', 1)[1].split('\n## ', 1)[0]
        blocks = re.findall(r'^```python\n(.*?)^```', usage, flags=re.MULTILINE | re.DOTALL)
        self.assertTrue(blocks, 'The README must contain a Python usage example')
        with tempfile.TemporaryDirectory(prefix='psd-example-') as directory:
            result = subprocess.run([sys.executable, '-'], input='\n'.join(blocks), cwd=directory,
                                    capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), 'PSD: 100.1920 grid cells')
        self.assertEqual(result.stderr, '')

    def test_github_manual_module_copy(self):
        source = Path(__file__).resolve().parents[1] / 'src/psd_python.py'
        with tempfile.TemporaryDirectory(prefix='psd-manual-copy-') as directory:
            shutil.copy2(source, Path(directory) / 'psd_python.py')
            code = ('from pathlib import Path\n'
                    'import numpy as np\n'
                    'import psd_python as psd\n'
                    'assert Path(psd.__file__).resolve() == Path("psd_python.py").resolve()\n'
                    'assert psd.calculate_PSD(np.ones((2, 3)), np.ones((2, 3))) == 0\n')
            result = subprocess.run([sys.executable, '-c', code], cwd=directory,
                                    capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, '')
        self.assertEqual(result.stderr, '')


if __name__ == '__main__':
    unittest.main(verbosity=2)
