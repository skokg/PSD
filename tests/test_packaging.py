"""Verify installed artifacts, and optionally inspect PSD_SDIST_PATH."""
import base64
import hashlib
from importlib import metadata
import os
from pathlib import Path, PurePosixPath
import tarfile
import unittest

import psd_python as psd


class PackagingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.distribution = metadata.distribution('psd_python')

    def test_metadata(self):
        info = self.distribution.metadata
        self.assertEqual(self.distribution.version, '2.0.0')
        self.assertEqual(info['Name'], 'psd_python')
        self.assertEqual(info['Requires-Python'], '>=3.10')
        self.assertEqual(info['License-Expression'], 'MIT')
        self.assertEqual(info['Description-Content-Type'], 'text/markdown')
        self.assertIn('# PSD Python:', info.get_payload())
        readme = Path(__file__).resolve().parents[1] / 'README.md'
        self.assertEqual(info.get_payload().strip(), readme.read_text(encoding='utf-8').strip())
        self.assertEqual(set(self.distribution.requires),
                         {'numpy>=1.23.5', 'scipy>=1.9.3', 'scikit-image>=0.20'})

    def test_universal_wheel(self):
        wheel = self.distribution.read_text('WHEEL')
        self.assertIn('Root-Is-Purelib: true', wheel)
        self.assertIn('Tag: py3-none-any', wheel)

    def test_installed_module_and_public_api(self):
        installed = Path(self.distribution.locate_file('psd_python.py')).resolve()
        self.assertEqual(Path(psd.__file__).resolve(), installed)
        source = Path(__file__).resolve().parents[1] / 'src/psd_python.py'
        self.assertEqual(installed.read_bytes(), source.read_bytes())
        self.assertEqual(set(psd.__all__),
                         {'calculate_PSD', 'calculate_PSS', 'construct_circular_kernel_field'})
        self.assertTrue(all(callable(getattr(psd, name)) for name in psd.__all__))
        files = {str(entry) for entry in self.distribution.files}
        self.assertNotIn('PSD_library.py', files)
        self.assertFalse(any(name.startswith(('tests/', 'examples/', 'src/')) for name in files))

    def test_installed_record(self):
        checked = set()
        for entry in self.distribution.files:
            if entry.hash is None:
                continue
            with self.subTest(file=str(entry)):
                content = Path(self.distribution.locate_file(entry)).read_bytes()
                digest = hashlib.new(entry.hash.mode, content).digest()
                encoded = base64.urlsafe_b64encode(digest).rstrip(b'=').decode('ascii')
                self.assertEqual(encoded, entry.hash.value)
                self.assertEqual(len(content), entry.size)
                checked.add(str(entry))
        self.assertIn('psd_python.py', checked)
        self.assertTrue(any(name.endswith('/METADATA') for name in checked))

    def test_license(self):
        entries = [p for p in self.distribution.files
                   if str(p).endswith('.dist-info/licenses/LICENSE')]
        self.assertEqual(len(entries), 1)
        text = Path(self.distribution.locate_file(entries[0])).read_text()
        self.assertIn('Copyright (c) 2022 Gregor Skok', text)
        self.assertIn('Permission is hereby granted', text)

    @unittest.skipUnless(os.environ.get('PSD_SDIST_PATH'), 'Set PSD_SDIST_PATH to inspect the sdist.')
    def test_source_archive(self):
        required = {'pyproject.toml', 'setup.py', 'MANIFEST.in', 'README.md', 'LICENSE',
                    'src/psd_python.py', 'tests/test_psd_python.py',
                    'tests/test_packaging.py', 'tests/test_examples.py'}
        with tarfile.open(os.environ['PSD_SDIST_PATH'], 'r:gz') as archive:
            members = archive.getmembers()
            roots = {PurePosixPath(member.name).parts[0] for member in members}
            self.assertEqual(roots, {'psd_python-2.0.0'})
            files = set()
            for member in members:
                path = PurePosixPath(member.name)
                self.assertNotIn('..', path.parts)
                self.assertFalse(path.is_absolute())
                self.assertFalse(member.issym() or member.islnk())
                if member.isfile():
                    files.add(str(PurePosixPath(*path.parts[1:])))
            self.assertTrue(required <= files, required - files)
            extras = files - required
            self.assertIn('PKG-INFO', extras)
            # Setuptools releases may omit regenerable egg-info metadata.
            self.assertLessEqual(extras, {'PKG-INFO', 'src/psd_python.egg-info/PKG-INFO',
                                     'src/psd_python.egg-info/SOURCES.txt',
                                     'src/psd_python.egg-info/dependency_links.txt',
                                     'src/psd_python.egg-info/requires.txt',
                                     'src/psd_python.egg-info/top_level.txt', 'setup.cfg'})


if __name__ == '__main__':
    unittest.main(verbosity=2)
