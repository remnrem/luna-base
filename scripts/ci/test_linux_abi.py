#!/usr/bin/env python3
"""Regression tests for Linux release compatibility checks."""
import unittest
from check_linux_abi import validate

DYNAMIC = '(NEEDED) Shared library: [libc.so.6]\n(NEEDED) Shared library: [libm.so.6]'


class CompatibilityTest(unittest.TestCase):
    def test_accepts_baseline(self):
        self.assertEqual(validate(DYNAMIC, 'Name: GLIBC_2.2.5 Name: GLIBC_2.28')[1], (2, 28))

    def test_rejects_reported_hpc_failures(self):
        for version in ('2.35', '2.38', '2.43'):
            with self.subTest(version=version), self.assertRaisesRegex(RuntimeError, 'exceed baseline'):
                validate(DYNAMIC, 'Name: GLIBC_' + version)

    def test_rejects_external_gcc_runtime(self):
        for library in ('libstdc++.so.6', 'libgcc_s.so.1', 'lib_lightgbm.so', 'libfftw3.so.3'):
            with self.subTest(library=library), self.assertRaisesRegex(RuntimeError, 'Unexpected shared'):
                validate(DYNAMIC + f'\n(NEEDED) Shared library: [{library}]', 'Name: GLIBC_2.28')

    def test_rejects_cpp_symbol_requirements(self):
        for symbol in ('CXXABI_1.3.15', 'GLIBCXX_3.4.33'):
            with self.subTest(symbol=symbol), self.assertRaisesRegex(RuntimeError, 'external C\\+\\+'):
                validate(DYNAMIC, 'Name: GLIBC_2.28 Name: ' + symbol)

    def test_rejects_hidden_build_paths(self):
        with self.assertRaisesRegex(RuntimeError, 'search path'):
            validate(DYNAMIC + '\n(RUNPATH) Library runpath: [/work/deps/lib]', 'Name: GLIBC_2.28')

    def test_rejects_named_or_private_glibc(self):
        for tag in ('GLIBC_PRIVATE', 'GLIBC_ABI_DT_RELR'):
            with self.subTest(tag=tag), self.assertRaises(RuntimeError):
                validate(DYNAMIC, 'Name: GLIBC_2.28 Name: ' + tag)

    def test_missing_requirements_do_not_pass(self):
        with self.assertRaises(RuntimeError):
            validate(DYNAMIC, '')


if __name__ == '__main__':
    unittest.main()
