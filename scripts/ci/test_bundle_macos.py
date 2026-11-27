#!/usr/bin/env python3
"""Exercise relocation with a real Mach-O dependency chain, without a DuckDB SDK."""
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tarfile
import tempfile
import unittest
from unittest.mock import patch

from bundle_macos import bundle, rpaths
from check_duckdb import check_macos_package, check_static_lightgbm


class StaticDependencyPolicyTest(unittest.TestCase):
    def test_rejects_shared_dependencies_on_both_platforms(self):
        for system, extension in (("Darwin", ".dylib"), ("Linux", ".so.3")):
            for name in ("lib_lightgbm", "libomp", "libgomp", "libiomp5", "libfftw3"):
                with self.subTest(system=system, library=name), \
                     patch("check_duckdb.platform.system", return_value=system), \
                     patch("check_duckdb.subprocess.check_output", return_value=name + extension):
                    with self.assertRaisesRegex(RuntimeError, "Unexpected shared"):
                        check_static_lightgbm(Path("package"))

    def test_allows_system_runtime(self):
        for system in ("Darwin", "Linux"):
            with patch("check_duckdb.platform.system", return_value=system), \
                 patch("check_duckdb.subprocess.check_output", return_value="/usr/lib/libc++.1.dylib libc.so.6 libstdc++.so.6"):
                check_static_lightgbm(Path("package"))


@unittest.skipUnless(platform.system() == "Darwin", "requires macOS Mach-O tools")
class RelocationTest(unittest.TestCase):
    def test_transitive_libraries_without_build_directory(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "runner build's libraries"
            source.mkdir()
            package = root / "luna-dev"
            package.mkdir()
            (source / "omp.c").write_text("int omp_fixture(void) { return 7; }\n", encoding="utf-8")
            (source / "lgbm.c").write_text("extern int omp_fixture(void); int lgbm_fixture(void) { return omp_fixture()+1; }\n", encoding="utf-8")
            (source / "main.c").write_text("extern int lgbm_fixture(void); int main(void) { return lgbm_fixture()==8 ? 0 : 1; }\n", encoding="utf-8")
            subprocess.run(["cc", "-dynamiclib", str(source / "omp.c"), "-o", str(source / "libomp.dylib"),
                            "-Wl,-install_name,@rpath/libomp.dylib"], check=True)
            subprocess.run(["cc", "-dynamiclib", str(source / "lgbm.c"), "-o", str(source / "lib_lightgbm.dylib"),
                            "-L" + str(source), "-lomp", "-Wl,-install_name,@rpath/lib_lightgbm.dylib",
                            "-Wl,-rpath," + str(source)], check=True)
            subprocess.run(["cc", str(source / "main.c"), "-o", str(package / "luna"),
                            "-L" + str(source), "-l_lightgbm", "-Wl,-rpath," + str(source)], check=True)
            for name in ("destrat", "destrat-sqlite", "behead"):
                shutil.copy2(package / "luna", package / name)
            with self.assertRaisesRegex(RuntimeError, "Non-relocatable"):
                check_macos_package(package)
            bundle(package)
            source.rename(root / "unavailable-build-libraries")
            archive = root / "package.tar.gz"
            with tarfile.open(archive, "w:gz") as out:
                out.add(package, arcname="luna-dev")
            relocated = root / "download"
            with tarfile.open(archive) as incoming:
                incoming.extractall(relocated, filter="data")
            binaries = relocated / "luna-dev"
            check_macos_package(binaries)
            env = {k: v for k, v in os.environ.items() if not k.startswith("DYLD_")}
            for name in ("luna", "destrat", "destrat-sqlite", "behead"):
                binary = binaries / name
                self.assertFalse(any(p.startswith("/") for p in rpaths(binary)))
                subprocess.run([str(binary)], cwd=relocated, env=env, check=True)
            self.assertTrue((binaries / "lib/lib_lightgbm.dylib").is_file())
            self.assertTrue((binaries / "lib/libomp.dylib").is_file())


if __name__ == "__main__":
    unittest.main()
