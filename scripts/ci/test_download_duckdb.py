#!/usr/bin/env python3
"""Offline regression checks for the local SDK downloader."""
import argparse
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("download_duckdb", Path(__file__).resolve().parents[1] / "download_duckdb.py")
sdk = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sdk)


class DownloadTest(unittest.TestCase):
    def test_all_asset_digests_are_pinned(self):
        for value in [*sdk.STATIC.values(), *sdk.SHARED.values(), sdk.LICENSE_SHA]:
            self.assertRegex(value, r"^[0-9a-f]{64}$")

    def test_platform_selection(self):
        for system, arch, libc, expected in (
            ("Darwin", "arm64", "", "osx-arm64"),
            ("Darwin", "x86_64", "", "osx-amd64"),
            ("Linux", "aarch64", "glibc", "linux-arm64"),
            ("Linux", "x86_64", "musl", "linux-amd64-musl"),
            ("Windows", "AMD64", "", "windows-mingw"),
        ):
            with patch.object(sdk.platform, "system", return_value=system), \
                 patch.object(sdk.platform, "machine", return_value=arch), \
                 patch.object(sdk.platform, "libc_ver", return_value=(libc, "")), \
                 patch.object(sdk.Path, "glob", return_value=[]), \
                 patch.dict(sdk.os.environ, {"MSYSTEM": "UCRT64"}):
                self.assertEqual(sdk.native_platform(), expected)

    def test_windows_requires_ucrt64(self):
        for environment in ("MINGW64", "MINGW32", "CLANG64", "MSYS", ""):
            with self.subTest(environment=environment), \
                 patch.object(sdk.platform, "system", return_value="Windows"), \
                 patch.object(sdk.platform, "machine", return_value="AMD64"), \
                 patch.dict(sdk.os.environ, {"MSYSTEM": environment}):
                with self.assertRaisesRegex(RuntimeError, "UCRT64"):
                    sdk.native_platform()

    def test_sdk_from_old_windows_environment_is_not_reused(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            names = {"include/" + name for name in sdk.HEADERS} | {
                "include/luna_duckdb_static.h", "lib/libduckdb_static.a", "DUCKDB-LICENSE"}
            hashes = {}
            for name in names:
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"fixture")
                hashes[name] = sdk.digest(path)
            (root / "manifest.json").write_text(json.dumps({
                "version": sdk.VERSION, "platform": "windows-mingw", "msystem": "MINGW64", "sha256": hashes,
            }), encoding="utf-8")
            with patch.object(sdk, "native_platform", return_value="windows-mingw"), \
                 patch.dict(sdk.os.environ, {"MSYSTEM": "UCRT64"}):
                with self.assertRaisesRegex(RuntimeError, "Refusing to overwrite"):
                    sdk.install(argparse.Namespace(prefix=root, cache_dir=None, with_shared=False))

    def test_verified_cache_is_reused(self):
        with tempfile.TemporaryDirectory() as temp:
            cache = Path(temp)
            (cache / "archive").write_bytes(b"verified")
            with patch.object(sdk, "run") as run:
                sdk.fetch("https://example.invalid", hashlib.sha256(b"verified").hexdigest(), cache, "archive")
                run.assert_not_called()

    def test_bad_download_does_not_replace_cache(self):
        with tempfile.TemporaryDirectory() as temp:
            cache = Path(temp)
            (cache / "archive").write_bytes(b"old")
            def download(command, **kwargs):
                Path(command[command.index("--output")+1]).write_bytes(b"bad download")
            with patch.object(sdk, "run", side_effect=download):
                with self.assertRaisesRegex(RuntimeError, "SHA-256 mismatch"):
                    sdk.fetch("https://example.invalid", "0"*64, cache, "archive")
            self.assertEqual((cache / "archive").read_bytes(), b"old")

    def test_existing_directory_is_preserved(self):
        with tempfile.TemporaryDirectory() as temp:
            prefix = Path(temp)
            keep = prefix / "user-file"
            keep.write_text("keep", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "Refusing to overwrite"):
                sdk.install(argparse.Namespace(prefix=prefix, cache_dir=None, with_shared=False))
            self.assertEqual(keep.read_text(encoding="utf-8"), "keep")

    def test_only_requested_regular_files_are_extracted(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            archive = root / "archive.tar.gz"
            with tarfile.open(archive, "w:gz") as out:
                for name in ("duckdb.h", "../outside"):
                    member = tarfile.TarInfo(name)
                    member.size = 4
                    out.addfile(member, io.BytesIO(b"data"))
            target = root / "sdk"
            target.mkdir()
            sdk.extract_files(archive, ["duckdb.h"], target)
            self.assertEqual((target / "duckdb.h").read_bytes(), b"data")
            self.assertFalse((root / "outside").exists())


if __name__ == "__main__":
    unittest.main()
