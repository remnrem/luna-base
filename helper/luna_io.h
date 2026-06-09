
//    --------------------------------------------------------------------
//
//    This file is part of Luna.
//
//    LUNA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Luna is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Luna. If not, see <http://www.gnu.org/licenses/>.
//
//    Please see LICENSE.txt for more details.
//
//    --------------------------------------------------------------------

#ifndef __LUNA_IO_H__
#define __LUNA_IO_H__

// UTF-8-safe file I/O wrappers (C++17, cross-platform).
//
// On Linux/macOS the filesystem is natively UTF-8, so std::ifstream and
// fopen() already accept UTF-8 byte strings unchanged.
//
// On Windows the MSVC/MinGW CRT interprets narrow-string paths using the
// process ANSI code page, which is typically not UTF-8.  The wrappers here
// convert UTF-8 → wchar_t and call the wide-character variants (_wfopen,
// std::filesystem::u8path) so that Chinese (and other non-ASCII) filenames
// work uniformly on all platforms.
//
// Usage:
//   Replace:  std::ifstream IN( path.c_str(), std::ios::in );
//   With:     std::ifstream IN = LunaIO::open_ifstream( path );
//
//   Replace:  FILE* f = fopen( path.c_str(), "rb" );
//   With:     FILE* f = LunaIO::fopen_utf8( path, "rb" );

#include <fstream>
#include <cstdio>
#include <string>
#include <filesystem>

#if defined(WINDOWS) || defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#ifdef IN
#undef IN
#endif
#ifdef OUT
#undef OUT
#endif
#ifdef small
#undef small
#endif
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif
#endif

namespace LunaIO {

namespace fs = std::filesystem;

inline std::ifstream open_ifstream(const std::string& utf8_path,
                                   std::ios_base::openmode mode = std::ios_base::in)
{
#if defined(WINDOWS) || defined(_WIN32)
    return std::ifstream(fs::u8path(utf8_path), mode);
#else
    return std::ifstream(utf8_path, mode);
#endif
}

inline std::ofstream open_ofstream(const std::string& utf8_path,
                                   std::ios_base::openmode mode = std::ios_base::out)
{
#if defined(WINDOWS) || defined(_WIN32)
    return std::ofstream(fs::u8path(utf8_path), mode);
#else
    return std::ofstream(utf8_path, mode);
#endif
}

// UTF-8-safe replacement for fopen().  Needed for EDF binary reads/writes
// because the EDF format requires a raw FILE* handle.
inline FILE* fopen_utf8(const std::string& utf8_path, const char* mode)
{
#if defined(WINDOWS) || defined(_WIN32)
    const int wpath_len = MultiByteToWideChar(CP_UTF8, 0, utf8_path.c_str(), -1, nullptr, 0);
    const int wmode_len = MultiByteToWideChar(CP_UTF8, 0, mode,            -1, nullptr, 0);
    if (wpath_len <= 0 || wmode_len <= 0) return nullptr;
    std::wstring wpath(wpath_len, L'\0');
    std::wstring wmode(wmode_len, L'\0');
    MultiByteToWideChar(CP_UTF8, 0, utf8_path.c_str(), -1, wpath.data(), wpath_len);
    MultiByteToWideChar(CP_UTF8, 0, mode,            -1, wmode.data(), wmode_len);
    return _wfopen(wpath.c_str(), wmode.c_str());
#else
    return std::fopen(utf8_path.c_str(), mode);
#endif
}

} // namespace LunaIO

#endif // __LUNA_IO_H__
