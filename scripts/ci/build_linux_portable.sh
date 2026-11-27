#!/usr/bin/env bash
# Runs INSIDE manylinux_2_28; checkout/cache/upload actions run on the host.
set -euo pipefail

export PATH="/opt/python/cp312-cp312/bin:${PATH}"
export PYTHONIOENCODING=utf-8 PYTHONUTF8=1
export CC=gcc CXX=g++
export CFLAGS="-O2 -fPIC -march=x86-64 -mtune=generic"
export CXXFLAGS="${CFLAGS}"
export LDFLAGS="-static-libgcc -static-libstdc++"
mkdir -p "${HOME}"
[[ "$(getconf GNU_LIBC_VERSION)" == 'glibc 2.28' ]]
[[ "$(uname -m)" == x86_64 ]]
python3 --version
g++ --version
cmake --version

: "${FFTW_VERSION:?}" "${LIGHTGBM_REF:?}" "${LINUX_BUILD_IMAGE:?}"
deps="${PWD}/.ci-deps/linux-manylinux228"
mkdir -p "${deps}"

if [[ "${NATIVE_DEPS_MODE:-luna}" == all || ! -f "${deps}/fftw/lib/libfftw3.a" || ! -f "${deps}/lightgbm/lib_lightgbm.a" ]]; then
  work="$(mktemp -d "${deps}/native-build.XXXXXX")"
  curl --fail --location --retry 3 "https://www.fftw.org/fftw-${FFTW_VERSION}.tar.gz" -o "${work}/fftw.tar.gz"
  tar -xzf "${work}/fftw.tar.gz" -C "${work}"
  (
    cd "${work}/fftw-${FFTW_VERSION}"
    ./configure --disable-shared --enable-static --with-pic --prefix="${deps}/fftw"
    make -j"$(nproc)"
    make install
    cp COPYING "${deps}/FFTW-LICENSE"
  )
  git clone --recursive --branch "${LIGHTGBM_REF}" --depth 1 https://github.com/microsoft/LightGBM "${work}/LightGBM"
  cmake -S "${work}/LightGBM" -B "${work}/lgbm-build" -DBUILD_STATIC_LIB=ON -DUSE_OPENMP=OFF -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_BUILD_TYPE=Release
  cmake --build "${work}/lgbm-build" --parallel "$(nproc)"
  mkdir -p "${deps}/lightgbm/include"
  cp -R "${work}/LightGBM/include/LightGBM" "${deps}/lightgbm/include/"
  cp "${work}/LightGBM/lib_lightgbm.a" "${deps}/lightgbm/"
  cp "${work}/LightGBM/LICENSE" "${deps}/LIGHTGBM-LICENSE"
  cp "${work}/LightGBM/VERSION.txt" "${deps}/lightgbm/"
  rm -rf "${work}"
fi

if [[ "${NATIVE_DEPS_MODE:-luna}" == all || ! -f "${deps}/duckdb/lib/libduckdb_static.a" ]]; then
  python3 scripts/ci/build_duckdb.py --prefix "${deps}/duckdb" --source "${deps}/duckdb-source" --build-dir "${deps}/duckdb-build" --generator 'Unix Makefiles' --jobs 2
  rm -rf "${deps}/duckdb-source" "${deps}/duckdb-build"
fi

# Never reuse objects compiled outside this baseline container.
make ARCH=LINUX DUCKDB=0 LGBM=0 clean
make -j"$(nproc)" ARCH=LINUX PORTABLE_LINUX=1 \
  CC='gcc -march=x86-64 -mtune=generic' CXX='g++ -march=x86-64 -mtune=generic' \
  DUCKDB=1 DUCKDB_PATH="${deps}/duckdb" FFTW="${deps}/fftw" \
  LGBM=1 LGBM_LINK=static LGBM_PATH="${deps}/lightgbm" ORT=0 \
  luna destrat destrat-sqlite behead

python3 scripts/ci/test_linux_abi.py
python3 scripts/ci/check_linux_abi.py .
./luna __LUNA_TESTS__ lgbm verbose
python3 scripts/ci/check_duckdb.py --static-lightgbm --linux-baseline --library "${deps}/duckdb/lib/libduckdb.so"

version="${LUNA_PACKAGE_VERSION:-dev}"
[[ "${version}" =~ ^[A-Za-z0-9._-]+$ ]]
package="$(mktemp -d)"
mkdir "${package}/luna-${version}"
cp luna destrat destrat-sqlite behead "${package}/luna-${version}/"
cp "${deps}/duckdb/DUCKDB-LICENSE" "${deps}/duckdb/DUCKDB-THIRD-PARTY-NOTICES.txt" "${deps}/duckdb/manifest.json" "${package}/luna-${version}/"
cp "${deps}/LIGHTGBM-LICENSE" "${deps}/FFTW-LICENSE" "${package}/luna-${version}/"
python3 - "${package}/luna-${version}/LINUX-BUILD.json" <<'PYTHON'
import json, os, platform, subprocess, sys
from pathlib import Path
Path(sys.argv[1]).write_text(json.dumps({
    'image': os.environ['LINUX_BUILD_IMAGE'], 'glibc_baseline': '2.28',
    'architecture': platform.machine(), 'static_gcc_runtime': True,
    'compiler': subprocess.check_output(['g++', '--version'], text=True),
    'lightgbm': os.environ['LIGHTGBM_REF'], 'fftw': os.environ['FFTW_VERSION'],
}, indent=2) + '\n', encoding='utf-8')
PYTHON
tar -czf linux-luna.tar.gz -C "${package}" "luna-${version}"
python3 scripts/ci/check_duckdb.py --static-lightgbm --linux-baseline --archive linux-luna.tar.gz --library "${deps}/duckdb/lib/libduckdb.so"
rm -rf "${package}"
