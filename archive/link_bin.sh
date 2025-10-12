#!/bin/bash
# Notes ( http://clarkkromenaker.com/post/library-dynamic-loading-mac/ )

myArray=("luna" "destrat" "behead" "fixrows")

rm -rf ../mac-arm-luna
mkdir ../mac-arm-luna
cd ../mac-arm-luna

fftw_lib='/opt/homebrew/opt/fftw/lib/libfftw3.3.dylib'
omp_lib='/opt/homebrew/opt/libomp/lib/libomp.dylib'

cp $fftw_lib .
cp $omp_lib .
cp ../luna .
cp ../destrat .
cp ../behead .
cp ../fixrows .

for str in ${myArray[@]}; do
 
    # Fixing LC_LOAD_DYLIB entries (change dynamic library loader path)
    # Please check the field LC_LOAD_DYLIB in "otool -l luna" output
    install_name_tool -change $fftw_lib @rpath/libfftw3.3.dylib ./$str
    install_name_tool -change $omp_lib @rpath/libomp.dylib ./$str

    # Instructs the loader to search a list of paths to find the dynamic library
    # Please check the field LC_RPATH in "otool -l luna" output
    # Because I used @rpath in the earlier step, the executable must specify a list of search path in LC_RPATH
    install_name_tool -add_rpath @executable_path/. ./$str

done

# Prints paths to all dynamic libraries used by the executable
otool -L ./luna
otool -L ./destrat
echo

# To view what libraries an executable in loading
DYLD_PRINT_LIBRARIES=1 ./luna

echo
echo "Please check the newly created directory named ../mac-arm-luna"
