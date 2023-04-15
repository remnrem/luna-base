#!/bin/bash

# Notes ( http://clarkkromenaker.com/post/library-dynamic-loading-mac/ )

myArray=("luna" "destrat" "behead" "fixrows")

mkdir ../mac-arm-luna
cd ../mac-arm-luna

cp /opt/homebrew/opt/fftw/lib/libfftw3.3.dylib .
cp /opt/homebrew/opt/libomp/lib/libomp.dylib .
cp ../luna .
cp ../destrat .
cp ../behead .
cp ../fixrows .

for str in ${myArray[@]}; do

    # Instructs the loader to search a list of paths to find the dynamic library
    # Please check the field LC_RPATH in "otool -l luna" output
    install_name_tool -add_rpath @executable_path/. ./$str

    # Change dynamic library loader path
    # Please check the field LC_LOAD_DYLIB in "otool -l luna" output
    install_name_tool -change /opt/homebrew/opt/fftw/lib/libfftw3.3.dylib @rpath/libfftw3.3.dylib ./$str
    install_name_tool -change /opt/homebrew/opt/libomp/lib/libomp.dylib @rpath/libomp.dylib ./$str

done

# Prints paths to all dynamic libraries used by the executable
otool -L ./luna
otool -L ./destrat

echo
echo "Please check the newly created directory named ../mac-arm-luna"
