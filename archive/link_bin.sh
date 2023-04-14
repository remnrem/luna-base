#!/bin/bash

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
    install_name_tool -add_rpath @executable_path/. ./$str
    install_name_tool -change /opt/homebrew/opt/fftw/lib/libfftw3.3.dylib @rpath/libfftw3.3.dylib ./$str
    install_name_tool -change /opt/homebrew/opt/libomp/lib/libomp.dylib @rpath/libomp.dylib ./$str
done

otool -L ./luna
otool -L ./destrat

echo
echo "Please check the newly created directory named ../mac-arm-luna"
