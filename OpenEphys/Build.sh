#!/bin/bash

Branch=${1,,}; Branch=${Branch^}
Branches="Development Testing Master"
Here=$(pwd)
OEPath=~/Software/Git/OpenEphys
OEInstallPath="$OEPath"${Branch:0:1}
CPUs=$(cat /proc/cpuinfo | grep processor | wc -l)

cd $OEPath
echo "Fetching upstream..."
git fetch upstream

if [[ " $Branches " =~ " $Branch " ]]; then
    echo "Merging upstream..." 
    git checkout ${Branch,,}
    git merge upstream/${Branch,,}
    
    echo ""
    echo "Compiling GUI..." 
    cd Builds/Linux
    make clean
    make -j$CPUs
    if [ $? -eq 0 ]; then
        echo "GUI compiled."
    else
        echo "GUI compilation failed."
        exit
    fi
    
    echo ""
    echo "Compiling plugins..." 
    make -f Makefile.plugins clean
    make -j$CPUs -f Makefile.plugins
    if [ $? -eq 0 ]; then
        echo "Plugins compiled."
    else
        echo "Plugins compilation failed."
        exit
    fi
    
    cd ../../
    cp Resources/DLLs/Linux-USB3/x64/libokFrontPanel.so Resources/Bitfiles/rhd2000_usb3.bit Builds/Linux/build/
    cp -R Builds/Linux/build/ $OEInstallPath
    
    echo ""
    echo "GUI installed at $OEInstallPath"
    echo ""

    cd $Here
    
else
    echo ""
    echo "Usage: $0 [development | testing | master]"
    echo ""
fi
