#!/bin/bash

if [ "${1,,}" == --help ]; then
    printf "
Usage: 
    OrganizeFileNames.sh [prefix] [extension of files to be renamed]

This will organize files by order in folder. Ex: 

    $ ls
    abcde.png    lfjkdls.txt    ndkfwiueh.txt    qrstu.png    wkjhfr.txt
    
    $ bash OrganizeFileNames.sh MyFiles txt
    Files that will be modified:
    lfjkdls.txt    ndkfwiueh.txt    wkjhfr.txt
    Are you sure you want to continue? [y/N] y

    Done renaming files. Log recorded to file ./RenameLog .
    
    $ ls
    abcde.png    MyFiles01.txt    MyFiles02.txt    MyFiles03.txt    qrstu.png

Note that only files with the selected extension will be modified.

"
    exit 0
fi

echo "Files that will be modified:"
ls *."$2"
echo "Are you sure you want to continue? [y/N]"

read Ans

if [ \( "${Ans,,}" == y \) -o \( "${Ans,,}" == yes \) ]; then
	Number=1
	for File in *."$2"; do 
		aa=$(printf "%02d" "$Number")
		mv "$File" "$1""$aa"."${File##*.}"; 
		echo ""$File" moved to "$1""$aa"."${File##*.}"" >> RenameLog;
		((Number++))
	done
    echo "Done renaming files. Log recorded to file ./RenameLog ."
else
	echo "Aborted."
fi

