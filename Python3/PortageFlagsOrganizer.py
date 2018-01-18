#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:16:38 2018

@author: malfatti
"""

import os
from argparse import ArgumentParser
from datetime import datetime

Parser = ArgumentParser()
Parser.add_argument('path', help='Path containing flags files, like "/etc/portage/package.use/"', default=None)
Args = Parser.parse_args()

Date = datetime.now().strftime("%Y%m%d%H%M%S")


if Args.path[-1] == '/': Path = Args.path
else: Path = Args.path + '/'

FlagFiles = os.listdir(Path)

print('Getting flags...')
Flags = []
for File in FlagFiles:
    with open(Path+File, 'r') as F:
        for Flag in F.readlines():
            if '#' not in Flag and Flag != '\n'  and Flag not in Flags:
                Flags.append(Flag)
Flags.sort()

print('Getting categories...')
Files = []
for Flag in Flags:
    Flag = Flag.split('/')[0]
    
    if Flag[0] in ['=', '<', '>']:
        if Flag[1] == '=': Flag = Flag[2:]
        else: Flag = Flag[1:]
    
    if Flag not in Files: Files.append(Flag)
Files.sort()

print('Backing up files to', Path+'.Bkp-'+Date+'...')
os.makedirs(Path+'.Bkp-'+Date)
for File in FlagFiles: os.rename(Path+File, Path+'.Bkp-'+Date+'/'+File)

print('Writing new files...')
for File in Files:
    with open(Path+File, 'w') as F: F.writelines(sorted([_ for _ in Flags if File in _]))

print('Done.')        