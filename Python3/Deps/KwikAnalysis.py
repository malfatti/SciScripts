# -*- coding: utf-8 -*-
"""
    Copyright (C) 2015  T. Malfatti
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

import glob
import Kwik
import matplotlib.pyplot as plt
import numpy as np
import os
import shelve
from scipy import signal

def RemoveDateFromFolderName():
    RenameFolders = input('Rename folders in KwikFiles/* (BE CAREFUL)? [y/N] ')
    if RenameFolders in ['y', 'Y', 'yes', 'Yes', 'YES']:
        DirList = glob.glob('KwikFiles/*'); DirList.sort()
        for FolderName in DirList:
            NewFolderName = ''.join([FolderName[:10], FolderName[21:]])
            NewFolderName = NewFolderName.replace("-", "")
            os.rename(FolderName, NewFolderName)
            print(FolderName, ' moved to ', NewFolderName)
        del(RenameFolders, DirList, FolderName, NewFolderName)


def ABR(FileName, ABRCh=[1, 16], ABRTimeBeforeTTL=0, ABRTimeAfterTTL=12, 
        ABRTTLCh=1):
    print('empty for now')


def GPIAS(FileName):
    with shelve.open(FileName) as Shelve:
        DataInfo = Shelve['DataInfo']
        Freqs = Shelve['Freqs']
        FreqSlot = Shelve['FreqSlot']
    
    for Key, Value in DataInfo.items():
        exec(str(Key) + '=' + 'Value')
    del(Key, Value)
    