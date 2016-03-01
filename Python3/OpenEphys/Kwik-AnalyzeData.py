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

Experiment: stimulation of the brainstem using light and sound, recording with
a silicon probe (16 channels) + 2 tungsten wires + reference screw.
"""
#%% Set experiment details

# Channels
LFPCh = [1, 16]
BrokenCh = []

# ABR
ABRLength = 12    # in ms
ABRTTLCh = 1
#==========#==========#==========#==========#

import glob
import Kwik
import os

BrokenCh = [_ - 1 for _ in BrokenCh]
Channels = [_ for _ in list(range(16)) if _ not in BrokenCh]


## Remove date from folder name
RenameFolders = input('Rename folders (BE CAREFUL)? [y/N] ')
if RenameFolders in ['y', 'Y', 'yes', 'Yes', 'YES']:
    DirList = glob.glob('KwikFiles/*'); DirList.sort()
    for FolderName in DirList:
        NewFolderName = ''.join([FolderName[:10], FolderName[21:]])
        NewFolderName = NewFolderName.replace("-", "")
        os.rename(FolderName, NewFolderName)
        print(FolderName, ' moved to ', NewFolderName)
    del(RenameFolders, DirList, FolderName, NewFolderName)


## Create folders
os.makedirs('Figs', exist_ok=True)    # Figs folder


## ABR traces
DirList = glob.glob('KwikFiles/*'); DirList.sort()
for RecFolder in DirList:
    FilesList = glob.glob(''.join([RecFolder, '/*'])); FilesList.sort()
    for File in FilesList:
        if '.kwd' in File:
            try:
                Raw = Kwik.load(File)
            except OSError:
                    print('File ', File, " is corrupted :'(")
            
        elif '.kwe' in File:
            try:
                Events = Kwik.load(File)
            except OSError:
                print('File ', File, " is corrupted :'(")
            
        elif '.kwik' in File:
            try:
                Events = Kwik.load(File)
            except OSError:
                print('File ', File, " is corrupted :'(")
            
        elif '.kwx' in File:
            try:
                Spks = Kwik.load(File)
            except OSError:
                print('File ', File, " is corrupted :'(")
    
    if 'Raw' not in globals():
        print('.kwd file is corrupted. Skipping dataset...')
        continue
    
    if 'Events' not in globals():
        print('.kwe/.kwik file is corrupted. Skipping dataset...')
        continue
    
    print('Data from ', RecFolder, ' loaded.')
    
    Rate = Raw['info']['sample_rate']
    NoOfSamples = int(ABRLength*Rate*(10**-3))
    