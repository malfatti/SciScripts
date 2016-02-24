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

LFPCh = [1, 16]
BrokenCh = []
FilterData = True
#==========#==========#==========#==========#

import glob
import Kwik
import os

BrokenCh = [_ - 1 for _ in BrokenCh]
Channels = [_ for _ in list(range(16)) if _ not in BrokenCh]


## Remove date from folder name
DirList = glob.glob('KwikFiles/*'); DirList.sort()
for FolderName in DirList:
    NewFolderName = ''.join([FolderName[11:]])
    NewFolderName = NewFolderName.replace("-", "")
    os.rename(FolderName, NewFolderName)
    print(FolderName, ' moved to ', NewFolderName)



DirList = glob.glob('OpenEphysFiles/*'); DirList.sort()

Raw = Kwik.load()
Rate = Data['info']['sample_rate']