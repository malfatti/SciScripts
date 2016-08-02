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

Functions for manipulating specific .mat files.
"""

import Hdf5F
import KwikAnalysis
import numpy as np
import os
from glob import glob
from scipy import io

def WriteDataToMMSS(FileName, StimType=['Sound'], Override={}):
    DirList = glob('KwikFiles/*'); DirList.sort()
    for Stim in StimType:
        if Override != {}: 
            if 'Stim' in Override.keys():
                Stim = Override['Stim']
        
        Exps = Hdf5F.LoadExpPerStim(Stim, DirList, FileName)
        
        for FInd, RecFolder in enumerate(Exps): 
            Path = os.getcwd() + '/' + RecFolder
            ClusterFile = Path + '/SpkClusters.hdf5'
            Clusters = Hdf5F.LoadClusters(ClusterFile)
            
            Path = os.getcwd() + '/' + RecFolder +'/SepRec/'
            os.makedirs(Path, exist_ok=True)
            
            for Rec in Clusters.keys():
                RecS = "{0:02d}".format(int(Rec))
                        
                print('Writing files for clustering... ', end='')
                WF = []; ST = []; ChList = []
                for Ch in Clusters[RecS].keys():
                    WF.append(Clusters[RecS][Ch]['Spikes'][:])
                    ST.append(Clusters[RecS][Ch]['Timestamps'][:])
                    ChList.append(Ch)
                
                data = {'waveforms': np.array(WF, dtype='object'),
                        'spiketimes': np.array(ST, dtype='object'),
                        'ChList': np.string_(ChList)}
                
                MatName = ''.join(['Exp_', RecS, '.mat'])
                io.savemat(Path+MatName, {'data': data})
                print('Done.')
    return(None)

