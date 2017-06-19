#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Klusta module

@author: malfatti
"""
import os

from DataAnalysis.DataAnalysis import Pairwise
from IO.IO import RunProcess
from IO.Txt import DictPrint, DictWrite


def PrbWrite(File, Channels=list(range(15,-1,-1)), Spacing=25):
    Prb = {'0': {}}
    Prb['0']['channels'] = Channels
    Prb['0']['graph'] = list(Pairwise(Prb['0']['channels']))
    Pos = list(range(0, len(Prb['0']['channels'])*Spacing, Spacing))
    Prb['0']['geometry'] = {str(Ch):(0,Pos[C]) for C,Ch in enumerate(Prb['0']['channels'])}
    DictWrite(File, 'channel_groups = '+DictPrint(Prb))
    
    return(None)


def PrmWrite(File, experiment_name, prb_file, raw_data_files, sample_rate, 
             n_channels, dtype, spikedetekt={}, klustakwik2={}):
    
    traces = {
        'raw_data_files': raw_data_files,
        'sample_rate': sample_rate,
        'n_channels': n_channels,
        'dtype': dtype,
    }
    
    with open(File, 'w') as F:
        F.write('experiment_name = "'+experiment_name+'"')
        F.write('\n\n')
        F.write('prb_file = "'+prb_file+'"')
        F.write('\n\n')
        F.write('traces = '+DictPrint(traces, breaklineat=50))
        F.write('\n\n')
        F.write('spikedetekt = '+DictPrint(spikedetekt))
        F.write('\n\n')
        F.write('klustakwik2 = '+DictPrint(klustakwik2))
        F.write('\n\n')
    
    return(None)


def Run(PrmFile, Path, Overwrite=False):
    Here = os.getcwd(); os.chdir(Path)
    Klusta = os.environ['HOME']+'/Software/Miniconda3/envs/klusta/bin/klusta'
    Cmd = [Klusta, PrmFile]
    if Overwrite: Cmd.append('--overwrite')
    
    print('Entering in directory', Path, '...')
    print('Clustering spikes...')
    ReturnCode = RunProcess(Cmd, PrmFile+'.log'); os.chdir(Here)
    print('Going back to', Here, '...')
    
    if ReturnCode: print('Error: ReturnCode', ReturnCode)
    else: print('Done clustering.')
    return(None)


def Phy_Run(KwikFile, Path):
    Here = os.getcwd(); os.chdir(Path)
    Phy = os.environ['HOME']+'/Software/Miniconda3/envs/klusta/bin/phy'
    Cmd = [Phy, 'kwik-gui', KwikFile, '--debug']
    
    print('Entering in directory', Path, '...')
    print('Loading phy...')
    ReturnCode = RunProcess(Cmd, KwikFile+'.log'); os.chdir(Here)
    print('Going back to', Here, '...')
    
    if ReturnCode: print('Error: ReturnCode', ReturnCode)
    else: print('Done.')
    return(None)

