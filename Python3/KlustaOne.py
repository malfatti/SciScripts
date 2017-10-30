#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 11:03:23 2017

@author: malfatti
"""
#%%
import os, subprocess
import numpy as np
from glob import glob
from imp import load_source
from itertools import tee


DataFolder = '/home/malfatti/NotSynced/Tmp/SergioData/juj006d03'
ExpName = DataFolder.split(os.sep)[-1]
KlustaFolder = os.sep.join([DataFolder, 'KlustaFiles'])
SpkTSFile = glob(os.sep.join([DataFolder, '*.spike']))[0]
SpkWFFile = glob(os.sep.join([DataFolder, '*.swave']))[0]
InfoFile = glob(os.sep.join([DataFolder, '*-info.txt']))[0]

# To parse from InfoFile
SpkChNo = 3
TrialNo = 40
SpkWFLen = 38
Rate = 32000


def DictPrint(value, htchar='    ', itemchar=' ', breaklineat='auto', lfchar='\n', indent=0):
    ''' Modified from y.petremann's code.
        Added options to set item separator for list or tuple and to set a number
        of items per line, or yet, to calculate items per line so it will not 
        have more than 80 chars per line.
        Source: https://stackoverflow.com/a/26209900 '''
    
    nlch = lfchar + htchar * (indent + 1)
    if type(value) is dict:
        items = [
            nlch + repr(key) + ': ' + DictPrint(value[key], htchar, itemchar, breaklineat, lfchar, indent + 1)
            for key in value
        ]
        
        return '{%s}' % (','.join(items) + lfchar + htchar * indent)
    
    elif type(value) is list or type(value) is tuple:
        items = [
            itemchar + DictPrint(item, htchar, itemchar, breaklineat, lfchar, indent + 1)
            for item in value
        ]
        
        if breaklineat == 'auto':
           bl = int((80 - (len(htchar)*(indent + 1)))/
                (int((sum([len(i)+4 for i in items])-len(itemchar)-1)/len(items))))
         
        else: bl = breaklineat
        
        if not bl: bl = 1
       
        if len(items) > bl:
            for i in list(range(bl, len(items), bl)):
                items[i] = lfchar + htchar*(indent+1) + '  ' + items[i]
        
        return '[%s]' % (','.join(items))
    
    elif type(value) is np.ndarray:
        value = value.tolist()
        items = DictPrint(value, htchar, itemchar, breaklineat, lfchar, indent)
        return items
    
    else:
        return repr(value)


def Pairwise(iterable):
    """ from https://docs.python.org/3.6/library/itertools.html#itertools-recipes
    s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def PrbWrite(File, Channels=list(range(3)), Spacing=100):
    Prb = {'0': {}}
    Prb['0']['channels'] = Channels
    Prb['0']['graph'] = list(Pairwise(Prb['0']['channels']))
    Pos = list(range(0, len(Prb['0']['channels'])*Spacing, Spacing))
    Prb['0']['geometry'] = {str(Ch):(0,Pos[C]) 
                            for C,Ch in enumerate(Prb['0']['channels'])}
    
    UserFile = os.sep.join(File.split(os.sep)[:-2]) + os.sep + ExpName + '-Klusta.prb'
    if os.path.isfile(UserFile):
        with open(UserFile, 'r') as F: UserPrb = load_source('UserPrb', '', F)
        Prb = {**Prb, **UserPrb.channel_groups}
    
    with open(File, 'w') as F: F.write('channel_groups = '+DictPrint(Prb))
    return(None)


def PrmWrite(File, experiment_name, prb_file, raw_data_files, sample_rate, 
             n_channels, dtype, SpkWFLen):
    
    traces = {
        'raw_data_files': raw_data_files,
        'sample_rate': sample_rate,
        'n_channels': n_channels,
        'dtype': dtype,
    }
    
    spikedetekt = {
        'filter_low': 300,
        'filter_high_factor': 0.5,  # will be multiplied by the sample rate
        'filter_butter_order': 3,
    
        # Data chunks.
        'chunk_size_seconds': 1.,
        'chunk_overlap_seconds': .015,
    
        # Threshold.
        'n_excerpts': 50,
        'excerpt_size_seconds': 1.,
        'use_single_threshold': True,
        'threshold_strong_std_factor': 4.5,
        'threshold_weak_std_factor': 2.,
        'detect_spikes': 'negative',
    
        # Connected components.
        'connected_component_join_size': 1,
    
        # Spike extractions.
        'extract_s_before': int(SpkWFLen/2),
        'extract_s_after': int(SpkWFLen/2),
        'weight_power': 2,
    
        # Features.
        'n_features_per_channel': 5,
        'pca_n_waveforms_max': 10000,
    
    }
    
    klustakwik2 = {
         'prior_point': 1,
         'mua_point': 2,
         'noise_point': 1,
         'points_for_cluster_mask': 100,
         'penalty_k': 0.0,
         'penalty_k_log_n': 1.0,
         'max_iterations': 1000,
         'num_starting_clusters': 500,
         'use_noise_cluster': True,
         'use_mua_cluster': True,
         'num_changed_threshold': 0.05,
         'dist_thresh': 9.210340371976184,
         'full_step_every': 1,
         'split_first': 20,
         'split_every': 40,
         'max_possible_clusters': 5*n_channels,
         'max_quick_step_candidates': 100000000, # this uses around 760 MB RAM
         'max_quick_step_candidates_fraction': 0.4,
         'always_split_bimodal': False,
         'subset_break_fraction': 0.01,
         'break_fraction': 0.0,
         'fast_split': False,
         'max_split_iterations': None,
         'consider_cluster_deletion': True,
         'num_cpus': 8,
    }
    
    UserFile = os.sep.join(File.split(os.sep)[:-2]) + os.sep + ExpName + '-Klusta.prm'
    if os.path.isfile(UserFile):
        with open(UserFile, 'r') as F: Prm = load_source('Prm', '', F)
        if 'experiment_name' in dir(Prm): experiment_name = Prm.experiment_name
        if 'experiment_name' in dir(Prm): experiment_name = Prm.experiment_name
        if 'traces' in dir(Prm): traces = {**traces, **Prm.traces}
        if 'spikedetekt' in dir(Prm): spikedetekt = {**spikedetekt, **Prm.spikedetekt}
        if 'klustakwik2' in dir(Prm): klustakwik2 = {**klustakwik2, **Prm.klustakwik2}
    
    with open(File, 'w') as F:
        F.write('experiment_name = "'+experiment_name+'"')
        F.write('\n\n')
        F.write('prb_file = "'+prb_file+'"')
        F.write('\n\n')
        # F.write('traces = '+DictPrint(traces, breaklineat=50))
        F.write('traces = '+DictPrint(traces))
        F.write('\n\n')
        F.write('spikedetekt = '+DictPrint(spikedetekt))
        F.write('\n\n')
        F.write('klustakwik2 = '+DictPrint(klustakwik2))
        F.write('\n\n')
    
    return(None)


def RunProcess(Cmd, LogFile=''):
    if LogFile == '': print('Logging disabled, outputting to STDOUT.')
    else: print('Check progress in file', LogFile)
    
    try:
        if LogFile == '': Log = subprocess.PIPE
        else:  Log = open(LogFile, 'w')
        
        P = subprocess.Popen(Cmd,
                             stdout=Log,
                             stderr=subprocess.STDOUT)
        
        print('Process id:', P.pid )
        P.communicate()[0]; ReturnCode = P.returncode
        if LogFile != '': Log.close()
    
    except Exception as e:
        ReturnCode = 1; print(e)
    
    return(ReturnCode)


def Run(PrmFile, Path, Overwrite=False):
    Here = os.getcwd(); os.chdir(Path)
    Klusta = os.environ['HOME']+'/Software/Miniconda3/envs/klusta/bin/klusta'
    # Klusta = KLUSTAPATH
    Cmd = [Klusta, PrmFile]
    if Overwrite: Cmd.append('--overwrite')
    
    print('Entering in directory', Path, '...')
    print('Clustering spikes...')
    ReturnCode = RunProcess(Cmd, PrmFile+'.log'); os.chdir(Here)
    print('Going back to', Here, '...')
    
    if ReturnCode: print('Error: ReturnCode', ReturnCode)
    else: print('Done clustering.')
    return(None)


SpkTSRaw = np.fromfile(SpkTSFile, dtype='>i4')
SpkTS = abs(np.diff(SpkTSRaw))
SpkTS = np.where(SpkTS > (SpkTS.mean() + (SpkTS.std() * 3)))[0] + 1
SpkTS = [S[1:] for S in np.split(SpkTSRaw, SpkTS)]

SpkTimestamps = [SpkTS[S:S+SpkChNo] for S in range(0, TrialNo*SpkChNo, SpkChNo)]
SpkTS = [np.concatenate(Trial) for Trial in SpkTimestamps]

SpkWFRaw = np.fromfile(SpkWFFile, dtype='>i2')
SpkWF = [1]
for S in range(TrialNo*SpkChNo):
    Index = SpkWF[-1] + (SpkWFRaw[SpkWF[-1]]*SpkWFLen) + 4
    SpkWF.append(Index)

SpkWF = [S[1:] for S in np.split(SpkWFRaw, SpkWF)]
SpkWF = [S for S in SpkWF if len(S) > SpkWFLen]
SpkWF = [SpkWF[S:S+SpkChNo] for S in range(0, TrialNo*SpkChNo, SpkChNo)]
SpkWF = [[np.split(Ch[:-(len(Ch)%SpkWFLen)], (len(Ch)-len(Ch)%SpkWFLen)/SpkWFLen) for Ch in Trial] for Trial in SpkWF]

DataLen = [max(Trial)+(2*SpkWFLen) for Trial in SpkTS]
SpkWFArrays = [np.random.randint(-3, 4, size=(DataLen[T],SpkChNo)).astype(np.int16) for T in range(TrialNo)]
for T, Trial in enumerate(SpkTimestamps):
    for C, Ch in enumerate(Trial):
        for S, Spk in enumerate(Ch):
            Start = int(Spk - (SpkWFLen/2))
            End = int(Spk + (SpkWFLen/2))
            SpkWFArrays[T][Start:End, C] = SpkWF[T][C][S]


os.makedirs(KlustaFolder, exist_ok=True)
StrF = "{0:0"+str(len(str(TrialNo)))+"d}"
for T, Trial in enumerate(SpkWFArrays):
    DataFile = ExpName+'_Trial'+StrF.format(T)
    with open(KlustaFolder + os.sep + DataFile + '.dat', 'wb') as File: 
        File.write(Trial.tobytes())

FilesPrefix = KlustaFolder + os.sep + ExpName
raw_data_files = sorted(glob(FilesPrefix + '*.dat'))
PrbWrite(FilesPrefix+'.prb', Channels=list(range(SpkChNo)))
PrmWrite(FilesPrefix+'.prm', ExpName,  FilesPrefix+'.prb', raw_data_files, 
         Rate, SpkChNo, '<i2', SpkWFLen)

Run(ExpName+'.prm', KlustaFolder, Overwrite=True)

File = FilesPrefix+'.prm'
with open(File, 'r') as F: Prb = F.read()


# SpkNo = 0
# for Trial in SpkTimestamps:
#     for Ch in Trial:
#         SpkNo += len(Ch)

    
# Params = {'backend': 'TkAgg'}
# from matplotlib import rcParams; rcParams.update(Params)
# from matplotlib import pyplot as plt

## Template - remove completely when finished {
# from IO import Hdf5, Txt
# from DataAnalysis.DataAnalysis import NestedClean

# TestFile = '../../../NotSynced/Tmp/SergioData/KlustaTest/hybrid_10sec.kwx'
# PrbFile = '/'.join(TestFile.split('/')[:-1]) + '/A16-25.prb'
# KwikA, KwikB = Hdf5.DataLoad('/', TestFile)
# KwikA, KwikB = NestedClean(KwikA), NestedClean(KwikB)

# del(KwikB['application_data'])
# del(KwikB['channel_groups'])
# del(KwikB['recordings']['0']['raw'])

# with open(PrbFile, 'r') as F: PrbStr = F.read()
# Prb = Txt.literal_eval(PrbStr.split(' = ')[1])
## }

# for Trial in SpkWF: 
#     for Ch in Trial:
#         for Spk in Ch: plt.plot(Spk)
#         plt.show()

# from itertools import accumulate
# Offsets = [max(_) for _ in SpkTS]
# Offsets = list(accumulate(Offsets))
# TimeSamples = [El+Offsets[T] for T, Trial in enumerate(SpkTS) for El in Trial]
# TimeFractional = np.concatenate(SpkTS)

# Recs = [np.repeat(T, len(Trial)) for T, Trial in enumerate(SpkTS)]
# Recs = np.concatenate(Recs)

# KwikData = {
#     'channel_groups': {'0': {'spikes': {
#                                  'recording': Recs, 
#                                  'time_fractional': TimeFractional, 
#                                  'time_samples': TimeSamples
# }}}}

# KwikAttrs = {
#     'creator_version': np.array([b'klusta 3.0.16'], dtype='|S13'), 
#     'kwik_version': 2, 
#     'name': np.array([os.sep.join([KlustaFolder, ExpName]).encode()]), 
#     'recordings': {str(Key): {
#                        'name': np.array([''.join('recording_',str(Key)).encode()]), 
#                        'sample_rate': Rate}
#                    for Key in range(TrialNo)}
# }

# KwxAttrs = {
#     'creator_version': np.array([b'klusta 3.0.16'], dtype='|S13'), 
#     'kwik_version': 2, 
#     'name': np.array([os.sep.join([KlustaFolder, ExpName]).encode()]), 
#     'recordings': {str(Key): {
#                        'name': np.array([''.join('recording_',str(Key)).encode()]), 
#                        'sample_rate': Rate}
#                    for Key in range(TrialNo)}
# }
