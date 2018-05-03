# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@year: 2015
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts

This is a script to generate a white noise sound at several frequencies and
intensities, play and record them at the same time.
Next, it calculates the intensity in RMS and dB for each frequency at each
amplification factor. In our setup, we use this script to calibrate the audio
equipment.
"""
#%% Import
import os, pandas
import numpy as np
import sounddevice as SD

from DataAnalysis.DataAnalysis import SignalIntensity
from IO import Hdf5, SigGen
from datetime import datetime
from time import sleep

Params = {'backend': 'Qt5Agg'}
from matplotlib import rcParams; rcParams.update(Params)
import matplotlib.pyplot as plt


# \usepackage{setspace}
# \usepackage{titlesec}
# \titleformat{\chapter}{\normalfont\LARGE\bfseries}{\thechapter.}{1em}{}
# \cleardoublepage
def WriteTexTable(SoundIntensity, Date, Folder):
    Freqs = list(SoundIntensity.keys())
    AmpFList = sorted(list(SoundIntensity[Freqs[0]].keys()), reverse=True)
    TexTable = pandas.DataFrame([[float(AmpF)] + [SoundIntensity[Freq][AmpF]['dB'] 
                                              for Freq in SoundIntensity]
                             for AmpF in AmpFList])
    
    with open(Folder+'/'+'IntensityTable_'+Date+'.tex', 'w') as File:
        File.write(
        r"""%% Configs =====
\documentclass[12pt,a4paper]{report}
\usepackage[left=0.5cm,right=0.5cm,top=0.5cm,bottom=0.5cm]{geometry}
\usepackage{longtable}
% Document ======
\begin{document}
\section{Sound measurements}
"""     )
        File.write(TexTable.to_latex(longtable=True))
        File.write(r"""
\end{document}
"""     )
    
    return(None)


def Intensity_AmpF(SoundIntensity, Date, Folder='.', Ext=['svg'], Save=True, Show=True):
    Colors = ['r', 'g', 'b', 'm', 'k', '#ffa500', '#00b2b2']
    Freqs = list(SoundIntensity.keys())
    Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
    
    plt.figure(figsize=(6, 6))
    for FKey in SoundIntensity:
        FInd = Freqs.index(FKey)
        AmpFList = sorted(list(SoundIntensity[FKey].keys()), reverse=True)
       
        FigTitle = 'Intensity curve per Freq per AmpF'
        YLabel = 'Intensity [dBSPL]'
        XLabel = 'Sound amplification factor'
        LineLabel = FKey + ' Hz'
        
        plt.semilogx(AmpFList, [SoundIntensity[FKey][AmpF]['dB'] for AmpF in AmpFList], 
                     label=LineLabel, color=Colors[FInd])
    
    plt.ylabel(YLabel); plt.xlabel(XLabel)
    plt.legend(loc='upper left')
    plt.locator_params(tight=True)
    plt.tick_params(direction='out')
    plt.axes().spines['right'].set_visible(False)
    plt.axes().spines['top'].set_visible(False)
    plt.axes().yaxis.set_ticks_position('left')
    plt.axes().xaxis.set_ticks_position('bottom')
    plt.title(FigTitle)
    
    if Save:
        for E in Ext:
            plt.savefig(Folder+'/'+'SoundMeasurement_'+Date+'.'+E, format=E)
    
    if Show: plt.show()
    else: plt.close()
    
    return(None)


def PxxSp_Freq(SoundIntensity, PSD, Date, Group, FileName, Folder='.', Ext=['svg'], Save=True, Show=True):
    Colormaps = [plt.get_cmap(Map) for Map in ['Reds', 'Greens', 'Blues', 
                                               'Purples', 'Greys', 'Oranges', 'GnBu']]
    
    Freqs = list(SoundIntensity.keys())
    Freqs.sort(key=lambda x: [int(y) for y in x.split('-')[-1]])
    
    Intensities = [80, 70, 60, 50, 40, 30, 0]
    AmpFList = SigGen.dBToAmpF(Intensities, Group, FileName)
    for AFKey in AmpFList:
        AmpFList[AFKey] = [str(_) for _ in AmpFList[AFKey]]
    
    Fig, Axes = plt.subplots(len(Freqs), sharex=True, figsize=(6, 12))
    for FKey in Freqs:
        Freq = Freqs.index(FKey)
        
        FigTitle = 'Power spectral density'
        AxTitle = FKey + ' Hz'
        YLabel = 'PSD [V²/Hz]'
        XLabel = 'Frequency [Hz]'
        
        # AmpFList = {FKey: sorted(list(SoundIntensity[FKey].keys()), reverse=True)} ## Override to delete
        for AKey in AmpFList[FKey]:
            AmpF = AmpFList[FKey].index(AKey)
            
            Colors = [Colormaps[Freq](255-(AmpF*255//len(AmpFList)))]
            LineLabel = str(round(SoundIntensity[FKey][AKey]['dB'])) + ' dB'
            # try:
            #     LineLabel = str(round(SoundIntensity[FKey][AKey]['dB'])) + ' dB'
            # except KeyError:
            #     AKey = '0.0'
            #     LineLabel = str(round(SoundIntensity[FKey][AKey]['dB'])) + ' dB'
            
            for IAmpF in SoundIntensity[FKey]:
                IdB = SoundIntensity[FKey][IAmpF]['dB']
                
                if SoundIntensity[FKey][AKey]['dB'] == IdB:
    #                print(IAmpF)
                    Axes[Freq].semilogy(PSD[FKey][IAmpF][0], 
                                        PSD[FKey][IAmpF][1], 
                                        label=LineLabel, color=Colors[0])
                else: continue
            
        Axes[Freq].set_xlim(left=5000, right=20000)
        Axes[Freq].set_ylabel(YLabel)
        Axes[Freq].spines['right'].set_visible(False)
        Axes[Freq].spines['top'].set_visible(False)
        Axes[Freq].yaxis.set_ticks_position('left')
        Axes[Freq].xaxis.set_ticks_position('bottom')
        Axes[Freq].set_title(AxTitle)
    #    if Freq < 2:
    #        Axes[Freq].legend(loc='lower right')
    #    else:
    #        Axes[Freq].legend(loc='lower left')
    Axes[-1].set_xlabel(XLabel)
    Fig.suptitle(FigTitle)
    Fig.tight_layout()
    Fig.subplots_adjust(top=0.93)
    
    if Save:
        for E in Ext:
            Fig.savefig(Folder+'/'+'PowerSpectrumDensity_'+Date+'.'+E, format=E)
    
    if Show: plt.show()
    else: plt.close()
    
    return(None)


## Set parameters of the experiment
SoundSystem = 'Jack-IntelOut-Marantz-IntelIn'
Setup = 'UnitRec'
SBOutAmpF = Hdf5.DataLoad(SoundSystem+'/SBOutAmpF', SigGen.CalibrationFile)[0]
OutMax = 1/SBOutAmpF

## Sound (Durations in sec)
Rate = 192000
SoundPulseDur = 2
NoiseFrequency = [[8000, 10000], [9000, 11000]] # Override
# NoiseFrequency = [[8000, 10000], [9000, 11000], [10000, 12000], [12000, 14000], 
#                   [14000, 16000], [16000, 18000], [8000, 18000]]

TTLAmpF = 0
# Mic sensitivity, from mic datasheet, in dB re V/Pa or in V/Pa
MicSens_dB = -47.46
#MicSens_VPa = 0.0042364

Folder = os.environ['DATAPATH']+'/Tests/SoundMeasurements'
FileName = Folder + '/' + 'SoundMeasurements.hdf5'
# FileName = 'Test.hdf5'
Group = '/'.join([SoundSystem, Setup])

os.makedirs(Folder, exist_ok=True)

SoundAmpF = [OutMax, 0.4, 0.3] # Override
# SoundAmpF = np.hstack((
#                 np.flipud(np.logspace(np.log10(1e-4), np.log10(OutMax), 299)),
#                 np.array(0.0)
#             ))
# SoundAmpF = np.hstack((
#                 np.flipud(np.logspace(np.log10(1e-6), np.log10(1e-4), 20)),
#                 np.array(0.0)
#             ))
SoundAmpF = np.array([round(_,6) for _ in SoundAmpF])

## Prepare dict w/ experimental setup
Date = datetime.now().strftime("%Y%m%d%H%M%S")
DataInfo = dict((Name, eval(Name)) 
                for Name in ['Rate', 'SoundPulseDur', 'Date','SoundAmpF', 
                             'NoiseFrequency', 'TTLAmpF', 'SoundSystem', 
                             'Setup', 'MicSens_dB', 'Folder'])

Freqs = [str(Freq[0]) + '-' + str(Freq[1]) for Freq in NoiseFrequency]
SoundAmpF = {Freq: SoundAmpF for Freq in Freqs}

# Prepare audio objects
SD.default.device = 'system'
SD.default.samplerate = Rate
SD.default.blocksize = 384
SD.default.channels = 2

# Warn user
FullTime = (len(SoundAmpF[Freqs[0]])*len(NoiseFrequency)*(SoundPulseDur))/60
FullTime = str(round(FullTime, 2))
print('')
print('Full test will take', FullTime, 'min to run.')
print('Current time: ', datetime.now().strftime("%H:%M:%S"))
input('Press any key to start... ')
print('This can be loud - cover your ears!!')
print('')
for i in range(5, 0, -1): print(i, end=' '); sleep(1)
print('')

print('Sound measurement running...')
SoundRec = {}
for Freq in NoiseFrequency:
    FKey = str(Freq[0]) + '-' + str(Freq[1])
    SoundRec[FKey] = {}
    Sound = SigGen.SoundStim(Rate, SoundPulseDur, SoundAmpF, 
                                        [Freq], TTLAmpF, SoundSystem, 
                                        TTLs=False, Map=[1,2])
    
    for AKey in Sound[FKey]:
        print(FKey, AKey)
        SoundRec[FKey][AKey] = SD.playrec(Sound[FKey][AKey], blocking=True, channels=1)
        SoundRec[FKey][AKey] = SoundRec[FKey][AKey][2000:] # Temp override
    
    print('Done playing/recording', FKey + '.')
    Hdf5.DataWrite(SoundRec[FKey], Group+'/SoundRec/'+FKey, FileName, Overwrite=True)
    del(Sound, SoundRec[FKey])

print('Finished recording \O/')


## Analysis
SoundRec = Hdf5.DataLoad(Group+'/SoundRec', FileName)[0]
SBInAmpF = Hdf5.DataLoad(SoundSystem+'/SBInAmpF', SigGen.CalibrationFile)[0]
# NoiseRMS = Hdf5.DataLoad(SoundSystem+'/NoiseLevel', SigGen.CalibrationFile)[0]
# Noise *= SBInAmpF

if 'MicSens_dB' in DataInfo:
    DataInfo['MicSens_VPa'] = 10**(DataInfo['MicSens_dB']/20)
    dB = True
#elif 'MicSens_VPa' not in DataInfo:
#    print('No mic sensitivity info. Impossible to calculate dBSPL values.')
#    dB = False

print('Calculating PSD, RMS and dBSLP...')
SoundIntensity = {}; PSD = {}
for F, Freq in SoundRec.items():
    SoundIntensity[F] = {}; PSD[F] = {}
    FreqBand = [int(_) for _ in F.split('-')]
    
    # NoiseRMS = SignalIntensity(Noise, DataInfo['Rate'], FreqBand, DataInfo['MicSens_VPa'])
    # NoiseRMS = NoiseRMS['RMS']
    
    for A, AmpF in Freq.items():
        AmpF = AmpF * SBInAmpF
        SoundIntensity[F][A], PSD[F][A] = SignalIntensity(AmpF[:,0], 
                                                          DataInfo['Rate'], 
                                                          FreqBand, 
                                                          DataInfo['MicSens_VPa'])#,
                                                     # NoiseRMS)
        SoundRec[F][A] = None

del(SoundRec)

## Save analyzed data
print('Saving analyzed data...')
os.makedirs(Folder, exist_ok=True)
Hdf5.DataWrite(SoundIntensity, Group+'/SoundIntensity', FileName, Overwrite=True)
Hdf5.DataWrite(PSD, Group+'/PSD', FileName, Overwrite=True)

WriteTexTable(SoundIntensity, DataInfo['Date'], Folder)
Intensity_AmpF(SoundIntensity, DataInfo['Date'], Folder)
PxxSp_Freq(SoundIntensity, PSD, DataInfo['Date'], Group, FileName, Folder)
print('Done.')