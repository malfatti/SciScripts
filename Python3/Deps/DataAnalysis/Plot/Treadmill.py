#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:06:47 2017

@author: malfatti
"""
from Plot.Plot import Set, Spectrogram


## Level 0
def AllChs(Treadmill, FigName=None, Ext='svg', Save=False, 
                     Visible=True):
    Params = Set(Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    if 'V' in Treadmill.keys(): ChNo = len(Treadmill) - 1
    else: ChNo = len(Treadmill)
    
    Fig, SxxAx = plt.subplots(ChNo,1,figsize=(8, 3*ChNo))
    for C, Ch in Treadmill.items():
        if C == 'V': continue
        
        T, F, Sxx, VMeans = Ch['T'], Ch['F'], Ch['Sxx'], Ch['VMeans']
        
        Spectrogram(SxxAx[int(C)-1], T, F, Sxx, HighFreqThr=100, 
                         Line=True, LineX=T, LineY=VMeans, LineColor='r', 
                         LineYLabel='Speed [m/s]')
    
    if not FigName: FigName = 'Spectrogram.' + Ext
    if Save: Fig.savefig(FigName, format=Ext)
    if Visible: plt.show()
    Fig.clf(); plt.close()


def AllPerCh(T, F, Sxx, SxxMaxs, SxxPerV, VMeans, VMeansSorted, 
                       VInds, FigName=None, Ext='svg', Save=False, Visible=True):
    Params = Set(Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    SAx = [0,0]
    Fig, SxxAx = plt.subplots(2,1,figsize=(10, 8))
    Fig.subplots_adjust(bottom=0.15, right=0.85)
    
    SxxAx[0].pcolormesh(T, F, Sxx, cmap='inferno'); SxxAx[0].set_ylim(0,100)
    SxxAx[0].yaxis.set_ticks_position('left')
    SxxAx[0].set_xlabel('Time [s]'); SxxAx[0].set_ylabel('Frequency [Hz]')
    SxxAx[0].spines['left'].set_position(('outward', 5))
    
    SAx[0] = SxxAx[0].twinx()
    SAx[0].plot(T, SxxMaxs, 'r', lw=2); 
    SAx[0].set_ylim(-max(SxxMaxs)/3, max(SxxMaxs))
    SAx[0].spines['right'].set_position(('outward', 55))
    SAx[0].yaxis.set_label_position('right'); SAx[0].yaxis.set_ticks_position('right')
    SAx[0].set_ylabel('ThetaPxx/DeltaPxx', color='red')
    
    VAx = Fig.add_axes(SxxAx[0].get_position()); VAx.patch.set_visible(False)
    VAx.yaxis.set_label_position('right'); VAx.yaxis.set_ticks_position('right')
    VAx.spines['bottom'].set_visible(False)
    
    VAx.plot(T, VMeans, 'g', lw=2)
    VAx.set_ylabel('Mean speed [m/s]', color='green')
    VAx.spines['right'].set_position(('outward', 5))
    
    SxxAx[1].pcolormesh(VMeansSorted, F, SxxPerV, cmap='inferno'); SxxAx[1].set_ylim(0,100)
    SxxAx[1].yaxis.set_ticks_position('left'); SxxAx[1].xaxis.set_ticks_position('bottom')
    SxxAx[1].set_xlabel('Velocity [m/s]', color='k')
    SxxAx[1].set_ylabel('Frequency [Hz]', color='k')
    
    SAx[1] = SxxAx[1].twinx()
    SAx[1].plot(VMeansSorted, SxxMaxs[VInds], 'r'); 
    SAx[1].set_ylim(-max(SxxMaxs)/3, max(SxxMaxs))
    SAx[1].spines['right'].set_position(('outward', 5))
    SAx[1].xaxis.set_ticks_position('bottom')
    SAx[1].set_ylabel('ThetaPxx/DeltaPxx', color='r')
    
    if not FigName: FigName = 'SpectrogramsPerCh.' + Ext
    if Save: Fig.savefig(FigName, format=Ext)
    if Visible: plt.show()
    Fig.clf(); plt.close()


def ChPair(Ch1, Ch2, Ch1Title='', Ch2Title='', FigName=None, 
                     Ext='svg', Save=False, Visible=True):
    Params = Set(Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    Titles = [Ch1Title, Ch2Title]
    Fig, SxxAx = plt.subplots(2,1,figsize=(8, 3*4))
    for C, Ch in enumerate([Ch1, Ch2]):
        T, F, Sxx, VMeans = Ch['T'], Ch['F'], Ch['Sxx'], Ch['VMeans']
        
        Spectrogram(SxxAx[C], T, F, Sxx, HighFreqThr=100, 
                         Line=True, LineX=T, LineY=VMeans, LineColor='r', 
                         LineYLabel='Speed [m/s]')
        
        SxxAx[C].set_title(Titles[C])
    
    if not FigName: FigName = 'Spectrogram.' + Ext
    if Save: Fig.savefig(FigName, format=Ext)
    if Visible: plt.show()
    Fig.clf(); plt.close()


def SepChs(Treadmill, FigName=None, Ext='svg', Save=False, 
                     Visible=True):
    Params = Set(Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    for C, Ch in Treadmill.items():
        if C == 'V': continue
        
        Fig, SxxAx = plt.subplots(1,1,figsize=(8, 3))
        
        T, F, Sxx, VMeans = Ch['T'], Ch['F'], Ch['Sxx'], Ch['VMeans']
        
        Spectrogram(SxxAx, T, F, Sxx, HighFreqThr=100, 
                         Line=True, LineX=T, LineY=VMeans, LineColor='r', 
                         LineYLabel='Speed [m/s]')
    
    if not FigName: FigName = 'Spectrogram.' + Ext
    if Save: Fig.savefig(FigName, format=Ext)
    if Visible: plt.show()
    else: plt.close()

