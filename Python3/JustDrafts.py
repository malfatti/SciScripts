# -*- coding: utf-8 -*-
"""
Just drafts
"""
#%%
for Key in UnitRec:    
    ClusterNo = len(UnitRec[Key]['Spks'])
    if ClusterNo == 0: print(Key, 'is lost'); continue
    
    Fig, Axes = plt.subplots(ClusterNo,2, figsize=(6, 3*ClusterNo))
    
    for Cluster in range(ClusterNo):
        SpkNo = len(UnitRec[Key]['Spks'][Cluster])
        print(str(SpkNo), 'Spks in cluster', str(Cluster))
        print('Max of', max(UnitRec[Key]['PSTH'][Cluster]),  'Spks in PSTH')
        
        if not SpkNo:
            print('No Spk data on cluster', str(Cluster) + '. Skipping...')
            Thrash[Key] = UnitRec[Key].copy()
            continue
        
        PSTHPeak = max(UnitRec[Key]['PSTH'][Cluster])
        PSTHMean = np.mean(UnitRec[Key]['PSTH'][Cluster])
#        if max(UnitRec[Key]['PSTH'][Cluster]) < 4: 
#            print('No peaks in PSTH. Skipping cluster', str(Cluster), '...')
#            continue
        
        if SpkNo > 100: SpkNo = 100
        
        for Spike in range(SpkNo):
            if ClusterNo == 1: Axes[0].plot(UnitRec[Key]['Spks'][Cluster][Spike], 'r')
            else: Axes[Cluster][0].plot(UnitRec[Key]['Spks'][Cluster][Spike], 'r')
        
        if ClusterNo == 1:
            Axes[0].set_title('Peak='+str(PSTHPeak)+' Mean='+str(PSTHMean))
            Axes[0].plot(np.mean(UnitRec[Key]['Spks'][Cluster], axis=0), 'k')
            Axes[1].bar(XValues, UnitRec[Key]['PSTH'][Cluster])
        else:
            Axes[Cluster][0].set_title('Peak='+str(PSTHPeak)+' Mean='+\
                                       str(PSTHMean)+' Std='+str(PSTHStd))
            Axes[Cluster][0].plot(np.mean(UnitRec[Key]['Spks'][Cluster]), 'k')
            Axes[Cluster][1].bar(XValues, UnitRec[Key]['PSTH'][Cluster])


#UnitRec = Units['Sound_NaCl']['00']['00']
##UnitRec = Units[Stim][FIndS][RecS]
#Thrash = {}; PSTHStd = []
#for Key in UnitRec:    
#    ClusterNo = len(UnitRec[Key]['Spks'])
#    if ClusterNo == 0: print(Key, 'is lost'); continue
#    
#    Fig, Axes = plt.subplots(ClusterNo,2, figsize=(6, 3*ClusterNo))
#    
#    for Cluster in range(ClusterNo):
#        SpkNo = len(UnitRec[Key]['Spks'][Cluster])
#        print(str(SpkNo), 'Spks in cluster', str(Cluster))
#        print('Max of', max(UnitRec[Key]['PSTH'][Cluster]),  'Spks in PSTH')
#        
#        if not SpkNo:
#            print('No Spk data on cluster', str(Cluster) + '. Skipping...')
#            Thrash[Key] = UnitRec[Key].copy()
#            continue
#        
#        if max(UnitRec[Key]['PSTH'][Cluster]) < 4: 
#%%
F = h5py.File(FileName)
for key in F['ExpInfo'].keys():
    N = "{0:02d}".format(int(key))
    F['ExpInfo'][N] = F['ExpInfo'][key]
    del(F['ExpInfo'][key])

for ind, key in enumerate(F['ExpInfo'].keys()):
    if ind <5: F['ExpInfo'][key].attrs['StimType'] = np.string_(['Sound_NaCl'])
    else: F['ExpInfo'][key].attrs['StimType'] = np.string_(['Sound_CNO'])
F.close()



#%%
#import h5py
#from glob import glob

Files = glob('*/*.hdf5'); Files.sort()
AnalysisFile = Files[0][:14] + '-Analysis.hdf5'

Analysis = h5py.File(AnalysisFile)

for Ind, File in enumerate(Files):
    F = h5py.File(File)
    for Key in F.keys():
        for Exp in ['ABRs-', 'GPIAS']:
            if Exp in Key:
                F.copy('/'+Key, Analysis)
                del(F[Key])
    F.close()

print(list(Analysis.keys()))
Analysis.close()


#%%
def GPIASAAA(FileName, GPIASTimeBeforeTTL=50, GPIASTimeAfterTTL=150, FilterLow=3, 
          FilterHigh=300, FilterOrder=4, GPIASTTLCh=2, PiezoCh=1):
    
    print('set paths...')
    os.makedirs('Figs', exist_ok=True)    # Figs folder
    DirList = glob.glob('KwikFiles/*'); DirList.sort()
    
    for RecFolder in DirList:
        print('Load files...')
        FilesList = glob.glob(''.join([RecFolder, '/*']))
        Files = {}
        for File in FilesList:
            if '.kwd' in File:
                try:
                    Raw = Kwik.load(File, 'all')
                    Files['kwd'] = File
                except OSError:
                        print('File ', File, " is corrupted :'(")
                
            elif '.kwe' in File:
                try:
                    Events = Kwik.load(File)
                    Files['kwe'] = File
                except OSError:
                    print('File ', File, " is corrupted :'(")
                
            elif '.kwik' in File:
                try:
                    Events = Kwik.load(File)
                    Files['kwik'] = File
                except OSError:
                    print('File ', File, " is corrupted :'(")
            
#            elif '.db' in File:
#                with shelve.open(File[:-3]) as Shelve: 
#                    DataInfo = Shelve['DataInfo']
        
        DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
        DataInfo['SoundBackgroundAmpF'] = Hdf5F.LoadDict(
                                             '/DataInfo/SoundBackgroundAmpF', 
                                             FileName, Attrs=False)
        DataInfo['SoundPulseAmpF'] = Hdf5F.LoadDict('/DataInfo/SoundPulseAmpF', 
                                                    FileName, Attrs=False)
        
        print('Check if files are ok...')
        if 'Raw' not in locals():
            print('.kwd file is corrupted. Skipping dataset...')
            continue
        
        if 'Events' not in locals():
            print('.kwe/.kwik file is corrupted. Skipping dataset...')
            continue
        
        if 'DataInfo' not in locals():
            print('No data info. Skipping dataset...')
            continue
        
        print('Data from ', RecFolder, ' loaded.')
        
        print('Preallocate memory...')
        GPIAS = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
        for Freq in range(len(DataInfo['NoiseFrequency'])):
            GPIAS[Freq] = [[0] for _ in range(round(DataInfo['NoOfTrials']*2))]
        
#        AllTTLs = GPIAS.copy()
        
        print('Get TTL data...')
        EventID = Events['TTLs']['user_data']['eventID']
        EventCh = Events['TTLs']['user_data']['event_channels']
        EventRec = Events['TTLs']['recording']
        EventSample = Events['TTLs']['time_samples']
        
        TTLChs = np.nonzero(np.bincount(EventCh))[0]
        TTLRecs = np.nonzero(np.bincount(EventRec))[0]
        TTLsPerRec = {_Rec: [EventSample[_] for _ in range(len(EventRec)) 
                             if EventRec[_] == _Rec 
                             and EventCh[_] ==  GPIASTTLCh-1 
                             and EventID[_] == 1]
                      for _Rec in TTLRecs}
        TTLsDPerRec = {_Rec: [EventSample[_] for _ in range(len(EventRec)) 
                             if EventRec[_] == _Rec 
                             and EventCh[_] ==  GPIASTTLCh-1 
                             and EventID[_] == 0]
                      for _Rec in TTLRecs}
        TTLRising = Kwik.get_rising_edge_times(Files['kwe'], GPIASTTLCh-1)
        
        Rate = Raw['info']['0']['sample_rate']
        NoOfSamplesBefore = int(round((GPIASTimeBeforeTTL*Rate)*10**-3))
        NoOfSamplesAfter = int(round((GPIASTimeAfterTTL*Rate)*10**-3))
        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
        XValues = (range(-NoOfSamplesBefore, 
                         NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
        
        print('Set filter...')
        passband = [FilterLow/(Rate/2), FilterHigh/(Rate/2)]
        f2, f1 = signal.butter(FilterOrder, passband, 'bandpass')
        
        for Rec in range(len(Raw['data'])):
            TTLNo = [0]
            for _ in range(1, len(TTLsPerRec)+1):
                TTLNo = TTLNo + [len(TTLsPerRec[_-1]) + TTLNo[-1]]

            TTLNo = [0] + [TTLNo[_]-1 for _ in range(1, len(TTLNo))]
            
            RawTime = list(Raw['timestamps'][str(Rec)])
            
            if Rec == 0:
                sTTLNo = 0
            else:
                sTTLNo = TTLNo[Rec] + 1
            
            print('Slicing and filtering Rec ', str(Rec), '...')
            for TTL in range(sTTLNo, TTLNo[Rec+1]):
                TTLLoc = RawTime.index(TTLRising[TTL])
                    
                Start = TTLLoc-NoOfSamplesBefore
                End = TTLLoc+NoOfSamplesAfter
                
                Index = TTL-sTTLNo
                
                gData = Raw['data'][str(Rec)][Start:End, PiezoCh-1]
#                gData = [float(gData[_]) for _ in range(NoOfSamples)]
                gData = abs(signal.hilbert(gData))
                gData = signal.filtfilt(f2, f1, gData, padtype='odd', padlen=0)
                
                Freq = DataInfo['FreqOrder'][Rec][0]; 
                Trial = DataInfo['FreqOrder'][Rec][1];
                GPIAS[Freq][Trial] = gData[:]
#                AllTTLs[Freq][Trial] = [TTLStart, TTLEnd]
                
                del(TTLLoc, Start, End, Freq, Trial, gData)
        
        for Freq in range(len(DataInfo['NoiseFrequency'])):
            gData = GPIAS[Freq][:]
            NoGapAll = [gData[_] for _ in range(len(gData)) if _%2 == 0]
            GapAll = [gData[_] for _ in range(len(gData)) if _%2 != 0]
            NoGapSum = list(map(sum, zip(*NoGapAll)))
            GapSum = list(map(sum, zip(*GapAll)))
            
            gData = [0, 0]
            gData[0] = [_/DataInfo['NoOfTrials'] for _ in NoGapSum]
            gData[1] = [_/DataInfo['NoOfTrials'] for _ in GapSum]
            gData[0] = signal.savgol_filter(gData[0], 5, 2, mode='nearest')
            gData[1] = signal.savgol_filter(gData[1], 5, 2, mode='nearest')
            GPIAS[Freq] = gData[:]
            
#            tData = AllTTLs[Freq][:]
#            TTLNoGapAll = [tData[_] for _ in range(len(tData)) if _%2 == 0]
#            TTLGapAll = [tData[_] for _ in range(len(tData)) if _%2 != 0]
#            TTLNoGapSum = list(map(sum, zip(*TTLNoGapAll)))
#            TTLGapSum = list(map(sum, zip(*TTLGapAll)))
#            
#            tData = [0, 0]
#            tData[0] = [int(round(_/DataInfo['NoOfTrials'])) for _ in TTLNoGapSum]
#            tData[1] = [int(round(_/DataInfo['NoOfTrials'])) for _ in TTLGapSum]
#            AllTTLs[Freq] = tData[:]
            
            del(NoGapAll, GapAll, NoGapSum, GapSum, gData)
#            del(TTLNoGapAll, TTLGapAll, TTLNoGapSum, TTLGapSum, tData)
            
        FileName = DataInfo['AnimalName'] + '-GPIAS-' + RecFolder[10:]
        print('Saving data to ' + FileName)
#        with shelve.open(FileName) as Shelve:
#            Shelve['GPIAS'] = GPIAS
#            Shelve['AllTTLs'] = AllTTLs
#            Shelve['XValues'] = XValues
#            Shelve['RecFolder'] = RecFolder
#            Shelve['DataInfo'] = DataInfo
        print('Done.')
    
    print('Finished.')
    return(None)

