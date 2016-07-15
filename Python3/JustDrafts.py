# -*- coding: utf-8 -*-
"""
Just drafts
"""
#%%
#Spk = Units['Sound_NaCl']['00']['00']['Ch05']['SpkWF'][:]
#Hist = Units['Sound_NaCl']['00']['00']['Ch05']['NoOfSpks'][:]
UnitRec = Test['Sound_NaCl']['02']['00']
#UnitRec = Units[Stim][FIndS][RecS]
Thrash = {}
for Key in UnitRec:    
    ClusterNo = len(UnitRec[Key]['SpkWF'])
    if ClusterNo == 0: print(Key, 'is lost'); continue
    
    Fig, Axes = plt.subplots(ClusterNo,2, figsize=(6, 3*ClusterNo))
    
    for Cluster in range(ClusterNo):
        SpkNo = len(UnitRec[Key]['SpkWF'][Cluster])
        print(str(SpkNo), 'Spks in cluster', str(Cluster))
        print('Max of', max(UnitRec[Key]['NoOfSpks'][Cluster]),  'Spks in PSTH')
        
        if max(UnitRec[Key]['NoOfSpks'][Cluster]) < 4: 
            print('No peaks in PSTH. Skipping cluster', str(Cluster), '...')
            continue
        
        if not SpkNo:
            print('No Spk data on cluster', str(Cluster) + '. Skipping...')
            Thrash[Key] = UnitRec[Key].copy()
            continue
        
        if SpkNo > 100: SpkNo = 100
        
        for Spike in range(SpkNo):
            if ClusterNo == 1: Axes[0].plot(UnitRec[Key]['SpkWF'][Cluster][Spike], 'r')
            else: Axes[Cluster][0].plot(UnitRec[Key]['SpkWF'][Cluster][Spike], 'r')
        
        if ClusterNo == 1:
            Axes[0].plot(np.mean(UnitRec[Key]['SpkWF'][Cluster], axis=0), 'k')
            Axes[1].bar(XValues, UnitRec[Key]['NoOfSpks'][Cluster])
        else:
            Axes[Cluster][0].plot(np.mean(UnitRec[Key]['SpkWF'][Cluster]), 'k')
            Axes[Cluster][1].bar(XValues, UnitRec[Key]['NoOfSpks'][Cluster])

for SKey in Units.keys():
            for FKey in Units[SKey].keys():
                for RKey in Units[SKey][FKey]:
                    for Ch in Units[SKey][FKey][RKey]:
                        for Key, Data in Units[SKey][FKey][RKey][Ch].items():
                            plt.figure()
                            
                            plt.plot(Data)

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



def UnitAnalysis(FileName, StimTTLCh=-1, PSTHTimeBeforeTTL=0, 
                 PSTHTimeAfterTTL=300, StimType=['Sound'], AnalogTTLs=False, 
                 Board='OE', OverrideRec=[]):
    print('Load DataInfo...')
    DirList = glob('KwikFiles/*'); DirList.sort()
    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
    
    CustomAdaptor = [5, 6, 7, 8, 9, 10 ,11, 12, 13, 14, 15, 16, 1, 2, 3, 4]
    A16 = {'ProbeTip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
           'ProbeHead': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}
    ChannelMap = GetProbeChOrder(A16['ProbeTip'], A16['ProbeHead'], CustomAdaptor)
    
    Units = {}
    for Stim in StimType:
        Exps = Hdf5F.LoadExpPerStim(Stim, DirList, FileName)
        Units[Stim] = {}
        
        if isinstance(StimTTLCh, dict): UnitTTLCh = StimTTLCh[Stim]
        else: UnitTTLCh = StimTTLCh
        
        for FInd, RecFolder in enumerate(Exps):        
            if AnalogTTLs: Raw, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
            else: Raw, Events, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
            
            OEProc = GetProc(Raw, Board)[0]
            
            Path = getcwd() + '/' + RecFolder +'/SepCh/'
            makedirs(Path, exist_ok=True)
            
            Rate = Raw[OEProc]['info']['0']['sample_rate']
            NoOfSamplesBefore = int(round((PSTHTimeBeforeTTL*Rate)*10**-3))
            NoOfSamplesAfter = int(round((PSTHTimeAfterTTL*Rate)*10**-3))
            NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
            XValues = (range(-NoOfSamplesBefore, 
                             NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
            
            FIndS = "{0:02d}".format(FInd)
            Units[Stim][FIndS] = {}
            for Rec in range(len(Raw[OEProc]['data'])):
                if OverrideRec != []: Rec = OverrideRec
                RecS = "{0:02d}".format(Rec)
                
                print('Separating channels according to ChannelMap...')
                Data = [Raw[OEProc]['data'][str(Rec)][:, _-1] * 
                        Raw[OEProc]['channel_bit_volts'][str(Rec)][_-1] * 1000
                        for _ in ChannelMap]
                
                TTLs = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, UnitTTLCh, 
                                          OEProc)
                
                print('Writing files for clustering... ', end='')
                FileList = []
                for Ind, Ch in enumerate(Data):
                    MatName = 'Exp' + Files['100_kwd'][-13:-8] + '_' + \
                              RecS + '-Ch' + "{0:02d}".format(Ind+1) + '.mat'
                    
                    FileList.append(MatName)
                    io.savemat(Path+MatName, {'data': Ch})
                
                TxtFile = open(Path+'Files.txt', 'w')
                for File in FileList: TxtFile.write(File+'\n')
                TxtFile.close()
                print('Done.')
                
                CallWaveClus(Rate, Path)
                ClusterList = glob(Path+'times_*'); ClusterList.sort()
                
                Units[Stim][FIndS][RecS] = {}
                print('Preparing histograms and spike waveforms...')
                for File in ClusterList:
                    Clusters = io.loadmat(File)
                    
                    ClusterClasses = np.unique(Clusters['cluster_class'][:,0])
                    Ch = File[-8:-4] 
                    
                    Units[Stim][FIndS][RecS][Ch] = {}
                    Units[Stim][FIndS][RecS][Ch]['NoOfSpks'] = [
                                        np.zeros(len(XValues)) 
                                        for _ in range(len(ClusterClasses))]
                    
                    Units[Stim][FIndS][RecS][Ch]['SpkWF'] = [
                                        [] for _ in range(len(ClusterClasses))]
                    
                    for Cluster in range(len(ClusterClasses)):
                        ClassIndex = Clusters['cluster_class'][:,0] == Cluster
                        if not len(Clusters['spikes'][ClassIndex,:]): continue
                        
                        for TTL in range(len(TTLs)):
                            Firing = Clusters['cluster_class'][ClassIndex, 1] \
                                     - TTLs[TTL]
                            Firing = Firing[(Firing >= XValues[0]) * 
                                            (Firing < XValues[-1])]
                            SpkCount = np.histogram(Firing, 
                                                    np.hstack((XValues, 300)))[0]
                            
                            Units[Stim][FIndS][RecS][Ch]['NoOfSpks'][Cluster] = \
                              Units[Stim][FIndS][RecS][Ch]['NoOfSpks'][Cluster] \
                              + SpkCount
                            
                            del(Firing, SpkCount)
                        
                        Units[Stim][FIndS][RecS][Ch]['SpkWF'][Cluster] = \
                                        Clusters['spikes'][ClassIndex,:]
                    
                    print(Ch, 'have', str(len(ClusterClasses)), 'clusters')
                    del(Clusters)
                
                del(Data, TTLs)
                print('Cleaning...')
                ToDelete = glob(Path+'*')
                for File in ToDelete: remove(File)
                if OverrideRec != []: break
            
            del(Raw)
            removedirs(Path)
    
    AnalysisFile = '../' + DataInfo['AnimalName'] + '-Analysis.hdf5'
    Hdf5F.WriteUnits(Units, XValues, AnalysisFile)
    return(None)


