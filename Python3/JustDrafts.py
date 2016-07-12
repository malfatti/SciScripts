# -*- coding: utf-8 -*-
"""
Just drafts
"""
F = h5py.File(FileName)
for key in F['ExpInfo'].keys():
    N = "{0:02d}".format(int(key))
    F['ExpInfo'][N] = F['ExpInfo'][key]
    del(F['ExpInfo'][key])

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





def ABRAnalogTTLs(FileName, ABRCh=[1, 16], ABRTimeBeforeTTL=0, ABRTimeAfterTTL=12, 
        ABRTTLCh=1, FilterLow=300, FilterHigh=3000, FilterOrder=4, 
        StimType='Sound'):
    """
    Analyze ABRs from data recorded with OpenEphys. A '*ABRs.hdf5' file will be 
    saved in cwd, containing:
        - ABRs dict, where data will be saved as 
          ABRs[Ear][Freq][AmpF][DVCoord][Trial], where:
              Ear = 0 (right) or 1 (left)
              Freq = index of DataInfo['NoiseFrequency']
              AmpF = index of DataInfo['SoundAmpF']['Freq']
              DVCoord = string with DV coordinate at the moment of recording
              Trial = Trial number - 1 (so if it is one trial, Trial=0)
              
        - XValues array, for x axis of the plot;
        
        - DataInfo dict, where all info will be saved.
    
    For this function to work:
        - The Kwik folders must be in 'KwikFiles/';
        - There must be a *Exp.hdf5 file containing all experimental settings 
          (see Python3/SoundBoardControl/SoundAndLaserStimulation.py, 1st cell);
    """
    print('set paths...')
    os.makedirs('Figs', exist_ok=True)    # Figs folder
    DirList = glob.glob('KwikFiles/*'); DirList.sort()
    
    print('Load DataInfo...')
    DataInfo, Exps = Hdf5F.ExpDataInfo(FileName, DirList, StimType, 
                                               Var='Exps')
    
    print('Preallocate memory...')
    ABRs = [[], []]
    ABRs[0] = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
    ABRs[1] = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
    for Freq in range(len(DataInfo['NoiseFrequency'])):
        Key = str(DataInfo['NoiseFrequency'][Freq][0]) + '-' \
              + str(DataInfo['NoiseFrequency'][Freq][1])
        ABRs[0][Freq] = [{} for _ in range(len(DataInfo['SoundAmpF'][Key]))]
        ABRs[1][Freq] = [{} for _ in range(len(DataInfo['SoundAmpF'][Key]))]
    del(Freq)
    
    for RecFolder in Exps:
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
        
        ExpInfo = Hdf5F.ExpExpInfo(FileName, RecFolder, DirList)
        
        print('Check if files are ok...')
        if 'Raw' not in locals():
            print('.kwd file is corrupted. Skipping dataset...')
            continue
        
        print('Data from ', RecFolder, ' loaded.')
        
        if '0' not in list(Raw['data'].keys()):
            print('Rec numbers are wrong. Fixing...')
            for iKey in Raw.keys():
                Recs = list(Raw[iKey].keys())
                Recs = [int(_) for _ in Recs]; Min = min(Recs)
                
                for Key in Recs:
                    Raw[iKey][str(Key-Min)] = Raw[iKey].pop(str(Key))
                
            print('Fixed.')
        
        Rate = Raw['info']['0']['sample_rate']

        NoOfSamplesBefore = ABRTimeBeforeTTL*int(Rate*10**-3)
        NoOfSamplesAfter = ABRTimeAfterTTL*int(Rate*10**-3)
        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
        XValues = list((range((NoOfSamplesBefore)*-1, 0)/Rate)*10**3) + \
                  list((range(NoOfSamplesAfter)/Rate)*10**3)
        
        Lf2, Lf1 = signal.butter(FilterOrder, FilterLow/(Rate/2), 'highpass')
        Hf2, Hf1 = signal.butter(FilterOrder, FilterHigh/(Rate/2), 'lowpass')
        
        for Rec in range(len(Raw['data'])):
#            RecChNo = Raw['data'][str(Rec)].shape[1]
#            TTLCh = Raw['data'][str(Rec)][:, ABRTTLCh + (RecChNo-9)]
            TTLCh = Raw['data'][str(Rec)][:, -1]
            Threshold = max(TTLCh)/5
            TTLs = []
            for _ in range(1, len(TTLCh)):
                if TTLCh[_] > Threshold:
                    if TTLCh[_-1] < Threshold:
                        TTLs.append(_)
            
            rABR = [[0 for _ in range(NoOfSamples)] 
                    for _ in range(len(TTLs))]
            lABR = [[0 for _ in range(NoOfSamples)] 
                    for _ in range(len(TTLs))]
            sTTL = [[0 for _ in range(NoOfSamples)] 
                    for _ in range(len(TTLs))]
            
            print('Slicing and filtering ABRs Rec ', str(Rec), '...')
            for TTL in range(len(TTLs)):
                TTLLoc = int(TTLs[TTL])
                Start = TTLLoc-NoOfSamplesBefore
                End = TTLLoc+NoOfSamplesAfter
                
                rABR[TTL] = Raw['data'][str(Rec)][Start:End, ABRCh[0]-1] * \
                            Raw['channel_bit_volts'][str(Rec)][ABRCh[0]-1]
                
#                lABR[TTL] = Raw['data'][str(Rec)][Start:End, ABRCh[1]-1] * \
#                            Raw['channel_bit_volts'][str(Rec)][ABRCh[1]-1]
                
                rABR[TTL] = signal.filtfilt(Lf2, Lf1, rABR[TTL], padlen=0)
#                lABR[TTL] = signal.filtfilt(Lf2, Lf1, lABR[TTL], padlen=0)
                
                del(TTLLoc, Start, End)
            
            # Mean
            try:
                rABR = np.mean(rABR, axis=0)
            except ValueError:
                print('Last TTL was lost :(')
                rABR[-1] = rABR[-2]
                rABR = np.mean(rABR, axis=0)
            lABR = np.mean(lABR, axis=0)
#            rData = np.mean(rABR, axis=0)
#            lData = np.mean(lABR, axis=0)
            
            rABR = signal.filtfilt(Hf2, Hf1, rABR, padlen=0)
#            lABR = signal.filtfilt(Hf2, Hf1, lABR, padlen=0)
#            rData = signal.filtfilt(Hf2, Hf1, rData, padlen=0)
#            lData = signal.filtfilt(Hf2, Hf1, lData, padlen=0)
            
            if ExpInfo['DVCoord'] not in ABRs[0][ExpInfo['Hz']][Rec]:
                ABRs[0][ExpInfo['Hz']][Rec][ExpInfo['DVCoord']] = [rABR]
                ABRs[1][ExpInfo['Hz']][Rec][ExpInfo['DVCoord']] = [lABR]
            else:
                ABRs[0][ExpInfo['Hz']][Rec][ExpInfo['DVCoord']].append(rABR)
                ABRs[1][ExpInfo['Hz']][Rec][ExpInfo['DVCoord']].append(lABR)
    
    if 'XValues' in locals():
        print('Saving data to ' + FileName)
        GroupName = 'ABRs-' + datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        with h5py.File(FileName) as F:
            if 'ABRs' in F.keys():
                F[GroupName] = F['ABRs']; del(F['ABRs'])
            
            F.create_group('ABRs')
            F['ABRs'].attrs['XValues'] = XValues
            F['ABRs'].create_group('Right'); F['ABRs'].create_group('Left')
                
            for Freq in range(len(ABRs[0])):
                for AmpF in range(len(ABRs[0][Freq])):
                    for DV in ABRs[0][Freq][AmpF].keys():
                        Path = str(Freq) + '/' + \
                               str(AmpF) + '/' + \
                               DV
                        
                        F['ABRs']['Right'].create_group(Path)
                        F['ABRs']['Left'].create_group(Path)
                        del(Path)
                        
                        for Trial in range(len(ABRs[0][Freq][AmpF][DV])):
                            F['ABRs']['Right'][str(Freq)][str(AmpF)][DV][str(Trial)] = \
                                ABRs[0][Freq][AmpF][DV][Trial][:]
                            F['ABRs']['Left'][str(Freq)][str(AmpF)][DV][str(Trial)] = \
                                ABRs[1][Freq][AmpF][DV][Trial][:]
    
    print('Done.')
    return(None)


