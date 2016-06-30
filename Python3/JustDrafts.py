# -*- coding: utf-8 -*-
"""
Just drafts
"""

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


