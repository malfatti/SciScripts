#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 19:33:26 2018

@author: malfatti
"""
#%% Reanalysis again
from DataAnalysis import DataAnalysis
from glob import glob
from IO import Bin, Intan, Klusta, Txt


def ClusterizeSpks(Animals, ProbeSpacing, Log=None):
    if not Log: Log = {K: [] for K in ['Done', 'Errors', 'ErrorsLog', 'Skipped']}
    
    for A, Animal in Animals.items():
        try:
            for Es, Exps in enumerate(Animal):
                if A+'-'+str(Es) in Log['Done']: continue
                
                for E, Exp in enumerate(Exps):
                    Data = Intan.Load(Exp); Rate = 25000
                    
                    ExpName = Exp.split('/')[-1][:-4]
                    ExpFolder = '/'.join(Exp.split('/')[:-2]) + '/KlustaFiles'
                    
                    DataFile = ''.join([ExpName, '_Exp', "{0:02d}".format(int(Es)), 
                                       '-', "{0:02d}".format(int(E))])
                    RawDataInfo = {'Rate': Rate}
                    Bin.Write(DataFile+'.dat', ExpFolder, Data, RawDataInfo)
                
                FilesPrefix = ExpFolder+'/'+ExpName+'_Exp' + "{0:02d}".format(int(Es))
                RawDataInfo = Txt.DictRead(ExpFolder+'/'+DataFile+'-Info.dict')
                raw_data_files = glob(FilesPrefix + '*.dat')
                raw_data_files.sort()
                
                PrbFile = FilesPrefix+'.prb'
                Map = DataAnalysis.RemapCh('A16', 'RHAHeadstage')
                Klusta.PrbWrite(PrbFile, Map, ProbeSpacing)
                
                Klusta.PrmWrite(FilesPrefix+'.prm', ExpName+'_Exp' + "{0:02d}".format(int(Es)), 
                                PrbFile, raw_data_files, RawDataInfo['Rate'],
                                RawDataInfo['Shape'][1], RawDataInfo['DType'])
                
                Klusta.Run(ExpName+'_Exp' + "{0:02d}".format(int(Es))+'.prm', 
                           ExpFolder, Overwrite=True)
            
                Log['Done'].append(A+'-'+str(Es))
        
        except Exception as e:
            Log['Errors'].append(A+'-'+str(Es))
            Log['ErrorsLog'].append(e)
            print(''); print(e); print('')


#%%
Animals = dict(
    CaMKIIaChR2n20 = [
        ['/home/cerebro/Malfatti/MScData/20140902-CaMKIIaChR2n20-UnitRec/Intan/s45Light_155245.int',
         #'/home/cerebro/Malfatti/MScData/20140902-CaMKIIaChR2n20-UnitRec/Intan/s45SoundLight_155903.int',
         '/home/cerebro/Malfatti/MScData/20140902-CaMKIIaChR2n20-UnitRec/Intan/s45Sound_155441.int'],
        ['/home/cerebro/Malfatti/MScData/20140902-CaMKIIaChR2n20-UnitRec/Intan/s4mmS_154056.int',
         #'/home/cerebro/Malfatti/MScData/20140902-CaMKIIaChR2n20-UnitRec/Intan/s4mm_154456.int',
         '/home/cerebro/Malfatti/MScData/20140902-CaMKIIaChR2n20-UnitRec/Intan/s4mmL_153644.int']
    ],
    
    CaMKIIaChR2n21 = [
        ['/home/cerebro/Malfatti/MScData/20140917-CaMKIIaChR2n21-UnitRec/Intan/n42L_161907.int',
         '/home/cerebro/Malfatti/MScData/20140917-CaMKIIaChR2n21-UnitRec/Intan/n42S_162158.int',
         '/home/cerebro/Malfatti/MScData/20140917-CaMKIIaChR2n21-UnitRec/Intan/n42SL_162727.int'],
        ['/home/cerebro/Malfatti/MScData/20140917-CaMKIIaChR2n21-UnitRec/Intan/n45L_154601.int',
         '/home/cerebro/Malfatti/MScData/20140917-CaMKIIaChR2n21-UnitRec/Intan/n45S_154908.int',
         '/home/cerebro/Malfatti/MScData/20140917-CaMKIIaChR2n21-UnitRec/Intan/n45SL_160259.int']
    ],
    
    CaMKIIaChR2n22 = [
        ['/home/cerebro/Malfatti/MScData/20141027-CaMKIIaChR2n22-UnitRec/Intan/n35L_141027_171433.int',
         '/home/cerebro/Malfatti/MScData/20141027-CaMKIIaChR2n22-UnitRec/Intan/n35LS_141027_174533.int',
         '/home/cerebro/Malfatti/MScData/20141027-CaMKIIaChR2n22-UnitRec/Intan/n35S_141027_171722.int'],
        ['/home/cerebro/Malfatti/MScData/20141027-CaMKIIaChR2n22-UnitRec/Intan/n40L_141027_175250.int',
         '/home/cerebro/Malfatti/MScData/20141027-CaMKIIaChR2n22-UnitRec/Intan/n40LS_141027_175934.int',
         '/home/cerebro/Malfatti/MScData/20141027-CaMKIIaChR2n22-UnitRec/Intan/n40S_141027_175540.int'],
        ['/home/cerebro/Malfatti/MScData/20141027-CaMKIIaChR2n22-UnitRec/Intan/n43L_141027_180804.int',
         '/home/cerebro/Malfatti/MScData/20141027-CaMKIIaChR2n22-UnitRec/Intan/n43LS_141027_181511.int',
         '/home/cerebro/Malfatti/MScData/20141027-CaMKIIaChR2n22-UnitRec/Intan/n43S_141027_181114.int']
    ],
    
    CaMKIIaChR2n23 = [
        ['/home/cerebro/Malfatti/MScData/20150529-CaMKIIaChR2n23-UnitRec/Intan/CaMKIIaChR2n23-40L_170540.int',
         '/home/cerebro/Malfatti/MScData/20150529-CaMKIIaChR2n23-UnitRec/Intan/CaMKIIaChR2n23-40S_170252.int',
         '/home/cerebro/Malfatti/MScData/20150529-CaMKIIaChR2n23-UnitRec/Intan/CaMKIIaChR2n23-40SL_173656.int'],
        ['/home/cerebro/Malfatti/MScData/20150529-CaMKIIaChR2n23-UnitRec/Intan/CaMKIIaChR2n23-42L_174424.int',
         '/home/cerebro/Malfatti/MScData/20150529-CaMKIIaChR2n23-UnitRec/Intan/CaMKIIaChR2n23-42S_174140.int',
         '/home/cerebro/Malfatti/MScData/20150529-CaMKIIaChR2n23-UnitRec/Intan/CaMKIIaChR2n23-42SL_175200.int'],
        ['/home/cerebro/Malfatti/MScData/20150529-CaMKIIaChR2n23-UnitRec/Intan/CaMKIIaChR2n23-45L_180359.int',
         '/home/cerebro/Malfatti/MScData/20150529-CaMKIIaChR2n23-UnitRec/Intan/CaMKIIaChR2n23-45S_180035.int',
         '/home/cerebro/Malfatti/MScData/20150529-CaMKIIaChR2n23-UnitRec/Intan/CaMKIIaChR2n23-45SL_180736.int']
    ],
    
    CaMKIIaChR2n24 = [
        ['/home/cerebro/Malfatti/MScData/20150602-CaMKIIaChR2n24-UnitRec/Intan/CaMKIIaChR2n24-415L_172612.int',
         '/home/cerebro/Malfatti/MScData/20150602-CaMKIIaChR2n24-UnitRec/Intan/CaMKIIaChR2n24-415S_172855.int',
         '/home/cerebro/Malfatti/MScData/20150602-CaMKIIaChR2n24-UnitRec/Intan/CaMKIIaChR2n24-415SL_173131.int'],
        ['/home/cerebro/Malfatti/MScData/20150602-CaMKIIaChR2n24-UnitRec/Intan/CaMKIIaChR2n24-45L_173705.int',
         '/home/cerebro/Malfatti/MScData/20150602-CaMKIIaChR2n24-UnitRec/Intan/CaMKIIaChR2n24-45S_174001.int',
         '/home/cerebro/Malfatti/MScData/20150602-CaMKIIaChR2n24-UnitRec/Intan/CaMKIIaChR2n24-45SL_174333.int'],
        ['/home/cerebro/Malfatti/MScData/20150923-CaMKIIaChR2n24-UnitRec/Intan/CChR2-4000Bo_142407.int',
         '/home/cerebro/Malfatti/MScData/20150923-CaMKIIaChR2n24-UnitRec/Intan/CChR2-4000L_143307.int',
         '/home/cerebro/Malfatti/MScData/20150923-CaMKIIaChR2n24-UnitRec/Intan/CChR2-4000S_134214.int'],
        ['/home/cerebro/Malfatti/MScData/20150923-CaMKIIaChR2n24-UnitRec/Intan/CChR2-4300Bo_144945.int',
         '/home/cerebro/Malfatti/MScData/20150923-CaMKIIaChR2n24-UnitRec/Intan/CChR2-4300L_144628.int',
         '/home/cerebro/Malfatti/MScData/20150923-CaMKIIaChR2n24-UnitRec/Intan/CChR2-4300S_144151.int'],
        ['/home/cerebro/Malfatti/MScData/20150923-CaMKIIaChR2n24-UnitRec/Intan/CChR2-4500Bo_150305.int',
         '/home/cerebro/Malfatti/MScData/20150923-CaMKIIaChR2n24-UnitRec/Intan/CChR2-4500L_145953.int',
         '/home/cerebro/Malfatti/MScData/20150923-CaMKIIaChR2n24-UnitRec/Intan/CChR2-4500S_145653.int'],
        ['/home/cerebro/Malfatti/MScData/20150923-CaMKIIaChR2n24-UnitRec/Intan/MS-4000Bo_162644.int',
         '/home/cerebro/Malfatti/MScData/20150923-CaMKIIaChR2n24-UnitRec/Intan/MS-4000L_162402.int',
         '/home/cerebro/Malfatti/MScData/20150923-CaMKIIaChR2n24-UnitRec/Intan/MS-4000S_162023.int'],
        ['/home/cerebro/Malfatti/MScData/20150923-CaMKIIaChR2n24-UnitRec/Intan/MS-4500Bo_165731.int',
         '/home/cerebro/Malfatti/MScData/20150923-CaMKIIaChR2n24-UnitRec/Intan/MS-4500L_165441.int',
         '/home/cerebro/Malfatti/MScData/20150923-CaMKIIaChR2n24-UnitRec/Intan/MS-4500S_165152.int']
    ],
    
    CaMKIIaArch3n01 = [
        ['/home/cerebro/Malfatti/MScData/20150319-CaMKIIaArch3n01-UnitRec/Intan/Arch367L_143837.int',
         '/home/cerebro/Malfatti/MScData/20150319-CaMKIIaArch3n01-UnitRec/Intan/Arch367S_143507.int',
         '/home/cerebro/Malfatti/MScData/20150319-CaMKIIaArch3n01-UnitRec/Intan/Arch367SL_144132.int'],
        ['/home/cerebro/Malfatti/MScData/20150319-CaMKIIaArch3n01-UnitRec/Intan/Arch40L_150233.int',
         '/home/cerebro/Malfatti/MScData/20150319-CaMKIIaArch3n01-UnitRec/Intan/Arch40S_145922.int',
         '/home/cerebro/Malfatti/MScData/20150319-CaMKIIaArch3n01-UnitRec/Intan/Arch40SL_150909.int'],
        ['/home/cerebro/Malfatti/MScData/20150319-CaMKIIaArch3n01-UnitRec/Intan/Arch45L_153548.int',
         '/home/cerebro/Malfatti/MScData/20150319-CaMKIIaArch3n01-UnitRec/Intan/Arch45S_153244.int',
         '/home/cerebro/Malfatti/MScData/20150319-CaMKIIaArch3n01-UnitRec/Intan/Arch45SL_153835.int']
    ],
    
    CaMKIIaArch3n02 = [
        ['/home/cerebro/Malfatti/MScData/20150416-CaMKIIaArch3n02-UnitRec/IntanLightOnly/Arch3n2-35L_180342.int',
         '/home/cerebro/Malfatti/MScData/20150416-CaMKIIaArch3n02-UnitRec/IntanSound/Arch3n2-35S_175902.int',
         '/home/cerebro/Malfatti/MScData/20150416-CaMKIIaArch3n02-UnitRec/IntanSound/Arch3n2-35SL_180751.int'],
        ['/home/cerebro/Malfatti/MScData/20150416-CaMKIIaArch3n02-UnitRec/IntanLightOnly/Arch3n2-40L_184008.int',
         '/home/cerebro/Malfatti/MScData/20150416-CaMKIIaArch3n02-UnitRec/IntanSound/Arch3n2-40S_184417.int',
         '/home/cerebro/Malfatti/MScData/20150416-CaMKIIaArch3n02-UnitRec/IntanSound/Arch3n2-40SL_183521.int'],
        ['/home/cerebro/Malfatti/MScData/20150416-CaMKIIaArch3n02-UnitRec/IntanLightOnly/Arch3n2-43L_185924.int',
         '/home/cerebro/Malfatti/MScData/20150416-CaMKIIaArch3n02-UnitRec/IntanSound/Arch3n2-43LS_185520.int',
         '/home/cerebro/Malfatti/MScData/20150416-CaMKIIaArch3n02-UnitRec/IntanSound/Arch3n2-43S_184905.int'],
        ['/home/cerebro/Malfatti/MScData/20150416-CaMKIIaArch3n02-UnitRec/IntanLightOnly/Arch3n2-45L_190432.int',
         '/home/cerebro/Malfatti/MScData/20150416-CaMKIIaArch3n02-UnitRec/IntanSound/Arch3n2-45S_191253.int',
         '/home/cerebro/Malfatti/MScData/20150416-CaMKIIaArch3n02-UnitRec/IntanSound/Arch3n2-45SL_190841.int']
    ],
    
    CaMKIIaArch3n03 = [
        ['/home/cerebro/Malfatti/MScData/20150417-CaMKIIaArch3n03-UnitRec/IntanTroubles/Arch3n3-x2-37L_143710.int',
         '/home/cerebro/Malfatti/MScData/20150417-CaMKIIaArch3n03-UnitRec/IntanTroubles/Arch3n3-x2-37LS_143226.int',
         '/home/cerebro/Malfatti/MScData/20150417-CaMKIIaArch3n03-UnitRec/IntanTroubles/Arch3n3-x2-37S_144131.int'],
        ['/home/cerebro/Malfatti/MScData/20150417-CaMKIIaArch3n03-UnitRec/IntanLightOnly/Arch3n3-40L_124123.int',
         '/home/cerebro/Malfatti/MScData/20150417-CaMKIIaArch3n03-UnitRec/IntanSound/Arch3n3-40LS_123724.int',
         '/home/cerebro/Malfatti/MScData/20150417-CaMKIIaArch3n03-UnitRec/IntanSound/Arch3n3-40S_123259.int'],
        ['/home/cerebro/Malfatti/MScData/20150417-CaMKIIaArch3n03-UnitRec/IntanLightOnly/Arch3n3-43L_124636.int',
         '/home/cerebro/Malfatti/MScData/20150417-CaMKIIaArch3n03-UnitRec/IntanSound/Arch3n3-43S_125706.int',
         '/home/cerebro/Malfatti/MScData/20150417-CaMKIIaArch3n03-UnitRec/IntanSound/Arch3n3-43SL_125126.int'],
        ['/home/cerebro/Malfatti/MScData/20150417-CaMKIIaArch3n03-UnitRec/IntanLightOnly/Arch3n3-45L_131728.int',
         '/home/cerebro/Malfatti/MScData/20150417-CaMKIIaArch3n03-UnitRec/IntanSound/Arch3n3-45LS_131210.int',
         '/home/cerebro/Malfatti/MScData/20150417-CaMKIIaArch3n03-UnitRec/IntanSound/Arch3n3-45S_130436.int']
    ],
    
    CaMKIIaArch3n04 = [
        ['/home/cerebro/Malfatti/MScData/20150604-CaMKIIaArch3n04-UnitRec/Intan/CaMKIIaArch3n4-35S_190501.int',
         '/home/cerebro/Malfatti/MScData/20150604-CaMKIIaArch3n04-UnitRec/Intan/CaMKIIaArch3n4-35SL_190800.int',
         '/home/cerebro/Malfatti/MScData/20150604-CaMKIIaArch3n04-UnitRec/Intan/CaMKIIaArch3n4-43S_191146.int'],
        ['/home/cerebro/Malfatti/MScData/20150604-CaMKIIaArch3n04-UnitRec/Intan/CaMKIIaArch3n4-43SL_191426.int',
         '/home/cerebro/Malfatti/MScData/20150604-CaMKIIaArch3n04-UnitRec/Intan/CaMKIIaArch3n4-45S_191823.int',
         '/home/cerebro/Malfatti/MScData/20150604-CaMKIIaArch3n04-UnitRec/Intan/CaMKIIaArch3n4-45SL_192118.int']
    ],
    
    DIOChR2n03 = [
        ['/home/cerebro/Malfatti/MScData/20141112-DIOChR2n03-UnitRec/Intan/n37_L_164402.int',
         '/home/cerebro/Malfatti/MScData/20141112-DIOChR2n03-UnitRec/Intan/n37_S_163002.int',
         '/home/cerebro/Malfatti/MScData/20141112-DIOChR2n03-UnitRec/Intan/n37_SL_164658.int'],
        ['/home/cerebro/Malfatti/MScData/20141112-DIOChR2n03-UnitRec/Intan/n40_L_165345.int',
         '/home/cerebro/Malfatti/MScData/20141112-DIOChR2n03-UnitRec/Intan/n40_S_165839.int',
         '/home/cerebro/Malfatti/MScData/20141112-DIOChR2n03-UnitRec/Intan/n40_SL_170356.int'],
        ['/home/cerebro/Malfatti/MScData/20141112-DIOChR2n03-UnitRec/Intan/n43_L_172757.int',
         '/home/cerebro/Malfatti/MScData/20141112-DIOChR2n03-UnitRec/Intan/n43_S_172341.int',
         '/home/cerebro/Malfatti/MScData/20141112-DIOChR2n03-UnitRec/Intan/n43_SL_173437.int']
    ],
    
    DIOChR2n04 = [
        ['/home/cerebro/Malfatti/MScData/20141121-DIOChR2n04-UnitRec/Intan/L3.5_161903.int',
         '/home/cerebro/Malfatti/MScData/20141121-DIOChR2n04-UnitRec/Intan/LS3.5_163031.int',
         '/home/cerebro/Malfatti/MScData/20141121-DIOChR2n04-UnitRec/Intan/S3.5_162410.int'],
        ['/home/cerebro/Malfatti/MScData/20141121-DIOChR2n04-UnitRec/Intan/L4_163724.int',
         '/home/cerebro/Malfatti/MScData/20141121-DIOChR2n04-UnitRec/Intan/LS4_164501.int',
         '/home/cerebro/Malfatti/MScData/20141121-DIOChR2n04-UnitRec/Intan/S4_164024.int'],
        ['/home/cerebro/Malfatti/MScData/20141121-DIOChR2n04-UnitRec/Intan/L4.5_165842.int',
         '/home/cerebro/Malfatti/MScData/20141121-DIOChR2n04-UnitRec/Intan/LS4.5_170614.int',
         '/home/cerebro/Malfatti/MScData/20141121-DIOChR2n04-UnitRec/Intan/S4.5_170146.int']
    ],
    
    DIOChR2n05 = [
        ['/home/cerebro/Malfatti/MScData/20150602-DIOChR2n05-UnitRec/Intan/DioChR2n5-315L_123316.int',
         '/home/cerebro/Malfatti/MScData/20150602-DIOChR2n05-UnitRec/Intan/DioChR2n5-315S_124102.int',
         '/home/cerebro/Malfatti/MScData/20150602-DIOChR2n05-UnitRec/Intan/DioChR2n5-315SL_124440.int'],
        ['/home/cerebro/Malfatti/MScData/20150602-DIOChR2n05-UnitRec/Intan/DioChR2n5-45L_130351.int',
         '/home/cerebro/Malfatti/MScData/20150602-DIOChR2n05-UnitRec/Intan/DioChR2n5-45S_130711.int',
         '/home/cerebro/Malfatti/MScData/20150602-DIOChR2n05-UnitRec/Intan/DioChR2n5-45SL_131034.int'],
        ['/home/cerebro/Malfatti/MScData/20150602-DIOChR2n05-UnitRec/Intan/DioChR2n5-45L_134651.int',
         '/home/cerebro/Malfatti/MScData/20150602-DIOChR2n05-UnitRec/Intan/DioChR2n5-45S_134957.int',
         '/home/cerebro/Malfatti/MScData/20150602-DIOChR2n05-UnitRec/Intan/DioChR2n5-45SL_135241.int']
    ],
)
