# -*- coding: utf-8 -*-
"""
Script for data analysis

Experiment: stimulation of the brainstem using light and sound, recording with
a silicon probe (16 channels) + 2 tungsten wires + reference screw.
"""
#%% Set experiment
Rate = 100000
LFPCh = [1, 16]
BrokenCh = []
FilterData = True

BrokenCh = [_ - 1 for _ in BrokenCh]
Channels = [_ for _ in list(range(16)) if _ not in BrokenCh]

#%% Import
import glob
import OpenEphys
import os

import matplotlib.pyplot as plt

#%% Rename files

DirList = glob.glob('OpenEphysFiles/*'); DirList.sort()

# Remove date from folder name
for FolderName in DirList:
    NewFolderName = ''.join([FolderName[11:]])
    NewFolderName = NewFolderName.replace("-", "")
    os.rename(FolderName, NewFolderName)
    print(FolderName, ' moved to ', NewFolderName)


DirList = glob.glob('OpenEphysFiles/*'); DirList.sort()

# Add 0 to ch numbers < 10
for RecFolder in DirList:
    DataFiles = glob.glob(''.join([RecFolder,'/*.continuous']))
    EventsFile = glob.glob(''.join([RecFolder,'/all_chann*.events']))
    DataFiles.sort()
    
    for FileName in DataFiles:
        if FileName[-13] == 'H':
            NewFileName = ''.join([FileName[:-12], '0', FileName[-12:]])
            os.rename(FileName, NewFileName)
            print(FileName, ' moved to ', NewFileName)

    DataFiles = glob.glob(''.join([RecFolder,'/*.continuous']))
    DataFiles.sort()


#%% Generate TTLs

os.makedirs('DataArrays', exist_ok=True)
DirList = glob.glob('OpenEphysFiles/*'); DirList.sort()

for RecFolder in DirList:
    DataFiles = glob.glob(''.join([RecFolder,'/*.continuous']))
    EventsFile = ''.join([RecFolder,'/all_channels.events'])
    DataFiles.sort()
    
    TempData = OpenEphys.loadContinuous(DataFiles[0]) # Any ch would work ;)
    Events = OpenEphys.loadEvents(EventsFile)
                                                                                
    # Draw TTL
    print('Drawing TTL channel...')                                              
    if len(Events['timestamps']) > 100: # arbitrary value, I never use < 100 stimulation pulses
        TTLSound = [0] * len(TempData['data'])
        for EvTS in range(len(Events['timestamps'])):
            if Events['eventType'][EvTS] == 3 # if event is a TTL
                if Events['eventId'][EvTS] == 1 # if TTL is on
                    TTLSound(find(TempTimestamps==EventsTimestamps(EvTS)): ...  
                        find(TempTimestamps==EventsTimestamps(EvTS+1))) = 1;    
                end                                                             
            end                                                                 
        end                                                                     
    else                                                                        
        TTLSound = zeros(size(TempData, 1), size(TempData, 2));                 
    end                                                                         
    disp('Done.')


#%% Filter raw data
OpenEphys.loadFolderToArray(RecFolder)
import glob
import numpy
import os
import Intan

SoundIntFiles = glob.glob('IntanSound/*.int')
LightOnlyIntFiles = glob.glob('IntanLightOnly/*.int')
N = numpy.size(SoundIntFiles, axis=0)
O = numpy.size(LightOnlyIntFiles, axis=0)

os.makedirs('SpkSound', exist_ok=True)
os.makedirs('SpkLightOnly', exist_ok=True)

LFPCh = numpy.array([16, 1]) # Larger number first!!!

for SoundIntFile in SoundIntFiles:
    MyData = Intan.ReadData(SoundIntFile)
    AuxNo = (MyData['aux'].mean(axis=0)>0).nonzero()
    aux=[aux,zeros(length(aux),1)]; % add 1 more empty row to aux
    print(SoundIntFile)

    if length(AuxNo)>0
        kk=1;
        [peaks locs] = findpeaks(double(aux(:,AuxNo(kk))));
        AuxLoc(kk).locs=locs;
        AuxLoc(kk).Aux = aux(:,AuxNo(kk));

        for jj=1:length(locs)
            aux(locs(jj):locs(jj)+length(FakeTTL)-1,7)=FakeTTL;
            kk=2;
            [peaks locs2] = findpeaks(double(aux(:,7)));
            AuxLoc(kk).locs=locs2;
            AuxLoc(kk).Aux = aux(:,7);
        end
    end

    % test TTL

    %plot(t,aux)
    %hold on
    %plot(t(AuxLoc(2).locs),ones(length(AuxLoc(2).locs)),'ro')
    %hold off
    %pause

    % the little balls should be aligned w/ the pulses

    Spk=zeros(size(y,1),size(y,2));

    for jj=1:size(y,1)
        Spk(jj,:) = eegfilt(double(y(jj,:)),25000,300,3000);
        jj
    end

    Spk(LFPCh,:)=[];
    FName=SoundFiles(ii).name;
    Fstr=findstr(FName,'.int');
    FName=FName(1:Fstr-1);
    
     if length(AuxNo)>0
    disp(['Saving as ', FName]);
    save(['SpkSound/',FName,'-Spk.mat'],'Spk');
    save(['SpkSound/',FName,'-AuxLoc.mat'],'AuxLoc');
    save(['SpkSound/',FName,'-AuxT.mat'],'aux','t');
    clear AuxNo
    else
        disp(['Saving as ', FName]);
        AuxLoc=0;
        save(['SpkSound/',FName,'-Spk.mat'],'Spk');
        save(['SpkSound/',FName,'-AuxLoc.mat'],'AuxLoc');
        save(['SpkSound/',FName,'-AuxT.mat'],'aux','t');
        clear AuxNo
        disp('No Aux file!!!')
    end
end

for ii=1:O
    [t,amps,y,aux] = read_intan_data_leao(['IntanLightOnly/',LightOnlyFiles(ii).name]);
    AuxNo = find(mean(aux)>0);
    disp(LightOnlyFiles(ii).name)
    if length(AuxNo)>0
    for kk=1:length(AuxNo)
        [peaks locs] = findpeaks(double(aux(:,AuxNo(kk))));
        AuxLoc(kk).locs=locs;
        AuxLoc(kk).Aux = aux(:,AuxNo(kk));
    end
    end
    Spk=zeros(size(y,1),size(y,2));
    for jj=1:size(y,1)
        Spk(jj,:) = eegfilt(double(y(jj,:)),25000,300,3000);
        jj
    end
    Spk(LFPCh,:)=[];
    FName=LightOnlyFiles(ii).name;
    Fstr=findstr(FName,'.int');
    FName=FName(1:Fstr-1);
    if length(AuxNo)>0
    disp(['Saving as ', FName]);
    save(['SpkLightOnly/',FName,'-Spk.mat'],'Spk');
    save(['SpkLightOnly/',FName,'-AuxLoc.mat'],'AuxLoc');
    clear AuxNo
    else
        disp(['Saving as ', FName]);
        AuxLoc=0;
        save(['SpkLightOnly/',FName,'-Spk.mat'],'Spk');
        save(['SpkLightOnly/',FName,'-AuxLoc.mat'],'AuxLoc');
        clear AuxNo
        disp('No Aux file!!!')
    end
end
