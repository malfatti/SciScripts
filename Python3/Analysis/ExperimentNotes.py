#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@date: 2017-12-18
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""
import numpy as np
import os

from glob import glob
from IO import Txt

from DataAnalysis.Plot import Plot
Params = Plot.Set(Params=True)
from matplotlib import rcParams; rcParams.update(Params)
import matplotlib.pyplot as plt


NotesPath = os.environ['HOME']+'/Malfatti/Nebula/Documents/PhD/Notes/Experiments'
Notes = [N for N in glob(NotesPath+'/*.py') if 'Prevention' in N or
                                               'Control' in N or
                                               'VTF' in N]
Notes.sort()

Dead = {K: [] for K in ['Date', 'Ket', 'Xyl', 'Weight', 'Age', 'Shots']}
Alive = {K: [] for K in ['Date', 'Ket', 'Xyl', 'Weight', 'Age', 'Shots']}

for Note in Notes:
    Info = Txt.DictRead(Note)
    
    for E, Exp in Info['Experiments'].items():
        if 'UnitRec' in E: continue
        
        if 'Comments' in Exp.keys():
            Comms = list(Exp['Comments'].values())
            Died = [True for C in Comms if 'ied' in C or 'dead' in C]
            
            if Died:
                Dead['Date'].append(E.split('-')[0])
                Dead['Ket'].append(Exp['Injections']['Ketamine']['00'][0])
                Dead['Xyl'].append(Exp['Injections']['Xylazine']['00'][0])
                Dead['Shots'].append(len(Exp['Injections']['Ketamine'].keys()) - 3)
                Dead['Weight'].append(Exp['Animal']['Weight'])
                Dead['Age'].append(Exp['Animal']['Age'])
            else:
                Alive['Date'].append(E.split('-')[0])
                Alive['Ket'].append(Exp['Injections']['Ketamine']['00'][0])
                Alive['Xyl'].append(Exp['Injections']['Xylazine']['00'][0])
                Alive['Shots'].append(len(Exp['Injections']['Ketamine'].keys()) - 3)
                Alive['Weight'].append(Exp['Animal']['Weight'])
                Alive['Age'].append(Exp['Animal']['Age'])

Colors = ['r', 'b']
for E, End in enumerate([('Dead', Dead), ('Alive', Alive)]):
    Ax = plt.axes()
    Ax.scatter(End[1]['Age'], End[1]['Weight'], c=Colors[E], label=End[0])
    # X = [I[0] for I in sorted(enumerate(End['Age']), key = lambda x:x[1])]
    # Y = [End['Weight'][I] for I in X]
    # X = [End['Age'][I] for I in X]
    
    # plt.plot(X, Y)
AxArgs = {'xlabel': 'Age [days]', 
          'ylabel': 'Weight [g]',}
Plot.Set(Ax=Ax, AxArgs=AxArgs)
plt.legend(loc='best')
plt.savefig('WeightxAge.pdf', format='pdf', dpi=300)
plt.show()