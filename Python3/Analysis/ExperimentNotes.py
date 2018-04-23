#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@date: 2017-12-18
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""
# import numpy as np
import os

from glob import glob
from IO import Txt
from IO.IO import RunProcess

from DataAnalysis.Plot import Plot
Params = Plot.Set(Params=True)
from matplotlib import rcParams; rcParams.update(Params)
import matplotlib.pyplot as plt


def TimeKey(Key):
    MM = 0 if int(Key[4:]) <30 else 1
    MM = "{0:02d}".format(int(Key[2:4])+MM)
    HH = Key[:2]
    
    if MM == 60: 
        HH = "{0:02d}".format(int(HH)+1)
        MM = '00'
    
    Time = HH + 'h' + MM
    return(Time)


def Note2Text(Note):
    Info = Txt.DictRead(Note)
    
    Text = ['~~~ ' + Info['Animal']['Name'] + ' ~~~', '', '']
    
    Text += ['### Animal info', '']
    for K,V in Info['Animal'].items():
        if K == 'Name': continue
        Text.append(K+': '+V)
    Text += ['', '']
    
    Text += ['### Experiments', '']
    for E,Exp in Info['Experiments'].items():
        Text += ['## '+E, '']
        
        for EB, ExpBlock in Exp.items():
            if EB in ['Animal', 'Setup']:
                Text.append('* '+EB)
                for K,V in ExpBlock.items(): Text.append(K+': '+str(V))
                Text.append('')
            
            elif EB == 'Injections':
                Text.append('* '+EB)
                for I,Inj in ExpBlock.items():
                    Text.append(I+':')
                    
                    if type(Inj['Location']) == str:
                        Loc = Inj['Location']
                    elif type(Inj['Location']) == dict:
                        Loc = 'AP: '+str(Inj['Location']['AP']) + '\n              ' + \
                              'ML: '+str(Inj['Location']['ML']) + '\n              ' + \
                              'DV: '+str(Inj['Location']['AP']) + '\n              '
                        
                    Text.append('    Location: '+Loc)
                    
                    ThisDrug = [
                            '    ' + \
                            TimeKey(Drug[2]) + ': ' + \
                            str(Drug[0])+''+Inj['Units'][0] + ', ' + \
                            str(Drug[1])+''+Inj['Units'][1]
                        for D, Drug in Inj.items()
                        if D not in ['Units', 'Location']
                    ]
                    ThisDrug.sort()
                    
                    Text += ThisDrug+['']
        
        for EB, ExpBlock in Exp.items():
            if EB == 'Comments':
                Text.append('* '+EB)
                for K,V in ExpBlock.items(): 
                    KK = TimeKey(K) if len(K) == 6 else K
                    Text.append(KK+': '+str(V))
                Text.append('')
        
        Text.append('')
    
    Text = '\n'.join(Text)
    return(Text)


def Note2Tex(Note):
    Text = Note2Text(Note)
    Replace = [('%', '\\%'), ('µl', '\\si{\\ul}'), ('###', '\\section'), 
               ('##', '\\subsection'), ('#')]
    
    for R in Replace:
        Text = Text.replace(R)
    
    return(Text)


def NoteTexHead():
    return(
    r'''%% Expnotes
\documentclass[11pt,a4paper,notitlepage]{article}
\usepackage[utf8]{inputenc} % Allow accents
\usepackage{siunitx} % Allow SI units
\usepackage[version=3]{mhchem} % Allow chemical formulas
\usepackage{enumitem} % Enhanced itemize
''')


def NoteWrite(Note, Format='txt'):
    if Format == 'txt':
        Text = Note2Text(Note)
        with open(Note.split('.')[0]+'.txt', 'w') as F: F.write(Text)
    
    elif Format in ['tex', 'pdf']:
        Path = '/'.join(Note.split('/')[:-1])+'/Tex'
        File = Note.split('/')[-1].split('.')[0]
        os.makedirs(Path, exist_ok=True)
        
        Head = NoteTexHead()
        Text = Note2Tex(Note)
        with open(Path+'/'+File+'.tex', 'w') as F:
            F.write(Head)
            F.write('\\begin{document}')
            F.write(Text)
            F.write('\\end{document}')
        
        if Format == 'pdf':
            Cmd = ['pdflatex', Path+'/'+File+'.tex']
            ReturnCode = RunProcess(Cmd, Path+'/'+File+'.log')
            
            if ReturnCode: print('An error ocurred: no pdf generated.')
            else:
                os.rename(Path+'/'+File+'.pdf', Path+'/../'+File+'.pdf')
                os.rmdir(Path)
        
        elif Format == 'tex':
            os.rename(Path+'/'+File+'.tex', Path+'/../'+File+'.tex')
            os.rmdir(Path)
            

NotesPath = os.environ['HOME']+'/Nebula/Documents/PhD/Notes/Experiments/Animals'
Notes = [N for N in glob(NotesPath+'/*.dict')]
Notes.sort()

for Note in Notes:
    Text = NotesWrite(Note)
    
    with open(Note.split('.')[0]+'.tex', 'w') as NoteFile: 
        NoteFile.write('\n'.join(Text))













#%% Report animals that died during surgeries
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