# -*- coding: utf-8 -*-
"""
    Copyright (C) 2015  T. Malfatti                                             
                                                                                
    This program is free software: you can redistribute it and/or modify        
    it under the terms of the GNU General Public License as published by        
    the Free Software Foundation, either version 3 of the License, or           
    (at your option) any later version.                                         
                                                                                
    This program is distributed in the hope that it will be useful,             
    but WITHOUT ANY WARRANTY; without even the implied warranty of              
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               
    GNU General Public License for more details.                                
                                                                                
    You should have received a copy of the GNU General Public License           
    along with this program.  If not, see <http://www.gnu.org/licenses/>.       

Generate chaotic fractal note sequences

In order to create chaotic fractal phase space trajectories, codes written by 
Dr. Sean Whalen (year?) and by Gribble and Kistemaker (2012) were modified and 
applied in a totally different context :)

Sources:
http://www.node99.org/tutorials/ar/
http://www.gribblelab.org/compneuro/2_Modelling_Dynamical_Systems.html

The code in "Giving names to numbers" section was adapted from Kroger (2012). 

Source:
http://pedrokroger.net/mfgan/

The code in "Making it readable and putting together" was written in 
collaboration with a great colleague and friend, J. Targino. 

"""
#%% Defining functions

import numpy
import matplotlib as plot

from mpl_toolkits.mplot3d.axes3d import Axes3D
from itertools import groupby
from operator import itemgetter

## Fourth-order Runge-Kutta integrator
def rk4(odes, LorenzCoord, LorenzParameters, dt=0.01):
    k1 = dt * odes(LorenzCoord, LorenzParameters)
    k2 = dt * odes(LorenzCoord + 0.5 * k1, LorenzParameters)
    k3 = dt * odes(LorenzCoord + 0.5 * k2, LorenzParameters)
    k4 = dt * odes(LorenzCoord + k3, LorenzParameters)
    return LorenzCoord + (k1 + 2 * k2 + 2 * k3 + k4) / 6


## Generate data
def generate(DataLength, odes, LorenzCoord, LorenzParameters):
    data = numpy.zeros([LorenzCoord.shape[0], DataLength])
    # Since the first iterates are transient, they can removed.
    for i in range(800):
        LorenzCoord = rk4(odes, LorenzCoord, LorenzParameters)

    for i in range(DataLength):
        LorenzCoord = rk4(odes, LorenzCoord, LorenzParameters)
        data[:, i] = LorenzCoord

    return data

    
## Define values for Lorenz equation
def lorenz_odes(LorenzCoord,LorenzParameters):
  # unpack the LorenzCoord vector
  x = LorenzCoord[0]
  y = LorenzCoord[1]
  z = LorenzCoord[2]

  # these are our constants
  rho = LorenzParameters[0]
  beta = LorenzParameters[1]  
  sigma = LorenzParameters[2]

  # compute LorenzCoord derivatives
  xd = sigma*(y - x)
  yd = (rho-z)*x - y
  zd = x*y - beta*z

  # return the LorenzCoord derivatives
  return numpy.array([xd, yd, zd])


## Organize data to generate Lorenz attractor
def lorenz_generate(DataLength):
    return generate(DataLength, lorenz_odes, LorenzCoord, LorenzParameters)


## Define values for Rossler equation
def rossler_odes(RosslerCoord,RosslerParameters):
  # unpack the RosslerCoord vector
  x = RosslerCoord[0]
  y = RosslerCoord[1]
  z = RosslerCoord[2]

  # these are our constants
  a = RosslerParameters[0]
  b = RosslerParameters[1]  
  c = RosslerParameters[2]

  # compute RosslerCoord derivatives
  xd = -y - z
  yd = x + a*y
  zd = b + z*(x - c)

  # return the RosslerCoord derivatives
  return numpy.array([xd, yd, zd])


## Organize data to generate Lorenz attractor
def rossler_generate(DataLength):
    return generate(DataLength, rossler_odes, RosslerCoord, RosslerParameters)


## Plot the data
def PlotData(Data):
    if (Data.all == LorenzData.all) or (Data.all == RosslerData.all):
        # 2D plot
        plot.pylab.plot(Data[0])
        #3D plot
        figure = plot.pylab.figure()
        axes = Axes3D(figure)
        axes.plot3D(Data[0], Data[1], Data[2])
        figure.add_axes(axes)
        plot.pylab.show()
    else:
        print('Are you bananas?')


## Giving names to numbers

def NoteCode(Note):
    """
    This will generate music that can contain all notes from C0 to B6, where 
    C3 is the central C. To change this behavior (for example limiting the 
    notes you want to a specific scale) remove the notes you don't want, 
    update the DataLength accordingly (see comment there) and change the 
    number of chunks to subdivide the notes list (see comment in "Mapping 
    notes and intensities" section). Or just call the musician :)
    """
    NoteCodes = "c,, d,, e,, ges,, aes,, bes,, \
                 c, d, e, ges, aes, bes, \
                 c d e ges aes bes \
                 c' d' e' ges' aes' bes' \
                 c'' d'' e'' ges'' aes'' bes'' \
                 c''' d''' e''' ges''' aes''' bes''' \
                 c'''' d'''' e'''' ges'''' aes'''' bes'''' \
                 ".split()
    return NoteCodes[Note]


def DurationCode(Dur):
    DurationCodes = "64 64. 32 32. 16 16. 8 8. 4 4. 2 2. 1".split()
    return DurationCodes[Dur]
    

def IntensityCode(Int):
    IntensityCodes = "pppp ppp pp p mp mf f ff fff ffff".split()
    return IntensityCodes[Int]
    

def ChunkList(l, n):
    n = max(1, n)
    return [l[i:i + n] for i in range(0, len(l), n)]


#%% Setting values and generating data
"""
I'm using the inicial coordinates as follows:
x = fundamental frequency of the note *10e-3
y = Tempo (BPM/60 = Hz; will be used to set duration of each note) *10e-3
z = intensity of the note *10e-3

see https://en.wikipedia.org/wiki/Piano_key_frequencies
and https://en.wikipedia.org/wiki/Tempo#Basic_tempo_markings
and http://ada.evergreen.edu/~arunc/intro_doc/node11.htm

Basically, you set the first note and let chaos do the rest :)
"""
LorenzX = 0.0979989 
LorenzY = 0.0011667
LorenzZ = 0.08422
LorenzCoord = numpy.array([LorenzX, LorenzY, LorenzZ])

RosslerX = 10
RosslerY = 0
RosslerZ = 0
RosslerCoord = numpy.array([RosslerX, RosslerY, RosslerZ])

"""
According to McGUlNNESS (1983), using parameters r = 40, sigma = 16 and 
b = 4 will generate a Lorenz attractor with Hausdorff dimension of 
2.06 ± 0.01
"""
LorenzParameters = numpy.array([40, 4, 16])
RosslerParameters = numpy.array([0.15, 0.2, 10.0])


"""
DataLength sets the "time", which is, the number of notes to generate.
For the notes to be correctly mapped, choose a number that is divisible by the 
number of notes used (84) and the number of intensities used (10)
"""
DataLength = 3360
        
# Generate data
LorenzData = lorenz_generate(DataLength)
RosslerData = rossler_generate(DataLength)

# Define which data to use
Data = LorenzData
    
# Plot it
PlotData(Data)


#%% Map notes, intensities and durations
"""
In order to transform this trajectory into a melody, the x axis should be 
mapped into notes, the y axis (yplot) should be mapped into durations, and 
the z axis (zplot) should be mapped into intensities.
"""

## Make all values > 0 (when it has neg values) by adding the minimal value to 
#  all values

if (Data[0].min() < 0):
    XCoordPos = [Data[0][XCoord] - Data[0].min() for XCoord in range(Data[0].size)]
else:
    XCoordPos = [Data[0][XCoord] for XCoord in range(Data[0].size)]


if (Data[1].min() < 0):
    YCoordPos = [Data[1][YCoord] - Data[1].min() for YCoord in range(Data[1].size)]
else:
    YCoordPos = [Data[1][YCoord] for YCoord in range(Data[1].size)]


if (Data[2].min() < 0):
    ZCoordPos = [Data[2][ZCoord] - Data[2].min() for ZCoord in range(Data[2].size)]
else:
    ZCoordPos = [Data[2][ZCoord] for ZCoord in range(Data[2].size)]


## Mapping notes and intensities
"""
Mapping notes and intensities to index code (to further apply functions to 
give names to notes).
"""
SortedXCoordPos = sorted(XCoordPos)
SortedZCoordPos = sorted(ZCoordPos)

# 84 and 10 are the number of notes and intensities possible, respectively. If 
# you change the notes or intensities, remember to change the numbers here.
ChunkedXCoordPos = numpy.array(ChunkList(SortedXCoordPos, int(len(XCoordPos)/42)))
ChunkedZCoordPos = numpy.array(ChunkList(SortedZCoordPos, int(len(ZCoordPos)/10)))

NotesIndexCode = [numpy.nonzero(ChunkedXCoordPos == RawNotes)[0][0] for RawNotes in XCoordPos]
IntensitiesIndexCode = [numpy.nonzero(ChunkedZCoordPos == RawIntensities)[0][0] for RawIntensities in ZCoordPos]


## Mapping durations
"""
Since the space between the values are logarithmic, the mapping should be done
considering that the sum of the space between the values is equals 1, thus the 
other spaces will be a percentage of it.

I considered 14 durations: semibreve, minim, crotchet, quaver, semiquaver, 
demisemiquaver, hemidemisemiquaver, and their dotted augmentations (except 
dotted semibreve).
See https://en.wikipedia.org/wiki/Note_value
"""

DurationsIndexCode = [0] * len(YCoordPos);
DurationRange = max(YCoordPos) - min(YCoordPos)

Hemidemisemiquaver = DurationRange * 0.004514672686230248;
DottedHemidemisemiquaver = DurationRange * 0.011286681715575621;
Demisemiquaver = DurationRange * 0.020316027088036117;
DottedDemisemiquaver = DurationRange * 0.033860045146726865;
Semiquaver = DurationRange * 0.05191873589164785;
DottedSemiquaver = DurationRange * 0.07900677200902935;
Quaver = DurationRange * 0.11512415349887133;
DottedQuaver = DurationRange * 0.16930022573363432;
Crotchet = DurationRange * 0.24153498871331827;
DottedCrotchet = DurationRange * 0.34988713318284426;
Minim = DurationRange * 0.49435665914221216;
DottedMinim = DurationRange * 0.7110609480812641;
Semibreve = DurationRange

for DurEl in YCoordPos:
    if DurEl <= Hemidemisemiquaver:
        DurationsIndexCode[YCoordPos.index(DurEl)] = int(0);
    elif Hemidemisemiquaver < DurEl <= DottedHemidemisemiquaver:
        DurationsIndexCode[YCoordPos.index(DurEl)] = int(1);
    elif DottedHemidemisemiquaver < DurEl <= Demisemiquaver:
        DurationsIndexCode[YCoordPos.index(DurEl)] = int(2);
    elif Demisemiquaver < DurEl <= DottedDemisemiquaver:
        DurationsIndexCode[YCoordPos.index(DurEl)] = int(3);
    elif DottedDemisemiquaver < DurEl <= Semiquaver:
        DurationsIndexCode[YCoordPos.index(DurEl)] = int(4);
    elif Semiquaver < DurEl <= DottedSemiquaver:
        DurationsIndexCode[YCoordPos.index(DurEl)] = int(5);
    elif DottedSemiquaver < DurEl <= Quaver:
        DurationsIndexCode[YCoordPos.index(DurEl)] = int(6);
    elif Quaver < DurEl <= DottedQuaver:
        DurationsIndexCode[YCoordPos.index(DurEl)] = int(7);
    elif DottedQuaver < DurEl <= Crotchet:
        DurationsIndexCode[YCoordPos.index(DurEl)] = int(8);
    elif Crotchet < DurEl <= DottedCrotchet:
        DurationsIndexCode[YCoordPos.index(DurEl)] = int(9);
    elif DottedCrotchet < DurEl <= Minim:
        DurationsIndexCode[YCoordPos.index(DurEl)] = int(10);
    elif Minim < DurEl <= DottedMinim:
        DurationsIndexCode[YCoordPos.index(DurEl)] = int(11);
    else:
        DurationsIndexCode[YCoordPos.index(DurEl)] = int(12);
    


#%% Making it readable and putting together
"""
Now the lists containing notes, intensities and durations in an index number 
code will be translated and concatenated into a list containing notes with the 
durations and intensities coded in lilypond format.
"""

NotesReadable = [NoteCode(int(NREl)) for NREl in NotesIndexCode]
DurationsReadable = [DurationCode(int(DREl)) for DREl in DurationsIndexCode]
IntensitiesReadable = [IntensityCode(int(IREl)) for IREl in IntensitiesIndexCode]

MelodyRaw = [NRVal+DurationsReadable[NoteIndex]+chr(92)+
             IntensitiesReadable[NoteIndex] for NoteIndex,NRVal in 
                                                enumerate(NotesReadable)]

## Writing the text file

# Breaking lines each 7 elements
ElementsPerLine = 7
MelodySize = len(MelodyRaw)
NumberOfLines = int(MelodySize/ElementsPerLine)

Melody = '\n'.join([' '.join(MelodyRaw[(i*ElementsPerLine):(ElementsPerLine+i*ElementsPerLine)]) for i in range(NumberOfLines)])

## Write the text file
TextFile = open('Melody.txt','w')
TextFile.write(Melody);
TextFile.close()

"""
Remember: this is not music. This is an arbitrary mapping of coordinates
of a phase space trajectory into a melody. By definition, there will
rarely be an interesting melody in it, since the notes will end up having
basically steps, and almost no skips. 

(see https://en.wikipedia.org/wiki/Steps_and_skips) 

Also, there will be a great number of note repetitions, so it will probably be 
veeery monotonous.

In summary, generate your notes and call a musician :)
"""
