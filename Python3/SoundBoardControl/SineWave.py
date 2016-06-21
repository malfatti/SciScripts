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

"""

import array
import math
import pyaudio

Rate = 128000
Freq = 200; Time = 10
    
print('Generating sound...')
#Pulse = [math.sin(2*math.pi*Freq*(_/Rate)) * 3.3/1.7 #+ 0.5 
#         for _ in range(round(Rate*Time))]
#Pulse[-1] = 0

Pulse = [1]*round((Rate*0.01)//2) + [-1]*round((Rate*0.01)//2) + [0]*round(Rate*0.08)
Pulse = Pulse*100
Pulse[-1] = 0; Pulse[0] = 0

print('Interleaving channels...')
List = [0]*(2*len(Pulse))
for _ in range(len(Pulse)):
    List[_ *2] = 0
    List[_ *2+1] = Pulse[_]

Pulse = array.array('f')
Pulse.fromlist(List)
Pulse = bytes(Pulse)

print('Generating sound objects...')
p = pyaudio.PyAudio()
Stimulation = p.open(format=pyaudio.paFloat32,
                channels=2,
                rate=Rate,
                output=True)

print('Playing...')
Stimulation.write(Pulse)