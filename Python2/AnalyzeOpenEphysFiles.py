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

Just sketches for now.
"""
#%%
import OpenEphys
import matplotlib.pyplot as plt

OpenEphysData = OpenEphys.load('100_CH6.continuous')
OpenEphysEvents = OpenEphys.load('all_channels.events')

RawData = OpenEphysData["data"].tolist()
Events = OpenEphysEvents["timestamps"].tolist()
EventChannels = OpenEphysEvents["channel"].tolist()
EventTypes = OpenEphysEvents["eventType"].tolist()

TTL0 = [0] * len(RawData)
TTL1 = [0] * len(RawData)

