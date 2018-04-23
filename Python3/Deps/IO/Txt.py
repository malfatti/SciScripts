#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:39:12 2017

@author: malfatti
"""
import numpy as np

from ast import literal_eval



def DictFlat(Var, UpKey='', KeySep='_', Flat={}):
    if type(Var) == dict:
        for K, V in Var.items():
            NewKey = UpKey + KeySep + K if UpKey else K
            Flat = {**Flat, **DictFlat(Var[K], NewKey, KeySep, Flat)}
        return(Flat)
    else:
        Flat[UpKey] = Var
        return(Flat)


def DictPrint(value, htchar='    ', itemchar=' ', breaklineat='auto', lfchar='\n', indent=0):
    ''' Modified from y.petremann's code.
        Added options to set item separator for list or tuple and to set a number
        of items per line, or yet, to calculate items per line so it will not 
        have more than 80 chars per line.
        Source: https://stackoverflow.com/a/26209900 '''
    
    nlch = lfchar + htchar * (indent + 1)
    if type(value) is dict:
        items = [
            nlch + repr(key) + ': ' + DictPrint(value[key], htchar, itemchar, breaklineat, lfchar, indent + 1)
            for key in value
        ]
        
        return '{%s}' % (','.join(items) + lfchar + htchar * indent)
    
    elif type(value) is list or type(value) is tuple:
        items = [
            itemchar + DictPrint(item, htchar, itemchar, breaklineat, lfchar, indent + 1)
            for item in value
        ]
        
        if breaklineat == 'auto':
           bl = int((80 - (len(htchar)*(indent + 1)))/
                (int((sum([len(i)+4 for i in items])-len(itemchar)-1)/len(items))))
         
        else: bl = breaklineat
        
        if not bl: bl = 1
       
        if len(items) > bl:
            for i in list(range(bl, len(items), bl)):
                items[i] = lfchar + htchar*(indent+1) + '  ' + items[i]
        
        return '[%s]' % (','.join(items))
    
    elif type(value) is np.ndarray:
        value = value.tolist()
        items = DictPrint(value, htchar, itemchar, breaklineat, lfchar, indent)
        return items
#        items = [
#            itemchar + DictPrint(item, htchar, itemchar, breaklineat, lfchar, indent + 1)
#            for item in value
#        ]
#        
#        if breaklineat == 'auto':
#           bl = int((80 - (len(htchar)*(indent + 1)))/
#                (int((sum([len(i)+4 for i in items])-len(itemchar)-1)/len(items))))
#         
#        else: bl = breaklineat
#        
#        if not bl: bl = 1
#       
#        if len(items) > bl:
#            for i in list(range(bl, len(items), bl)):
#                items[i] = lfchar + htchar*(indent+1) + '  ' + items[i]
#        
#        return '[%s]' % (','.join(items))
    
#     elif type(value) is tuple:
#         items = [
#             itemchar + DictPrint(item, htchar, itemchar, breaklineat, lfchar, indent + 1)
#             for item in value
#         ]
#         
#         if breaklineat == 'auto':
#            bl = (80 - (len(htchar)*(indent + 1)))// \
#                 ((sum([len(i)+4 for i in items])-2)//len(items))
#         else: bl = breaklineat
#         
#         if not bl: bl = 1
#        
#         if len(items) > bl:
#             for i in list(range(bl, len(items), bl)):
#                 items[i] = lfchar + htchar*(indent+1) + '  ' + items[i]
#         
#         return '(%s)' % (','.join(items))
    
    else:
        return repr(value)


def DictRead(File):
    Dict = literal_eval(open(File).read())
    return(Dict)


def DictWrite(File, Var):
    with open(File, 'w') as F: F.write(DictPrint(Var))
    return(None)

