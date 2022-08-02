#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 08:57:45 2022

@author: publico
"""
import numpy as np
import matplotlib.pyplot as plt
import pyvisa
import time
from funcGen import AFG
from osciloscope import Osciloscope

rm = pyvisa.ResourceManager('@py')
instrumentos = rm.list_resources()

osc = Osciloscope(1)

# osci=rm.list_resources()[1]
# gen=rm.list_resources()[0]

t, V = osc.getWindow(1)

plt.plot(t,V)

#%%
from funcGen import AFG
gen=AFG(0)
#%%
gen.barrido(0.5, 10, 100, 5, 5, 1)
#%%
osc.guardar(1)

#%%
data=np.loadtxt('/home/publico/Desktop/Grupo 6/'+'medicion_tiempo.txt')
t=data[0]
V=data[1]
plt.plot(t,V)

#%%
t=[]
V=[]
def barrido(amplitud, freqini, freqfin, step, Time, channel):
    gen.setVoltage(amplitud, channel)
    Freq=freqini
    gen.setFrequency(Freq, channel)
    while Freq < freqfin:
        v1,v2,v3=osc.getMeasValues()
        t.append(v1)
        V.append(v2)
        Freq=Freq+step
        gen.setFrequency(Freq, channel) 
        time.sleep(Time)
barrido(1,10,100,10,2,1)
data_to_save = np.vstack((t, V))
folderpath = '/home/publico/Desktop/Grupo 6/'
np.savetxt(folderpath  + 'medicion_barrido.txt', data_to_save)

            
#%%

