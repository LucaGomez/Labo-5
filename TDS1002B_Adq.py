# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 10:44:21 2021

@author: User
"""


import pyvisa as visa
import time
import numpy as np
from matplotlib import pyplot as plt


rm = visa.ResourceManager()

instrumentos = rm.list_resources()  
print(instrumentos)
# Esto lista todos los dispositivos disponibles, uno de los cuales
# deberia ser "USB0::0x0699::0x0368::C017044::INSTR", donde los terminos
# indican "puerto::marca::modelo::nro_de_serie" del instrumento.
#%%
# Elijo el elemento que corresponde en instrumentos
#Con ese nombre abro el vinculo con el osciloscopio

osci=rm.open_resource(instrumentos[0])
#osc=rm.open_resource('USB0::0x0699::0x0363::C065093::INSTR')
#chequeo la comunicación
print(osci.query('*IDN?'))

#MEDIMOS
#para leer las curvas del canal 1 y 2 necesito tomar los datos de la configuración
xze,xin=osci.query_ascii_values('WFMPRE:XZE?;XIN?',separator=';') #conf. base de tiempo
yze1,ymu1,yoff1=osci.query_ascii_values('WFMPRE:CH1:YZE?;YMU?;YOFF?',separator=';') #conf. vertical canal 1
yze2,ymu2,yoff2=osci.query_ascii_values('WFMPRE:CH2:YZE?;YMU?;YOFF?',separator=';') #conf. vertical canal 2

## Modo de transmision: Binario (osea, la onda digitalizada)
osci.write('DAT:ENC RPB') 
osci.write('DAT:WID 1') 

#leo las curvas como datos binarios
osci.write('DAT:SOU CH1' )
data1=osci.query_binary_values('CURV?', datatype='B',container=np.array)
osci.write('DAT:SOU CH2')    
data2=osci.query_binary_values('CURV?', datatype='B',container=np.array)

#transformo los datos 
tiempo = xze + np.arange(len(data1)) * xin #tiempos en s
data1V=(data1-yoff1)*ymu1+yze1 #tensión canal 1 en V
data2V=(data2-yoff2)*ymu2+yze2 #tensión canal 2 en V 

osci.close()

#graficamos los datos
plt.plot(tiempo,data1V,label='Canal 1')
plt.plot(tiempo,data2V,label='Canal 2')
plt.legend()

#%%
#guardamos los datos
import pandas as pd

mediciones=np.zeros([3,2500])
mediciones=[tiempo,data1V,data2V]
mediciones=np.transpose(mediciones)
df=pd.DataFrame(mediciones)
#print(time.localtime())
df.to_csv('Mediciones'+str(time.localtime()[0])+'-'+str(time.localtime()[1])+'-'+str(time.localtime()[2])+'-'+str(time.localtime()[3])+'-'+str(time.localtime()[4])+'-'+str(time.localtime()[5])+'.csv')

