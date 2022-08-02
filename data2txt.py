# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 14:47:38 2022

@author: Labo5 v2022
"""

import numpy as np
folderpath = 'C:\\Users\\Martin\\Nextcloud\\Labo 5\\Codigos intstrumentacion\\writing\\'


#%% open file and writes
f = open(folderpath  + 'readme.txt', 'a')
line = 'something \n'
f.write(line)
f.close()    
#%% open file and writes with context manager
with open(folderpath  + 'readme.txt', 'a') as f:
    line = 'something \n'
    f.write(line)
#%% reads file
with open(folderpath  + 'readme.txt', 'r') as f:
    line = f.readlines()
#%% save a numpy array to a file
a = np.array([1,2,3, 5])
np.savetxt(folderpath  + 'readme.txt', a, fmt = '%.2f')
#%% appends a numpy array to a file as a row
with open(folderpath  + 'readme.txt', 'a') as f:
    a = np.array([1,2,3])
    np.savetxt(f, a, fmt = '%.2f', newline = ' ')
    f.write('\n')
