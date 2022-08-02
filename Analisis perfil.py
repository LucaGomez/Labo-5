# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 17:37:25 2022

@author: luo
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage,signal
from IPython import get_ipython
from scipy.io import loadmat
from scipy.optimize import curve_fit
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from scipy import ndimage

from IPython import get_ipython

get_ipython().run_line_magic('matplotlib', 'qt5')
plt.rcParams.update({'font.size': 14})

#%%

imagen_escala = plt.imread(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 3\calibracion dosvortex.jpg')
#imagen_escala = plt.imread(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\escala perladoasd2.jpg')


plt.figure()
plt.imshow(imagen_escala)
plt.xlabel('Posicion X (px)')
plt.ylabel('Posicion Y (px)')

#%%

ptos= plt.ginput(2)
ptos=np.array(ptos)

#la roto para que quede mas o menos bien 

p_x=ptos[:,0]
p_y=ptos[:,1]
deltax=abs(p_x[1]-p_x[0])
deltay=abs(p_y[1]-p_y[0])

#%%
ang=abs(abs(180*np.arctan(deltay/deltax)/np.pi))

esc=ndimage.rotate(imagen_escala, ang)

matriz_esc=esc
#height width (ver los datos en prop de la foto) and intensidad o algo asi, dejarlo 1 por default, 
plt.figure()
plt.imshow(matriz_esc)

#%%

puntos_esc = plt.ginput(9)
puntos_esc = np.array(puntos_esc)

coord_esc=puntos_esc[:,0]

# La función np.diff() calcula la diferencia entre elmentos consecutivos
diff_esc= np.diff(coord_esc)

# Promedio de los valores
escala=np.mean( diff_esc)  # tomo el valor medio
escala_error=np.std(  diff_esc )/np.sqrt(len(diff_esc)) # tomo el error de la media

print(f'Escala: ( {escala}  ± {escala_error}  ) px/mm')

#%%
'''
#DA FEO CON EL KERNEL

s = [[1, 2, 1],  
     [0, 0, 0], 
     [-1, -2, -1]]

# Calcula la convolucion de la imagen con el kernel (H aplica el filtro en la dirección x y V en la dirección y)
H= signal.convolve2d(matriz_esc, s)
V= signal.convolve2d(matriz_esc, np.transpose(s))
img_bordes = (H**2 + V**2)**0.5

#grafica dicha convolución
plt.figure()
plt.imshow(img_bordes)
plt.title('Filtro de detección de bordes')

'''

#%%

#ANALIZAR ES POR ACA

imagenlento= plt.imread(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 3\IMG_5624_escala.png')

plt.figure()
plt.imshow(imagenlento)
plt.xlabel('Posicion X (px)')
plt.ylabel('Posicion Y (px)')

#%%

pts= plt.ginput(2)
pts=np.array(pts)

#la roto para que quede mas o menos bien 

ptx,pty=pts[:,0],pts[:,1]
deltaptx,deltapty=abs(ptx[1]-ptx[0]),abs(pty[1]-pty[0])

#%%
angpt=abs(abs(180*np.arctan(deltapty/deltaptx)/np.pi))

matriz=ndimage.rotate(imagenlento, angpt)

plt.figure()
plt.imshow(matriz)

#%%

eje_y = np.arange(matriz.shape[0])/escala
eje_x = np.arange(matriz.shape[1])/escala

#eje_y = np.arange(matriz.shape[0])/escala-x0
#eje_x = np.arange(matriz.shape[1])/escala-y0

plt.figure() 
plt.imshow(matriz, extent=[ eje_x.min(), eje_x.max(), eje_y.min(), eje_y.max()]  )
plt.xlabel('X [cm]')
plt.ylabel('Y [cm]')

centros=plt.ginput(1)
x0,y0=np.array(centros)

#%%

puntos= plt.ginput(20)
puntos=np.array(puntos)

xp,yp=puntos[:,0],puntos[:,1]

plt.figure()
plt.plot(xp,yp,label='Agitación lenta',color='b')
#plt.plot(xp,yp,label='Agitación rápida')
plt.xlabel('Radio (cm)',fontsize=14)
plt.ylabel('Altura (cm)',fontsize=14)
plt.title('Perfil del vórtice - agua')
plt.legend(loc='upper right')
plt.show()


#%%

imagenrapida= plt.imread(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 3\IMG_5624_escala.png')

plt.figure()
plt.imshow(imagenrapida)
plt.xlabel('Posicion X (px)')
plt.ylabel('Posicion Y (px)')

pts= plt.ginput(2)
pts=np.array(pts)

ptx,pty=pts[:,0],pts[:,1]
deltaptx,deltapty=abs(ptx[1]-ptx[0]),abs(pty[1]-pty[0])

#%%
angpt=abs(abs(180*np.arctan(deltapty/deltaptx)/np.pi))

matriz1=ndimage.rotate(imagenrapida, angpt)

plt.figure()
plt.imshow(matriz1)

#%%


#eje_y = np.arange(matriz1.shape[0])/escala-x1
#eje_x = np.arange(matriz1.shape[1])/escala-y1

eje_y1 = np.arange(matriz1.shape[0])/escala
eje_x1 = np.arange(matriz1.shape[1])/escala

plt.figure() 
plt.imshow(matriz1, extent=[ eje_x1.min(), eje_x1.max(), eje_y1.min(), eje_y1.max()]  )
plt.xlabel('X [cm]')
plt.ylabel('Y [cm]')

centros1=plt.ginput(1)
x1,y1=np.array(centros1)

#%%

puntos1= plt.ginput(20)
puntos1=np.array(puntos1)

xp1,yp1=puntos1[:,0],puntos1[:,1]

plt.figure()
plt.plot(xp,yp,label='Agitación lenta',color='b')
plt.plot(xp1,yp1,label='Agitación rápida',color='orange')
plt.xlabel('Radios [cm]',fontsize=14)
plt.ylabel('Altura [cm]',fontsize=14)
plt.title('Perfil del vórtice - agua')
plt.legend(loc='upper right')
plt.show()

#%%

def z_est(x,A,B):
    h=1
    return h-(A**2*B**4)/(98000*2*x)

def z_nest(x,A,B):
    return (A**2/(2*98000))*(x**2-x**4/B**4+x**6/(3*B**6))

#%%

n,n1=len(yp),len(yp1)
errory,errory1=np.full((n),escala_error),np.full((n1),escala_error)

popte, pcove = curve_fit(z_est, xp,yp,sigma=errory,absolute_sigma=True)
rdata= np.linspace(min(xp),max(yp),len(xp))
perre=np.diag(pcove)

poptne, pcovne = curve_fit(z_nest, xp1,yp1,sigma=errory1,absolute_sigma=True)
rdata1= np.linspace(min(xp1),max(yp1),len(xp1))
perrne=np.diag(pcovne)

plt.figure()
plt.title('z(r) - agua')
#plt.plot(matrizProm_Rp_Vtp[:,0], matrizProm_Rp_Vtp[:,1], '.',label='datos')
plt.errorbar(xp, yp, yerr =errory,fmt='.',label='Datos lento',zorder=1,color='b')
plt.errorbar(xp1, yp1, yerr =errory1,fmt='.',label='Datos rápido',zorder=1,color='orange')
plt.plot(rdata,z_est(rdata,*popte), label='Ajuste lento',zorder=2,color='r')
plt.plot(rdata1,z_nest(rdata1,*poptne), label='Ajuste rápido',zorder=2,color='g')
plt.xlabel("Radios [cm]", fontsize=14)
plt.ylabel("Altura [cm]", fontsize = 14)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.legend()

fitted_function,fitted_function1=z_est(xp,*popte),z_nest(xp1,*poptne)
dof,dof1=(len(yp)-len(popte)),(len(yp1)-len(poptne))
chisquare=sum((yp-fitted_function)**2/(errory**2))
chisquare1=sum((yp1-fitted_function1)**2/(errory1**2))
chi2red=chisquare/dof
chi2red1=chisquare1/dof1


print('Chi-square reducido es', chi2red, chi2red1)
print('wm est es', popte[0],'+-', perre[0])
print('wm no est es', poptne[0],'+-', perrne[0])
print('rm est es',  popte[1],'+-', perre[1])
print('rm no est es',  poptne[1],'+-', perrne[1])



#%%

fsalida=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 3\Datos perfil.txt','a')


fsalida.write('Escala px/cm,')
fsalida.write('Error de escala px/cm,')
fsalida.write('centro x0,')
fsalida.write('centro y0,')
fsalida.write('centro x1,')
fsalida.write('centro y1,')
fsalida.write('Chi-squared red est,')
fsalida.write('Chi-squared red no est,')
fsalida.write('wm est,')
fsalida.write('error wm est,')
fsalida.write('wm no est,')
fsalida.write('error wm no est,')
fsalida.write('rm est')
fsalida.write('error rm est')
fsalida.write('rm no est')
fsalida.write('error rm no est'+'\n')


fsalida.write('%10.10f' % (escala)+',' )
fsalida.write('%10.10f' % (escala_error)+',')
fsalida.write('%10.10f' % (x0)+',')
fsalida.write('%10.10f' % (y0)+',')
fsalida.write('%10.10f' % (x1)+',')
fsalida.write('%10.10f' % (y1)+',')
fsalida.write('%10.10f' % (chi2red)+',')
fsalida.write('%10.10f' % (chi2red1)+',')
fsalida.write('%10.10f' % (popte[0])+',')
fsalida.write('%10.10f' % (perre[0])+',')
fsalida.write('%10.10f' % (popte[1])+',')
fsalida.write('%10.10f' % (perre[1])+',')
fsalida.write('%10.10f' % (poptne[0])+',')
fsalida.write('%10.10f' % (perrne[0])+',')
fsalida.write('%10.10f' % (poptne[1])+',')
fsalida.write('%10.10f' % (perrne[1])+'\n')
fsalida.close()


