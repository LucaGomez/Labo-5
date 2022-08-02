# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 19:55:39 2022

@author: luo
"""
import numpy as np
import matplotlib.pyplot as plt
from IPython import get_ipython
from scipy.io import loadmat
from scipy.optimize import curve_fit
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from scipy import ndimage

from IPython import get_ipython

get_ipython().run_line_magic('matplotlib', 'qt5')
plt.rcParams.update({'font.size': 12})

#%%

#COMENTE TODO ESTO PORQUE ESTAN LAS CALIBRACIONES DE CADA DIA ABAJO

'''
#Calibracion escala espacial px vs cm

imagen_escala = plt.imread(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 3\calibracion sumidero L72 T112 W1136 H800.jpg')
#imagen_escala = plt.imread(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\escala perladoasd2.jpg')

plt.figure()
plt.imshow(imagen_escala)
plt.xlabel('Posicion X (px)')
plt.ylabel('Posicion Y (px)')

ptos= plt.ginput(2)
ptos=np.array(ptos)

#la roto para que quede mas o menos bien 

p_x=ptos[:,0]
p_y=ptos[:,1]
deltax=abs(p_x[1]-p_x[0])
deltay=abs(p_y[1]-p_y[0])
ang=abs(90-abs(180*np.arctan(deltay/deltax)/np.pi))

esc=ndimage.rotate(imagen_escala, ang)

#%%

#grafico la figura para tomar los puntos
mesc=esc[600:1000,400:900]
plt.figure()
plt.imshow(mesc)
plt.xlabel('X (pixel)',fontsize=14)
plt.ylabel('Y (pixel)',fontsize=14)
plt.title('Calibración espacial')

#%%

col=mesc[10:600,211]
x=np.linspace(1,len(col),len(col))
plt.figure()
plt.plot(x,col, 'r',label='Datos obtenidos')
plt.xlabel("Posicion en eje Y (mm)")
plt.ylabel("Intensidad")
plt.title("Interfranja")
plt.legend(loc='upper right')
plt.grid()
plt.show()

#%%

puntos_esc = plt.ginput(9)
puntos_esc = np.array(puntos_esc)
puntos_esc1 = plt.ginput(2)
puntos_esc1 = np.array(puntos_esc1)

coord_esc=np.hstack((puntos_esc[:,0],puntos_esc1[:,0]))
coord_escx=np.hstack((puntos_esc[:,1],puntos_esc1[:,1]))
diff_esc= np.diff(coord_esc)

# Promedio de los valores
escala=np.mean(diff_esc)
escala_error=np.std(diff_esc)/np.sqrt(len(diff_esc)) 

print('la escala de px/5mm es ',escala ,'+-',escala_error)

#%%

coord_esc=np.array([136.57291667, 166.61666667, 195.35416667, 225.39791667, 256.3125, 286.79166667, 315.52916667, 344.70208333,375.61666667, 406.96666667, 435.70416667, 465.74791667])
coord_escx=np.array([118.67057394, 112.53132351, 112.75870315, 103.20875804,102.07185981,  98.43378548,  99.11592442,  97.97902619,102.07185981, 106.61945272, 103.66351733,  99.57068371])

plt.figure()
plt.plot(x,col, 'r',label='Perfil',zorder=1)
plt.scatter(coord_esc,coord_escx,label='Mínimos',zorder=2,color='b')
plt.xlabel("Posicion en eje Y (pixel)",fontsize=14)
plt.ylabel("Intensidad",fontsize=14)
plt.title("Calibración espacial")
plt.legend(loc='upper right')
plt.show()

'''





#%%

#Calibracion dia 1
#escp=31.1183248154386002/5*10 #px por cm
#escp_error=0.3961946214313535*10

#Calibracion dia 2
escp=29.93713202275223/5*10 #px por cm
escp_error=0.13440475595279747*10

#Calibracion dia 3: dos vortices
#escp=21.55473505638159
#escp_error=0.23189202751466623

#Calibracion dia 3: sumidero
#escp=31.34482526881721
#escp_error=0.30160302831451574

#%%

#Funciones corro una vez y listo

#convertir en polares

def converTo_polares(MXY,eMXY):
    MRT= np.empty((len(MXY), 4), float)
    eMRT= np.empty((len(eMXY), 4), float)
    for i in range(len(MXY)):
        x = MXY[i][0]
        y = MXY[i][1]
        u = MXY[i][2]
        v = MXY[i][3]
        ex = eMXY[i][0]
        ey = eMXY[i][1]
        eu = eMXY[i][2]
        ev = eMXY[i][3]
        
        #el modulo de (x,y)
        r = np.sqrt(x**2 + y**2)
        er = np.sqrt((2*x*ex/r)**2+(2*y*ey/r)**2)
        
        #el angulo en radianes de (x,y)
        t = np.angle(x + y*1j)
        et = np.sqrt((2*x*ex/t)**2+(2*y*ey/t)**2)
        
        MRT[i][0] = r
        MRT[i][1] = t
        eMRT[i][0] = er
        eMRT[i][1] = et
        
        #defino las velocidades como 0 inicialmente
        vrad = 0
        vtan = 0
        evrad = 0
        evtan = 0
        
        #hay que tener cuidado con la posicion de la matriz que esta en el origen, pues ahi r = 0
        if r != 0 :
            #vel radial
            vrad = (1/r)*(x*u + y*v) 
            evrad = np.sqrt((vrad*eu/u)**2+(vrad*ev/v)**2)
            #vel tangencial r*(titapunto)
            vtan = (1/r)*(x*v - y*u)
            evtan = np.sqrt((vtan*eu/u)**2+(vtan*ev/v)**2)
            
        MRT[i][2] = vrad
        MRT[i][3] = vtan
        eMRT[i][2] = evrad
        eMRT[i][3] = evtan
    return MRT,eMRT


#promedio y desviacion estandar de la dispersion 

def homogVtan_R(MRT, p):
    heightMRT = len(MRT)
    #MRT es la matriz que contiene en la primer y cuarta columna los datos a promediar
    #p es la cantidad de datos que quiero obtener del promediado, si p->inf se reobtiene los datos originales
    
    #cuantos puntos quiero para la curva?, se define al definir la cantidad de bloques en los que se promedia
    blockBorders = np.linspace(min(MRT[:,0]), max(MRT[:,0]), p)
    #se devuelve esta lista transformada a array
    r_prom_Vtan = []
    r_prom_Vtan_ds=[]
    
    for j in range(p):
        rs_toprom = []
        vt_toprom = []
        for i in range(heightMRT):
            if j!= p - 1 and blockBorders[j] < MRT[i][0] < blockBorders[j+1]:     
                vt_toprom.append(MRT[i][3])
                rs_toprom.append(MRT[i][0])     
            #print(blockBorders[j])
            if j == p - 1 and blockBorders[j-2] < MRT[i][0] :
                vt_toprom.append(MRT[i][3])
                rs_toprom.append(MRT[i][0])               
        
        if len(rs_toprom) != 0:
            r_prom_Vtan.append([np.mean(rs_toprom),np.mean(vt_toprom)])
            r_prom_Vtan_ds.append([np.std(rs_toprom),np.std(vt_toprom)])
        
        matrizRp_Vtp = np.array(r_prom_Vtan)
        matrizRp_Vtp_ds = np.array(r_prom_Vtan_ds)
    return matrizRp_Vtp, matrizRp_Vtp_ds


def homogVr_R(MRT, p):
    heightMRT = len(MRT)
    #MRT es la matriz que contiene en la primer y cuarta columna los datos a promediar
    #p es la cantidad de datos que quiero obtener del promediado, si p->inf se reobtiene los datos originales
    
    #cuantos puntos quiero para la curva?, se define al definir la cantidad de bloques en los que se promedia
    blockBorders = np.linspace(min(MRT[:,0]), max(MRT[:,0]), p)
    #se devuelve esta lista transformada a array
    r_prom_Vr = []
    r_prom_Vr_ds=[]
    
    for j in range(p):
        rs_toprom = []
        vr_toprom = []
        for i in range(heightMRT):
            if j!= p - 1 and blockBorders[j] < MRT[i][0] < blockBorders[j+1]:     
                vr_toprom.append(MRT[i][2])
                rs_toprom.append(MRT[i][0])     
            #print(blockBorders[j])
            if j == p - 1 and blockBorders[j-2] < MRT[i][0] :
                vr_toprom.append(MRT[i][2])
                rs_toprom.append(MRT[i][0])               
        
        if len(rs_toprom) != 0:
            r_prom_Vr.append([np.mean(rs_toprom),np.mean(vr_toprom)])
            r_prom_Vr_ds.append([np.std(rs_toprom)/np.sqrt(len(rs_toprom)),np.std(vr_toprom)/np.sqrt(len(vr_toprom))])
        
        matrizRp_Vrp = np.array(r_prom_Vr)
        matrizRp_Vrp_ds = np.array(r_prom_Vr_ds)
    return matrizRp_Vrp, matrizRp_Vrp_ds


#%%

#Abro archivos a analizar

#promedio0=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\PIVlab promedio agua.txt')
promedio0=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\PIVlab.27.2.1.txt')
archivo='27.2 % - vortice 1'

matrizXY0 = np.loadtxt(promedio0,delimiter= ',', skiprows= 3, usecols = (0,1,2,3))

px0,py0=matrizXY0[:,0], matrizXY0[:,1] #posicion de vx y vy
pvx0,pvy0=matrizXY0[:,2], matrizXY0[:,3] #valores de vx y vy

#std=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\PIVlab.std.27.2.1.txt')
#matrizSTD = np.loadtxt(std,delimiter= ',', skiprows= 3, usecols = (0,1,2,3))
#epx0,epy0=matrizSTD[:,0], matrizSTD[:,1] #posicion de vx y vy
#epvx0,epvy0=matrizSTD[:,2], matrizSTD[:,3] #valores de vx y vy

plt.figure()
plt.quiver(px0,py0,pvx0,pvy0)
plt.xlabel('X [pixel]',fontsize=14)
plt.ylabel('Y [pixel]',fontsize=14)
plt.title('{}'.format(archivo))


#%%

centro = plt.ginput(1)
x0,y0 = np.array(centro[0])
ex0,ey0=0,0

#%%

esct=34.4321 #cambiar esto con la calibracion tiempo vs frame despues, por ahora tengo la automatica de 53 fps
esct_error=34.4321-60/(5.16-2.85)

epx0,epy0,epvx0,epvy0=0,0,0,0

x,y,vx,vy=(px0-x0)/escp, (py0-y0)/escp, pvx0*esct/escp, pvy0*esct/escp

ex,ey=np.sqrt((epx0/escp)**2+(ex0/escp)**2+(escp_error*x/escp)**2),np.sqrt((epy0/escp)**2+(ey0/escp)**2+(escp_error*y/escp)**2)
evx,evy=np.sqrt((epvx0*vx/pvx0)**2+(esct_error*vx/esct)**2+(escp_error*vx/escp)**2),np.sqrt((epvy0*vy/pvy0)**2+(esct_error*vy/esct)**2+(escp_error*vy/escp)**2)


plt.figure()
plt.quiver(x,y,vx,vy)
plt.ylabel('Posicion Y (cm)',fontsize=14)
plt.xlabel('Posicion X (cm)',fontsize=14)
plt.title('{}'.format(archivo))


#%%

#armo la matriz centrada y calibrada
Mxy,eMxy=np.array([x,y,vx,vy]),np.array([ex,ey,evx,evy])
M_xy,eM_xy=Mxy.T, eMxy.T

#obtengo las dos matrices a ajustar, la OG en polares y la promediada de la OG
M_rt,eM_rt = converTo_polares(M_xy,eM_xy)

matrizProm_Rp_Vtp,matrizProm_Rp_Vtp_sd = homogVtan_R(M_rt, 110)
matrizProm_Rp_Vtp_ev,matrizProm_Rp_Vtp_ev_sd = homogVtan_R(eM_rt, 110)
matrizProm_Rp_Vrp,matrizProm_Rp_Vrp_sd = homogVr_R(M_rt, 110)

#%%

#Estan negativos por el sentido de los vectores de v, no afecta en los resultados excepto el signo de Gamma (el vortice)
#para dejarlo bonito lo multiplico por -1

M_rt[:,3]=-1*M_rt[:,3]
matrizProm_Rp_Vtp[:,1]=-1*matrizProm_Rp_Vtp[:,1]

a_nan = np.where(np.isnan(M_rt[:,3])) 
M_rt[a_nan,3]=0
b_nan = np.where(np.isnan(matrizProm_Rp_Vtp[:,1])) 
matrizProm_Rp_Vtp[b_nan,1]=0
c_nan = np.where(np.isnan(matrizProm_Rp_Vtp_sd[:,1])) 
matrizProm_Rp_Vtp_sd[c_nan,1]=0
d_nan = np.where(np.isnan(matrizProm_Rp_Vrp[:,1])) 
matrizProm_Rp_Vrp[d_nan,1]=0

#%%

'''
#burgerking
def func(x,a,b,c,d):
    return ((a/(x-d))*(1 - np.exp(-((x-d)**2/b**2))))+c

'''

def func(x,a,b,c):
    return ((a/(x-c))*(1 - np.exp(-((x-c)**2/b**2))))

#%%

#grafico la OG y la promediada

plt.figure()
plt.plot(M_rt[:,0], M_rt[:,3], '.',label='Velocidad tangencial',color='b')
plt.plot(matrizProm_Rp_Vtp[:,0],matrizProm_Rp_Vtp[:,1], '.', label='Promedio-velocidad tangencial',color='r')
plt.xlabel("Radios [cm]", fontsize=15)
plt.ylabel("Velocidad tangencial [cm/s]", fontsize = 15)
plt.title('$V_θ$(r) - {}'.format(archivo))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc='upper right')


plt.figure()
plt.plot(M_rt[:,0], M_rt[:,2], '.',label='Velocidad radial',color='b')
plt.plot(matrizProm_Rp_Vrp[:,0],matrizProm_Rp_Vrp[:,1], '.', label='Promedio-velocidad radial',color='r')
plt.xlabel("Radios [cm]", fontsize=15)
plt.ylabel("Velocidad radial [cm/s]", fontsize = 15)
plt.title('$V_r$(r) - {}'.format(archivo))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc='upper right')

#%%
#Grafico de Vel. Tangencial en funcion del radio PERO PARA EL PROMEDIADO ESPACIAL
#fijo rmin para deshacerme del hueco alrededor del vortice (si es que hay, sino rmin=0)
#fijo rmax para deshacerme de efecto de borde

rmin= 0.0 #max o menos a ojo 
rmax= 6.5

#filtro y tomo el rango entre rmin y rmax para el ajuste
subMatrix0 = matrizProm_Rp_Vtp[matrizProm_Rp_Vtp[:,0] >rmin]
subMatrix = subMatrix0[subMatrix0[:,0] <rmax]
subMatrix_sd0=matrizProm_Rp_Vtp_sd[matrizProm_Rp_Vtp[:,0] >rmin]
subMatrix_sd = subMatrix_sd0[subMatrix0[:,0] <rmax]

#ajuste burger al promedio de la dispersion con error en Y=desviacion estandar del promedio
popt_prom, pcov_prom = curve_fit(func, subMatrix[:,0], subMatrix[:,1],sigma= subMatrix_sd[:,1],absolute_sigma=True)
rD_prom = np.linspace(min(subMatrix[:,0]),max(subMatrix[:,0]),len(subMatrix[:,0]))
perr_prom=np.diag(pcov_prom)

plt.figure()
plt.title('$V_θ$(r) promedio - {}'.format(archivo))
#plt.plot(matrizProm_Rp_Vtp[:,0], matrizProm_Rp_Vtp[:,1], '.',label='datos')
plt.errorbar(matrizProm_Rp_Vtp[:,0], matrizProm_Rp_Vtp[:,1], xerr = matrizProm_Rp_Vtp_sd[:,0], yerr =matrizProm_Rp_Vtp_sd[:,1],fmt='.',label='Datos',zorder=1,color='b')
plt.plot(rD_prom,func(rD_prom,*popt_prom), label='Ajuste',zorder=2,color='r')
plt.xlabel("Radios [cm]", fontsize=15)
plt.ylabel("Velocidad tangencial [cm/s]", fontsize = 15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()

from scipy import stats

fitted_function=func(subMatrix[:,0],*popt_prom)
dof=len(subMatrix[:,1])-len(popt_prom)
subMatrix_sd[subMatrix_sd == 0] = 0.1
chisquare=sum((subMatrix[:,1]-fitted_function)**2/(subMatrix_sd[:,1]**2))
chi2red=chisquare/dof
p=1-stats.chi2.cdf(chisquare,dof)


gamma=popt_prom[0]*2*np.pi
egamma=perr_prom[0]*2*np.pi

print('Chi-square es', chisquare)
print('Chi-square reducido es', chi2red)
print('Gamma es', gamma,'+-', egamma)
print('rc o sqrt(4nu/alpha) es', popt_prom[1],'+-', perr_prom[1])
print('r0 es', popt_prom[2],'+-', perr_prom[2])
#print('CTE es', popt_prom[3],'+-', perr_prom[3])



#%%

#Graficar con bandas de prediccion y confianza, un poco demasiado capaz

import uncertainties as unc
import uncertainties.unumpy as unp

# calculate parameter confidence interval
a,b,c,d=unc.correlated_values(popt_prom, pcov_prom)

# calculate regression confidence interval
px = np.linspace(matrizProm_Rp_Vtp[0,0],matrizProm_Rp_Vtp[len(matrizProm_Rp_Vtp[:,0])-1,0], 250)
py=((a/(px-d))*(1 - unp.exp(-((px-d)**2/b**2))))+c
nom = unp.nominal_values(py)
std = unp.std_devs(py)


def predband(x, xd, yd, p, func, conf=0.95):
    # x = requested points
    # xd = x data
    # yd = y data
    # p = parameters
    # func = function name
    alpha = 1.0 - conf    # significance
    N = xd.size          # data sample size
    var_n = len(p)  # number of parameters
    # Quantile of Student's t distribution for p=(1-alpha/2)
    q = stats.t.ppf(1.0 - alpha / 2.0, N - var_n)
    # Stdev of an individual measurement
    se = np.sqrt(1. / (N - var_n) * \
                 np.sum((yd - func(xd,*p)) ** 2))
    # Auxiliary definitions
    sx = (x - xd.mean()) ** 2
    sxd = np.sum((xd - xd.mean()) ** 2)
    # Predicted values (best-fit model)
    yp = func(x,*p)
    # Prediction band
    dy = q * se * np.sqrt(1.0+ (1.0/N) + (sx/sxd))
    # Upper & lower prediction bands.
    lpb, upb = yp - dy, yp + dy
    return lpb, upb

lpb, upb = predband(px,matrizProm_Rp_Vtp[:,0], matrizProm_Rp_Vtp[:,1],popt_prom, func, conf=0.95)

plt.figure()
plt.errorbar(matrizProm_Rp_Vtp[:,0], matrizProm_Rp_Vtp[:,1], xerr = matrizProm_Rp_Vtp_sd[:,0], yerr =matrizProm_Rp_Vtp_sd[:,1],fmt='.',label='Datos',zorder=1,color='b')
plt.plot(rD_prom,func(rD_prom,*popt_prom), label='Ajuste',zorder=2,color='r')
# uncertainty lines (95% confidence)
plt.plot(px, nom - 1.96 * std, '--', color='purple',lw=2,label='Banda de confianza del 95%',zorder=2)
plt.plot(px, nom + 1.96 * std, '--', color='purple',zorder=2)
# prediction band (95% confidence)
plt.plot(px, lpb, 'g--', lw=1,label='Banda de predicción del 95%',zorder=2)
plt.plot(px, upb, 'g--',zorder=2)
plt.title('$V_θ$(r) promedio - {}'.format(archivo))
plt.xlabel("Radios [cm]", fontsize=15)
plt.ylabel("Velocidad tangencial [cm/s]", fontsize = 15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc='upper right')
plt.show()

#%%

#algunas dimenciones

vol27_2=16*np.pi*((26.5/(2*np.pi))**2)

vol50=11.2*np.pi*((32.5/(2*np.pi))**2)

a_circular=np.pi*((76.5/(2*np.pi))**2)

h_27_2=vol27_2/a_circular
h_50=vol50/a_circular

print('altura de 27.2% es', h_27_2)
print('altura de 50% es', h_50)
print('radio del recipiente grande es',(76.5/(2*np.pi)))

#%%
#guardo los promedios

data_to_save = np.vstack((matrizProm_Rp_Vtp[:,0], matrizProm_Rp_Vtp[:,1],matrizProm_Rp_Vtp_sd[:,0],matrizProm_Rp_Vtp_sd[:,1] ))
np.savetxt(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\promedios agua bis.txt', data_to_save)

#%%
#guardo los datos de ajuste en un archivo txt
fsalida=open(r'C:\Users\luo\OneDrive\文档\Labo 5\Fluidos\Dia 2\27.2\Dato ajustes nuevo.txt','a')
'''
fsalida.write('Archivo ,')
fsalida.write('Escala px/cm,')
fsalida.write('Error de escala px/cm,')
fsalida.write('Escala fps,')
fsalida.write('Error de escala fps,')
fsalida.write('centro x0,')
fsalida.write('centro y0,')
fsalida.write('rmax,')
fsalida.write('rmin,')
fsalida.write('Ajuste promedio Chi-squared,')
fsalida.write('Ajuste promedio Chi-squared red,')
fsalida.write('Ajuste promedio pvalue,')
fsalida.write('Ajuste promedio Gamma/2pi,')
fsalida.write('Ajuste promedio error Gamma/2pi,')
fsalida.write('Ajuste promedio sqrt(4nu/alpha),')
fsalida.write('Ajuste promedio error sqrt(4nu/alpha),')
fsalida.write('Ajuste promedio CTE,')
fsalida.write('Ajuste promedio error CTE,')
fsalida.write('Ajuste promedio r0,')
fsalida.write('Ajuste promedio error r0'+'\n')
'''

fsalida.write('%10.10f' % (escp)+',' )
fsalida.write('%10.10f' % (escp_error)+',')
fsalida.write('%10.10f' % (esct)+',')
fsalida.write('%10.10f' % (esct_error)+',')
fsalida.write('%10.10f' % (x0)+',')
fsalida.write('%10.10f' % (y0)+',')
fsalida.write('%10.10f' % (rmax)+',')
fsalida.write('%10.10f' % (rmin)+',')
fsalida.write('%10.10f' % (chisquare)+',')
fsalida.write('%10.10f' % (chi2red)+',')
fsalida.write('%10.10f' % (p)+',')
fsalida.write('%10.10f' % (gamma)+',')
fsalida.write('%10.10f' % (egamma)+',')
fsalida.write('%10.10f' % (popt_prom[1])+',')
fsalida.write('%10.10f' % (perr_prom[1])+',')
fsalida.write('%10.10f' % (popt_prom[2])+',')
fsalida.write('%10.10f' % (perr_prom[2])+',')
fsalida.write('%10.10f' % (popt_prom[3])+',')
fsalida.write('%10.10f' % (perr_prom[3])+'\n')
fsalida.close()


#%%

#Rankine, creo que no llegamos a escribir esto ademas no aporta mucho

def Rankine1(x,a,b):
    c=np.sqrt(2.7)    
    return (a/(c**2))*x+b
def Rankine2(x,a,b):
    b=0.5
    return a/x

#rc=popt_prom[1]
rc=np.sqrt(2.65)
#Tomo dif array para el ajuste de rankine
rankMatrix1 = subMatrix[subMatrix[:,0]<rc**2]
rankMatrix_sd1= subMatrix_sd[subMatrix[:,0] <rc**2]
rankMatrix2 = subMatrix[subMatrix[:,0] >rc**2]
rankMatrix_sd2= subMatrix_sd[subMatrix[:,0] >rc**2]

popt1,pcov1 = curve_fit(Rankine1,rankMatrix1[:,0], rankMatrix1[:,1],sigma= rankMatrix_sd1[:,1],absolute_sigma=True)
popt2,pcov2 = curve_fit(Rankine2,rankMatrix2[:,0], rankMatrix2[:,1],sigma= rankMatrix_sd2[:,1],absolute_sigma=True)
perr1, perr2=np.diag(pcov1), np.diag(pcov2)

r_rank1 = np.linspace(min(rankMatrix1[:,0]),max(rankMatrix1[:,0]),len(rankMatrix1[:,0]))
r_rank2 = np.linspace(min(rankMatrix2[:,0]),max(rankMatrix2[:,0]),len(rankMatrix2[:,0]))

fitted_values1=Rankine1(rankMatrix1[:,0],*popt1)
fitted_values2=Rankine2(rankMatrix2[:,0],*popt2)

plt.figure()
plt.title('Rankine: $V_θ$(r) promedio - {}'.format(archivo))
plt.errorbar(matrizProm_Rp_Vtp[:,0], matrizProm_Rp_Vtp[:,1], xerr = matrizProm_Rp_Vtp_sd[:,0], yerr =matrizProm_Rp_Vtp_sd[:,1],fmt='.',label='Datos',zorder=1,color='b')
plt.plot(r_rank1,Rankine1(r_rank1,*popt1), label='Ajuste $r<r_c$',color='r',zorder=2)
plt.plot(r_rank2,Rankine2(r_rank2,*popt2), label='Ajuste $r>r_c$',color='g',zorder=2)
plt.xlabel("Radios [cm]", fontsize=15)
plt.ylabel("Velocidad tangencial [cm/s]", fontsize = 15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()


dof1=len(rankMatrix1[:,1])-len(popt1)
dof2=len(rankMatrix2[:,1])-len(popt2)
rankMatrix_sd1[rankMatrix_sd1 == 0] = 0.1
rankMatrix_sd2[rankMatrix_sd2 == 0] = 0.1
chisquare1=sum((rankMatrix1[:,1]-fitted_values1)**2/(rankMatrix_sd1[:,1]**2))
chisquare2=sum((rankMatrix2[:,1]-fitted_values2)**2/(rankMatrix_sd2[:,1]**2))
chi2red1=chisquare1/dof1
chi2red2=chisquare2/dof2

print('Chi-square 1 es', chisquare1)
print('Chi-square 2 es', chisquare2)
print('Chi-square reducido 1 es', chi2red1)
print('Chi-square reducido 2 es', chi2red2)
print('Gamma1/2pi es', popt1[0],'+-', perr1[0])
print('Gamma2/2pi es', popt2[0],'+-', perr2[0])
print('CTE 1 es', popt1[1],'+-', perr1[1])
#print('CTE 2 es', popt2[1],'+-', perr2[1])
