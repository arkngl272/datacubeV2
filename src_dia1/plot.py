#!/usr/bin/env python
# -*- coding: utf8 -*-

### Importando modulos  ###################################################

import os 

### Limpiando la pantalla de la terminal  #################################

os.system ("clear")

print ''
print ''
print 'Importando modulos...'
import numpy as np
import astropy.io.fits as fits
import astropy.io.ascii as ascii
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
rc('text', usetex=True)

### Lleyendo los datos de los archivos necesarios  ########################

print 'Leyendo datos...'

Dist = 10.

fitsimg = "../data/oa30_i0-zip/t5_p0_q0_oa30_R10_Mcl0.97_i0_total.fits"
hdulist = fits.open(fitsimg)
header = (hdulist[0].header)
pix = header["CDELT1"]

sedfile = "../data/oa30_i0-zip/t5_p0_q0_oa30_R10_Mcl0.97_i0_sed.dat"
tdata = ascii.read(sedfile)
lamb = tdata['col1']
tflx = tdata['col2']
agn = tdata['col3']
dflx = tflx - agn
conv = pix/(Dist*1e6/206265.)

regfile = "sed.dat"
regdata = ascii.read(regfile)
reg1 = regdata[0][:]*conv**2
reg2 = regdata[1][:]*conv**2
reg3 = regdata[2][:]*conv**2
reg4 = regdata[3][:]*conv**2
tot = regdata[4][:]*conv**2

yma = max(tflx)*1.2
ymi = yma*0.001
xma = max(lamb)*1.01
xmi = min(lamb)*0.97

### Graficando datos  ####################################################

print 'Graficando...'

plt.figure(figsize=(15,10))
plt.loglog(lamb,tflx,linestyle='-',color='black',label= "Total Flx")
plt.loglog(lamb,agn,linestyle='-',color='blue',label= "AGN")
plt.loglog(lamb,dflx,linestyle='-',color='red',label="Total Flx - AGN")
plt.loglog(lamb,tot,linestyle='--',color='black',linewidth=4,label= "Suma total")
plt.loglog(lamb,reg1,linestyle='-',color='green',label= "Region 1")
plt.loglog(lamb,reg2,linestyle='-',color='magenta',label= "Region 2")
plt.loglog(lamb,reg3,linestyle='-',color='purple',label= "Region 3")
plt.loglog(lamb,reg4,linestyle='-',color='gray',label= "Region 4")
plt.xlabel(r'$\lambda$    [$\mu$m]',fontsize=20)
plt.ylabel(r'$\lambda$F$_\lambda$   [W m$^{-2}$]',fontsize=20)
plt.tick_params(labelsize=20)
plt.xlim(xmi,xma)
plt.ylim(ymi,yma)
fig = plt.gcf()


plt.legend()
plt.show()
