#!/usr/bin/env python
# -*- coding: utf8 -*-

import os
import time
import pyregion 
import numpy as np
import astropy.io.fits as fits
import pyregion._region_filter as filter

### Limpiando la pantalla de la terminal  #################################

os.system ("clear")

### Lectura de la imagen fits #############################################

modelname="../data/oa30_i0-zip/t5_p0_q0_oa30_R10_Mcl0.97_i0_total.fits"
hdulist =fits.open(modelname)
header = hdulist[0].header
cube = hdulist[0].data[:,:,:]
nx = header ["NAXIS1"]
ny = header ["NAXIS2"]
nl = header ["NAXIS3"]

### Datos de las regiones #################################################

source =pyregion.open("../data/oa30_i0-zip/regions.reg")
nreg = len(source)
size = (ny, nx)
tot_sed = np.zeros((nl,nreg))

### Calculando Total Spectral Energy Distribution #########################

print("\n\nProcesando Regiones: \n")

start = time.time()
for reg in range(nreg):
    print ("Trabajando en region: %s"%reg)    
    shapes = pyregion.ShapeList([source[reg]])
    mask = shapes.get_mask(shape=size)

    for wln in range(nl):
        for x in range(nx):
            for y in range(ny):              
                if mask[y,x]:
                    tot_sed[wln,reg] = tot_sed[wln,reg] + cube[wln,y,x]

### Calculando Total SED x Region #########################################

seds=[]
for i in range(nreg-1):
    seds.append(tot_sed[:,i]-tot_sed[:,i+1])

end = time.time()


### Escribiendo datos #####################################################

fname = "sed.dat"
print("Escriendo datos en %s !!\n"%fname)
fmt = "%12.4e %12.4e %12.4e %12.4e %12.4e \n"
ff = open(fname,'w')
ff.write("#reg1 reg2 reg3 reg4 tot \n")
for v in zip(seds[0], seds[1],seds[2],seds[3], tot_sed[:,0]):
    ff.write(fmt%v)
ff.close()

print("\nFin del proceso, Tiempo total: %s\n"%(end - start))

