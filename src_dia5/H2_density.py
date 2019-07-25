#!/usr/bin/env python
# -*- coding: utf8 -*-

####################################################
#                                                  #
#    This routine calculates molecular Gas mass    #
#        (H2) in a M51 region from CO data         #
#                                                  #
####################################################

import pyregion 
import numpy as np
import astropy.io.fits as fits
import pyregion._region_filter as filter

pc2cm = 3.08568e18			#Conversion factor from pc to cm
img="../src/fit_CO.fits"
hdulist =fits.open(img)
header = hdulist[0].header
cube = hdulist[0].data[:,:,:]
nx = header ["NAXIS1"]
ny = header ["NAXIS2"]
nl = header ["NAXIS3"]
size = (ny, nx)

source =pyregion.open("../data/dust_gas_CO.reg")
tot_sed = 0.0
mask = source.get_mask(shape=size)

for i in range(nx):
	for n in range(ny):
		if mask[n,i] == True and cube[3,n,i]>0:
			tot_sed = tot_sed + cube[3,n,i]

dnmol=2e20*tot_sed

nmol=(11.635528347*pc2cm)**2*dnmol	#Pixel 0.3 arcseg

H2mass=2*1.6726219e-27*nmol/2e30

Gasmoltotal=1.36*H2mass

print "Molecular Gas Mass:","{:.2e}".format(Gasmoltotal),"[M"u'\u2609'"]"

