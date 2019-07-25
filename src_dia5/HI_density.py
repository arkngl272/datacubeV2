#!/usr/bin/env python
# -*- coding: utf8 -*-

####################################################
#                                                  #
#    This routine calculates molecular Gas mass    #
#        	(HI) in a M51 region               #
#                                                  #
####################################################

import pyregion 
import numpy as np
import astropy.io.fits as fits
import pyregion._region_filter as filter

pc2cm = 3.08568e18			#Conversion factor from pc to cm
img="../src/fit_HI.fits"
hdulist =fits.open(img)
header = hdulist[0].header
cube = hdulist[0].data[:,:,:]
nx = header ["NAXIS1"]
ny = header ["NAXIS2"]
nl = header ["NAXIS3"]
size = (ny, nx)

source =pyregion.open("../data/dust_gas_HI.reg")
tot_sed = 0.0
mask = source.get_mask(shape=size)

for i in range(nx):
	for n in range(ny):
		if mask[n,i] == True and cube[2,n,i]>0:
			tot_sed = tot_sed + cube[2,n,i]

dnat=1.82e18*tot_sed

nat=dnat*(58.182295*pc2cm)**2		#Pixel 1.5 arcseg

Hmass=1.6726219e-27*nat/2e30

print "Atomic Mass:","{:.2e}".format(Hmass),"[M"u'\u2609'"]"

