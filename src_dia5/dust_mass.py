#!/usr/bin/env python
# -*- coding: utf8 -*-

#################################################################
#                                                               #
#  This program calculates Dust Mass in a region of the Galaxy  #
#                                                               #
#################################################################

import pyregion 
import numpy as np
import astropy.io.fits as fits
import pyregion._region_filter as filter

pc2cm = 3.08568e18			#Conversion factor from pc to cm
img="../src/results_map.fits"
hdulist =fits.open(img)
header = hdulist[0].header
cube = hdulist[0].data[:,:,:]
nx = header ["NAXIS1"]
ny = header ["NAXIS2"]
nl = header ["NAXIS3"]
size = (ny, nx)

source =pyregion.open("../data/dust_gas_results.reg")
tot_sed = np.zeros(nl)
mask = source.get_mask(shape=size)

for m in (1, nl-1):
     for i in range(nx):
         for n in range(ny):
		if mask[n,i] == True and cube[m,n,i]>0:
			tot_sed[m] = tot_sed[m] + cube[m,n,i]

Dmass=tot_sed[1]

print "Dust Mass:","{:.2e}".format(Dmass),"[M"u'\u2609'"]"
