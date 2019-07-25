#!/usr/bin/env python
# -*- coding: utf8 -*-

#################################################################
#                                                               #
#      This program calculates Dust Mass in the Galaxy          #
#                                                               #
#################################################################

import numpy as np
import astropy.io.fits as fits

img = "results_map.fits"
hdulist = fits.open(img)
dmasspx=hdulist[0].data[1,:,:]
header=hdulist[0].header
nx=header["NAXIS1"]
ny=header["NAXIS2"]

DMASS=0.0
for i in range(nx):
	for j in range(ny):
		if dmasspx[j,i] >0.:
			DMASS=DMASS+dmasspx[j,i]

print "M51 Dust Total Mass:","{:.2e}".format(DMASS),"[M"u'\u2609'"]"

