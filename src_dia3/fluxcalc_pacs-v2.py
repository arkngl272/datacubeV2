#!/usr/bin/env python
# -*- coding: utf8 -*-
print ''
print ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print ' %                                                                  % '
print ' %         This routine measures the flux of a source from          % '
print ' %          Herschel data, imaged with PACS and calcuates           % '
print ' %    the flux uncertainty taking all error sources into account    % '
print ' %                                                                  % '
print ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% '
print ''
print ''
print 'Importing modules...'
import numpy
from astropy.io import fits
import pyregion
import pyregion._region_filter as filter
print 'DONE!'
#	first open the fits file, reads the header and the planes
img = "../../data/data-hershel/m51_pacs_160.fits"
########################################################################################
#                                                                                      #
#   The region file must contain firstly the region where the flux is measured, then   #
#   a number of regions where the sky is computed. They must be in image coordinates.  #
#                                                                                      #
########################################################################################
reg_name = '../../data/data-hershel/regiones/pacR.reg'
hdulist = fits.open(img)
flxmap = hdulist[1].data
errmap = hdulist[3].data
header = hdulist[1].header
nx = header["NAXIS1"]
ny = header["NAXIS2"]
size = (ny, nx)		#	this is the size of the image
#	first of all measure the on-source flux
source =  pyregion.open(reg_name)
nsky = len(source) - 1	#	total apertures to measure the sky level
del source[1:]
smask = source.get_mask(shape=size)
flux = 0.0
errsq = 0.0
on_s_pix = 0
for i in range(nx):
	for n in range(ny):
		if smask[n,i] == True:
			flux = flux + flxmap[n,i]
			errsq = errsq + errmap[n,i]**2
			on_s_pix = on_s_pix + 1  
#	now does the sky measurements
sky_ap = numpy.zeros(nsky)
pix_ap = numpy.zeros(nsky)
sky_pix = []		# This is the array which will contain all the pixels' values in the sky apertures
for k in range(nsky):
	sky =  pyregion.open(reg_name)
	del sky[0:k+1]
	del sky[1:nsky]
	bmask = sky.get_mask(shape=size)
	sflx = 0.0
	pix = 0
	for i in range(nx):
		for n in range(ny):
			if bmask[n,i] == True:
				sflx = sflx + flxmap[n,i]
				pix = pix + 1  
				sky_pix.append(flxmap[n,i])
	sky_ap[k] = sflx
	pix_ap[k] = pix
sk_p_pix = sum(sky_ap)/sum(pix_ap)	#	this is the average sky per pixel
skyave = sky_ap/pix_ap				#	average sky per pixel
tflux = flux-on_s_pix*sk_p_pix		#	this is the background subtracted flux 
#	Now calculating the various uncertainties components
#  ######################################################
#	this is the calibration uncertainty
err1 = 0.05*tflux	
err2 = numpy.sqrt(errsq)
#	this is the uncertainty on the sky measurement. It includes correlated and uncorrelated noise
sig_unc = (numpy.std(sky_pix))**2 * on_s_pix
sig_corr = (numpy.std(skyave) * on_s_pix)**2
#print 'Uncorrelated noise: ',numpy.sqrt(sig_unc),' Jy'
#print 'Correlated noise: ',numpy.sqrt(sig_corr),' Jy'
err3 = numpy.sqrt(sig_unc+sig_corr)
tot_err = numpy.sqrt(err1**2 + err2**2 + err3**2)
rel_err = tot_err/tflux*100
print 'err1=', err1
print 'err2=', err2
print 'err3=', err3
print 'Background-subtracted flux:', tflux,'±',tot_err,' Jy  (',rel_err,' %)'

print 'Done!!!'


[10.81,5.27,7.00,7.00,7.05]