#!/usr/bin/env python
# -*- coding: utf8 -*-
#import ds9   #	doesn't bloody work
import numpy as np
from astropy.io import fits
import pyregion
import pyregion._region_filter as filter

# define some SPIRE-typical constants
beam = [465, 822, 1769]						# beam area, in arcsec^2
ebar = [1.0432705,1.0521442,1.1025926] 		# effective beam area ratio (beta=2 SPIRE UM, table 5.8, page 96; see also M. Smith email, 16/04/14)
cnoise = [0.0058,0.0063,0.0068]				# confusion noise, in Jy/beam, as from Nguyen et al. 2010
kk = 2.*np.pi/(360.*3600.)					# 1 arcsec in rad
img ="../../data/data-hershel/m51_spi_500.fits" #raw_input('Name of the SPIRE data: ')
reg_name ="../../data/data-hershel/regiones/regionspire500.reg" #raw_input('Name of the DS9 file containing the regions for flux and sky measurements: ')
hdulist = fits.open(img)
wavel = 350#input('Wavelength of the image: ')
if wavel == 250:
	wind = 0
if wavel == 350:
	wind = 1
if wavel == 500:
	wind = 2
header = hdulist[1].header
fluxmap = hdulist[1].data
errmap = hdulist[2].data
nx = header["NAXIS1"]
ny = header["NAXIS2"]
size = (ny, nx)		#	this is the size of the image
pixsz = np.abs(header["CDELT1"])*3600.
## converts from MJy/sr to Jy/pix
fluxmap = fluxmap*1e6*kk**2*pixsz**2
errmap =errmap *1e6 * (kk * pixsz) **2
##	first of all measure the on-source flux
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
			if str(fluxmap[n,i]) != 'nan':
				flux = flux + fluxmap[n,i]
				errsq = errsq + errmap[n,i]**2
				on_s_pix = on_s_pix + 1  	#on-source pixels
sky_ap = np.zeros(nsky)
pix_ap = np.zeros(nsky)
sky_pix = []		# 	This is the array which will contain all the pixels' values in the sky apertures
for k in range(nsky):
	sky = pyregion.open(reg_name)
	del sky[0:k+1]
	del sky[1:nsky]
	bmask = sky.get_mask(shape=size)
	sflx = 0.0
	pix = 0
	for i in range(nx):
		for n in range(ny):
			if bmask[n,i] == True:
				if str(fluxmap[n,i]) != 'nan':
					sflx = sflx + fluxmap[n,i]
					pix = pix + 1  
					sky_pix.append(fluxmap[n,i])
	sky_ap[k] = sflx
	pix_ap[k] = pix
sk_p_pix = sum(sky_ap)/sum(pix_ap)	#	this is the average sky per pixel
skyave = sky_ap/pix_ap				#	average sky per pixel
tflux = flux-on_s_pix*sk_p_pix		#	this is the background subtracted flux 
#	Now calculating the various uncertainties components
#  ######################################################
#	this is the calibration uncertainty
err1 = 0.07*tflux	
#	this is the instrumental uncertainty
err2 = np.sqrt(errsq)
#	this is the confusion noise
err3 = cnoise[wind]*np.sqrt(on_s_pix*pixsz**2/beam[wind])
#	this is the background uncertainty
#		instrumental background
eflx = 0.0
pix_b = 0
for k in range(nsky):
	sky = pyregion.open(reg_name)
	del sky[0:k+1]
	del sky[1:nsky]
	bmask = sky.get_mask(shape=size)
	for i in range(nx):
		for n in range(ny):
			if bmask[n,i] == True:
				if str(fluxmap[n,i]) != 'nan':
					eflx = eflx + errmap[n,i]**2	# sums pixels in the errormap, in the "sky"regions
					pix_b = pix_b + 1
err4_1 = np.sqrt(eflx)*on_s_pix/pix_b
print "Instrumental background:",err4_1," Jy"
#		confusion background
err4_2 = cnoise[wind]*np.sqrt(pixsz**2/beam[wind])*on_s_pix/np.sqrt(pix_b)
print "Confusion background:",err4_2," Jy"
#		large-scale background
err4_3 = np.std(sky_ap)/np.sqrt(nsky)
print "Large-scale background:",err4_3," Jy"
err4 = np.sqrt(err4_1**2+err4_2**2+err4_3**2)
print "Calibration error:", err1," Jy"
print "Instrumental error:", err2," Jy"
print "Confusion noise:", err3," Jy"
print "Background error:", err4," Jy"
tot_err = np.sqrt(err1**2 + err2**2 + err3**2 + err4**2)
rel_err = tot_err/tflux*100
print 'Background-subtracted flux:', tflux,'±',tot_err,' Jy  (',rel_err,' %)'
print 'Average sky per pixel value:', sk_p_pix,' Jy'


"""
Instrumental background: 0.06909057823157644  Jy
Confusion background: 0.213286303697584  Jy
Large-scale background: 0.22837460087525382  Jy
Calibration error: 14.331282862647516  Jy
Instrumental error: 0.0432166700647272  Jy
Confusion noise: 0.2118646067962301  Jy
Background error: 0.3200304886573826  Jy
Background-subtracted flux: 204.73261232353593 ± 14.33648639991923  Jy  ( 7.002541625983598  %)
Average sky per pixel value: 0.0009287600937766091  Jy

"""

"""
Instrumental background: 0.053564005540125297  Jy
Confusion background: 0.17406191251517708  Jy
Large-scale background: 0.12574949347004608  Jy
Calibration error: 5.9209109730993585  Jy
Instrumental error: 0.03297025259795511  Jy
Confusion noise: 0.17309141713596368  Jy
Background error: 0.22131332356168007  Jy
Background-subtracted flux: 84.58444247284797 ± 5.92766514023834  Jy  ( 7.007985117524597  %)
Average sky per pixel value: 0.0025877744631146995  Jy
"""

"""
Instrumental background: 0.18744122516450384  Jy
Confusion background: 0.9123620603873758  Jy
Large-scale background: 2.5040469278591866  Jy
Calibration error: -5.097964543386519  Jy
Instrumental error: 2.3778496221471594  Jy
Confusion noise: 0.47169252170933246  Jy
Background error: 2.671664230221773  Jy
Background-subtracted flux: -72.8280649055217 ± 6.245293820570375  Jy  ( -8.575394428881589  %)
Average sky per pixel value: 0.007249913457169627  Jy
"""

"""

Instrumental background: 0.03700488097723203  Jy
Confusion background: 0.18078065644344227  Jy
Large-scale background: 0.03206584898395391  Jy
Calibration error: 1.9669285137685721  Jy
Instrumental error: 0.019096932794854044  Jy
Confusion noise: 0.15659085047523097  Jy
Background error: 0.18729448905753407  Jy
Background-subtracted flux: 28.098978768122457 ± 1.9821131126153961  Jy  ( 7.0540396822679226  %)
Average sky per pixel value: 0.0026859347968290763  Jy
"""