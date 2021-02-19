# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 16:33:02 2020

@author: Angelo Hollett angelokh26@gmail.com A00419149

Perform calibration for astrophotography. In particular, calculate the mean, 
median, and standard deviation of a master image.

"""
import os
import astropy.io.fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

master = []

def process_images():
	for dirs, subdirs, files in os.walk('./'):
		for file in files:
			if file.startswith("Bias_139"):
				
					
				hdulist = astropy.io.fits.open(file)
				header = hdulist[0].header
				data = hdulist[0].data
				
				print ("FITS File Name:", file, "\n")
				
				EXPOSURE = header['EXPTIME']
				EGAIN = header['EGAIN']
				print ("EXPOSURE:", EXPOSURE)
				print ("GAIN:", EGAIN)
				
				SCALED = data/16   #Scale to ADUs
				
				#print (SCALED)
				mean = np.mean(SCALED)
				sd = np.std(SCALED)
				median = np.mean(SCALED)
				print ("MEAN:", mean, "\nST DEV:", sd, "\nMEDIAN", median, "\n")
				
				#print (data)
				master.append(SCALED)
				#y.append(SCALED[1])
				
				#plt.figure()
				#plt.imshow(data, norm = LogNorm())
				hdulist.close()


def process_image(imagefile):				
	hdulist = astropy.io.fits.open(imagefile)
	header = hdulist[0].header
	data = hdulist[0].data
	
	print ("FITS File Name:", imagefile, "\n")
	
	EXPOSURE = header['EXPTIME']
	EGAIN = header['EGAIN']
	print ("EXPOSURE:", EXPOSURE)
	print ("GAIN:", EGAIN)
	
	SCALED = data/16   #Scale to ADUs
	
	mean = np.mean(SCALED)
	sd = np.std(SCALED)
	median = np.mean(SCALED)
	print ("MEAN:", mean, "\nST DEV:", sd, "\nMEDIAN", median, "\n")
	
	plt.figure()
	plt.imshow(data, norm = LogNorm())
	plt.gca().invert_yaxis()
	hdulist.close()

	return;

#path = './bias_at_gain_0/'
#imagename = path + 'Bias_0-0001.fit'
#process_image(imagename)

process_images()

mean = np.mean(master, axis=0)
median = np.median(master, axis=0)
stnddev = np.std(master, axis=0)

# plot the master images
# =============================================================================
# plt.figure()
# plt.title("Master Bias Image (Median)")
# plt.imshow(median, cmap = 'gist_yarg', norm = LogNorm())
# plt.gca().invert_yaxis()
# 
# plt.figure()
# plt.title("Master Bias Image (Mean)")
# plt.imshow(mean, cmap = 'gist_yarg', norm = LogNorm())
# plt.gca().invert_yaxis()
# =============================================================================

mastermean = np.mean(mean)
mastermedian = np.median(median)
masterstd = np.std(stnddev) 

print("The mean of the master image is:", mastermean)
print("The median of the master image is", mastermedian)
print("The standard deviation of the master image is", masterstd)
