# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 12:35:53 2020

A program to perform photometry

The main function of this code is to extract data from fits file images
of a target star and to use a reference and check star to compute magnitudes
for a lightcruve. The lightcurve data are stored to a file named 'mags.csv'
and will contain the HJD along with the relevant magnitudes.

If the images are taken over many nights, then data taken on the same night 
should be averaged, which is not incorporated into this program. 

@author: Angelo Hollett angelokh26@gmail.com
"""
import os
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from photutils import DAOStarFinder
from astropy.stats import mad_std
from photutils import datasets
from photutils import aperture_photometry, CircularAperture, SkyCircularAperture, SkyCircularAnnulus, CircularAnnulus
from astropy.wcs import WCS
from astropy import wcs
from astropy.table import Table, hstack


# Empty arrays to be filled after photometry is performed
target_mags = []
reference_mags = []

MJD_dates = []

for subdir, dirs, files in os.walk('./'):
    for file in files:
        if file.endswith(".fit"):

			
			#loop through and display the images for visible checking artifacts
# =============================================================================
# 			print (file)
# 			image_file = file
# 			fits.info(image_file)
# 			image_data = fits.getdata(image_file, ext=0)
# 			plt.figure()
# 			plt.imshow(image_data, cmap='gist_gray_r', norm = LogNorm())
# 			plt.colorbar()
# 			plt.title(image_file)
#           break
# =============================================================================


            hdulist = fits.open(file)
            data = hdulist[0].data
            header = hdulist[0].header

            w = wcs.WCS(header)
            image_file = fits.open(file)
			
            print("The MJD date of this observation is ", hdulist[0].header['HJD-MID'])
            MJD_dates.append(hdulist[0].header['HJD-MID'])
            
            print('')
            print('')
            print('')
			
            r_source = 7.0
            r_in = 9.0
            r_out = 12.0
			
			# Give the ra,dec of a target and reference star
            world = np.array([[308.22621,46.60128]])
            targetpix = w.wcs_world2pix(world,1)
            world = np.array([[308.15274048,46.58061218]])
            refpix = w.wcs_world2pix(world,1)
			
			
            positions = [(targetpix[0,0], targetpix[0,1]), (refpix[0,0], refpix[0,1])]
			
            aperture = CircularAperture(positions, r=r_source)
            annulus = CircularAnnulus(positions, r_in, r_out)

            # define the apertures

            apertures = CircularAperture(positions, r=7)
            annulus_apertures = CircularAnnulus(positions, r_in=10, r_out=13)

            # call aperture functions for aperture and annulus
            rawflux_table = aperture_photometry(data, apertures)
            bkgflux_table = aperture_photometry(data, annulus_apertures)

            # hstack just combines two arrays
            phot_table = hstack([rawflux_table, bkgflux_table],table_names=['raw', 'bkg'])

            print("================================================================================")
            print(phot_table)
            print("================================================================================")


            #store the flux values from the phot_table
# =============================================================================
#             target_flux = (phot_table[0][4])
#             target_annulus_flux = (phot_table[0][5])
# 
#             reference_flux = (phot_table[1][4])
#             reference_annulus_flux = (phot_table[1][5])
# 
#             # Collect the x,y, coorinates of target and reference, in pixels		
#             x_target = (phot_table[0][1]).value
#             x_ref = (phot_table[1][1]).value
# 
# 
#             y_target = (phot_table[0][2]).value
#             y_ref = (phot_table[1][2]).value
# =============================================================================

            # Create the apertures and annuli
            #apperture_target = CircularAperture([x_target,y_target], r=r_source)
            #apperture_ref = CircularAperture([x_ref,y_ref], r=5.)

            #annulus_target = CircularAperture([x_target,y_target], r=r_out)
            #annulus_ref = CircularAperture([x_ref,y_ref], r=r_out)

			# used for putting apertures on all stars
            #bkg_sigma = mad_std(image)  
            #daofind = DAOStarFinder(fwhm=4., threshold=10.*bkg_sigma)  
            #sources = daofind(image)
            #for col in sources.colnames:  
            #	sources[col].info.format = '%.8g'  # for consistent table output

            #print(sources)         


            #positions = np.transpose((sources['xcentroid'], sources['ycentroid']))  
            #apertures = CircularAperture(positions, r=r_source)
            #annuli =  SkyCircularAperture(positions, r=r_in)
            #annulus_target = SkyCircularAnnulus()


            #phot_table = aperture_photometry(image, apertures)  
            #for col in phot_table.colnames:  
            #	phot_table[col].info.format = '%.8g'  # for consistent table output

            #PLOTTING the appertures
            # =============================================================================
            # 			plt.imshow(image, cmap='gray_r', origin='lower', norm = LogNorm())
            # 			plt.title(image_fit)
            # 			apertures.plot(color='blue', lw=1.5, alpha=0.5)
            # =============================================================================

            #annuli.plot(color='blue', lw=1.5, alpha=0.5)

            # PLOTTING annulus
            # =============================================================================
            # 			apperture_target.plot(color='red', lw=1.5, alpha=1.0)
            # 			
            # 			
            # 			apperture_ref.plot(color='orange', lw=1.5, alpha=1.0)
            # 			
            # 			
            # 			annulus_target.plot(color='red', lw=1.5, alpha=1.0)
            # 			
            # 			
            # 			annulus_ref.plot(color='orange', lw=1.5, alpha=1.0)
            # =============================================================================



            # PLOTTING the stars
            # =============================================================================
            # 			
            #           image_data = fits.getdata(image_fit, ext=0)
            #	        plt.figure()
            # 		
            # 			plt.imshow(image_data, cmap='gist_gray_r', norm = LogNorm())
            # 			plt.colorbar()
            # 			plt.title(image_fit)
            # 			plt.gca().invert_yaxis()
            # =============================================================================

	
			
			# Litereature magnitude of the reference star
			
            lit_ref_mag = 9.558
            
			
            # Background treatment

            bkg_mean = phot_table['aperture_sum_bkg'] / annulus_apertures.area

            bkg_sum = bkg_mean * apertures.area
            print (bkg_sum)

           # Background subtraction

            final_sum = phot_table['aperture_sum_raw'] - bkg_sum
            phot_table['residual_aperture_sum'] = final_sum

           # Calculating magnitudes for target and reference star

            target_mag = lit_ref_mag-(2.5*np.log10(final_sum[0]/final_sum[1]));
            ref_mag = lit_ref_mag-(2.5*np.log10(final_sum[1]/final_sum[0]));
			
			
            target_mags.append(target_mag)
            reference_mags.append(ref_mag)
            
            
            image_file.close()
            


# =============================================================================
# out_name = "mags.csv"
# 
# np.savetxt(out_name, np.column_stack((MJD_dates,target_mags,reference_mags)),
#            delimiter=",", fmt='%s',
#            header='MJD_dates, target_mags, reference_mags')
# 
# plt.plot(target_mags,MJD_dates,ls='none', marker='o', ms=10)
# =============================================================================
