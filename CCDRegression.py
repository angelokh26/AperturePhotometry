"""
Created on Wed Apr 1 09:30:16 2020

This program is used to characterize a CCD camera for astrophotography. In 
particular, to examine the collected events as a function of exposure time
at a certain gain level.

@author: Angelo Hollett - angelokh26@gmail.com
"""

import pandas as pd

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np


#================data access===================
#datafile = ascii.read("data.txt")                         # Read data from text file
#x = np.array(datafile["x"])                              # x data
#y = np.array(datafilefile["y"])                              # y data
#y_err = np.array(data["y_err"])                       # y_err
#==============================================

#d = pd.DataFrame(data={'x':[round(x,2) for x in x],
#                       'y':[round(x,2) for x in y],
#                       'y_err':[round(x,2) for x in y_err]})

#fig, (ax1, ax2) = plt.subplots(2, 1, figsize=[10,20])
#fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)


d = pd.read_csv("Gain_139.txt")


x = d['x']
y = d['y']

y_errf = np.sqrt(y)
y_err = np.sqrt(y[:14])

# a simple linear fit

# set up the model
def model(x,m,b):
    return m*x+b

fig1, (ax1) = plt.subplots(1, 1, figsize=[10,3])

# LINEAR fit
init_guess = [1,1]
fit = curve_fit(model, x[:14],y[:14], sigma=y_err, p0=init_guess, absolute_sigma=True)

# unpack the results LINEAR
ans,cov = fit
fit_m,fit_b = ans
fit_sm,fit_sb = np.sqrt(np.diag(np.absolute(cov)))

# plot the LINEAR data and fit results
ax1.errorbar(x,y,y_errf,fmt='k.', label="data")
ax1.set_title("Regression for counts as a function of exposure time at 139 gain.")
ax1.set_ylabel("Counts (ADU)")
ax1.set_xlabel("Exposure (s)")
print("covariance:")
print(cov)
print ('Variance:', (np.std(y))**2)

# print the LINEAR fit results:
print("The fited values for the linear fit model are:")
print("Linear slope m: %.2f +/- %.2f"%(fit_m,fit_sm))
print("Y-axis intercept / vertical translation b: %.2f +/- %.2f"%(fit_b,fit_sb))
print

curve = (fit_m)*x[:14] + fit_b
ax1.plot(x[:14], curve, '-r', label='Model')
#ax1.plot(x, y,ls='none', marker='o', color='k')
#ax1.errorbar(x,y,yerr=y_err, xerr=0)