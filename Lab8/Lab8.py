

# # Lab 8 : Star Formation 



import numpy as np
from astropy import units as u
from astropy import constants as const

import matplotlib
import matplotlib.pyplot as plt


# # Part A
# 
# Create a function that returns the SFR for a given luminosity (NUV, FUV, TIR, Halpha)
# 
# $Log( {\rm SFR} (M_\odot/year)) = Log(Lx (erg/s)) - Log(Cx)$ 
# 
# Including corrections for dust absorption 
# 
# Kennicutt & Evans 2012 ARA&A Equation 12 and Table 1, 2


def StarFormationRate(L, Type, TIR=0):
    '''
    Function that computes the star formation rate of a galaxy following
    Kennicutt & Evans 2012 Eq 12 (ARA&A 50).

    Parameters
    ----------
    L : 'float'
        Luminosity of the galaxy in a particular waveband (erg/s)
    Type : 'string'
        The wavelength: 'FUV', 'NUV', 'TIR', 'Halpha'
    TIR : 'float'
        Total Infrared Luminosity is erg/s (default=0)

    Returns
    -------
        sfr: 'float'
            Log of the star formation rate (Msun/yr)

    '''

    if (Type == 'FUV'):
        logCx = 43.35 # Calibration from the Table 1 (K&E 2012)
        TIRc = 0.46 # Correction for dust absorption from Table 2 (K&E 2012)
        
    elif (Type == 'NUV'):
        logCx = 43.17 # Calibration from the Table 1 (K&E 2012)
        TIRc = 0.27 # Correction for dust absorption from Table 2 (K&E 2012)
        
    elif (Type == 'Halpha'):
        logCx = 41.27 # Calibration from the Table 1 (K&E 2012)
        TIRc = 0.0024# Correction for dust absorption from Table 2 (K&E 2012)
        
    elif (Type == 'TIR'):
        logCx = 43.41 # Calibration from the Table 1 (K&E 2012)
        TIRc = 0 # Correction for dust absorption from Table 2 (K&E 2012)
        
    else:
        print('Missing wavelength: FUV, NUV, TIR, Halpha')
        
    # Correct the luminosity for dust using TIR
    Lnew = L + TIRc*TIR
        
    # Star formation rate
    SFR = np.log10(Lnew) - logCx
    
    return SFR
    
    
# Let's try to reproduce SFRs derived for galaxies from UV luminosities measured with Galex. 
# 
# Using Table 1 from Lee et al. 2009
# https://ui.adsabs.harvard.edu/abs/2009ApJ...706..599L/abstract
# 
# We will use galaxy properties from NED:
# https://ned.ipac.caltech.edu/




#  WLM Dwarf Irregular Galaxy

# From NED, WLM NUV luminosity is 1.71E+07 Lsun. 
# From NED: WLM NIR luminosity is 	2.48E+06 LSun. 
# From NEDL WLM FIR luminosity is 7.84E+05 Lsun.

LsunErgS = const.L_sun.to(u.erg/u.s).value

NUV_WLM = 1.71e7*LsunErgS

TIR_WLM = 2.48e6*LsunErgS + 7.84e5*LsunErgS

print(StarFormationRate(NUV_WLM, 'NUV', TIR_WLM))


#  NGC24 Sc galaxy

NUV_NGC= 2.96E+08 *LsunErgS

TIR_NGC = 8.34E+08*LsunErgS + 3.09E+08*LsunErgS

print(StarFormationRate(NUV_NGC, 'NUV', TIR_NGC))



# # Part B Star formation main sequence
# 
# 1) Write a function that returns the average SFR of a galaxy at a given redshift. 
# 
# 2) What is the average SFR of a MW mass galaxy today? at z=1?
# 
# 3) Plot the SFR main sequence for a few different redshifts from 1e9 to 1e12 Msun.
# 
# 
# From Whitaker 2012:
# 
# log(SFR) = $\alpha(z)({\rm log}M_\ast - 10.5) + \beta(z)$
# 
# $\alpha(z) = 0.7 - 0.13z$
# 
# $\beta(z) = 0.38 + 1.14z - 0.19z^2$


# Step 1

# Step 2

# MW at z=0

def SFRMainSequence(Mstar, z):
    '''
    This is a function that calculates the avergae SFR of a galaxy
    as a function of stellar mass. 

    Parameters
    ----------
    Mstar : 'float'
        Stellar mass of the galaxy in solar masses. 
    z : 'float'
        Redshift

    Returns
    -------
        logSFR: 'flaot'
            log(SFR(Msun/yr))

    '''

    alpha = 0.7 - 0.13*z
    
    beta = 0.38 + 1.14*z - 0.19*z**2
    
    logSFR = alpha*(np.log10(Mstar) - 10.5)+beta
    
    return logSFR


# MW at z = 1

MW_disk = 8e10

print(10**SFRMainSequence(MW_disk, 0))
print(10**SFRMainSequence(MW_disk, 1))

# Step 3


# create an array of stellar masses


Mass = np.linspace(1e9, 1e12)


fig = plt.figure(figsize=(8,8), dpi=500)
ax = plt.subplot(111)

# add log log plots

plt.loglog(Mass, 10**SFRMainSequence(Mass, 0), color = 'b', linewidth = 3, label = 'z=0')

plt.loglog(Mass, 10**SFRMainSequence(Mass, 1), color = 'r', linewidth = 3, label = 'z=1')


plt.loglog(Mass, 10**SFRMainSequence(Mass, 2), color = 'g', linewidth = 3, label = 'z=2')


plt.loglog(Mass, 10**SFRMainSequence(Mass, 3), color = 'm', linewidth = 3, label = 'z=3')

plt.loglog(Mass, 10**SFRMainSequence(Mass, 4), color = 'y', linewidth = 3, label = 'z=4')

plt.loglog(Mass, 10**SFRMainSequence(Mass, 5), color = 'k', linewidth = 3, label = 'z=5')

plt.loglog(Mass, 10**SFRMainSequence(Mass, 6), color = 'cyan', linewidth = 3, label = 'z=6')

#plt.loglog(Mass, 10**SFRMainSequence(Mass, 12), color = 'hotpink', linewidth = 3, label = 'z=12')

# Add axis labels
plt.xlabel('Log (Mstar (M$_\odot$))', fontsize=12)
plt.ylabel('Log(SFR (M$_\odot$/year))', fontsize=12)


#adjust tick label font size
label_size = 12
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

plt.show()

# # Part C  Starbursts

#Use your `StarFormationRate` code to determine the typical star formation rates for the following systems with the listed Total Infrared Luminosities (TIR): 

# Normal Galaxies: $10^{10}$ L$_\odot$
#
# LIRG: $10^{11}$ L$_\odot$
# 
# ULIRG: $10^{12} $ L$_\odot$
# 
# HLIRG: $10^{13} $ L$_\odot$



# normal galaxies 

TIR_Normal = 1e10*LsunErgS

print(10**StarFormationRate(TIR_Normal, 'TIR'))


# LIRGs  

TIR_LIRG= 1e11*LsunErgS

print(10**StarFormationRate(TIR_LIRG, 'TIR'))


# ULIRGs

TIR_ULIRGS= 1e12*LsunErgS

print(10**StarFormationRate(TIR_ULIRGS, 'TIR'))



# HLIRGs

TIR_HLIRGS= 1e13*LsunErgS

print(10**StarFormationRate(TIR_HLIRGS, 'TIR'))




