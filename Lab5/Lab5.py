
# # Lab 5 ASTR 400B 

# Import Modules 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy import constants as const # import astropy constants
import astropy.units as u


# # Part A :  Mass to Light Ratios 
# 
# Wolf et al. 2010 
# 
# $M(<R_{half}) = \frac {4}{G}\sigma^2 R_e$
# 
# Where $R_{half}$ = 3D half mass radius 
# and $R_e$ is the 2D half mass radius of stars (observed)
# 
# Determine which of the following two systems are galaxies:
# 
# The system 47 Tuc is observed with:  $\sigma = 17.3$ km/s, $R_e = 0.5$ pc, $L_v \sim 10^5 L_\odot$ 
# 
# The system Willman I is observed with: $\sigma = 4.3$ km/s, $R_e = 25$ pc, $L_v = 10^3 L_\odot$


# Gravitational Constant in the desired units
# kpc^3/Gyr^2/Msun
Grav = const.G.to(u.kpc**3/u.Gyr**2/u.Msun)

def WolfMass(sigma, re):
    """ Function that defines the Wolf mass estimator from Wolf+ 2010
    PARAMETERS
    ----------
        sigma: astropy quantity
            1D line of sight velocity dispersion in km/s
        re: astropy quantity
            Effective radius, 2D radius enclosing half the
            stellar mass in kpc
    OUTPUTS
    -------
        mWolf: Returns the dynamical mass within the 
            half light radius in Msun
    """
    sigmaKpcGyr = sigma.to(u.kpc/u.Gyr)
    mWolf = 4/Grav*sigmaKpcGyr**2*re
    return mWolf

# # Part B :  Stellar to Halo Mass Relation
# 
# Following the work of [Moster et al. 2013 (MNRAS, 428, 3121)](https://ui.adsabs.harvard.edu/abs/2013MNRAS.428.3121M/abstract)
# 
# 
# `Equation 2:`                  $ \frac{m}{M} = 2N \left [ \left ( \frac{M}{M_1} \right)^{-\beta} + \left (\frac{M}{M_1} \right)^{\gamma} \right]$ 
# 
# $m$ = stellar mass, $M$ = halo mass
# 
# `Equation 11:`        log $M_1(z) = M_{10} + M_{11} \frac{z}{z+1} $ 
# 
# `Equation 12:`        $N(z) = N_{10} + N_{11} \frac{z}{z+1} $
# 
# `Equation 13:`         $\beta(z) = \beta_{10} + \beta_{11} \frac{z}{z+1} $
# 
# `Equation 14:`         $\gamma(z) = \gamma_{10} + \gamma_{11} \frac{z}{z+1} $

# # Q1 
# 
# Modify the class below by adding a function called `StellarMass` that uses the `SHMratio` function and returns the stellar mass.

# In[ ]:
    
mass_tuc = WolfMass(17.3*(u.km/u.s), 0.0005*(u.kpc))
print('Mass Tuc ' + f'{mass_tuc:.2e}')
mass_will = WolfMass(4.3*(u.km/u.s), 0.025*(u.kpc))
print('Mass Will ' + f'{mass_will:.2e}')
# If mass to light ratio is much larger it is a galaxy
print('Mass light ratio Tuc ' + f'{(mass_tuc/1e5).value:.5}')
print('Tuc is not galaxy')
print('Mass light ratio Will ' + f'{(mass_will/1e3).value:.5}')
print('Will is galaxy')



class AbundanceMatching:
    """ Class to define the abundance matching relations from 
    Moster et al. 2013"""
    
    
    def __init__(self, mhalo, z):
        """ Initialize the class
        PARAMETERS
        ----------
            mhalo: float
                Halo mass in Msun
            z: float
                redshift
        """
        
        #initializing the parameters:
        self.mhalo = mhalo # Halo Mass in Msun
        self.z = z  # Redshift
        
        
    def logM1(self):
        """eq. 11 of Moster 2013
        OUTPUT: 
            M1: float 
                characteristic mass in log(Msun)
        """
        M10      = 11.59
        M11      = 1.195 
        return M10 + M11*(self.z/(1+self.z))  
    
    
    def N(self):
        """eq. 12 of Moster 2013
        OUTPUT: 
            Normalization for eq. 2
        """
        N10      = 0.0351
        N11      = -0.0247
    
        return N10 + N11*(self.z/(1+self.z))
    
    
    def Beta(self):
        """eq. 13 of Moster 2013
        OUTPUT:  power of the low mass slope"""
        beta10      = 1.376
        beta11      = -0.826
    
        return beta10 + beta11*(self.z/(1+self.z))
    
    def Gamma(self):
        """eq. 14 of Moster 2013
        OUTPUT: power of the high mass slope """
        gamma10      = 0.608
        gamma11      = 0.329
    
        return gamma10 + gamma11*(self.z/(1+self.z))
    
    
    def SHMratio(self):
        """ 
        eq. 2 of Moster + 2013
        OUTPUT: 
            SHMratio float
                Stellar mass to halo mass ratio
        """
        M1 = 10**self.logM1() # Converting characteristic mass 
        # to Msun from Log(Msun)
        A = (self.mhalo/M1)**(-self.Beta())  # Low mass end
        B = (self.mhalo/M1)**(self.Gamma())   # High mass end
        Norm = 2*self.N() # Normalization
    
        SHMratio = Norm*(A+B)**(-1)
    
        return SHMratio
    
    def stellarMass(self):
        '''
        Method to compute the stellar mass of a galaxy using eq. 2 from Moster 13'
        
        Output:
            M1: float, stellar mass in Msun
        '''

        return self.mhalo*self.SHMratio()
    
# Q1: add a function to the class that takes the SHM ratio and returns 
# The stellar mass 
# # Part C : Plot the Moster Relation
# 
# Reproduce the below figure from Moster + 2013 
# Plot this for z=0, 0.5, 1, 2
# 
# ![mos](./MosterFig.png)



mh = np.logspace(10,15,1000) # Logarithmically spaced array


# Define Instances of the Class for each redshift
MosterZ0 = AbundanceMatching(mh,0)
MosterZ1 = AbundanceMatching(mh,1)
MosterZ2 = AbundanceMatching(mh,2)
MosterZ5 = AbundanceMatching(mh,5)
MosterZ9_5 = AbundanceMatching(mh,9.5)

fig,ax = plt.subplots(figsize=(10,8))


#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# Plot z = 0
plt.plot(np.log10(mh), np.log10(MosterZ0.stellarMass()),
         linewidth = 5, label='z=0')
plt.plot(np.log10(mh), np.log10(MosterZ1.stellarMass()),
         linewidth = 5, label='z=1')
plt.plot(np.log10(mh), np.log10(MosterZ2.stellarMass()),
         linewidth = 5, label='z=2')
plt.plot(np.log10(mh), np.log10(MosterZ5.stellarMass()),
         linewidth = 5, label='z=5')
plt.plot(np.log10(mh), np.log10(MosterZ9_5.stellarMass()),
         linewidth = 5, label='z=9.5')

# Continue plotting for the other redshifts here




# Axes labels 
plt.xlabel('log (M$_h$/M$_\odot$)',fontsize=22) 
plt.ylabel('log (m$_\star$/M$_\odot$)', fontsize=22)

# Legend
plt.legend(loc='lower right',fontsize='x-large')

# save the file 
plt.savefig('AbundanceMatching_Lab5.png')


# # Part D
# 
# # Q1
# 
# In traditional models of the Magellanic Clouds (prior to 2010), the LMC is thought to have a halo mass of order $3 \times 10^{10}$ M$_\odot$.  According to LCDM theory, what should be the stellar mass of such a halo?  
# 
# How does this compare against the actual observed stellar mass of the LMC at the present day of $3 \times 10^9$ M$_\odot$ ? 
# 
# What is the $\Lambda$CDM expected halo mass? What is the origin of any discrepancy? 

haloLMC1 = 3e10

LMC1 = AbundanceMatching(haloLMC1, 0)

LMC1star = LMC1.stellarMass()

print(LMC1star)

print(3e9/LMC1star)

haloLMC2 = 1.647e11

LMC2 = AbundanceMatching(haloLMC2, 0)

LMC2star = LMC2.stellarMass()

print(3e9/LMC2star)

# # Q2
# 
# 
# What is the expected stellar mass of an L* galaxy at z=0? 
# 
# What is the expected stellar mass of an L* galaxy at z = 2?  

print(f'Log M!, Characterisitic halo mass at z=0: {MosterZ0.logM1()}')

Mstarz0 = AbundanceMatching(10**(MosterZ0.logM1()), 0)

print(f'stellar Mass of a L* galalxy at z = 0: {np.round(Mstarz0.stellarMass()/1e10, 2)} x 1e10 Msun')


print(f'Log M!, Characterisitic halo mass at z=1: {MosterZ1.logM1()}')

Mstarz1 = AbundanceMatching(10**(MosterZ1.logM1()), 1)

print(f'stellar Mass of a L* galalxy at z = 1: {np.round(Mstarz1.stellarMass()/1e10, 2)} x 1e10 Msun')
