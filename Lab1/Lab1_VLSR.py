
# # In Class Lab 1
# Must be uploaded to your Github repository under a "Labs/Lab1" folder by 5 PM on Jan 31st 2023

# ## Part A:  The Local Standard of Rest
# Proper motion of Sgr A* from Reid & Brunthaler 2004
# $\mu = 6.379$ mas/yr 
# 
# Peculiar motion of the sun, $v_\odot$ = 12.24 km/s  (Schonrich 2010)
# 
# 
# $v_{tan} = 4.74 \frac{\mu}{\rm mas/yr} \frac{R_o}{\rm kpc} = V_{LSR} + v_\odot$
# 
# 
# ### a)
# 
# Create a function called VLSR to compute the local standard of res (V$_{LSR}$).
# 
# The function should take as input: the solar radius (R$_o$), the proper motion (mu)
# and the peculiar motion of the sun in the $v_\odot$ direction.
# 
# Compute V$_{LSR}$ using three different values R$_o$: 
# 1. Water Maser Distance for the Sun :  R$_o$ = 8.34 kpc   (Reid 2014 ApJ 783) 
# 2. GRAVITY Collaboration Distance for the Sun:  R$_o$ = 8.178 kpc   (Abuter+2019 A&A 625)
# 3. Value for Distance to Sun listed in Sparke & Gallagher : R$_o$ = 7.9 kpc 
# 



# Import Modules 
import numpy as np # import numpy
import astropy.units as u # import astropy units
from astropy import constants as const # import astropy constants



def VLSR(Ro, mu = 6.379, vsun = 12.24*u.km/u.s):
    '''This function will compute the velocity at the local standard of rest. 
                VLSR = 4.74 *mu*Ro - vsun
        
        Inputs:
            Ro: 'astropy quantity'
                The distance from the sun to the galactic center in kpc. 
            mu: 'float'
                The porper motion of Sag A* in mas/yr. 
                Default is from Reid and Brunthaler 2004. 
            vsun: 'astropu quantity'
                The peculiar motion of the sun in the v direction in km/s.
                Default is from Schonrich + 2010
                
        Outputs: 
            VLSR: 'astropy quantity'
                The velocity of the local standard of rest in km/s.
    '''
    
    return 4.74*mu*(Ro/u.kpc)*u.km/u.s - vsun

# Define our distances to the Galactic center from the Sun. 

RoReid = 8.34*u.kpc # Distance from Reid et al. 2014 in kpc.
RoGravity = 8.178*u.kpc # Distance from the Gravity Collab Abuter+2019 in kpc. 
RoSG = 7.9*u.kpc # Distance from the textbook Sparke and Gallagher. 
    

# Compute VLSR using Ro from Reid 2014.

VLSR_Reid = VLSR(RoReid)
print(np.round(VLSR_Reid))

VLSR_Gravity = VLSR(RoGravity)
print(np.round(VLSR_Gravity))

VLSR_SG = VLSR(RoSG)
print(np.round(VLSR_SG))

# ### b)
# 
# compute the orbital period of the sun in Gyr using R$_o$ from the GRAVITY Collaboration (assume circular orbit)
# 
# Note that 1 km/s $\sim$ 1kpc/Gyr

def TorbSun(R, V):
    '''This Function will compute the orbital period of the Sun.
    
                T = 2*pi*R / V
                
        Inputs:
            R: 'astropy quantity'
                Distance in kpc (Distance to the galactic center.)
            V: 'astropy quantity'
                Velocity in km/s (Velocity of the Sun.)
                
        Output: 
            'astropy quantity'
            orbital period in Gyr.
    '''

    VkpcGyr = V.to(u.kpc/u.Gyr) # converting v from km/s to kpc/Gyr
    T = 2*np.pi*R / VkpcGyr # Orbital Period
    
    return T
    




# ### c)
# 
# Compute the number of rotations about the GC over the age of the universe (13.8 Gyr)


VsunPeculiar = 12.24 * u.km/u.s
VSun = VLSR_Gravity + VsunPeculiar


T_Grav = TorbSun(RoGravity, VSun)
print(T_Grav)


# Age of the univers / Orbital Period
Age = 13.8*u.Gyr

print(Age / T_Grav)

# ## Part B  Dark Matter Density Profiles
# 
# ### a)
# Try out Fitting Rotation Curves 
# [here](http://wittman.physics.ucdavis.edu/Animations/RotationCurve/GalacticRotation.html)
# 
# 
# ### b)
# 
# 
# In the Isothermal Sphere model, what is the mass enclosed within the solar radius (R$_o$) in units of M$_\odot$? 
# 
# Recall that for the Isothermal sphere :
# $\rho(r) = \frac{V_{LSR}^2}{4\pi G r^2}$
# 
# Where $G$ = 4.4985e-6 kpc$^3$/Gyr$^2$/M$_\odot$, r is in kpc and $V_{LSR}$ is in km/s
# 
# What about at 260 kpc (in units of  M$_\odot$) ? 


# density profile rho = VLSR^2 / (4*pi*G*R^2)
# Mass = Integrate rho dV
#      = rho  4*pi*r**2 dr
#      = VLSR**2 /G/ (4*pi*r**2)  *  (4*pi*r**2) dr 
#      = VLSR**2 / G * r

Grav = const.G.to(u.kpc**3 / u.Gyr**2 / u.Msun)

def MassIso(r, VLSR):
    '''This Function will compute the Dark Matter mass enclosed within a given distance assuming an isothermal
    Sphere Mode for hte Dark Matter.
    M = VLSR**2 / G * r
    
    Inputs: 
        r: 'astropy quantity
            Distance to the Galactic Center in kpc.
        VLSR: 'astropy quantity
            Velocity of the Local Standard of Rest (km/s)
            
    Outputs: 
        M: Mass enclosed within r in units of Msun.'''

    VLSRkpcgyr = VLSR.to(u.kpc/u.Gyr) # Converting km/s to kpc / gyr
    
    M = (VLSRkpcgyr**2 / Grav )* r # Mass for isothermal sphere
    
    return M

MIsoSolar = MassIso(RoGravity, VLSR_Gravity)
print(f'{MIsoSolar:.3e}')

# Compute the mass within 260 kpc

MIso260 = MassIso(260*u.kpc, VLSR_Gravity)

print(f'{MIso260:.3e}')

    
# ## c) 
# 
# The Leo I satellite is one of the fastest moving satellite galaxies we know. 
# 
# 
# It is moving with 3D velocity of magnitude: Vtot = 196 km/s at a distance of 260 kpc (Sohn 2013 ApJ 768)
# 
# If we assume that Leo I is moving at the escape speed:
# 
# $v_{esc}^2 = 2|\Phi| = 2 \int G \frac{\rho(r)}{r}dV $ 
# 
# and assuming the Milky Way is well modeled by a Hernquist Sphere with a scale radius of $a$= 30 kpc, what is the minimum mass of the Milky Way (in units of M$_\odot$) ?  
# 
# How does this compare to estimates of the mass assuming the Isothermal Sphere model at 260 kpc (from your answer above)

# Potential for a Hernquist profile 
# phi(r) = -G*M / (r+a)

# USing the potential for a Hernquist profile, the equation for the escape speed becomes:
# vesc^2 = 2 * G*M / (r+a)

# Rearrange the escape speed equation for M: 
# M = vesc**2 / 2 / G * (r+a)

def MassFromVesc(vesc, r, a):
    '''This function determines the total mass needed for a given excape speed assuming a Hernquist profile for 
    the dark matter Halo.
    
        M = vesc**2*(r+a) / 2 / G
        
    Inputs: 
        vesc: 'astropy quantity' 
            The escape speed in km/s (or the speed of the satellite)
            
        r: 'astropy quanityt ' 
            The ddistance from the Galactic Center (kpc)
        
        a: 'astropy quantity' 
            The Hernquist scale length (kpc)
            
    Outputs: 
        M: 'astropy quantity' 
        Total mass within r in Msun'''
        
    vescKpcGyr = vesc.to(u.kpc/u.Gyr) # Converting velocity units to kpc.gyr
    
    M = vescKpcGyr**2 / 2 / Grav*(r+a) # Required mass
    
    return M

VLeoI = 196*u.km/u.s # Spedd of LEoI from Sohn 2013 ApJ 768
      
a = 30 * u.kpc # Scale Radius for the Hernquist Halo

r = 260 * u.kpc # Galacitc Center distance of Leo I

# Compute the Massneeded to keep LeoI bound

MLeoI = MassFromVesc(VLeoI, r, a)

print(f'{MLeoI:.2e}')

print(MIso260/MLeoI)


