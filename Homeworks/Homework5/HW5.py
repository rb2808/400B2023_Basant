# Author: Ritvik Basant
# Course: ASTRONOMY 400 B, Homework 4
# This python file contains the class which is used to calculate the mass profiles of different 
# galaxies and plotes them. In addition to this, it also computes the velocity curves
# for each galaxy and plots them.

# Importing modules
import numpy as np
import astropy.units as u
import astropy.constants as uc
import matplotlib.pyplot as plt

# Importing different functions and classes from previous homeworks. 
from ReadFile import Read
from CenterOfMass import CenterOfMass
from GalaxyMass import ComponentMass

class MassProfile: 
# This class computes the mass profiles, Hernquist Theoretical Mass profile, velocity curves, 
# and hernquist velocity curve for different galaxies. 

    def __init__(self, galaxy, snap):
        ''' Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            
            PARAMETERS
            ----------
            galaxy : `str`
                It is the galaxy name: MW, M31, or M33
            snap : `int: 1, 2, 3, ...`
                Time snapshot number.
        '''
        
        # Writing the filename using the snap number and the galaxy name.
        ilbl = '000' + str(snap)
        ilbl = ilbl[-3:]
        self.filename = "%s_"%(galaxy) + ilbl + '.txt'
     
        # read data in the given file using Read
        self.time, self.total, self.data = Read(self.filename)                                                                                             

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m']
        self.x = self.data['x']
        self.y = self.data['y']
        self.z = self.data['z']
        
        # Storing the name of the galaxy as a global property. 
        self.gname = galaxy
        
        
    def MassEnclosed(self, radius, ptype):
        ''' This function calculates the mass enclosed within a specific radius of a 
        particular particle type. 
        
            PARAMETERS
            ----
            
            radius: 'array'
                It is the array containing different radii for which we 
                need to calculate mass enclosed in kpc. 
            
            ptype: 'integer'
                Particle type
                
            Returns
            -------
                enclosedmass: 'array'
                    has the mass enclosed for the given particle time corresponding
                    to the radii in the radius array (M_sun). 
        '''
        
        # Calling the Center of Mass class from previous homework. 
        COM = CenterOfMass(self.filename, ptype)
        # Calculating the center of mass position for the given galaxy. 
        x_com, y_com, z_com = COM.COM_P(0.1)[0], COM.COM_P(0.1)[1], COM.COM_P(0.1)[2]

        # Calculating the indices for the specified particle type. 
        index = np.where(self.data['type'] == ptype)
        
        # Converting the x, y, z positions of the particles to the COM frame. 
        x_new= self.data['x'][index] - x_com/u.kpc
        y_new = self.data['y'][index] - y_com/u.kpc
        z_new = self.data['z'][index] - z_com/u.kpc
        
        # Calculating the magnitude of the particle's position vector. 
        r = np.sqrt(x_new**2 + y_new**2 + z_new**2)
        
        # Storing the masses for the specific particle type. 
        m_new = self.data['m'][index]
        
        # This array will contain the total mass enlcosed as a function of radius. 
        sums = np.zeros(len(radius))
        
        # This loop iterates over all the values in the radius array and computes the 
        # total mass enclosed within that radius. 
        for i in range(len(radius)):
            # This index defines the particles inside of the radius in the radii array. 
            index2 = np.where(r <= radius[i]/u.kpc)
            # Calculating the masses for all the enclosed particles. 
            masses = m_new[index2]
            # Summing up the masses and then storing it. 
            sums[i] = np.sum(masses) 
            
        return sums * 1e+10 * u.Msun
    

    def MassEnclosedTotal(self, rad):
        ''' This function will calculate the total mass enclosed (including all the particle types) 
        within a given radius and then return it.
        
        PARAMETERS
        ------
        rad: 'an array'
            It is an array that has different radii within which we need to calculate the total mass 
            enclosed in kpc.  
            
        Returns
        -------
            TME: 'array'
                This is the total mass enclosed at a given raidus from the rad array (M_sun). 
            
        
        '''
        # This array will store the total mass enclosed within a radius for different radii. 
        TME = np.zeros(len(rad))
        
        # These variables are for holding the arrays for masses of different components of the 
        # galaxy as a function of radius. 
        halo = 0
        bulge = 0
        disk = 0
        
        # This if statement computes the arrays for different components of the Milky Way 
        # galaxy using the MassEnclosed function from above. 
        if self.gname == 'MW':
            halo = (self.MassEnclosed(rad, 1)) / u.Msun
            bulge = (self.MassEnclosed(rad, 3)) / u.Msun
            disk = (self.MassEnclosed(rad, 2)) / u.Msun    
            
        # This if statement computes the arrays for different components of the M331 
        # galaxy using the MassEnclosed function from above. 
        elif self.gname == 'M31':
            halo = (self.MassEnclosed(rad, 1)) / u.Msun
            bulge = (self.MassEnclosed(rad, 3)) / u.Msun
            disk = (self.MassEnclosed(rad, 2)) / u.Msun  
            
        # This if statement computes the arrays for different components of the M33 
        # galaxy using the MassEnclosed function from above. 
        elif self.gname == 'M33':
            halo = (self.MassEnclosed(rad, 1)) / u.Msun
            disk = (self.MassEnclosed(rad, 2)) / u.Msun  
            bulge = np.zeros(len(rad))
            
        # Computing the total mass of all the components as a function of radius and 
        # storing it in the array. 
        for i in range(len(rad)):
            TME[i] = halo[i] + disk[i] + bulge[i] 
    
        return TME * u.Msun
    
    
    def HernquistMass(self, radius, Mhalo, scale_factor):
        
        '''This function computes the mass enclosed within a given radius
        using the theoretical Hernquist profile. 
        
        PARAMETERS
        -----
        radius: 'astropy quantity'
            This is the radius within which we need to compute the mass in kpc. 
             
        Mhalo: 'int'
            Total Halo Mass in M_sun.  
        
        scale_factor: 'int'
            It is the scale factor for the theoretical profile. 
            
        OUTPUTS
        -----
            m: 'astropy quantity'
                It is the halo mass enclosed within a radius. (M_sun)
        '''
        
        # These are the numerator and denominator from the formula for the Hernquist mass. 
        numerator = Mhalo * ((radius/u.kpc)**2)
        denominator = (scale_factor + (radius/u.kpc))**2
        
        # Total enclosed Hernquist Halo mass. 
        m = numerator / denominator
        
        return m
    
    
    def CircularVelocity(self, radii, ptype):
        '''This function takes in an array of radii and a particle type to calculate
        the circular velocity at each of the given radius. 
        

        Parameters
        ----------
        radii : 'array'
            This array contains various radii in kpc. 
        ptype : 'int'
            This is the particle type.

        Returns
        -------
        vel: 'array'
            This is array has the article specific velocity as a function of raidus for 
            a galaxy. (km /s)

        '''
        
        # Calculating the mass enclosed using the function above. 
        mass_encl = self.MassEnclosed(radii, ptype)
        
        # This array will store velocities for given particle type as a function of radius. 
        vel = np.zeros(len(radii))
        
        # This loop iterates to calculate the circular veocity for the given radius in the radii
        # array. The formula is simple and can be derived by equating gravitaional force with 
        # centrifugal force. 
        for i in range(len(radii)):
            # Adjusting the units of G so that the velocity is in km/s
            G = (uc.G).to(u.kpc*(u.km**2)/u.s**2/u.Msun)
            v = np.sqrt(G * mass_encl[i] / (radii[i]))
            vel[i] = v / (u.km / u.s)
            
        return vel * (u.km / u.s)
    
    def CircularVelocityTotal(self, radii):
        '''This function takes in an array of radii and then
        calculates the total circular velocity 

        Parameters
        ----------
        radii : 'array'
            This contains various radii for which we need to calculate the 
            total circular velocity in kpc. 

        Returns
        -------
        total_vel: 'array'
            Contains total circular velocity (in km/s) for every radius in the radii array. 

        '''
        # Calculating the total mass enclosed for the given radii in the radii array. 
        # Total mass means for all the components. 
        tot_mass_encl = self.MassEnclosedTotal(radii)
        
        # This array will store the total velocity as a function of radius. 
        tot_vel = np.zeros(len(radii))
        
        # This loop iterates over to calculate the total velocity using the same formula from 
        # CircularVelocity function. 
        for i in range(len(radii)):
            # Adjusting the untis. 
            G = (uc.G).to(u.kpc*(u.km**2)/u.s**2/u.Msun)
            v = np.sqrt(G * tot_mass_encl[i] / (radii[i]))
            tot_vel[i] = v / (u.km / u.s)
            
        return tot_vel * (u.km / u.s)
    
    def HerquinstVCirc(self, radius, scale, Mhalo):
        '''This function calculates the theoretical circular velocity using the
        Hernquist mass profile. 

        Parameters
        ----------
        radius : 'astropy quantity'
            The distance at which we need to calculate the velocity in kpc. 
            
        scale : 'int'
            Scale Factor
            
        Mhalo : 'astropy quantity'
            Total mass of the halo in M_sun.  

        Returns
        -------
            v: 'astropy quantity'
                This is the theoretical hernquist velocity (assuming spherical symmetry).
                (km/s)
        '''
        
        # Calculating the Theoretical Hernquist Mass for the given radius. 
        hern_mass = self.HernquistMass(radius, Mhalo, scale)
        
        # Adjusting the units of G. 
        G = (uc.G).to(u.kpc*u.km**2/u.s**2/u.Msun)
        
        # Calculating the velocity using teh same formula. 
        v = np.sqrt(G * hern_mass / (radius))
        
        return v
        
# Answering Question.
# I will comment for just one galaxy and all other would just follow similarly. 
if __name__ == '__main__' : 
    
    '''This part computes and plots the mass profiles and velocity curves for 
    Milky Way galaxy. '''
    
    # Calling the Mass profile class for Milky Way. 
    MW = MassProfile('MW', 0)
    
    # Defining the array for distance from the center of each galaxy. 
    r = np.arange(0.1, 30.5, 0.5) * u.kpc
    
    # These variables store the mass enclosed as a function of distance for every particle as well
    # as the total mass as a function of radius for Milky Way galaxy. 
    halo_mw = MW.MassEnclosed(r, 1)
    disk_mw = MW.MassEnclosed(r, 2)
    bulge_mw = MW.MassEnclosed(r, 3)
    total_mass_mw = MW.MassEnclosedTotal(r)
    
    # This variable uses the componentmass function from previous homeworks to calculate the total
    # mass of halo of the milky way. 
    halo_mass_mw = ComponentMass('MW_000.txt', 'Halo')
    
    # The scale factor for mass profile of milky way. 
    scale_mw_m = 62.5
    
    # This variables stores the Hernquinst theoretical mass for milky way. 
    hern_m_mw = MW.HernquistMass(r, halo_mass_mw, scale_mw_m)
    
    # Here we plot all these things. 
    plt.figure(figsize = (15, 15))
    plt.plot(r, halo_mw, label = 'Halo Mass', color = 'b')
    plt.plot(r, disk_mw, label = 'Disk Mass',  color = 'r')
    plt.plot(r, bulge_mw, label = 'Bulge Mass',  color = 'g')
    plt.plot(r, total_mass_mw, label = 'Total Mass', color = 'm')
    plt.scatter(r, hern_m_mw, color = 'k', label = 'Hernquist Mass Profile, a = ' + str(scale_mw_m), marker = '*', s = 60)
    plt.xlabel('Distance (kpc)', fontsize = 35)
    plt.ylabel('Mass ($10^{11}$ $M_{\odot}$)', fontsize = 35)
    plt.xticks(size = 20)
    plt.yticks(size = 20)
    plt.legend(fontsize = 20)
    plt.title('Milky Way Mass Distribution', fontsize = 35)
    plt.semilogy()
    plt.show()
    
    # These variables store the velocity as a fucntion of distance for all different particles
    # as well as the total velocity for milky way. 
    vel_halo_mw = MW.CircularVelocity(r, 1)
    vel_disk_mw = MW.CircularVelocity(r, 2)
    vel_bulge_mw = MW.CircularVelocity(r, 3)
    total_vel_mw = MW.CircularVelocityTotal(r)
    
    # Scale factor for the hernquist velocity 
    scale_mw_v = 62.5
    
    # Hernquist velocity for the given distance and the halo mass for milky way. 
    hern_vel_mw = MW.HerquinstVCirc(r, scale_mw_v, halo_mass_mw)
    
    # Here we are just plotting the velocity curves as a function of distance for milky way. 
    plt.figure(figsize = (15, 15))
    plt.plot(r, vel_halo_mw, label = 'Halo Velocity', color = 'b')
    plt.plot(r, vel_disk_mw, label = 'Disk Velocity',  color = 'r')
    plt.plot(r, vel_bulge_mw, label = 'Bulge Velocity',  color = 'g')
    plt.plot(r, total_vel_mw, label = 'Total Velocity', color = 'm')
    plt.scatter(r, hern_vel_mw, color = 'k', label = 'Hernquist Velocity Profile, a = ' + str(scale_mw_v), marker = '*', s = 60)
    plt.xlabel('Distance (kpc)', fontsize = 35)
    plt.ylabel('Velocity (km/s)', fontsize = 35)
    plt.xticks(size = 20)
    plt.yticks(size = 20)
    plt.legend(fontsize = 20)
    plt.title('Milky Way Velocity Distribution', fontsize = 35)
    plt.semilogy()
    plt.show()
    
    
    '''Repeating the same thing as above for M31 galaxy.'''
    M31 = MassProfile('M31', 0)
    
    halo_m31 = M31.MassEnclosed(r, 1)
    disk_m31 = M31.MassEnclosed(r, 2)
    bulge_m31 = M31.MassEnclosed(r, 3)
    total_mass_m31 = M31.MassEnclosedTotal(r)
    
    halo_mass_m31 = ComponentMass('M31_000.txt', 'Halo')
    
    scale_m31 = 61.5
    
    hern_m_m31 = M31.HernquistMass(r, halo_mass_m31, scale_m31)
    
    plt.figure(figsize = (15, 15))
    plt.plot(r, halo_m31, label = 'Halo Mass', color = 'b')
    plt.plot(r, disk_m31, label = 'Disk Mass',  color = 'r')
    plt.plot(r, bulge_m31, label = 'Bulge Mass',  color = 'g')
    plt.plot(r, total_mass_m31, label = 'Total Mass', color = 'm')
    plt.scatter(r, hern_m_m31, color = 'k', label = 'Hernquist Mass Profile, a = ' + str(scale_m31), marker = '*', s = 60)
    plt.xlabel('Distance (kpc)', fontsize = 35)
    plt.ylabel('Mass ($10^{11}$ $M_{\odot}$)', fontsize = 35)
    plt.xticks(size = 20)
    plt.yticks(size = 20)
    plt.legend(fontsize = 20)
    plt.title('M31 Mass Distribution', fontsize = 35)
    plt.semilogy()
    plt.show()
    
    vel_halo_m31 = M31.CircularVelocity(r, 1)
    vel_disk_m31 = M31.CircularVelocity(r, 2)
    vel_bulge_m31 = M31.CircularVelocity(r, 3)
    total_vel_m31 = M31.CircularVelocityTotal(r)
    
    scale_m31_v = 61.5
    
    hern_vel_m31 = M31.HerquinstVCirc(r, scale_m31_v, halo_mass_m31)
    
    plt.figure(figsize = (15, 15))
    plt.plot(r, vel_halo_m31, label = 'Halo Velocity', color = 'b')
    plt.plot(r, vel_disk_m31, label = 'Disk Velocity',  color = 'r')
    plt.plot(r, vel_bulge_m31, label = 'Bulge Velocity',  color = 'g')
    plt.plot(r, total_vel_m31, label = 'Total Velocity', color = 'm')
    plt.scatter(r, hern_vel_m31, color = 'k', label = 'Hernquist Velocity Profile, a = ' + str(scale_m31_v), marker = '*', s = 60)
    plt.xlabel('Distance (kpc)', fontsize = 35)
    plt.ylabel('Velocity (km/s)', fontsize = 35)
    plt.xticks(size = 20)
    plt.yticks(size = 20)
    plt.legend(fontsize = 20)
    plt.title('M31 Velocity Distribution', fontsize = 35)
    plt.semilogy()
    plt.show()
    
    '''Repeating the same thing for M33 galaxy. Only difference is that this galaxy does not have 
    bulge mass. '''
    M33 = MassProfile('M33', 0)

    halo_m33 = M33.MassEnclosed(r, 1)
    disk_m33 = M33.MassEnclosed(r, 2)
    total_mass_m33 = M33.MassEnclosedTotal(r)
    
    halo_mass_m33 = ComponentMass('M33_000.txt', 'Halo')
    
    scale_m_m33 = 27
    
    hern_m_m33 = M33.HernquistMass(r, halo_mass_m33, scale_m_m33)
    
    plt.figure(figsize = (15, 15))
    plt.plot(r, halo_m33, label = 'Halo Mass', color = 'b')
    plt.plot(r, disk_m33, label = 'Disk Mass',  color = 'r')
    plt.plot(r, total_mass_m33, label = 'Total Mass', color = 'm')
    plt.scatter(r, hern_m_m33, color = 'k', label = 'Hernquist Mass Profile, a = ' + str(scale_m_m33), marker = '*', s = 60)
    plt.xlabel('Distance (kpc)', fontsize = 35)
    plt.ylabel('Mass ($10^{11}$ $M_{\odot}$)', fontsize = 35)
    plt.xticks(size = 20)
    plt.yticks(size = 20)
    plt.legend(fontsize = 20)
    plt.title('M33 Mass Distribution', fontsize = 35)
    plt.semilogy()
    plt.show()
    
    vel_halo_m33 = M33.CircularVelocity(r, 1)
    vel_disk_m33 = M33.CircularVelocity(r, 2)
    total_vel_m33 = M33.CircularVelocityTotal(r)
    
    scale_m33_v = 27
    
    hern_vel_m33 = M33.HerquinstVCirc(r, scale_m33_v, halo_mass_m33)
    
    plt.figure(figsize = (15, 15))
    plt.plot(r, vel_halo_m33, label = 'Halo Velocity', color = 'b')
    plt.plot(r, vel_disk_m33, label = 'Disk Velocity',  color = 'r')
    plt.plot(r, total_vel_m33, label = 'Total Velocity', color = 'm')
    plt.scatter(r, hern_vel_m33, color = 'k', label = 'Hernquist Velocity Profile, a = ' + str(scale_m33_v), marker = '*', s = 60)
    plt.xlabel('Distance (kpc)', fontsize = 35)
    plt.ylabel('Velocity (km/s)', fontsize = 35)
    plt.xticks(size = 20)
    plt.yticks(size = 20)
    plt.legend(fontsize = 20)
    plt.title('M33 Velocity Distribution', fontsize = 35)
    plt.semilogy()
    plt.show()
    
    
    
    
    
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
    