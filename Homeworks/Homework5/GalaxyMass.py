# Author: Ritvik Basant
# Course: ASTRONOMY 400 B, Homework 3
# This python file contains the ComponentMass function that takes the file containing particle data for a specific galaxy. Then, 
# this function (ComponentMass) uses the input data file and the input galaxy component name to calculate the total mass of that
# component in the galaxy. For example, the two inputs: 'MW_000.txt' and 'Halo' would give you the total dark matter mass present
# in the Milky Way galaxy. 

# Importing the Read function for ReadFile.py package. 
from ReadFile import Read

# Importing libraries. 
import numpy as np
import astropy.units as u

def ComponentMass(file_name, par_type):
    '''This function will calculate the total mass of a specified componen in a specified galaxy by summing over all
    the available particles of a specific type in the data file. 
    
    Input: 
        file_name: string
                    This is the data file that contains the particle data for a galaxy. If the location of this file is not in the same
                    directory, add the path to this string. 
        par_type: string
                    This is the name of the component whose total mass we want to find.
                    
    Output: This function returns the total mass of the specifed component in the specified galaxy.'''
    
    ''' This part of the code uses the second input (par_type) to convert it to a number appropriate for the original 
    data file'''
    if par_type == 'Halo':
        particle_type = 1.0 # This variable stores the particle's type as a number if particle is Dark Matter. 
    elif par_type == 'Disk':
        particle_type = 2.0 # This variable stores the particle's type as a number if particle is in Disk.
    elif par_type == 'Bulge':
        particle_type = 3.0 # This variable stores the particle's type as a number if particle is in Halo.
        
    time, no_of_particles, data = Read(file_name)
    # This line runs the Read function from the ReadFile.py and then stores the time (as time), number of 
    # particles (as no_of_particles), and the particle's data array (as the data array)
        
    '''This part of the code uses the particle's type number (particle_type) to extract data of all those particles
    that have the input data type. Then it stores this data in separate arrays each specified for a specific quantity. 
    For example, for this homework, we only need to store the mass of the particles so we just store that.'''
    index = np.where(data['type'] == float(particle_type)) 
    # This variable stores all the indices for the input particle type from the actual data array.
    mass = data['m'][index]
    # This array stores the mass for all particles with the given particle type.
    
    total_mass = (np.sum(mass))
    # This sums over all the specified type particles and the variable total_mass stores it. 
    
    return np.round((total_mass * 1e+10) * u.M_sun, 3)
    # Returning the total mass in the units of 1e+12 solMass rounded of to 3 decimal places. 

def main():    
    '''This is the main function. Here we call the ComponentMass function multiple times to calculate the total mass of the three
    components in each of the three galaxies.'''
    
    # Here we calculate the total mass of the three components in the M33 galaxy using ComponentMass function.
    m33_1 = ComponentMass('M33_000.txt', 'Halo')/u.M_sun  # Total dark matter mass
    m33_2 = ComponentMass('M33_000.txt', 'Disk')/u.M_sun  # Total disk star mass  
    m33_3 = ComponentMass('M33_000.txt', 'Bulge')/u.M_sun # Total bulge star mass
    
    # Here we calculate the total mass of the three components in the Milky Way galaxy using ComponentMass function.
    mw_1 = ComponentMass('MW_000.txt', 'Halo')/u.M_sun  # Total dark matter mass
    mw_2 = ComponentMass('MW_000.txt', 'Disk')/u.M_sun  # Total disk star mass  
    mw_3 = ComponentMass('MW_000.txt', 'Bulge')/u.M_sun # Total bulge star mass

    # Here we calculate the total mass of the three components in the M31 galaxy using ComponentMass function.    
    m31_1 = ComponentMass('M31_000.txt', 'Halo')/u.M_sun  # Total dark matter mass
    m31_2 = ComponentMass('M31_000.txt', 'Disk')/u.M_sun  # Total disk star mass  
    m31_3 = ComponentMass('M31_000.txt', 'Bulge')/u.M_sun # Total bulge star mass
    
    # Here, we calculate the total star mass for each of the three galaxy by summing both the disk and bulge stellar masses. 
    m33_star_mass = np.round(m33_2+m33_3, 3) # Total star mass, rounded to 3 decimal places, for M33 galaxy.
    m31_star_mass = np.round(m31_2+m31_3, 3) # Total star mass, rounded to 3 decimal places, for M31 galaxy.
    mw_star_mass = np.round(mw_2+mw_3, 3)    # Total star mass, rounded to 3 decimal places, for Milky Way galaxy.
  
    # Here, we calculate the total mass for each of the three galaxy by summing all components: the disk, bulge stellar masses, and the DM mass.     
    m33_totalmass = np.round(m33_1 + m33_2 + m33_3, 3) # Total mass, rounded to 3 decimal places, for M33 galaxy.
    m31_totalmass = np.round(m31_1 + m31_2 + m31_3, 3) # Total mass, rounded to 3 decimal places, for M31 galaxy.
    mw_totalmass = np.round(mw_1 + mw_2 + mw_3, 3)     # Total mass, rounded to 3 decimal places, for Milky Way galaxy.
    
    # Here we calculate the baryonic fraction for each galaxy (f_baryonic)
    f_m33 = np.round(m33_star_mass / m33_totalmass, 3) # Baryonic fraction = total_star_mass / total_mass, for M33 galaxy.
    f_m31 = np.round(m31_star_mass / m31_totalmass, 3) # Baryonic fraction = total_star_mass / total_mass, for M31 galaxy.
    f_mw = np.round(mw_star_mass / mw_totalmass, 3)    # Baryonic fraction = total_star_mass / total_mass, for Milky Way galaxy.
    
    # Here, we calculate the total mass components for the local group.
    lg_star_mass = np.round(m33_star_mass + m31_star_mass + mw_star_mass, 3) 
    # Caclulating the total stellar mass in the LG by summing over all galaxies. 
    lg_dark_mass = np.round(m33_1 + m31_1 + mw_1, 3)  # Calculating the total dark matter mass in the LG.
    lg_disk_mass = np.round(m33_2 + m31_2 + mw_2, 3)  # Calculating the total disk mass in the LG. 
    lg_bulge_mass = np.round(m33_3 + m31_3 + mw_3, 3) # Calculating the total bulge mass in the LG. 
    
    
    lg_total_mass = np.round(lg_star_mass + lg_dark_mass, 3)  # Calculating the total mass in the LG. 
    f_local_group = np.round(lg_star_mass / lg_total_mass, 3) # Calculating the baryonic fraction of the LG. 


    '''Here, we write all the important quantities in the a text file.'''
    with open('Galaxy_Component_Masses.txt', 'w') as f:
        # Writing the header for the columns of the txt table file. 
        f.write('Galaxy_Name    Halo_Mass(1e+12 solMass)    Disk_Mass(1e+12 solMass)    Bulge_Mass(1e+12 solMass)    Total_Mass(1e+12 solMass)    f_bar\n')
        # Writing the calculated quantities for Milky Way galaxy.
        f.write('MW    '+str(mw_1)+'    '+str(mw_2)+'    '+str(mw_3)+'    '+str(mw_totalmass)+'    '
                +str(f_mw)+'\n')
        # Writing the calculated quantities for M31 galaxy.
        f.write('M31    '+str(m31_1)+'    '+str(m31_2)+'    '+str(m31_3)+'    '+str(m31_totalmass)+'    '
                +str(f_m31)+'\n')
        # Writing the calculated quantities for M33 galaxy.
        f.write('M33    '+str(m33_1)+'    '+str(m33_2)+'    '+str(m33_3)+'    '+str(m33_totalmass)+'    '
                +str(f_m33)+'\n')
        # Writing the calculated quantities for the Local Group.
        f.write('LG    '+str(lg_dark_mass)+'    '+str(lg_disk_mass)+'    '+str(lg_bulge_mass)+'    '+str(lg_total_mass)+'    '
                +str(f_local_group)+'\n')
    f.close()
