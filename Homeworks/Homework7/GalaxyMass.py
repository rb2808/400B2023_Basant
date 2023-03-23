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
    '''This function will calculate the total mass of a specified component in a specified galaxy by summing over all
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
    
    particle_type = 0
    
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
    
    return np.round((total_mass/100) * u.M_sun, 3)
    # Returning the total mass in the units of 1e+12 solMass rounded of to 3 decimal places. 

#print(ComponentMass('M31_000.txt', 1))