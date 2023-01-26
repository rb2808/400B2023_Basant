# Author: Ritvik Basant
# Course: ASTRONOMY 400 B, Homework 2
# This python file contains the ParticleInfo function that uses the data (x, y, z, vx, vy, vz, m, and type) 
# to calculate the distance (in kpc), speed (in km / s), and mass of the selected object. It also uses the Read function
# from ReadFile python file to unpack the MW_000.txt file. 

# Importing the Read function for ReadFile.py package. 
from ReadFile import Read

# Importing libraries. 
import numpy as np
import astropy.units as u

# This is the function that first checks which particle type is asked for and then converts it to the type number 
# appropriate for the original data file and then stores it. It then calls for the Read function from the ReadFile.py
# and using the input file, it runs the Read function to extract the total number of particles, the time, and the data
# of particles as an array. Then using the np.where function, this code selects the type of particle and then extracts 
# its data (x, y, z, vx, vy, vz, and m) from the data array using the specific particle number provided. Then, it uses 
# these values to calculate the 3-D distance and velocity of the selected object. Finally, it prints the total distance 
# (in kpc and lyrs), 3-d velocity (in km/s), and the mass (in units of solar masses). 
# Inputs: file_name (This is the name of file that contains the data we need. Add the actual path to the file as well. ), 
# par_type (This is the name of the particle: Dark Matter, Disk, or Halo), and particle_number (This is the particle 
# Number: 1, 2, 3, ....)
# Returns: This function will return the 3-D distance to the object (in kpc and lyrs), 3-D velocity of the object (in 
# km / s), and the mass of the selected object (in solar masses).
def ParticleInfo(file_name, par_type, particle_number):
    
    ''' This part of the code uses the second input (par_type) to convert it to a number appropriate for the original 
    data file'''
    if par_type == 'Dark Matter':
        particle_type = 1.0 # This variable stores the particle's type as a number if particle is Dark Matter. 
    elif par_type == 'Disk':
        particle_type = 2.0 # This variable stores the particle's type as a number if particle is in Disk.
    elif par_type == 'Halo':
        particle_type = 3.0 # This variable stores the particle's type as a number if particle is in Halo.
    
    time, no_of_particles, data = Read(file_name)
    # This line runs the Read function from the ReadFile.py and then stores the time (as time), number of 
    # particles (as no_of_particles), and the particle's data array (as the data array)
    
    kmps = u.km / u.s
    # Creating new unit (km /s) using the units in astropy 
    
    '''This part of the code uses the particle's type number (particle_type) to extract data of all those particles
    that have the input data type. Then it stores this data in separate arrays each specified for a specific quantity. 
    For example, x array will store the given particle type's x-position for all the particles.'''
    index = np.where(data['type'] == float(particle_type)) 
    # This variable stores all the indices for the input particle type from the actual data array.
    x = data['x'][index]
    # This array stores the x position for all particles with the given particle type.
    y = data['y'][index]
    # This array stores the y position for all particles with the given particle type.
    z = data['z'][index]
    # This array stores the z position for all particles with the given particle type.
    mass = data['m'][index]
    # This array stores the mass for all particles with the given particle type.
    vx = data['vx'][index]
    # This array stores the velocities in x direction for all particles with the given particle type.
    vy = data['vy'][index]
    # This array stores the velocities in y direction for all particles with the given particle type.
    vz = data['vz'][index]
    # This array stores the velocities in z direction for all particles with the given particle type.
    
    par_no = particle_number - 1
    # This variable converts the particle number into the appropriate index for the data file and then stores it.
    
    distance_to_object = np.round((np.sqrt((x[par_no]**2) + (y[par_no**2]) + (z[par_no]**2))), 3) * u.kpc
    # Calculates the total 3-D distance of the given particle by adding the squares of distances in x, y, and z and 
    # then taking the square root. It then stores it as in units of kpc using astropy. It is roundedto 3 decimal 
    # places. 
    speed_of_object = np.round((np.sqrt((vx[par_no]**2) + (vy[par_no]**2) + (vz[par_no]**2))), 3) * kmps
    # Calculates the total 3-D velocity of the given particle  by adding the squares of velocities in x, y, and z and 
    # then taking the square root. It then stores it as in units of km / s by multiplying it by kmps. It is rounded 
    # to 3 decimal places. 
    mass_of_object = mass[par_no] * 1e+10 * u.M_sun
    # Extracts the mass of the given particle from the data array. It then stores it in units of solar masses using 
    # astropy. 
    
    d_in_lyrs = np.round(distance_to_object.to(u.lyr), 3)
    # This variable converts the distance of the object from kpc into light years using the .to 
    # function and astropy, and then stores it. It is rounded to three decimal places. 
    
    '''This three lines print out the results'''
    print(f'The distance to the selected object is: {distance_to_object: .3e} or {d_in_lyrs: .3e}')
    print(f'The speed of the selected object is: {speed_of_object: .3e}')
    print(f'The mass of the selected object is: {mass_of_object: .3e}')

# This is the variable that stores the file name. Also, keep the data file in the same directory or else, change this 
# variable to account for the new path of the file. 
file_nm = 'MW_000.txt'

# Calls the ParticleInfo function and prints out the results. As of now, we have selected the 100th particle in the disk.
# Change it to anything else if want. 
ParticleInfo(file_nm, 'Disk', 100)