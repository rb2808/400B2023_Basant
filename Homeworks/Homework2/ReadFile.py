# Author: Ritvik Basant
# Course: ASTRONOMY 400 B, Homework 2
# This python file contains the Read function that reads in the data from the MW_000.txt file 
# containing the particle masses, type, distance (x, y, z), and velocities (x, y, z).

# Importing the Libraries
import numpy as np
import astropy.units as u

# This is the function that reads in the data from the MW_000.txt file containing the time, particle numbers, masses, type,
# distance (x, y, z), and velocities (x, y, z).
# Inputs: file_name (This is the name of file that contains the data we need. Add the actual path to the file as well. )
# Returns: This function will return the time (in units if Myears), total number of particles in file both as variables, 
# and will also return the particle type, particle's mass (in units of 1e+10 solar masses), particle's x, y, z positions 
# (all in units of kpc) and particle's velocities in x, y, z directions (in units of km /s).
def Read(file_name):
    file = open(file_name, 'r') 
    # Opening up the file in read mode.
    time = float((file.readline().split('    ')[1].strip('\n')).strip(' ')) * u.Myr 
    # Extracting the time from the first line of the input file and converting it into Million years and storing it.                                                                                     
    total_particles = float((file.readline().split('    ')[1].strip('\n')).strip(' '))   
    # Extracting the total number of particles from the second line of the input file and storing it.
    (file.readline())
    # This just iterates over the heading line (line 3) of the data file. 
    particle_data = np.genfromtxt(file, dtype = None, names = True)
    # Extracting the data (x, y, z, vx, vy, vz, m, and type) of particles from remaining lines and storing them in an
    # array. 
    return time, total_particles, particle_data
    # Returning everything (time, total number of particles, and particle data) back as variables. 


    