# Author: Ritvik Basant
# Course: ASTRONOMY 400 B, Homework 6
# This python file contains OrbitCOM function that iterates over multiple snapshot files for different galaxies
# and calculates the orbit parameters (x,y,z,vx,vy,vz) for various galaxies. 

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
import sys

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# Importing other important modules.
from ReadFile1 import Read1
from CenterOfMass import CenterOfMass

def OrbitCOM(galaxy, start = 0, end = 801, stepsize = 5):
    """Function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
        galaxy: string
            Identifies galaxy name (e.g., MW, M33, or M31)
        
        start: integer
            The first snapshot to be read in.
            
        end: integer
            The last snapshot to be read in. 
            
        stepsize: integer
            The interval over which the function returns the COM. 
          
    outputs: 
        Orbit_galaxy: txt file
            This txt file contains orbital parametrs from the simulation of the galaxy.
    """
    
    # compose the filename for output
    fileout = 'Orbit_'+str(galaxy)+'.txt'
    
    # set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta = 0.1
    volDec = 0
    
    # This is the bonus part: It check whether 'start', 'end', and 'stepsize' are all integers. If they are
    # not all integers, then the step_ids array would give an error. Thus to see if the step_ids array is empty
    # or not, this was the way to check. 
    try:         
        
        # This changes the volDec parameter based on the galaxy.
        if galaxy =='MW':
            volDec = 2.0
        elif galaxy =='M31':
            volDec = 2.0
        elif galaxy =='M33':
            volDec = 4.0
            
        # generate the snapshot id sequence 
        snap_ids = np.arange(start, end, stepsize)
        
        # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
        orbit = np.zeros([len(snap_ids), 7])
        
        # a for loop over files
        for i, snap_id in enumerate(snap_ids):
            
            # compose the data filename (be careful about the folder)
            ilbl = '000' + str(snap_id)
            ilbl = ilbl[-3:]
            file_name = "%s_"%(galaxy+'/'+galaxy) + ilbl + '.txt'
            
            # Initialize an instance of CenterOfMass class, using disk particles
            COM = CenterOfMass(file_name, 2, volDec)
    
            # Store the COM pos and vel. Remember that now COM_P required VolDec
            com_pos = COM.COM_P(delta)
            com_vel = COM.COM_V(com_pos[0], com_pos[1], com_pos[2])
            
            # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
            # note that you can store 
            # a[i] = var1, *tuple(array1)
            orbit[i] = (COM.time).value/1000, com_pos[0].value, com_pos[1].value, com_pos[2].value, com_vel[0].value, com_vel[1].value , com_vel[2].value  
            
            # print snap_id to see the progress
            print('Simulation in progress. Snap_id #'+str(snap_id))
            
        # write the data to a file
        # we do this because we don't want to have to repeat this process 
        # this code should only have to be called once per galaxy.
        np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                          .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
    
    # If start, end, or n is not an integers, then the code will throw an error. 
    except TypeError:
        print("'start', 'end', and 'n' have to be integers")
        
# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 
# Uncomment if you want to run the simulation. 

#OrbitCOM('MW')
#OrbitCOM('M31')
#OrbitCOM('M33')

# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  

def dist_vel(array1, array2):
    '''
    This function takes in two arrays (orbital parameters for two different galaxies) and then 
    returns the time, distance between two galaxies, and the relative velocities of the two 
    galaxies. 

    Parameters
    ----------
    array1 : 'array'
        This is the orbital parameters for one galaxy. 
    array2 : 'array'
        This is the orbital parameters for other galaxy. 

    Returns
    -------
    time : 'array'
        This is the time for each of the iteration of the simulation. 
    distance : 'array'
        The array containing relative distance between two galaxies as a function of time. 
    velocity : 'array'
        The array containing relative velocities between two galaxies as a function of time. 
    
    '''
    # Defining array for storing time.
    time = np.zeros([800, 1])
    
    # Defining array for storing relative distance.
    distance = np.zeros([800, 1])
    
    # Defining array for storing relative velocity.
    velocity = np.zeros([800, 1])
    
    # For loop to iterate over all the data points and then storing the magnitude of relaive distance and velocity. 
    # Here, time is also stored. 
    for i in range(len(array1)):
        time[i] = array1[i][0]
        distance[i] = np.sqrt((array1[i][1] - array2[i][1])**2 + (array1[i][2]-array2[i][2])**2 + (array1[i][3] - array2[i][3])**2)
        velocity[i] = np.sqrt((array1[i][4] - array2[i][4])**2 + (array1[i][5]-array2[i][5])**2 + (array1[i][6] - array2[i][6])**2)

    return time, distance, velocity
