# Author: Ritvik Basant
# Course: ASTRONOMY 400 B, Homework 7
# This python file contains a class that has OrbitIntegration function that can 
# integrate the orbit of two galaxies from their initial position and velocity vectors. 
# This integration uses Leap Frog method to integrate the orbits. 

# Importing essential modules
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const

# Importing other necessary functions from previous homeworks.
from CenterOfMass2 import CenterOfMass
from GalaxyMass import ComponentMass
from ReadFile1 import Read1
from OrbitCOM import dist_vel

# Defining the class M33AnalyticOrbit that will calculate the integrated orbits. 
class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, name_file): # **** add inputs
        '''
        Class to calculate the orbit evolution of M33 and M31 using Leap Frog integration. 
        ------
        Parameters:
            name_file: 'string'
                This is the txt file in which the integrated orbit will be stored. 
        '''

        # Get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        # Store the output file name
        self.filename = 'Integrated_Orbit_'+(name_file)+'.txt'
        
        # Get the current pos & vel of M33 using COM Class created in previous homeworks. 
        # Creating an instance of the CenterOfMass class for M33 
        COM_M33 = CenterOfMass('M33_000.txt', 2)
        # Storing the position VECTOR of the M33 COM
        r_m33 = (COM_M33.COM_P(0.1, 4.0))
        # Storing the velocity VECTOR of the M33 COM. 
        v_m33 = COM_M33.COM_V(r_m33[0], r_m33[1], r_m33[2]).value
        # Getting rid of the units in position vectors. 
        r_m33 = r_m33.value
        
        # Getting the current pos & vel of M31. 
        # Creating an instance of the  CenterOfMass class for M31 
        COM_M31 = CenterOfMass('M31_000.txt', 2)
        # Storing the position VECTOR of the M31 COM.
        r_m31 = COM_M31.COM_P(0.1, 2)
        # Storing the velocity VECTOR of the M31 COM. 
        v_m31 = COM_M31.COM_V(r_m31[0], r_m31[1], r_m31[2]).value
        # Getting rid of units in position vector. 
        r_m31 = r_m31.value
        
        # Storing the DIFFERENCE between the vectors posM33 - posM31
        # Creating two VECTORs self.r and self.v and have them be the
        # relative position and velocity VECTORS of M33
        self.r = np.array([r_m33[0] - r_m31[0], r_m33[1] - r_m31[1], r_m33[2] - r_m31[2]])
        self.v = np.array([v_m33[0] - v_m31[0], v_m33[1] - v_m31[1], v_m33[2] - v_m31[2]])
        
        # Gettting the mass of each component in M31 
        # Disk
        # Self.rdisk = scale length (no units)
        self.rdisk = 5
        # Self.Mdisk: set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk = ComponentMass('M31_000.txt', 'Disk').value*1e12
        # Bulge
        # Self.rbulge = set scale length (no units)
        self.rbulge = 1
        # Self.Mbulge: set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = ComponentMass('M31_000.txt', 'Bulge').value*1e12
        # Halo
        # Self.rhalo: set scale length from HW5 (no units)
        self.rhalo = 61.5
        # Self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = ComponentMass('M31_000.txt', 'Halo').value*1e12
    
    
    def HernquistAccel(self, M, r_a, R): 
        '''
        This function calculates the acceleration of the galaxy using Hernquist mass profile
        (for bulge and halo).

        Parameters
        ----------
        M : 'dimensionless quantity'
            It is the mass of a specfici component of the galaxy.
        r_a : 'dimensionless quantity'
            It is the scale length used in Hernquist mass profile. 
        R : 'np.array'
            This is a position vector for the relative position of the galaxy.

        Returns
        -------
        Hern : 'np.array'
            This is the acceleration vector of the galaxy (specific component).

        '''
        
        # Storing the magnitude of the position vector
        rmag = np.sqrt(R[0]**2 + R[1]**2 + R[2]**2)
        
        # Storing the Acceleration
        # Formula for acceleration = np.array([(-self.G*M/(rmag*((rmag + r_a)**2)))*R[0], 
        # (-self.G*M/(rmag*((rmag + r_a)**2)))*R[1],(-self.G*M/(rmag*((rmag + r_a)**2)))*R[2]])  
        # follow the formula in the HW instructions
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # Using -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        
        Hern = (-self.G*M/(rmag*((rmag + r_a)**2)))*R
        
        return Hern
    
    
    
    def MiyamotoNagaiAccel(self, M, r_d, pos_vec):
        '''
        This function uses Miyamoto Nagai method to calculate the acceleration of disk component
        of a galaxy. 

        Parameters
        ----------
        M : 'dimensionless quantity'
            This is the mass of a specific component of a galaxy. 
        r_d : 'dimensionless quantity'
            This is the scale length for the Miyamoto Nagai. 
        pos_vec : 'np.array'
            This is the position vector of the galaxy. 

        Returns
        -------
        Miya : 'np.array'
            This is the acceleration of galaxy using Miyamoto Nagai acceleration. 

        '''

        
        # Acceleration follow the formula in the HW instructions
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        #  multiply the whle thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        
        # Parameter for Miyamoto Nagai Acceleration profile. 
        z_d = self.rdisk / 5.0
        
        # Parameter for Miyamoto Nagai Acceleration profile.
        B = r_d + np.sqrt(pos_vec[2]**2 + z_d**2)
        
        # This is the (x**2 + y**2)**0.5
        R = np.sqrt(pos_vec[0]**2 + pos_vec[1]**2)
        
        # Formula for Miyamoto Nagai acceleration.
        Miya = (-self.G*M/((R**2 + B**2)**1.5))*np.array([pos_vec[0], pos_vec[1], 
                B*pos_vec[2] / (np.sqrt(pos_vec[2]**2 + z_d**2))])

        return Miya     
    
    def M31Accel(self, r_vec): 
        '''
        This function calculates the total acceleration on a galaxy. 

        Parameters
        ----------
        r_vec : 'np.array'
            This is the position vector of a galaxy. 

        Returns
        -------
        acc_total : 'np.array'
            This is the acceleration vector of the galaxy. 

        '''
        
        # Calculating the acceleration of halo component using Hernquist profile. 
        halo_acc = self.HernquistAccel(self.Mhalo, self.rhalo, r_vec)
        
        # Calculating the acceleration of bulge component using Hernquist profile. 
        bulge_acc = self.HernquistAccel(self.Mbulge, self.rbulge, r_vec)

        # Calculating the acceleration of bulge component using Hernquist profile.
        disk_acc = self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, r_vec)


            
        # Summing of the output of the acceleration functions - this will return a VECTOR 
        acc_total = np.array([halo_acc[0] + bulge_acc[0] + disk_acc[0], halo_acc[1] + bulge_acc[1] + disk_acc[1], 
                              halo_acc[2] + bulge_acc[2] + disk_acc[2]])
        
        # Returning the acceleration vector. 
        return acc_total
    
    
    
    def LeapFrog(self, r, v, dt): 
        '''
        This function will integrate the orbit of a galaxy using Leap Frog method. It works both
        ways, forward and backward in time.

        Parameters
        ----------
        r : 'np.array'
            This is the current position of a galaxy. 
        v : 'np.array'
            This is the current velocity of a galaxy. 
        dt : 'dimensionless quantity'
            This is the time step on which integration happens. 

        Returns
        -------
        dist_and_vel: 'np.array'
            This is the new position and velocity vector of the galaxy. 

        '''
        
        # Predicting the position at the next half timestep. rhalf is the new predicted position. 
        rhalf = r + v*dt/2
        
        # Calculating the total acceleration of the galaxy using M31Accel function. 
        acc = self.M31Accel(rhalf)
        
        # Predicting the final velocity at the next timestep using the acceleration field 
        # at the rhalf position.  
        vnew = v + acc*dt
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        # rnew stores the new location of the galaxy. 
        rnew =  r + (dt/2)*(v + vnew)
        
        # This dist_and_vel is a two dimensional array that stores the individual new position
        # and velocity of the galaxy. 
        dist_and_vel = np.array([rnew, vnew])
        
        return dist_and_vel
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        '''
        This function will perform the actual integration of the orbit of a galaxy using
        bove created functions. 

        Parameters
        ----------
        t0 : 'dimensionless quantity'
            This is the starting time, i.e., the lower limit for the orbit integration. 
        dt : 'dimensionless quantity'
            It is the time step for the integration. 
        tmax : 'dimensionless quantity'
            It is the upper limit for the orbit integration. 

        Returns
        -------
        None.

        '''

        # Initializing the time to the input starting time
        t = t0
        
        # Initializing an empty array of size: rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros((int(tmax/dt) + 2, 7))
        
        # Initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r), *tuple(self.v)
        # This above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        
        # Initializing a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # Starting the integration (advancing in time steps and computing LeapFrog at each step)
        while (t < tmax): # as long as t has not exceeded the maximal time 
            print('Integrating Orbit. Running Counter:', i)
            
            # Advancing the time by one timestep, dt
            t = t + dt
           
            # Storing the new time in the first column of the ith row
            orbit[i][0] = t
            
            # Advancing the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like: a,b,c = function(input)
            
            # This is the position vector that takes the last row of the orbit file. 
            position_vec = np.array([orbit[i-1][1], orbit[i-1][2], orbit[i-1][3]])
            
            # This is the velocity vector that takes the last row of the orbit file. 
            velocity_vec = np.array([orbit[i-1][4], orbit[i-1][5], orbit[i-1][6]])
            
            # This stores the newly integrated position and velocity of the orbit. 
            new_position, new_velocity  = self.LeapFrog(position_vec, velocity_vec, dt)
    
            # Storing the new position vector into the columns with indexes 1,2,3 of the ith 
            # row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            
            
            # Storing the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            orbit[i][1:4] = new_position
            orbit[i][4:7] = new_velocity            
            # Update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i = i + 1

        
        # Writing the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # There is no return function

# ANSWERING QUESTIONS
#######################
if __name__ == '__main__' : 
    # Creating the object for integration of orbits of M33 adnM31
    orbm31m33 = M33AnalyticOrbit('Orbit M33-M31')
    
    # Calling the orbit integration for both galaxies' orbit. 
    orbm31m33.OrbitIntegration(0, 0.01, 10)
    
    # Reading the integrated orbit file. 
    data = Read1('Integrated_Orbit_Orbit M33-M31.txt')
    
    # Distance array. This will store position. 
    distance = np.zeros([len(data)])
    
    # Velocity array. This will store the velocity. 
    vel = np.zeros([len(data)])
    
    # Time array. This will store the time. 
    time = np.zeros([len(data)])
    
    # This loop will create the magnitude of distance and velocity for the integrated orbit. 
    for i in range(len(data)):
        # Updating distance, vel, time arrays.
        distance[i] = np.sqrt(data[i][1]**2 + data[i][2]**2 + data[i][3]**2)
        vel[i] = np.sqrt(data[i][4]**2 + data[i][5]**2 + data[i][6]**2)
        time[i] = data[i][0]
    
    # Plotting the distance and velocity from HW 7 (this one) and HW6.
    fig, axs = plt.subplots(2, figsize = (15, 10))
    fig.suptitle('Relative Distance & Velocity b/w M31 and M33', fontsize = 35)
    axs[0].scatter(time, distance, s = 5, label = 'HW7', color = 'b')
    axs[0].tick_params(axis="x", labelsize=15)
    axs[0].tick_params(axis="y", labelsize=15)
    axs[0].set_ylabel('Distance (Kpc)', fontsize = 25)
    axs[1].scatter(time, vel, s = 5, label = 'HW7', color = 'b')
    axs[1].set_xlabel('Time (Gyrs)', fontsize = 25)
    axs[1].tick_params(axis="x", labelsize=15)
    axs[1].tick_params(axis="y", labelsize=15)
    axs[1].set_ylabel('Velocity (km/s)', fontsize = 25)

    # Reading in the Orbits from HW 6.
    M31_data = Read1('Orbit_M31.txt')
    M33_data = Read1('Orbit_M33.txt')
    
    # Calling dist_vel function to calculate the relative distance and velocity 
    # of previously calculated orbits. 
    time_m33_m31, dist_m33_m31, vel_m33_m31 = dist_vel(M33_data, M31_data)
    axs[0].scatter(time_m33_m31, dist_m33_m31, s = 5, label = 'HW6', color = 'r')
    axs[1].scatter(time_m33_m31, vel_m33_m31, s = 5, label = 'HW6', color = 'r')
    axs[0].legend(fontsize = 20)
    axs[1].legend(fontsize = 20)