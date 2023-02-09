# Author: Ritvik Basant
# Course: ASTRONOMY 400 B, Homework 4
# This python file contains the CenterofMass class that includes multiple other functions necessary
# to calculate the center of mass position and velocity for any galaxy. 

# Importing modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

# Importing Read function from ReadFile.py
from ReadFile import Read

class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot

    def __init__(self, filename, ptype):
        ''' Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            
            PARAMETERS
            ----------
            filename : `str`
                snapshot file
            ptype : `int; 1, 2, or 3`
                particle type to use for COM calculations
        '''
     
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m'][self.index]
        
        # write your own code to complete this for positions and velocities
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index] 
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]


    def COMdefine(self,a,b,c,m):
        ''' Method to compute the COM of a generic vector quantity by direct weighted averaging.
        
        PARAMETERS
        ----------
        a : `float or np.ndarray of floats`
            first vector component
        b : `float or np.ndarray of floats`
            second vector component
        c : `float or np.ndarray of floats`
            third vector component
        m : `float or np.ndarray of floats`
            particle masses
        
        RETURNS
        -------
        a_com : `float`
            first component on the COM vector
        b_com : `float`
            second component on the COM vector
        c_com : `float`
            third component on the COM vector
        '''
        
        # Here I am just printing the sum of m and the m array to see them.
        # Comment out or not, upto you. 
        # print('Here is sum of m:', np.sum(m))
        # print('m: ', m)
        
        # x - component of Center of mass
        a_com = np.dot(a, m) / np.sum(m)

        # Y - component of Center of mass
        b_com = np.dot(b, m) / np.sum(m)
        
        # Z - component of Center of mass
        c_com = np.dot(c, m) / np.sum(m)
        
        # Returns the 3 components separately
        return a_com, b_com, c_com
    
    
    def COM_P(self, delta = 0.1):
        '''Method to compute the position of the center of mass of the galaxy 
        using the shrinking-sphere method.

        PARAMETERS
        ----------
        delta : `float, optional`
            error tolerance in kpc. Default is 0.1 kpc
        
        RETURNS
        ----------
        p_COM : `np.ndarray of astropy.Quantity'
            3-D position of the center of mass in kpc
        '''                                                                     

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        x_COM, y_COM, z_COM = self.COMdefine(self.x, self.y, self.z, self.m)
        
        # Computing the magnitude of the COM position vector.
        r_COM = np.sqrt(x_COM**2 + y_COM**2 + z_COM**2)

        
        # iterative process to determine the center of mass                                                            

        # Changing the reference frame to COM frame                                                                          
        # Computing the difference between particle coordinates                                                          
        # and the first guess at COM position.
        x_new = self.x - x_COM
        y_new = self.y - y_COM
        z_new = self.z - z_COM
        
        # Calculating the magnitude of position vector.
        r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)

        # Finding the max 3D distance of all particles from the guessed COM                                               
        # and then calculating the reduced r_max. 
        r_max = max(r_new)/2.0
        
        # Picking an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        change = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (change > delta):
            # Selecting all particles within the reduced radius (starting from original x,y,z, m)
            self.index2 = np.where(r_new <= r_max) 
            x2 = self.x[self.index2]
            y2 = self.y[self.index2]
            z2 = self.z[self.index2]
            m2 = self.m[self.index2]

            # Refined COM position:                                                                                    
            # Computing the center of mass position using                                                                
            # the particles in the reduced radius
            x_COM2, y_COM2, z_COM2 = self.COMdefine(x2, y2, z2, m2)
            
            # Computing the new 3D COM position
            r_COM2 = np.sqrt(x_COM2**2 + y_COM2**2 + z_COM2**2)

            # Determining the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            change = np.abs(r_COM - r_COM2)
            
            # uncomment the following line if you want to check this                                                                                               
            #print ("CHANGE = ", change)                                                                                     

            # Before loop continues, reset : r_max, particle separations and COM                                            

            # reduce the volume by a factor of 2 again                                                                 
            r_max /= 2.0
            
            # check this.                                                                                              
            #print ("maxR", r_max)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            x_new = self.x - x_COM2
            y_new = self.y - y_COM2
            z_new = self.z - z_COM2
            # Calculating the position vector magnitude
            r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)

            # set the center of mass positions to the refined values                                                   
            x_COM = x_COM2
            y_COM = y_COM2
            z_COM = z_COM2
            r_COM = r_COM2

        # Creating an array (np.array) to store the COM position                                                                                                                                                       
        p_COM = np.array([round(x_COM, 2), round(y_COM, 2), round(z_COM, 2)])

        # set the correct units using astropy and round all values
        # and then return the COM positon vector
        return p_COM*u.kpc
        
        
    def COM_V(self, x_COM, y_COM, z_COM):
        ''' Method to compute the center of mass velocity based on the center of mass
        position.

        PARAMETERS
        ----------
        x_COM : 'astropy quantity'
            The x component of the center of mass in kpc
        y_COM : 'astropy quantity'
            The y component of the center of mass in kpc
        z_COM : 'astropy quantity'
            The z component of the center of mass in kpc
            
        RETURNS
        -------
        v_COM : `np.ndarray of astropy.Quantity'
            3-D velocity of the center of mass in km/s
        '''
        
        # the max distance from the center that we will use to determine 
        # the center of mass velocity                   
        rv_max = 15.0*u.kpc

        # Determining the position of all particles relative to the center of mass position (x_COM, y_COM, z_COM)
        xV = self.x - x_COM/u.kpc
        yV = self.y - y_COM/u.kpc
        zV = self.z - z_COM/u.kpc
        
        # Determining the magnitude of 3-D velocity
        rV = np.sqrt(xV**2 + yV**2 + zV**2)
        
        # Determining the index for those particles within the max radius
        indexV = np.where(rV <= (rv_max/u.kpc))
        
        # Determining the velocity and mass of those particles within the mas radius
        vx_new = self.vx[indexV]
        vy_new = self.vy[indexV]
        vz_new = self.vz[indexV]
        m_new =  self.m[indexV]
        
        # Computing the center of mass velocity using those particles
        vx_COM, vy_COM, vz_COM = self.COMdefine(vx_new, vy_new, vz_new, m_new)
        
        # Creating an array to store the COM velocity
        v_COM = np.array([round(vx_COM, 2), round(vy_COM, 2), round(vz_COM, 2)])

        # return the COM vector
        # set the correct units usint astropy
        # round all values                                                                                        
        return v_COM*(u.km/u.s)
    

# ANSWERING QUESTIONS
#######################
if __name__ == '__main__' : 

    # Create a Center of mass object for the MW, M31 and M33
    # below is an example of using the class for MW
    MW_COM = CenterOfMass("MW_000.txt", 2)

    # below gives you an example of calling the class's functions
    # MW:   store the position and velocity COM
    MW_COM_p = MW_COM.COM_P(0.1)
    print('This is final COM (Milky Way):', MW_COM_p, '; Magnitude:', np.round(np.sqrt(np.dot(MW_COM_p, MW_COM_p)), 3))
    MW_COM_v = MW_COM.COM_V(MW_COM_p[0], MW_COM_p[1], MW_COM_p[2])
    print('This is fincal COMV (Milky Way):',MW_COM_v, '; Magnitude:', np.round(np.sqrt(np.dot(MW_COM_v, MW_COM_v)), 3))

    # now write your own code to answer questions
    
    # For M31:
    M31_COM = CenterOfMass("M31_000.txt", 2)
    M31_COM_p = M31_COM.COM_P(0.1)
    print('This is final COM (M 31):', M31_COM_p, '; Magnitude:', np.round(np.sqrt(np.dot(M31_COM_p, M31_COM_p)), 3))
    M31_COM_v = M31_COM.COM_V(M31_COM_p[0], M31_COM_p[1], M31_COM_p[2])
    print('This is fincal COMV (M 31):',M31_COM_v, '; Magnitude:', np.round(np.sqrt(np.dot(M31_COM_v, M31_COM_v)), 3))
    
    # For M33:
    M33_COM = CenterOfMass("M33_000.txt", 2)
    M33_COM_p = M33_COM.COM_P(0.1)
    print('This is final COM (M 33):', M33_COM_p, '; Magnitude:', np.round(np.sqrt(np.dot(M33_COM_p, M33_COM_p)), 3))
    M33_COM_v = M33_COM.COM_V(M33_COM_p[0], M33_COM_p[1], M33_COM_p[2])
    print('This is fincal COMV (M 33):',M33_COM_v, '; Magnitude:', np.round(np.sqrt(np.dot(M33_COM_v, M33_COM_v)), 3))
    
    # Distance and Velocity for MW and M31
    dmwm31 = np.round(np.sqrt(np.dot(np.abs(MW_COM_p - M31_COM_p), np.abs(MW_COM_p - M31_COM_p))), 3) 
    # This is the distance between both galaxies. 
    vmwm31 = np.round(np.sqrt(np.dot(np.abs(MW_COM_v - M31_COM_v), np.abs(MW_COM_v - M31_COM_v))), 3)
    # This is the relative velocity between the two galaxies. 
    print('Separation between MW and M31:', dmwm31)
    print('Relative Velocity between MW and M31:', vmwm31)
    
    
    # Distance and Velocity for MW and M31
    dm33m31 = np.round(np.sqrt(np.dot(np.abs(M33_COM_p - M31_COM_p), np.abs(M33_COM_p - M31_COM_p))), 3)
    # This is the distance between both galaxies. 
    vm33m31 = np.round(np.sqrt(np.dot(np.abs(M33_COM_v - M31_COM_v), np.abs(M33_COM_v - M31_COM_v))), 3)
    # This is the relative velocity between the two galaxies. 
    print('Separation between M33 and M31:', dm33m31)
    print('Relative Velocity between M33 and M31:', vm33m31)
