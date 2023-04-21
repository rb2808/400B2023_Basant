### Course: ASTR 400 B Semester Project
### Author: Ritvik Basant
### Goals: This code is designed to calculate the evolution of Spin Angular Momentum of a galaxy. 
### It uses the galaxy merger simulation data files (given by Prof. Besla). Now, each of the numerical simulation
### of the galaxy uses three types of particles - bulge, halo, and disk. Suppose that we want to calculate the 
### spin of the halo particles. In order to do so, we put in the particle type 1. Then, this code will first select
### the data from the simulation files (x,y,z,vx,vy,vz) for the given particle type. Then, the code will calculate 
### the angular momentum of each particle by multiplying the mass of the particle with a space coordinate and a 
### velocity component. For eg, the x of component of angular momentum L_x = m*v_x*x. The code then stores the 
### angular momentum of each particle in three different arrays -- each individual array storing the angular 
### momentum in each direction -- x,y,z. In the end, then code sums up all the elements in the L_i array and 
### divides it by total mass of the particles. Thus, by the code completely runs, it gives us three quantities -- 
### SAM_x, SAM_y, and SAM_z at a given time, where SAM stands for Spin angular (avg) momentum.  
### This whole loop is then iterated for all the 800 files for each galaxy simulation. The evolution of SAM in all
### three directions is stored in form of txt files. Once these files are computed, the code then re-reads them and 
### plots the SAM x, y, z for all three galaxies as a function of time. I have also changed the particle type to 2
### to calculate how the disk angular momentum evolves during a merger.

# Importing Libraries and other essential modules.
import numpy as np
import matplotlib.pyplot as plt
from ParticleProperties import ParticleInfo

def specific_ang_mom(filename, ptype):
    '''
    This function calculates the specific angular momentum (SAM) of a galaxy for any given 
    particle type (halo, disk, or bulge) based on the equation:
        
        SAM = [summation(m_i * v_i * r_i)] / [summation(m_i)], where 'i' denotes the i'th particle. 
        
    The summation is on all the particles for a given type. 
    
    As the SAM is a vector, we use the above equation to calculate the individual SAM in x, y, and 
    z directions and then store them as three separte variables. 

    Parameters
    ----------
    filename : 'string'
        This is the file that contains the particle simulation data for a specific time for a specific
        galaxy. 
    ptype : 'integer'
        This is the type of particle (1: Halo, 2: Disk, 3: Bulge) whose SAM we need to calculate.

    Returns
    -------
    Total_L : 'np.array'
        This array contains the net x, y, and z components of the SAM.
    time : 'float'
        This is the 'time snapshot' that has been extracted from the simulation file. 

    '''
    
    # Using the ParticleInfo function (created in previous homework), we store the data (m, x, y, z, vx, vy, vz)
    # of all the particles of the given type (ptype) from the input file (filename) as an array ('data). Additionally,
    #  we also store the time snapshot from the input file as a variable called 'time'.
    data, time = ParticleInfo(ptype, filename)
    
    # These np.arrays will store the SAM for individual directions (x,y,z). The size of the arrays has been created using
    # using the size of the data array created above. 
    ang_x = np.zeros(len(data[0]))
    ang_y = np.zeros(len(data[0]))
    ang_z = np.zeros(len(data[0]))
    
    # This array stores the mass of the all the required particles. 
    mass = data[0] 
    
    # This is the iteration loop where we iterate over all the particles in the data array and calculate the
    # Angular Momentum in each direction for every particle. 
    for i in range(len(data[0])):
        ang_x[i] = data[1][i] * data[4][i] * mass[i]
        ang_y[i] = data[2][i] * data[5][i] * mass[i]
        ang_z[i] = data[3][i] * data[6][i] * mass[i]
    
    # These variables store the total (sum of all particles' individual AM) Angular Momentum in each direction. 
    x_comp = sum(ang_x)
    y_comp = sum(ang_y)
    z_comp = sum(ang_z)

    # This variable stores the SAM for the input file for the given particle type as an array. 
    # To do this, we just divide the AM in each direction by the total mass. 
    Total_L = np.array([x_comp, y_comp, z_comp]) / sum(mass)    
        
    return Total_L, time


def SAM_evolution(galaxy):
    '''
    This function iterates over all simulation files for a given galaxy and stores the
    evolution of SAM as a function of time in a txt file. 

    Parameters
    ----------
    galaxy :'string'
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    # 4 rows are for time, x, y, z components of L.
    SAM_evol = np.zeros((802, 4))
    
    for i in range(802):
        ilbl = '000' + str(i)
        ilbl = ilbl[-3:]
        file_name = "%s_"%(galaxy+'/'+galaxy) + ilbl + '.txt'
        
        L, time = specific_ang_mom(file_name, 1)
        
        SAM_evol[i][0] = time.value
        SAM_evol[i][1], SAM_evol[i][2], SAM_evol[i][3] = L[0], L[1], L[2]
        
        print(i)
        
    np.savetxt(galaxy+'_Angular_Momentum.txt', SAM_evol, fmt = "%11.3f"*4, comments='#',
               header="{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('T (Gyrs)', 'Lx', 'Ly', 'Lz'))

### Uncomment this part if you want to generate the SAM evolution files. If this happens to be the case
### then please include the simulation files in the same directory and change the location in above
### lines accordingly. 

#d1 = SAM_evolution('M33')
#d2 = SAM_evolution('M31')
#d3 = SAM_evolution('MW')

def Read(filename):
    """ Function to read in our data file
    
    Input:  
        filename: str
            e.g. "MW_000.txt"
        
    Outputs: 
        time: astropy quantity
            Time of snapshot in Myr
        total: float
            Total number of particles 
        data: array of floats
            An array with the particle data, including position 
            vectors, velocity vectors and mass
            
    Example usage:  time, total, data = Read("filename")
    """
    
    
    # open the file 
    file = open(filename,'r')
    
    #read header info line by line (line will be a string)
    # read first two lines FIRST and store as variable
    
    file.readline()
    
    # close file
    file.close()

    # read the remainder of the file, 
    # "dtype=None" specifies data type. None is default float
    # default delimiter is line is split using white spaces
    # "skip_header=3"  skipping the first 3 lines 
    # the flag "names=True" creates arrays to store the date
    #       with the column headers given in line 4 like "m", "x"
    
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    
    # Note, another option is loadtxt, skipping the first 3 rows.  
    # data = np.loadtxt(filename,skiprows=3)
    # But this loses the information in the headers
    # this will return the time of the snapshot, 
    #total number of particles 
    #and an array that stores the remainder of the data. 
    return data
    

# Reading Angular momentum files created above to plot the evolution of SAM
# for MW and M31. 
MW_data = Read("MW_Angular_Momentum.txt")
M31_data = Read("M31_Angular_Momentum.txt")

# Lx tells the direction and mw tells the galaxy. These arrays will store
# SAM for Milky Way galaxy.
Lxmw = np.zeros((len(MW_data), 1))
Lymw = np.zeros((len(MW_data), 1))
Lzmw = np.zeros((len(MW_data), 1))

# Will store the total magnitude
L_mag_mw = np.zeros((len(MW_data), 1))


# Lx tells the direction and mw tells the galaxy. These arrays will store
# SAM for M31 galaxy.
Lxm31 = np.zeros((len(M31_data), 1))
Lym31 = np.zeros((len(M31_data), 1))
Lzm31 = np.zeros((len(M31_data), 1))

# Will store the total magnitude
L_mag_m31 = np.zeros((len(M31_data), 1))

# These three arrays will store x, y, and z components of SAM
# for the total remnant. (They star from the beginning but we only plot
# after the final merger).
rem_x = np.zeros(len(MW_data))
rem_y = np.zeros(len(MW_data))
rem_z = np.zeros(len(MW_data))

# Will store the magnitude of SAM of the remnant.
rem_mag = np.zeros(len(MW_data))

# Stores time. 
t = np.zeros((len(MW_data), 1))

# Halo Particle Mass MW/M31 (in e+10 but not necessary)
hpmmw = 0.007899703
hpmm31 = 0.007683523

# This for loops converts the data from the txt files into data arrays.
# It also calculates the magnitude of SAM. 
for i in range(len(MW_data)):
    # Storing time. 
    t[i] = MW_data[i][0]
    
    # Storing SAM of MW. 
    Lxmw[i] = MW_data[i][1]
    Lymw[i] = MW_data[i][2]
    Lzmw[i] = MW_data[i][3]
    L_mag_mw[i] = np.sqrt(MW_data[i][1]**2 + MW_data[i][2]**2 + MW_data[i][3]**2)
    
    # Storing SAM of M31
    Lxm31[i] = M31_data[i][1]
    Lym31[i] = M31_data[i][2]
    Lzm31[i] = M31_data[i][3]
    L_mag_m31[i] = np.sqrt(M31_data[i][1]**2 + M31_data[i][2]**2 + M31_data[i][3]**2)
    
    # Adding SAM of both galxies in the ratio of the masses of the halo particles. 
    x_remnant = (MW_data[i][1]*hpmmw + M31_data[i][1]*hpmm31)/(hpmmw + hpmm31)
    y_remnant = (MW_data[i][3]*hpmmw + M31_data[i][2]*hpmm31)/(hpmmw + hpmm31)
    z_remnant = (MW_data[i][3]*hpmmw + M31_data[i][3]*hpmm31)/(hpmmw + hpmm31)
    
    # Storing remnant SAM. 
    rem_x[i] = x_remnant
    rem_y[i] = y_remnant
    rem_z[i] = z_remnant
    
    # Storing remnant SAM Magnitude. 
    rem_mag[i] = np.sqrt(x_remnant**2 + y_remnant**2 + z_remnant**2)
    
# These lists are storing values to plot the x axis on the plots. 
xax_x = [min(t), max(t)]
xax_y = [0, 0]

# These three values are the time of the closest approaches as estimated 
# from Homework 6's plot of relative distance between M31 and MW. 
t1 = 3.957 # First closest approach (in Gyrs)
t2 = 5.857 # Second closest approach (in Gyrs)
t3 = 6.486 # Final Merger (in Gyrs)


# Plotting the evolution of SAM for every galaxy (x,y,z separate plots). 
'''
fig, axs = plt.subplots(2, figsize = (18, 12))
fig.suptitle('Specific Angular Momentum - X Direction', fontsize = 35)
axs[0].set_ylabel('MW [kpc-km-$s^{-1}$]', fontsize = 20)
axs[1].set_ylabel('M31 [kpc-km-$s^{-1}$]', fontsize = 20)
axs[0].tick_params(axis='both', labelsize=15)
axs[1].tick_params(axis='both', labelsize=15)
axs[0].plot(t, Lxmw, color = 'b')
axs[1].plot(t, Lxm31, color = 'b')
axs[0].plot(xax_x, xax_y, color = 'k', linewidth = 1)
axs[1].plot(xax_x, xax_y, color = 'k', linewidth = 1)
axs[1].set_xlabel('Time (Gyrs)', fontsize = 15)

### Plotting the merge times
mwlx = [min(Lxmw), max(Lxmw)]
m31lx = [min(Lxm31), max(Lxm31)]
lx1time = [t1, t1]
lx2time = [t2. t2]
lx3time = [t3, t3]
axs[0].plot(lx1time, mwlx, color = 'm', linewidth = 2)
axs[0].plot(lx2time, mwlx, color = 'm', linewidth = 2)
axs[0].plot(lx3time, mwlx, color = 'm', linewidth = 2)
axs[1].plot(lx1time, m31lx, color = 'm', linewidth = 2)
axs[1].plot(lx2time, m31lx, color = 'm', linewidth = 2)
axs[1].plot(lx3time, m31lx, color = 'm', linewidth = 2)

fig, axs = plt.subplots(2, figsize = (18, 12))
fig.suptitle('Specific Angular Momentum - Y Direction', fontsize = 35)
axs[0].set_ylabel('MW [kpc-km-$s^{-1}$]', fontsize = 20)
axs[1].set_ylabel('M31 [kpc-km-$s^{-1}$]', fontsize = 20)
axs[0].tick_params(axis='both', labelsize=15)
axs[1].tick_params(axis='both', labelsize=15)
axs[0].plot(t, Lymw, color = 'r')
axs[1].plot(t, Lym31, color = 'r')
axs[0].plot(xax_x, xax_y, color = 'k', linewidth = 1)
axs[1].plot(xax_x, xax_y, color = 'k', linewidth = 1)
axs[1].set_xlabel('Time (Gyrs)', fontsize = 15)

### Plotting the merge times
mwly = [min(Lymw), max(Lymw)]
m31ly = [min(Lym31), max(Lym31)]
lx1time = [t1, t1]
lx2time = [t2. t2]
lx3time = [t3, t3]
axs[0].plot(ly1time, mwly, color = 'm', linewidth = 2)
axs[0].plot(ly2time, mwly, color = 'm', linewidth = 2)
axs[0].plot(ly3time, mwly, color = 'm', linewidth = 2)
axs[1].plot(ly1time, m31ly, color = 'm', linewidth = 2)
axs[1].plot(ly2time, m31ly, color = 'm', linewidth = 2)
axs[1].plot(ly3time, m31ly, color = 'm', linewidth = 2)

fig, axs = plt.subplots(2, figsize = (18, 12))
fig.suptitle('Specific Angular Momentum - Z Direction', fontsize = 35)
axs[0].set_ylabel('MW [kpc-km-$s^{-1}$]', fontsize = 20)
axs[1].set_ylabel('M31 [kpc-km-$s^{-1}$]', fontsize = 20)
axs[0].tick_params(axis='both', labelsize=15)
axs[1].tick_params(axis='both', labelsize=15)
axs[0].plot(t, Lzmw, color = 'g')
axs[1].plot(t, Lzm31, color = 'g')
axs[0].plot(xax_x, xax_y, color = 'k', linewidth = 1)
axs[1].plot(xax_x, xax_y, color = 'k', linewidth = 1)
axs[1].set_xlabel('Time (Gyrs)', fontsize = 15)

### Plotting the merge times
mwlz = [min(Lzmw), max(Lzmw)]
m31lz = [min(Lzm31), max(Lzm31)]
lx1time = [t1, t1]
lx2time = [t2. t2]
lx3time = [t3, t3]
axs[0].plot(lz1time, mwlz, color = 'm', linewidth = 2)
axs[0].plot(lz2time, mwlz, color = 'm', linewidth = 2)
axs[0].plot(lz3time, mwlz, color = 'm', linewidth = 2)
axs[1].plot(lz1time, m31lz, color = 'm', linewidth = 2)
axs[1].plot(lz2time, m31lz, color = 'm', linewidth = 2)
axs[1].plot(lz3time, m31lz, color = 'm', linewidth = 2)'''


# Plotting the overall evolution of both the galxies from start to end. 
fig, axs = plt.subplots(2, figsize = (18, 12))
fig.suptitle("Specific Angular Momentum Evolution\n of Milky Way "
             + 'and M31 Halos', fontsize = 35)
axs[0].plot(t, Lxmw, color = 'b', label='$\hat{x}$')
axs[0].plot(t, Lymw, color = 'r', label='$\hat{y}$')
axs[0].plot(t, Lzmw, color = 'g', label='$\hat{z}$')
axs[0].plot(t, L_mag_mw, color = 'navy', label = 'Mag')
axs[0].plot(xax_x, xax_y, color = 'k', linewidth = 1)
axs[0].set_ylabel('MW [kpc-km-$s^{-1}$]', fontsize = 20)
axs[0].tick_params(axis='both', labelsize=15)
axs[0].text(3.6, 2.25, '3.957 \n Gyrs', fontsize = 14)
axs[0].text(5.55, 2.25, '5.857 \n Gyrs', fontsize = 14)
axs[0].text(6.25, 2.25, '6.486 \n Gyrs', fontsize = 14)
axs[0].text(3.6, 2.25, '3.957 \n Gyrs', fontsize = 14)
axs[0].text(5.55, 2.25, '5.857 \n Gyrs', fontsize = 14)
axs[0].text(6.25, 2.25, '6.486 \n Gyrs', fontsize = 14)
axs[0].legend(fontsize = 15)
axs[1].plot(t, Lxm31, color = 'b', label='$\hat{x}$')
axs[1].plot(t, Lym31, color = 'r', label='$\hat{y}$')
axs[1].plot(t, Lzm31, color = 'g', label= '$\hat{z}$')
axs[1].plot(t, rem_mag, color = 'navy', label = 'Mag')
axs[1].plot(xax_x, xax_y, color = 'k', linewidth = 1)
axs[1].set_xlabel('Time (Gyrs)', fontsize = 15)
axs[1].set_ylabel('M31 [kpc-km-$s^{-1}$]', fontsize = 20)
axs[1].tick_params(axis='both', labelsize=15)
axs[1].text(3.6, 2.25, '3.957 \n Gyrs', fontsize = 14)
axs[1].text(5.55, 2.25, '5.857 \n Gyrs', fontsize = 14)
axs[1].text(6.25, 2.25, '6.486 \n Gyrs', fontsize = 14)
axs[1].text(3.6, 2.25, '3.957 \n Gyrs', fontsize = 14)
axs[1].text(5.55, 2.25, '5.857 \n Gyrs', fontsize = 14)
axs[1].text(6.25, 2.25, '6.486 \n Gyrs', fontsize = 14)

### Plotting the lines for merge times
mwlx = [min(Lxmw), max(Lxmw)]
m31lx = [min(Lxm31), 2]
lx1time = [3.957, 3.957]
lx2time = [5.857, 5.857]
lx3time = [6.486, 6.486]
axs[0].plot(lx1time, mwlx, '--', color = 'm', linewidth = 2)
axs[0].plot(lx2time, mwlx, '--', color = 'm', linewidth = 2)
axs[0].plot(lx3time, mwlx, '--', color = 'm', linewidth = 2)
axs[1].plot(lx1time, m31lx, '--', color = 'm', linewidth = 2)
axs[1].plot(lx2time, m31lx, '--', color = 'm', linewidth = 2)
axs[1].plot(lx3time, m31lx, '--', color = 'm', linewidth = 2)


###############
### Afer Merger
###############
    
# Plotting the pre-merger galaxies (before t1)
# and post merger (after t3) remnant SAM evolution. 
fig, axs = plt.subplots(3, figsize = (18, 12))
fig.suptitle("Evolution of Halo's Specific Angular Momentum\n [Pre-Me"
             +'rger Galaxies (MW & M31) and Post-Merger Remnant]', fontsize = 33)
axs[0].plot(t, Lxmw, color = 'b', label= '$\hat{x}$')
axs[0].plot(t, Lymw, color = 'r', label= '$\hat{y}$')
axs[0].plot(t, Lzmw, color = 'g', label= '$\hat{z}$')
axs[0].plot(t, L_mag_mw, color = 'navy', label = 'Mag')
axs[0].plot(xax_x, xax_y, color = 'k', linewidth = 1)
axs[0].set_xlim(0, t1)
axs[0].set_ylabel('MW [kpc-km-$s^{-1}$]', fontsize = 20)
axs[0].tick_params(axis='both', labelsize=15)
axs[0].legend(fontsize = 15)
axs[1].plot(t, Lxm31, color = 'b', label= '$\hat{x}$')
axs[1].plot(t, Lym31, color = 'r', label= '$\hat{y}$')
axs[1].plot(t, Lzm31, color = 'g', label= '$\hat{z}$')
axs[1].plot(t, L_mag_m31, color = 'navy', label = 'Mag')
axs[1].plot(xax_x, xax_y, color = 'k', linewidth = 1)
axs[1].set_xlim(0, t1)
axs[1].set_ylabel('M31 [kpc-km-$s^{-1}$]', fontsize = 20)
axs[1].tick_params(axis='both', labelsize=15)
axs[2].plot(t, rem_x, color = 'b', label= '$\hat{x}$')
axs[2].plot(t, rem_y, color = 'r', label= '$\hat{y}$')
axs[2].plot(t, rem_z, color = 'g', label= '$\hat{z}}$')
axs[2].plot(t, L_mag_m31, color = 'navy', label = 'Mag')
axs[2].plot(xax_x, xax_y, color = 'k', linewidth = 1)
axs[2].set_ylabel('Remnant [kpc-km-$s^{-1}$]', fontsize = 20)
axs[2].tick_params(axis='both', labelsize=15)
axs[2].set_xlabel('Time (Gyrs)', fontsize = 15)
axs[2].set_xlim(t3, max(t))

# Plotting the post-merger galaxies (after t3)
# and post merger (after t3) remnant SAM evolution. 
fig, axs = plt.subplots(3, figsize = (18, 12))
fig.suptitle("Evolution of Halo's Specific Angular Momentum\n [Post-Merger "
             +'Galaxies (MW & M31) and Post-Merger Remnant]', fontsize = 33)
axs[0].plot(t, Lxmw, color = 'b', label= '$\hat{x}$')
axs[0].plot(t, Lymw, color = 'r', label= '$\hat{y}$')
axs[0].plot(t, Lzmw, color = 'g', label= '$\hat{z}$')
axs[0].plot(t, L_mag_mw, color = 'navy', label = 'Mag')
axs[0].plot(xax_x, xax_y, color = 'k', linewidth = 1)
axs[0].set_ylabel('MW [kpc-km-$s^{-1}$]', fontsize = 20)
axs[0].tick_params(axis='both', labelsize=15)
axs[0].set_xlim(t3, max(t))
axs[0].legend(fontsize = 15)
axs[1].plot(t, Lxm31, color = 'b', label= '$\hat{x}$')
axs[1].plot(t, Lym31, color = 'r', label= '$\hat{y}$')
axs[1].plot(t, Lzm31, color = 'g', label= '$\hat{z}$')
axs[1].plot(t, L_mag_mw, color = 'navy', label = 'Mag')
axs[1].plot(xax_x, xax_y, color = 'k', linewidth = 1)
axs[1].set_ylabel('M31 [kpc-km-$s^{-1}$]', fontsize = 20)
axs[1].tick_params(axis='both', labelsize=15)
axs[1].set_xlim(t3, max(t))
axs[2].plot(t, rem_x, color = 'b', label= '$\hat{x}$')
axs[2].plot(t, rem_y, color = 'r', label= '$\hat{y}$')
axs[2].plot(t, rem_z, color = 'g', label= '$\hat{z}$')
axs[2].plot(t, L_mag_mw, color = 'navy', label = 'Mag')
axs[2].plot(xax_x, xax_y, color = 'k', linewidth = 1)
axs[2].set_ylabel('Remnant [kpc-km-$s^{-1}$]', fontsize = 20)
axs[2].tick_params(axis='both', labelsize=15)
axs[2].set_xlabel('Time (Gyrs)', fontsize = 15)
axs[2].set_xlim(t3, max(t))

##################################
### Plotting Magnitudes Separately
##################################

# Plotting the pre-merger galaxies (before t1)
# and post merger (after t3) remnant SAM evolution. 
fig, axs = plt.subplots(3, figsize = (18, 12))
fig.suptitle("Evolution of Halo's Specific Angular Momentum Magnitude\n [Pre-Merger "
             +'Galaxies (MW & M31) and Post-Merger Remnant]', fontsize = 33)
axs[0].set_ylabel('MW [kpc-km-$s^{-1}$]', fontsize = 20)
axs[0].tick_params(axis='both', labelsize=15)
axs[0].plot(t, L_mag_mw, color = 'navy')
axs[0].set_xlim(0, t1)
axs[1].set_ylabel('M31 [kpc-km-$s^{-1}$]', fontsize = 20)
axs[1].plot(t, L_mag_m31, color = 'navy')
axs[1].tick_params(axis='both', labelsize=15)
axs[1].set_xlim(0, t1)
axs[2].tick_params(axis='both', labelsize=15)
axs[2].set_ylabel('Remnant [kpc-km-$s^{-1}$]', fontsize = 20)
axs[2].plot(t, rem_mag, color = 'navy')
axs[2].set_xlabel('Time (Gyrs)', fontsize = 15)
axs[2].set_xlim(t3, max(t))

# Plotting the post-merger galaxies (after t3)
# and post merger (after t3) remnant SAM evolution. 
fig, axs = plt.subplots(3, figsize = (18, 12))
fig.suptitle("Evolution of Halo's Specific Angular Momentum Magnitude\n [Post-Merger "
             +'Galaxies (MW & M31) and Post-Merger Remnant]', fontsize = 33)
axs[0].plot(t, L_mag_mw, color = 'navy')
axs[0].set_ylabel('MW [kpc-km-$s^{-1}$]', fontsize = 20)
axs[0].tick_params(axis='both', labelsize=15)
axs[0].set_xlim(t3, max(t))
axs[1].plot(t, L_mag_m31, color = 'navy')
axs[1].set_ylabel('M31 [kpc-km-$s^{-1}$]', fontsize = 20)
axs[1].tick_params(axis='both', labelsize=15)
axs[1].set_xlim(t3, max(t))
axs[2].plot(t, rem_mag, color = 'navy')
axs[2].tick_params(axis='both', labelsize=15)
axs[2].set_ylabel('Remnant [kpc-km-$s^{-1}$]', fontsize = 20)
axs[2].set_xlabel('Time (Gyrs)', fontsize = 15)
axs[2].set_xlim(t3, max(t))


