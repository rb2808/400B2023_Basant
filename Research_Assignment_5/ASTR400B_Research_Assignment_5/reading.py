
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt


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
    
    # read and store time
    line1 = file.readline()
    
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


d = Read("MW_Angular_Momentum.txt")
Lxmw = []
Lymw = []
Lzmw = []
t = []

for i in range(len(d)):
    t.append(d[i][0])
    Lxmw.append(d[i][1])
    Lymw.append(d[i][2])
    Lzmw.append(d[i][3])

d1 = Read("M33_Angular_Momentum.txt")
Lxm33 = []
Lym33 = []
Lzm33 = []

for i in range(len(d1)):
    Lxm33.append(d1[i][1])
    Lym33.append(d1[i][2])
    Lzm33.append(d1[i][3])

d2 = Read("M31_Angular_Momentum.txt")
Lxm31 = []
Lym31 = []
Lzm31 = []

for i in range(len(d2)):
    Lxm31.append(d2[i][1])
    Lym31.append(d2[i][2])
    Lzm31.append(d2[i][3])

    
fig, axs = plt.subplots(3, figsize = (18, 12))
fig.suptitle('Specific Angular Momentum - X Direction', fontsize = 35)
axs[0].set_ylabel('MW [kpc-km-$s^{-1}$]', fontsize = 20)
axs[1].set_ylabel('M31 [kpc-km-$s^{-1}$]', fontsize = 20)
axs[2].set_ylabel('M33 [kpc-km-$s^{-1}$]', fontsize = 20)

axs[2].set_xlabel('Time (Gyrs)', fontsize = 15)
axs[0].tick_params(axis='both', labelsize=15)
axs[1].tick_params(axis='both', labelsize=15)
axs[2].tick_params(axis='both', labelsize=15)
axs[0].plot(t, Lxmw, color = 'b')
axs[1].plot(t, Lxm31, color = 'b')
axs[2].plot(t, Lxm33, color = 'b')

### Plotting the merge times
mwlx = [min(Lxmw), max(Lxmw)]
m31lx = [min(Lxm31), max(Lxm31)]
m33lx = [min(Lxm33), max(Lxm33)]
lx1time = [3.957, 3.957]
lx2time = [5.857, 5.857]
lx3time = [6.486, 6.486]
axs[0].plot(lx1time, mwlx, color = 'k', linewidth = 3)
axs[0].plot(lx2time, mwlx, color = 'k', linewidth = 3)
axs[0].plot(lx3time, mwlx, color = 'k', linewidth = 3)
axs[1].plot(lx1time, m31lx, color = 'k', linewidth = 3)
axs[1].plot(lx2time, m31lx, color = 'k', linewidth = 3)
axs[1].plot(lx3time, m31lx, color = 'k', linewidth = 3)
axs[2].plot(lx1time, m33lx, color = 'k', linewidth = 3)
axs[2].plot(lx2time, m33lx, color = 'k', linewidth = 3)
axs[2].plot(lx3time, m33lx, color = 'k', linewidth = 3)


fig, axs = plt.subplots(3, figsize = (18, 12))
fig.suptitle('Specific Angular Momentum - Y Direction', fontsize = 35)
axs[0].set_ylabel('MW [kpc-km-$s^{-1}$]', fontsize = 20)
axs[1].set_ylabel('M31 [kpc-km-$s^{-1}$]', fontsize = 20)
axs[2].set_ylabel('M33 [kpc-km-$s^{-1}$]', fontsize = 20)
axs[2].set_xlabel('Time (Gyrs)', fontsize = 15)
axs[0].tick_params(axis='both', labelsize=15)
axs[1].tick_params(axis='both', labelsize=15)
axs[2].tick_params(axis='both', labelsize=15)
axs[0].plot(t, Lymw, color = 'r')
axs[1].plot(t, Lym31, color = 'r')
axs[2].plot(t, Lym33, color = 'r')

### Plotting the merge times
mwly = [min(Lymw), max(Lymw)]
m31ly = [min(Lym31), max(Lym31)]
m33ly = [min(Lym33), max(Lym33)]
ly1time = [3.957, 3.957]
ly2time = [5.857, 5.857]
ly3time = [6.486, 6.486]
axs[0].plot(ly1time, mwly, color = 'k', linewidth = 3)
axs[0].plot(ly2time, mwly, color = 'k', linewidth = 3)
axs[0].plot(ly3time, mwly, color = 'k', linewidth = 3)
axs[1].plot(ly1time, m31ly, color = 'k', linewidth = 3)
axs[1].plot(ly2time, m31ly, color = 'k', linewidth = 3)
axs[1].plot(ly3time, m31ly, color = 'k', linewidth = 3)
axs[2].plot(ly1time, m33ly, color = 'k', linewidth = 3)
axs[2].plot(ly2time, m33ly, color = 'k', linewidth = 3)
axs[2].plot(ly3time, m33ly, color = 'k', linewidth = 3)



fig, axs = plt.subplots(3, figsize = (18, 12))
fig.suptitle('Specific Angular Momentum - Z Direction', fontsize = 35)
axs[0].set_ylabel('MW [kpc-km-$s^{-1}$]', fontsize = 20)
axs[1].set_ylabel('M31 [kpc-km-$s^{-1}$]', fontsize = 20)
axs[2].set_ylabel('M33 [kpc-km-$s^{-1}$]', fontsize = 20)
axs[2].set_xlabel('Time (Gyrs)', fontsize = 15)
axs[0].tick_params(axis='both', labelsize=15)
axs[1].tick_params(axis='both', labelsize=15)
axs[2].tick_params(axis='both', labelsize=15)
axs[0].plot(t, Lzmw, color = 'g')
axs[1].plot(t, Lzm31, color = 'g')
axs[2].plot(t, Lzm33, color = 'g')

### Plotting the merge times
mwlz = [min(Lzmw), max(Lzmw)]
m31lz = [min(Lzm31), max(Lzm31)]
m33lz = [min(Lzm33), max(Lzm33)]
lz1time = [3.957, 3.957]
lz2time = [5.857, 5.857]
lz3time = [6.486, 6.486]
axs[0].plot(lz1time, mwlz, color = 'k', linewidth = 3)
axs[0].plot(lz2time, mwlz, color = 'k', linewidth = 3)
axs[0].plot(lz3time, mwlz, color = 'k', linewidth = 3)
axs[1].plot(lz1time, m31lz, color = 'k', linewidth = 3)
axs[1].plot(lz2time, m31lz, color = 'k', linewidth = 3)
axs[1].plot(lz3time, m31lz, color = 'k', linewidth = 3)
axs[2].plot(lz1time, m33lz, color = 'k', linewidth = 3)
axs[2].plot(lz2time, m33lz, color = 'k', linewidth = 3)
axs[2].plot(lz3time, m33lz, color = 'k', linewidth = 3)
