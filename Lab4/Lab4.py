
#In Class Lab 4 Template
# G Besla ASTR 400B

# Load Modules
import numpy as np
import astropy.units as u

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib


# The Figure illustrates the color magnitude diagram (CMD) for the Carina Dwarf along with the interpreted 
# star formation history from isochrone fitting to the CMD.
# The image is from Tolstoy+2009 ARA&A 47 review paper about dwarf galaxies
# 
# ![Iso](./Lab4_Isochrones.png)
# 

# # This Lab:
# 
# Modify the template file of your choice to plot isochrones that correspond to the inferred star formation episodes (right panel of Figure 1) to recreate the dominant features of the CMD of Carina (left panel of Figure 1). 



# Some Notes about the Isochrone Data
# DATA From   http://stellar.dartmouth.edu/models/isolf_new.html
# files have been modified from download.  ( M/Mo --> M;   Log L/Lo --> L)
# removed #'s from all lines except column heading
# NOTE SETTINGS USED:  Y = 0.245 default   [Fe/H] = -2.0  alpha/Fe = -0.2
# These could all be changed and it would generate a different isochrone




# Filename for data with Isochrone fit for 1 Gyr
# These files are located in the folder IsochroneData
filename1 = "./IsochroneData/Isochrone1.txt"
filename2 = "./IsochroneData/Isochrone2.txt"
filename3 = "./IsochroneData/Isochrone3.txt"
filename4 = "./IsochroneData/Isochrone4.txt"
filename5 = "./IsochroneData/Isochrone5.txt"
filename6 = "./IsochroneData/Isochrone6.txt"
filename7 = "./IsochroneData/Isochrone7.txt"


# READ IN DATA
# "dtype=None" means line is split using white spaces
# "skip_header=8"  skipping the first 8 lines 
# the flag "names=True" creates arrays to store the date
#       with the column headers given in line 8 

# Read in data for an isochrone corresponding to 1 Gyr
data1 = np.genfromtxt(filename1,dtype=None,names=True,skip_header=8)
data2 = np.genfromtxt(filename2,dtype=None,names=True,skip_header=8)
data3 = np.genfromtxt(filename3,dtype=None,names=True,skip_header=8)
data4 = np.genfromtxt(filename4,dtype=None,names=True,skip_header=8)
data5 = np.genfromtxt(filename5,dtype=None,names=True,skip_header=8)
data6 = np.genfromtxt(filename6,dtype=None,names=True,skip_header=8)
data7 = np.genfromtxt(filename7,dtype=None,names=True,skip_header=8)
# Plot Isochrones 
# For Carina

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# Plot Isochrones

# Isochrone for 1 Gyr
# Plotting Color vs. Difference in Color 
plt.plot(data1['B']-data1['R'], data1['R'], color='red', linewidth=3, label='1 Gyr')
plt.plot(data2['B']-data2['R'], data2['R'], color='orange', linewidth=3, label='2 Gyr')
plt.plot(data3['B']-data3['R'], data3['R'], color='yellow', linewidth=3, label='3 Gyr')
plt.plot(data4['B']-data4['R'], data4['R'], color='green', linewidth=3, label='4 Gyr')
plt.plot(data5['B']-data5['R'], data5['R'], color='blue', linewidth=3, label='5 Gyr')
plt.plot(data6['B']-data6['R'], data6['R'], color='indigo', linewidth=3, label='6 Gyr')
plt.plot(data7['B']-data7['R'], data7['R'], color='violet', linewidth=3, label='7 Gyr')



# Add axis labels
plt.xlabel('B-R', fontsize=22)
plt.ylabel('M$_R$', fontsize=22)

#set axis limits
plt.xlim(-0.5,2)
plt.ylim(5,-2.5)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

#add figure text
plt.figtext(0.6, 0.15, 'Isochrone Carina', fontsize=22)

plt.savefig('IsochroneLab4.png')


# Question-2: It does not seem like there are younger ages. The peak star formation
#     occurs at the same the time period with only little variation. 
# Question-3: At this redshift, the burst star formation maybe due to high metallicities. 




