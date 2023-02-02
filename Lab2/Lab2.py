

# ASTR 400 B 
# In Class Lab 2

# Import Modules 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import quad # For integration
# Documentation and examples for quad : 
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html
# https://www.tutorialspoint.com/scipy/scipy_integrate.htm


# ## Part A:  Schechter Fxn
# 
# The galaxy luminosity function in the nearby universe is well described by a Schechter Function:
# 
# \begin{equation}
# \Phi(M)dM = ( 0.4 \, ln10 ) \, \phi_\ast \, 10^{0.4(M_\ast - M)(\alpha +1)} e^{-10^{0.4(M_\ast - M)}} dM
# \end{equation}
# 
# With the following parameters from Smith+2009 for Field Galaxies in SDSS at z$\sim$0.1 in the Kband:
# 
# 
#  $\phi_\ast$ =1.66 $  \times 10^{-2}$  $h^3$ Mpc$^{-3}$
# 
#  $\alpha$ =  -0.81 
# 
# 
#   M$_\ast$ =  M$_k^\ast$= -23.19  - 5*log($h$)
#   
#  $h$ = the Hubble constant in units of 100 km/s/Mpc . At z=0 this is 0.7. 
# But we are going to se $h$=1 here. Units will then be in "comoving" coordinates.
#   
#   This function is defined for you below:



def schechter_M(m,phi_star=0.0166,m_star=-23.19,alpha=-0.81): # always define function with small initial
    """ Function that computes the Schechter Luminosity Function for a given magnitude, 
    assuming default parameters for field galaxies in SDSS at z~0.1 in the Kband (Smith+2009)
    
    Inputs
        m : an array of floats
            an array of Kband magnitudes  (assumes -5*log(h) implicitly)
        phi_star:  float
            normalization of Schechter fxn (h^3 Mpc^-3)
        m_star:  float 
            knee of the Schechter fxn (K-band magnitude, assumes -5*log(h) implicitly)
        alpha:  float
            faint end slope of the Schechter fxn
    
    Output:
        schechterM: float
            number density of galaxies (comoving units) at the given magnitude m - 5*log(h)
            

    """

# You should divide up long functions instead of writing them out as one long set
    a = 0.4*np.log(10)*phi_star # Grouping all constants together
    b = 10**(0.4*(m_star-m)*(alpha+1.0)) # The Power Law, controlling the faint end slope
    c = np.exp(-10**(0.4*(m_star-m))) # The Exponential controlling the high mass end behavior
    
    schechterM = a*b*c # schechter function for the given magnitude
# i.e. don't do the below
#    return 0.4*np.log(10)*phistar*10**(0.4*(Mstar - M)*(alpha +1.0))*np.exp(-10**(0.4*(Mstar - M)))

    return schechterM


# # Q1 
# 
# Utilizing the defined function, plot the Schechter Function using the above parameter values over a magnitude range of -17 to -26. 
# Try to reproduce the black solid line in Smith+2009 MNRAS 397,868 [UKIDSS Survey] Figure below.
# 
# 
# ![Smith](./Smith09.png)



# # Q2 
# 
# Galaxies in the Virgo Cluster have different parameters, like $\alpha$=-1.35  (Ferrarese+2016 ApJ 824).
# 
# Overplot the Schechter Function with this new value of $\alpha$.  
# 
# Try a smaller value of $\alpha = -0.6$.
# 
# How does the function change?  What does this mean? 
# 




# Create an array to store Kband Magnitudes from -26 to -17


mK = np.arange(-26,-17,0.1)
print(mK[2])


# Plot the Schechter Function

fig = plt.figure(figsize=(10,10))  # sets the scale of the figure
ax = plt.subplot(111) 

# Plot the default values (y axis log)
# ADD HERE

ax.semilogy(mK,schechter_M(mK),color ='blue', linewidth=5, label ='Smith+09')                          # as only y is in log not x
# Q2 solutions: change alpha
# ADD HERE

ax.semilogy(mK,schechter_M(mK,alpha=-1.35), color= 'red', linestyle=':',linewidth=3 ,
            label =r'high $\alpha$') # chnaging default alpha value
# by 'r' we are telling them to enter latex mode

ax.semilogy(mK,schechter_M(mK,alpha=-0.6), color= 'green', linestyle='--',linewidth=3 ,
            label =r'low $\alpha$')

# Add labels
plt.xlabel(r'M$_k$ + 5Log($h$)', fontsize=22)
plt.ylabel(r'$\Phi$ (Mpc$^{-3}h^3$/mag)', fontsize=22)

#set axis limits
plt.xlim(-17,-26)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')

# Save to a file
plt.savefig('Schechter_M.png')




# # Q3
# 
# Build a function to compute the Schechter Function in terms of luminosity.
# 
# Use `quad` to numerically integrate the function to compute the fraction of the luminosity that lies above L* in the following three cases:  
# 
# $\alpha$=-0.7 (default), $\alpha$=-0.6, $\alpha$=1.85. 
# 
# 
# Schecheter Function $\Phi(L) = \frac{n_\ast}{L_\ast} (\frac{L}{L_\ast})  ^{\alpha}  e^{-L/L_\ast}$
# 
# $n_\ast$ = 0.008  $h^3$ Mpc$^{-3}$
# 
# $L_\star = 1.4 \times 10^{10} L_\odot$

def schechter_L(lum,n_star=8e-3,l_star = 1.4e10,alpha= -0.7):
    """ this is the function that compute the scheter luminosity function
    for a given luminosity
    defuatls are from Sparke & galaghar
    
    Inputs :
        lum = array of floats 
            Array of Luminosities (Lsun)
            
        n_star : float
        Normalization of the schelter function fxn ( h^3 Mpc ^-3)
        
        l_star : float
          Characteristic velocity luminosity ( knee of the schedter function) in units of L_sun
        
        alpha : float
        Faint end slope
        
        Outputs :
            
            scheter_l:float
            number density of galaxies for a given lunimosyty ( h^3 * Mpc^-3/Lsun)
    
    """
    
    # Break down equation into parts
    
    a = (lum/l_star)**alpha # faint end
    b= np.exp(-lum/l_star) # bright end
    c= n_star/l_star        # constants

    schechter_l = a*b*c    
    return schechter_l
    








# Understanding lambda functions
# Short cut -- defines and evaluates a function in one line ! 

# lambda says that a function follows, where the variables are a and b, and the function to be evaluated is a*b
x = lambda a, b : a * b
print(x(5, 6))  





# Example Usage of quad and lambda

print(quad(np.sin, 0, np.pi)) # integrating sin function from 0 to pi.


f = lambda x: np.sin(x)
print(quad(f, 0, np.pi)) # as f is sin x , so bascially we made it more complicated

# first element quad is the integral, second element is the error


def ex(x):
    return np.sin(x) 

print(quad(lambda x: ex(x), 0, np.pi))


# What fraction of the integrated luminoisty density lies above L*
# alpha= -0.7

l_upper = quad(lambda L: L*schechter_L(L), 1.4e10 , 1e14)
print(l_upper)

l_total = quad(lambda L: L*schechter_L(L) , 0.1 , 1e14)

print("Flux ratio (>L*)/LTotal" , np.round(l_upper[0]/l_total[0],3))

l_upper1 = quad(lambda L: L*schechter_L(L, alpha = -1.0), 1.4e10 , 1e14)
print(l_upper)

l_total1 = quad(lambda L: L*schechter_L(L , alpha = -1.0) , 0.1 , 1e14)

print("Flux ratio (>L*)/LTotal" , np.round(l_upper1[0]/l_total1[0],3))

l_upper2 = quad(lambda L: L*schechter_L(L, alpha = -1.85), 1.4e10 , 1e14)
print(l_upper)

l_total2 = quad(lambda L: L*schechter_L(L , alpha = -1.85) , 0.1 , 1e14)

print("Flux ratio (>L*)/LTotal" , np.round(l_upper2[0]/l_total2[0],3))





# ## Part B: IMF 
# 
# Create a function called `Salpeter` that defines the Salpeter IMF: 
# 
# \begin{equation}
# \xi(M) = \xi_0 (M/M_\odot)^{-\alpha}
# \end{equation}
# 
# $\alpha = 2.35$
# The function should take as input an array of stellar masses, M. 
# You will need to determine the normalization, $\xi_0$, by integrating this equation over mass from 0.1 to 120 M$_\odot$
# and setting the value to 1.  The function should then return $\xi(M)$, which will now represent the fractional number of stars. 
# 
# Integration:
# 
# `quad(lambda x:  fxn(x),xmin,xmax)`
# 
# quad returns an array with 2 values. you want the first value. 
# Note I've used a "lambda" expression.   Python's lambda expressions allow a 
# function to be created and passed around all in one line of code


def salpeter(m, m_min = 0.1, m_max = 120, alpha = 2.35):
    
    '''Function that defines the SAlpeter IMF. The function is nirmalied such that
    it returns the fraction of stars expected  (assuming stars range in mass
                                                              from m_min to m_max)
    
    Inputs: 
        m: array of floats
            Array of stellar masses (Msun)
            
        m_min: float
            minimum mass (Msun)
            
        m_max: float
            maximal mass (Msun)
            
        alpha: float
            power law for the salpeter function
            
    Outputs: 
        normalized_salpeter: a float
            normalized fraction of stars at a given mass.'''
            
            
            
    # determine the magnitude of the integral
    to_normalize = quad(lambda m: m**(-alpha), m_min, m_max)

    # Determine the normalization factor

    norm = 1 / to_normalize[0]    
    
    # return the normalized the Salpeter IMF
    norm_salpeter = norm*m**(-alpha)
    
    return norm_salpeter




# ## Q1: 
# Double Check: if you integrate your function from 0.1 to 120 you should return 1.0 
# 

Test = quad(lambda m: salpeter(m), 0.1, 120)

print(np.round(Test[0]))





# ## Q2: 
# Integrate your normalized function to compute the fraction of stars with stellar masses greater than the sun and less 
# than 120 M$_\odot$.


frac = quad(lambda m: salpeter(m), 1.0, 120)

print(np.round(frac[0], 3))



# ## Q3:
# 
# How might you modify the above to return the fraction of MASS ? instead of fraction of the total numbers of stars.

print(5000 * np.round(frac[0], 3))




def salpeter_mass(m, m_min = 0.1, m_max = 120, alpha = 2.35):
    
    '''Function that defines the SAlpeter IMF. The function is nirmalied such that
    it returns the fraction of mass (assuming stars range in mass from m_min to m_max)
                                                              
    
    Inputs: 
        m: array of floats
            Array of stellar masses (Msun)
            
        m_min: float
            minimum mass (Msun)
            
        m_max: float
            maximal mass (Msun)
            
        alpha: float
            power law for the salpeter function
            
    Outputs: 
        normalized_salpeter: a float
            normalized fraction of fraction of mass at a given m range.'''
            
            
            
    # determine the magnitude of the integral
    to_normalize = quad(lambda m: m*m**(-alpha), m_min, m_max)

    # Determine the normalization factor

    norm = 1 / to_normalize[0]    
    
    # return the normalized the Salpeter IMF
    norm_salpeter = norm*m*m**(-alpha)
    
    return norm_salpeter


# Determine the fraction of mass in stars that are more massive than the sun

frac2 = quad(lambda m: salpeter_mass(m), 1, 120)

print(np.round(frac2[0], 3))


# 500 Msun cluster

print(5000*np.round(frac2[0], 3))

# 100 Msun cluster
print(100*np.round(frac2[0], 3))




