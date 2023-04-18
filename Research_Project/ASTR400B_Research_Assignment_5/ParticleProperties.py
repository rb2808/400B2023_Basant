# GB : ASTR 400B Solutions to PSet 2 
# ParticleInfo

# load modules
import numpy as np
import astropy.units as u
from ReadFile import Read




def ParticleInfo(PType, filename):
  
    """ Function to return properties of a particle of a given type
    
    Input: 
        PType: int
            particle type, e.g. Halo: 1, Disk: 2, Bulge: 3
        PNum: int 
            particle number, e.g. 100)
        filename: str
            e.g. "MW_000.txt")
        
    Output: 
        R3D: astropy quantity
            Magnitude of 3D Pos vector (kpc)
        V3D: astropy quantity
            Magnitude of 3D Velocity vector (km/s)
        Mass: astropy quantity
            Mass of the Particle (Msun)
    """


    # read in the file 
    time, total, data = Read(filename)

    
    #create an array to store indexes of particles of desired Ptype
    index = np.where(data['type'] == PType)

    # create new arrays with the m, x, y, z, 
    # vx, vy, vz of just the Ptype desired
    # Add units using Astropy
    # Recall Mass was stored in units of Msun/1e10
    mnew = data['m'][index]*1e10*u.Msun
    xnew = data['x'][index]*u.kpc
    ynew = data['y'][index]*u.kpc
    znew = data['z'][index]*u.kpc
    vxnew = data['vx'][index]*u.km/u.s
    vynew = data['vy'][index]*u.km/u.s
    vznew = data['vz'][index]*u.km/u.s
    
    return np.array([mnew, xnew, ynew, znew, vxnew, vynew, vznew]), time