import numpy as np
import matplotlib
import matplotlib.pyplot as plt

AUTOEV = 27.2114
NMTOEV = 1239.84193
EVTOJ = 1.60218e-19
AGTOM = 1e-10
AUTOAG = 5.29177210903e-1
AUTOS = 2.4188843265857e-17
HBAR = 1.054571817e-34  # J*s
FermiEnergy = -0.146 * AUTOEV
BETA = 0.07756 # eV at 300K

def delta_to_Gaussian(center, gaussian_std, x):
    # center, gaussian_std must be scalar. 
    y = 1/(gaussian_std * np.sqrt(2*np.pi)) * np.exp(-(x - center)**2 / (2*gaussian_std**2))
    return y

def nonuniform_broadening_Abs(E_range, sqrtOS_range, lowest_gaussian_std, bulkBandGap):
    x_range = np.arange(min_E, max_E, lowest_gaussian_std/20)
    y_range = np.zeros(len(x_range))
    for i in range(len(E_range)):
        this_energy = E_range[i]
        this_sqrtOS = sqrtOS_range[i]
        this_gaussian_std = lowest_gaussian_std * (E_range[i] - bulkBandGap) / (E_range[0] - bulkBandGap)  
        y_range += delta_to_Gaussian(this_energy, this_gaussian_std, x_range) * this_sqrtOS * this_sqrtOS * this_energy
    return (x_range, y_range)

def nonuniform_broadening_Em(E_range, sqrtOS_range, lowest_gaussian_std, bulkBandGap):   # at 300K
    x_range = np.arange(min_E, max_E, lowest_gaussian_std/20)
    y_range = np.zeros(len(x_range))
    for i in range(len(E_range)):
        this_energy = E_range[i]
        this_sqrtOS = sqrtOS_range[i]
        this_gaussian_std = lowest_gaussian_std * (E_range[i] - bulkBandGap) / (E_range[0] - bulkBandGap)  
        BoltzmannFactor = np.exp(-(this_energy-E_range[0]) / (0.0257))   # kbT at 298K is 25.7 meV
        y_range += delta_to_Gaussian(this_energy, this_gaussian_std, x_range) * this_sqrtOS * this_sqrtOS * this_energy * BoltzmannFactor
    return (x_range, y_range)

def firstPeak(spectrum):
    if not isinstance(spectrum, np.ndarray) or spectrum.ndim != 2:
        raise ValueError("Input must be a 2D numpy array.")
    
    intensity = spectrum[:, 1]
    decreasing_index = np.where(np.diff(intensity) < 0)[0]

    if decreasing_index.size > 0:
        index = decreasing_index[0]
        energy, intensity = spectrum[index]
        return energy, intensity
    else:
        return None
