#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 16:40:04 2020

@author: pablo
"""


from Basics import Basics
from SSPs.SSPs import Pipe3Dssp as SSP
import numpy as np

class RecoverSFH(SSP, Basics):
    
    def __init__(self, **kwargs):
    
        Basics.__init__(self, **kwargs)
        SSP.__init__(self)
        
# =============================================================================
#   SFH FILES        
# =============================================================================
    def get_SFHweights(self, mode='individual'):
        """ 
        This method provides the relative weights of each single stellar 
        population within the considered SSP-library to the flux intensity at
        5635AA for each individual spaxel within the original MaNGA cube.
        In addition, the relative weights for each considered individual SSP 
        we include the corresponding weights for the 39 ages (averaged by 
        metallicity) and the 4 metallicities (averaged by age).
        """
        
        modes = {'individual':np.arange(0,156), 'age':np.arange(156, 195),
                 'metallicity':np.arange(195, 199)}
        
        return self.sfh_data.data[modes[mode], :, :]   

# =============================================================================
#   COMPUTATION METHODS
# =============================================================================
    def compute_ssp_masses(self, mode='individual', wl_ref=5635):
        """
        There are two possible modes to compute the masses:
            - 'individual' mode uses the 156 ssps and their corresponding 
               weights.
            - 'age' mode uses metallicity-averaged weights (39).
        """
        
        self.flux = self.get_median_flux(unit='Lsun')
        self.redshift = self.get_redshift()
        
        self.angular_diameter_distance = self.compute_angular_diameter_distance()
        self.spaxel_area = self.get_spaxel_area(unit='Kpc')
        
        self.luminosity = 4* np.pi* self.compute_distance(unit='cm')**2 * self.flux
        
                                
        ssp_mass_to_lum = self.ssp_present_mass_lum_ratio(mode)
        ssp_alive_stellar_mass = self.ssp_alive_stellar_mass(mode)      
        ssp_weights = self.get_SFHweights()        
                
        self.ssp_masses = (
               self.luminosity[np.newaxis, :, :] * self.median_flux_units
                    *ssp_weights *ssp_mass_to_lum[:, np.newaxis, np.newaxis]
                    *ssp_alive_stellar_mass[:, np.newaxis, np.newaxis]  
                            )

    def compute_luminosity_weighted_age(self):
        ages = self.ssp_ages(mode='individual')
        ssp_weights = self.get_SFHweights()        
        lum_weighted_age = np.sum(np.log10(ages[:, np.newaxis, np.newaxis]) * ssp_weights, axis=0)
        return lum_weighted_age
    
    
    def compute_mass_weighted_age(self):
        ages = self.ssp_ages(mode='individual')
        ssp_mass_to_lum = self.ssp_present_mass_lum_ratio(mode='individual')
        ssp_weights = self.get_SFHweights()        
        
        all_weights = np.sum(ssp_weights, axis=0)
        mask = all_weights==0               
        
        mass_weighted_age = np.sum(
        np.log10(ages[:, np.newaxis, np.newaxis])*ssp_weights*ssp_mass_to_lum[:, np.newaxis, np.newaxis], 
        axis=0)/np.sum(ssp_weights*ssp_mass_to_lum[:, np.newaxis, np.newaxis], axis=0)
        mass_weighted_age[mask] = 0
        
        return mass_weighted_age
                
    def compute_SFH(self, mode='individual', today=14e9):        
        
        try: 
            self.ssp_masses                    
        except:
            self.compute_ssp_masses(mode)
        
        ages = self.ssp_ages(mode='individual')            
        
        sort_ages = np.argsort(ages)
        
        sorted_ssp_masses = self.ssp_masses[sort_ages, :, :]
        
        sorted_ssp_masses = sorted_ssp_masses.reshape(39, 4, 
                                      self.ssp_masses.shape[1], 
                                      self.ssp_masses.shape[2])
        
        # Sum over metallicities
        self.total_ssp_mass = np.sum(sorted_ssp_masses, axis=1)
                            
        ages =  np.unique(ages[sort_ages])[::-1] #from old to young        
        self.total_ssp_mass = self.total_ssp_mass[::-1, :, :]
        
                
        self.stellar_mass_history = np.cumsum(
            self.total_ssp_mass, axis=0) 
        
        self.time_bins = ages[0] - ages 
        
        self.time = (self.time_bins[1:]+self.time_bins[:-1])/2
        
        self.star_formation_history = (
            self.stellar_mass_history[1:]-self.stellar_mass_history[:-1])/(
            self.time_bins[1:, np.newaxis, np.newaxis]-\
            self.time_bins[:-1, np.newaxis, np.newaxis])
                
        self.stellar_mass_history = (
                                     self.stellar_mass_history[1:, :, :]+\
                                     self.stellar_mass_history[0:-1, : , :]
                                     )/2
        
        self.specific_star_formation_history = self.star_formation_history/(
                                                    self.stellar_mass_history)        
            
    def compute_binned_SFH(self, mode='individual', today=14e9):        
        self.get_continuum_segmentation()
        
        try: 
            self.ssp_masses                    
        except:
            self.compute_ssp_masses(mode)
        
        ages = self.ssp_ages(mode='individual')            
        
        sort_ages = np.argsort(ages)
        
        sorted_ssp_masses = self.ssp_masses[sort_ages, :, :]
        
        sorted_ssp_masses = sorted_ssp_masses.reshape(39, 4, 
                                      self.ssp_masses.shape[1], 
                                      self.ssp_masses.shape[2])
        
        self.total_ssp_mass = np.sum(sorted_ssp_masses, axis=1)
            
        ages =  np.unique(ages[sort_ages])[::-1] #from old to young
        self.time_bins = ages[0] - ages 
        
        ages =  ages[::-1] #from old to young
                        
        self.total_ssp_mass = self.total_ssp_mass[::-1, :, :]
        
        area = self.get_spaxel_area()
        
        self.stellar_mass_history = np.cumsum(
            self.total_ssp_mass, axis=0) 
        
        self.binned_mass_history = np.zeros((self.segmentation_cube.shape[0],
                                                 ages.size))
        self.bin_area = np.zeros(self.segmentation_cube.shape[0])
        
        for i in range(self.segmentation_cube.shape[0]):
            n_spaxels = len( self.segmentation_cube[i, :, :][
                self.segmentation_cube[i, :, :] == True])
            
            self.binned_mass_history[i, :] = np.sum(
                self.stellar_mass_history[:, self.segmentation_cube[i, :, :]], 
                axis=(1))
            self.bin_area[i] = area*n_spaxels
            
            
        self.binned_sfh = np.diff(self.binned_mass_history)/np.diff(self.time_bins)
            
    def integrated_mass_history(self):
        return np.sum(self.stellar_mass_history, axis=(1,2))
    
    def integrated_star_formation_history(self):
        return np.sum(self.star_formation_history, axis=(1,2))
    
    def mass_to_density(self, mass_array):
        return mass_array/self.spaxel_area/4/np.pi        
    
# ----------------------------------------------------------------------------            
if __name__=='__main__':
    from matplotlib import pyplot as plt
    galaxy = RecoverSFH(survey='CALIFA', galaxy_id='NGC0001')    
    galaxy.compute_binned_SFH()
    MH = galaxy.binned_mass_history
    CMF = MH/MH[:, -1, np.newaxis]
    lookbacktime = (max(galaxy.time_bins) - galaxy.time_bins)/1e9
    
    SFR_fraction = (lookbacktime[0]/lookbacktime[np.newaxis, :])*(1-CMF)
    
    plt.figure()    
    plt.plot(lookbacktime, np.nanmean(SFR_fraction, axis=0), '.-')
    plt.fill_between(lookbacktime, np.nanmean(SFR_fraction, axis=0)- np.nanstd(SFR_fraction, axis=0),
                     np.nanmean(SFR_fraction, axis=0)+ np.nanstd(SFR_fraction, axis=0),
                     alpha=0.2)
    plt.xscale('log')
    
    # maps = Maps(survey='MANGA', galaxy_id='manga-9510-12705')
    # maps.plot_segmentation_regions()
    # plt.figure()
    # plt.imshow(maps.get_velocity_map(), vmin=-400, vmax=400, cmap='jet')
    # plt.colorbar()
    # plt.figure()
    # plt.plot(maps.get_surface_density().flatten(), maps.get_velocity_dispersion_map().flatten(), 'o')
    # plt.xlim(5,9)
    # maps = Maps(survey='MANGA', galaxy_id='manga-9510-12703')
    # plt.figure()
    # plt.imshow(maps.get_velocity_map(), vmin=-400, vmax=400, cmap='jet')
    # plt.colorbar()


# =============================================================================
# Mr. Krtxo       
# =============================================================================
