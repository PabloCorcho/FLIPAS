#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:32:46 2020

@author: pablo
"""


from Builder import Builder
import numpy as np

class Global_properties(Builder):
    def __init__(self, **kwargs):
        Builder.__init__(self, **kwargs)
        self.galaxy = self.build_data_bundle() 


    def get_stellar_mass(self):
        
        return self.galaxy.catalog[
                                    self.galaxy.own_keywords('stellar_mass')
                                   ][self.galaxy.catalog_entry]
    
    def get_redshift(self):
        
        return self.galaxy.catalog[
                                    self.galaxy.own_keywords('redshift')
                                   ][self.galaxy.catalog_entry]

    def compute_distance(self,  unit='pc'):
        units = {'cm':3.086e+24,'pc':1e6, 'Kpc':1e3, 'Mpc':1, 'Gpc':1e-3}
        speed_of_light = 3e5 # km/s
        Hubble_constant = 70        # km/s/Mpc
        distance = self.get_redshift()*speed_of_light/Hubble_constant
        return distance*units[unit]

    def compute_angular_diameter_distance(self, unit='Kpc'):                
        return self.compute_distance(unit)/(1+self.get_redshift())**2
        
    def get_spaxel_size(self, unit='Kpc'):
        self.spaxel_diameter = self.compute_angular_diameter_distance(unit)*\
                                            self.galaxy.fiber_angle_diameter
        return self.spaxel_diameter
    
    def get_spaxel_area(self, unit='Kpc'):
        self.spaxel_area=np.pi * self.get_spaxel_size(unit)**2
        return self.spaxel_area
        
class Maps(Builder):
    def __init__(self, **kwargs):
        Builder.__init__(self, **kwargs)
        self.galaxy = self.build_data_bundle() 
        
        self.ssp_data = self.galaxy.file_SSP()        
        self.sfh_data = self.galaxy.file_SFH()        
        
        
# =============================================================================
#   SSP FILES      
# =============================================================================
    
    def get_Vcontinuum_flux(self, unit='erg'):
        units = {'erg':1e-16, 'Lsun':1e-16/3.828e33}        
        self.v_band_units = units[unit]
        return self.ssp_data.data[0]
    
    def get_continuum_segmentation(self):
        self.segmentation = self.ssp_data.data[1]
        
        n_bins = len(np.unique(self.segmentation))
        bins = np.arange(0, n_bins) 
                
        self.segmentation_cube = np.zeros((n_bins, self.segmentation.shape[0], 
                                                   self.segmentation.shape[1]),
                                          dtype=bool)        
        for i in range(n_bins):
            self.segmentation_cube[i, :, :] = self.segmentation == bins[i]
            
    def get_dezonification(self):
        return self.ssp_data.data[2]
    
    def get_median_flux(self, unit='erg'):
        units = {'erg':1e-16, 'Lsun':1e-16/3.828e33}        
        self.median_flux_units = units[unit]
        return self.ssp_data.data[3]
    
    def get_sigma_median_flux(self, unit='erg'):
        units = {'erg':1e-16, 'Lsun':1e-16/3.828e33}        
        self.median_flux_units = units[unit]
        return self.ssp_data.data[4]
    
    
    def get_luminosity_weighted_age(self):
        return self.ssp_data.data[5]    
    def get_mass_weighted_age(self):
        return self.ssp_data.data[6]   
            
    def get_age_estimation_error(self):    
        return self.ssp_data.data[7]    
    
    def get_velocity_map(self):    
        #TODO: Substract recessional velocity 
        vmap = self.ssp_data.data[13]    
        central_v = vmap[vmap.shape[0]//2, vmap.shape[1]//2]    
        return vmap - central_v
    def get_velocity_sigma_map(self):
        return self.ssp_data.data[14]    
    
    def get_velocity_dispersion_map(self):    
        return self.ssp_data.data[15]    
    def get_velocity_dispersion_sigma_map(self):    
        return self.ssp_data.data[16]    
    
    def get_surface_density(self, dust_corr=True):
        """
        Stellar Mass density per pixel without dust correction
        --> log10(Msun/spaxels^2) 	
        
        *CAVEAT*: This map has not considered mass loss due to stellar death.        
        """        
        if dust_corr:
            return self.ssp_data.data[19]
        else:
            return self.ssp_data.data[18]
        
    
# =============================================================================
#   PLOTS  
# =============================================================================

    def plot_segmentation_regions(self, show=True, save=False):
        self.get_continuum_segmentation()
        self.segmentation[self.segmentation==0] = np.nan
        plt.figure()
        plt.imshow(self.segmentation, cmap='flag')
        plt.annotate(int(np.nanmax(self.segmentation)), xy=(.1,.9), xycoords='axes fraction')        
        if save:
            plt.savefig('output/segmentation/segmentation_regions_'+self.galaxy_id+'.png')            
        if not show:
            plt.close()
        # plt.colorbar()
    
    def plot_single_segmentation_region(self, which, show=True, save=False):
        self.get_continuum_segmentation()        
        
        voxel = self.segmentation_cube[which, :, :]
        n_spaxels = voxel[voxel==True].size
        
        plt.figure()
        plt.imshow(np.array(voxel, dtype=int), cmap='Greys')
        plt.annotate('Total regions:' + str(int(np.nanmax(self.segmentation))),
                     xy=(.01,.95), xycoords='axes fraction')        
        plt.annotate('Region:' + str(which), xy=(.01,.9), 
                     xycoords='axes fraction')        
        plt.annotate('Spaxels:' + str(n_spaxels), xy=(.01,.85), 
                     xycoords='axes fraction')        
        if save:
            plt.savefig('output/segmentation/segmentation_regions_'+self.galaxy_id+'.png')            
        if not show:
            plt.close()
        # plt.colorbar()


class Basics(Global_properties, Maps):
    def __init__(self, **kwargs):
        Global_properties.__init__(self, **kwargs)
        Maps.__init__(self, **kwargs)
        
        
# ----------------------------------------------------------------------------            
if __name__=='__main__':
    from matplotlib import pyplot as plt
    properties = Global_properties(survey='CALIFA', galaxy_id='NGC0001')
    maps = Maps(survey='CALIFA', galaxy_id='NGC0001')
    maps.plot_segmentation_regions(show=1)
    gal = Basics(survey='MANGA', galaxy_id='manga-9510-12705')
    maps.plot_segmentation_regions()
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


