#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 08:50:06 2020

@author: pablo
"""

from astropy.io import fits
import numpy as np 
from pandas import read_csv 


"""
This script contains the basic instrucctions to read the Pipe3D outputs 
for every given survey in a common format.
--> CURRENT VERSION CONTAINS METHODS FOR READING CALIFA AND MANGA OUTPUTS

For additional surveys a *survey*_starter class must be created.
"""
# =============================================================================
# MANGA
# =============================================================================
class MANGA_starter(object):    
   
    fiber_angle_diameter = 0.5/3600 *np.pi/180 # rad                    
    
    def __init__(self, **kwargs):
        self.catalog_path = '/home/pablo/obs_data/MANGA/Pipe3D/manga.Pipe3D-v2_4_3.fits'
        self.cubes_path = '/home/pablo/obs_data/MANGA/Pipe3D/cubes'        
                
        galaxy_id = kwargs['galaxy_id']
        
        self.cube_path = self.cubes_path+'/'+galaxy_id+'.Pipe3D.cube.fits'                                                            
        self.cube = fits.open(self.cube_path)
        
        self.catalog = fits.open(self.catalog_path)[1].data      
        self.catalog_entry = np.where(
                                   self.catalog['mangaid'] == galaxy_id
                                     )[0][0]
                                              
    def own_keywords(self, keyword):
        keyword_converter = { 'redshift':'redshift',
                             'stellar_mass': 'log_mass'
                             }
        return keyword_converter[keyword]
    
    def file_SSP(self):
        return self.cube[1]
    
    def file_SFH(self):
        return self.cube[2]
    
    def file_FLUX_ELINES(self):
        return self.cube[3]
    
    def file_INDICES(self):
        return self.cube[4]
                
    def close_all(self):
        self.cube.close()
        self.catalog.close()
# =============================================================================
# CALIFA
# =============================================================================
class CALIFA_starter(object):    
    
    fiber_angle_diameter = 2.7/3600 *np.pi/180 # rad                    
    
    def __init__(self, **kwargs):
        self.catalog_path = '/home/pablo/obs_data/CALIFA/Pipe3D/DR2_Pipe3D_obj.tab.csv'
        self.cubes_path = '/home/pablo/obs_data/CALIFA/Pipe3D/cubes'        
        
        galaxy_id = kwargs['galaxy_id']
        
        self.cubes_extension = self.cubes_path+'/'+galaxy_id                                                            
                    
        self.catalog = read_csv(self.catalog_path)
        
        self.catalog_keys = self.catalog.keys()
        
        self.catalog_entry = np.where(
                                    self.catalog[self.catalog_keys[0]] == galaxy_id
                                      )[0][0]
    def own_keywords(self, keyword):
        keyword_converter = { 'redshift':'# (5) Pipe3D redshift',
                             'stellar_mass':'# (6) log(Mass/Msun)'
                             }
        return keyword_converter[keyword]
        
    def file_SSP(self):
        
        hdul = fits.open(self.cubes_extension+'.SSP.cube.fits')
        ssp = hdul[0]
        # hdul.close()
        return ssp
    
    def file_SFH(self):
        hdul = fits.open(self.cubes_extension+'.SFH.cube.fits')        
        sfh = hdul[0]
        # hdul.close()
        return sfh
    
    def file_FLUX_ELINES(self):
        hdul = fits.open(self.cubes_extension+'.ELINES.cube.fits')        
        elines = hdul[0]
        # hdul.close()
        return elines
    
    def file_INDICES(self):
        raise NameError('Method not implemented')
                
    def close_all(self):
        raise NameError('Method not implemented')

# =============================================================================
# Others
# =============================================================================
class ANOTHER_SURVEY_starter(object):
        pass            
    
# =============================================================================
# BUILDER FUNCTION
# =============================================================================

class Builder(object):
    def __init__(self, **kwargs):
        self.survey = kwargs['survey']
        self.galaxy_id = kwargs['galaxy_id']
        
        if self.survey == 'MANGA':            
            self.starter = MANGA_starter
        
        elif self.survey == 'CALIFA':
            self.starter = CALIFA_starter
        
        
    def build_data_bundle(self):
        return self.starter(galaxy_id = self.galaxy_id)        
    
# =============================================================================
#         Example
# =============================================================================

if __name__=='__main__':
    galaxy = Builder(survey='CALIFA', galaxy_id='NGC0001').build_data_bundle()
    galaxy2 = Builder(survey='MANGA', galaxy_id='manga-9510-12703').build_data_bundle()
    
# By Mr. Krtxo       