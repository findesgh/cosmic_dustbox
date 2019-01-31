import astropy.units as _u
import numpy as _np
import warnings as _warnings
import os as _os
import miepython as miepy
from scipy import interpolate as _interp
import sys
Path =  _os.path.dirname(_os.path.dirname(_os.path.realpath(__file__))) + \
                 '//Solver//'  
sys.path.insert(0, Path)
import Solver as sv
###############################################################################
###############################################################################


class GrainSpecies(object):

    def __init__(self):
        return


class SpAsil(Grainspecies): 
    
        
        def getSilica_props(en):
            
            en = en.to(u.micron, equivalencies=u.spectral())
            if en.unit != _u.eV or en.unit != en.J :
                raise('The unit of Energy should be in eV or Joule')
            elif en.unit == en.J:
                en = en * _u.J.to(_u.eV)
                
            lamda = en.unit * _u.eV.to(_u.micron, \
                equivalencies=_u.spectral())         
            path = cls.loadFiles['Silica'] 
            silica = loadOprops(path)
            assert lamda.unit == _u.micron, 'The wave length is not in micron'        
            Silica_d = Find_Interpolate(lamda,silica)
            return Silica_d         


        def Sigma_abs(en,sizes):
            m = getSilica_props(en,sizes)
            QABS = []
            sigma_abs = []

            for i, size in enumerate(sizes):
                QABS append.(miepy.mie(m,size))
                sigma_abs.append.(_np.pi * size**2 * QABS[i])
            return _np.array(sigma_abs) 


        def Sigma_sca(en,sizes):
            m = getSilica_props(en,sizes)
            QSCA = []
            sigma_sca = []

            for i, size in enumerate(sizes):
                QSCA append.(miepy.mie(m,size))
                sigma_sca.append.(_np.pi * size**2 * QSCA[i])
            return _np.array(sigma_sca) 
        
    
    
    class SpGraphite(Grainspecies): 
                """
                """

        def getCarbon_props(en):
            
            if en.unit != _u.eV or en.unit != en.J :
                raise('The unit of Energy should be in eV or Joule')
            elif en.unit == en.J:
                en = en * _u.J.to(_u.eV)

            lamda = en.to(_u.micron, \
                          equivalencies=_u.spectral())         
        path1 = cls.loadFiles['carbon_par_0.01'] 
        c_par_0.01_ = loadOprops(path1)
        path2 = cls.loadFiles['carbon_par_0.1'] 
        c_par_0.1_ = loadOprops(path1)
        path3 = cls.loadFiles['carbon_per_0.01'] 
        c_per_0.01_d = loadOprops(path3)
        path4 = cls.loadFiles['carbon_per_0.1'] 
        c_per_0.1_d = loadOprops(path4)
        assert lamda.unit == _u.micron, 'The wave length is not in micron'
#        
        c_par_0.01_d = Find_Interpolate(lamda,c_par_0.01_)
        c_par_0.1_d = Find_Interpolate(lamda,c_par_0.1_)
        c_per_0.01_d = Find_Interpolate(lamda,c_per_0.01_)
        c_par_0.1_d = Find_Interpolate(lamda,c_per_0.1_)
#        
        return c_par_0.01_d , c_par_0.1_d , c_per_0.01_d , c_per_0.1_d  
        

        
        def Sigma_abs(en,sizes):
            
            c_par_0.01_, c_par_0.1_, c_per_0.01_, c_per_0.1_ =
                getCarbon_props(en,sizes)

            QABS_par_0.01_ = [] ; QABS_par_0.1_ = [] ; QABS_per_0.01_ = [] 
            QABS_per_0.1_ = []
            sigma_par_0.01_ = [] ; sigma_par_0.1_ = [] ; sigma_per_0.01_ = []
            sigma_per_0.01_ = []
            
            for i, size in enumerate(sizes):                
                QABS_par_0.01_ = miepy.mie(c_par_0.01_ ,size)
                QABS_par_0.1_  = miepy.mie(c_par_0.1_ ,size)
                QABS_per_0.01_ = miepy.mie(c_per_0.01_ ,size)
                QABS_per_0.1_  = miepy.mie(c_per_0.1_ ,size)
#
                sigma_par_0.01_.append.(_np.pi * size**2 * QABS_par_0.01_[i])
                sigma_par_0.1_.append.(_np.pi * size**2 * QABS_par_0.1_[i])
                sigma_per_0.01_.append.(_np.pi * size**2 * QABS_per_0.01_[i])
                sigma_per_0.1_.append.(_np.pi * size**2 * QABS_per_0.1_[i])
#                
            sigma_abs_0.01_ = ( 2 * _np.array(sigma_par_0.01_) + \
                           _np.array(sigma_per_0.01_) ) / 3 

            sigma_abs_0.1_ = ( 2 * _np.array(sigma_par_0.1_) + \
                           _np.array(sigma_per_0.1_) ) / 3 
#
            return sigma_abs_0.01_ , sigma_abs_0.1_

#
        def Sigma_sca(en,sizes):
#            
            c_par_0.01_, c_par_0.1_, c_per_0.01_, c_per_0.1_ =
                getCarbon_props(en,sizes)

            QSCA_par_0.01_ = [] ; QSCA_par_0.1_ = [] ; QSCA_per_0.01_ = [] 
            QSCA_per_0.1_ = []
            sigma_par_0.01_ = [] ; sigma_par_0.1_ = [] ; sigma_per_0.01_ = []
            sigma_per_0.01_ = []
            
            for i, size in enumerate(sizes):                
                QSCA_par_0.01_ = miepy.mie(c_par_0.01_ ,size)
                QSCA_par_0.1_  = miepy.mie(c_par_0.1_ ,size)
                QSCA_per_0.01_ = miepy.mie(c_per_0.01_ ,size)
                QSCA_per_0.1_  = miepy.mie(c_per_0.1_ ,size)
#
                sigma_par_0.01_.append.(_np.pi * size**2 * QSCA_par_0.01_[i])
                sigma_par_0.1_.append.(_np.pi * size**2 * QSCA_par_0.1_[i])
                sigma_per_0.01_.append.(_np.pi * size**2 * QSCA_per_0.01_[i])
                sigma_per_0.1_.append.(_np.pi * size**2 * QSCA_per_0.1_[i])
#                
            sigma_sca_0.01_ = ( 2 * _np.array(sigma_par_0.01_) + \
                           _np.array(sigma_per_0.01_) ) / 3 

            sigma_sca_0.1_ = ( 2 * _np.array(sigma_par_0.1_) + \
                           _np.array(sigma_per_0.1_) ) / 3 
#
            return sigma_sca_0.01_ , sigma_sca_0.1_
