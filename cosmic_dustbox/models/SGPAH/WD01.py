from cosmic_dustbox import sdist
from cosmic_dustbox import Gspcs as gs
import astropy.units as _u
import os as _os
import numpy as _np
import sys
import pandas


Path = _os.path.dirname(_os.path.dirname(_os.path.realpath(__file__))) + \
    '\\lib\\paramsWD01a.csv'  
sys.path.insert(0, Path)

class WD01(object): 
    with open(Path)as f:
        R_V = [] ; bc=[] ; alpha_g = [] ; beta_g = [] ; a_tg = [] ; a_cg = []
        C_g=[] ; alpha_s = []; beta_s = []; a_ts = [] ; a_cs = [] ; C_s =[]
        df = pandas.read_csv(f)
        R_V.append(df['R_V'])
        bc.append(df['b_C'])
        alpha_g.append(df['\\alpha_g'])
        beta_g.append(df['\\beta_g'])
        a_tg.append(df['a_{t,g} [m]'])
        a_cg.append(df['a_{c,g} [m]'])
        C_g.append(df['C_g'])
        alpha_s.append(df['\\alpha_s'])
        beta_s.append(df['\\beta_s'])
        a_ts.append(df['a_{t,s} [m]'])
        a_cs.append(df['a_{c,s} [m]'])
        C_s.append(df['C_s'])

        rho_g = 2.24 *_u.g/_u.cm**3
        b = _np.array([0.75 * bc[1],0.25 * bc[2]])
#    
    def __init__(self, Indx):
        
        def load_file(Indx, st):
        sdists = {   
    'gra': sdist.Log_Normal(3.5*_u.Angstrom, _np.inf, rho_g,0.4,bc,b,a0) + \
        	sdist.WD01_dst(3.5*_u.Angstrom, _np.inf , a_tg[Indx],
                        beta_g[Indx],a_cg[Indx], alpha_g[Indx],C_g[Indx]) , \
    'sil': sdist.WD01(  
        	sdist.WD01_dst(3.5*_u.Angstrom, _np.inf,a_ts[Indx], beta_s[Indx], \
            a_cs[Indx],alpha_s[Indx],C_s[Indx]) 
            }        
        return sdist[st]
        
        self.sdist_gr = load_file(Indx,'gra')
        self.sdist_sil = load_file(Indx,'sil')
    return     
#     