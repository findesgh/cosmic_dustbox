    loadFiles = { 'carbon_par_0.01': 'callindex.out_CpaD03_0.01.txt', 
                    'carbon_par_0.1':'callindex.out_CpaD03_0.10.txt',
                    'carbon_per_0.01':'callindex.out_CpeD03_0.01.txt' ,
                    'carbon_per_0.1': 'callindex.out_CpeD03_0.10.txt',
                    'silica': 'callindex.out_silD03.txt' }


    def Find_Interpolate(lamda,array):
        w = array[0]
        if lamda in w:
            ind = _np.where(w == lamda)
            Re_n_out = array[1][ind]
            Im_n_out = array[2][ind]
        else:    
            ind1 = _np.where(w >= lamda).min()
            ind2 = _np.where(w <= lamda).max()
            wght1 = (w[ind1] - lamda)/lamda
            wght2= (lamda - w[ind2]) / lamda         
            Re_n_out = (wght1* array[1][ind1] + wght2* array[1][ind2]) / 2 
            Im_n_out = (wght1* array[2][ind1] + wght2* array[2][ind2]) / 2
            Ref_Indx = Re_n_out + j * Im_n_out
        return Ref_Indx


class Crefin(object):

    def __init__(self, f):
        self.f = f
        return

    def __call__(self, a, en):
        return self.f(a, en)


class SGPAHCrefin(Crefin):

    @classmethod
    def loadOprops(path):
        """
        Load the optical properties files in CSV format.

        'path': path to Refractive Index and Material properties file

        returns:
            energy (1D) [eV] at which Q values where computed
            Real Refractive index and Imaginary Refractive index
        """
            rel_path = '\\lib\\' + path
            path = os.path.dirname(os.path.realpath(__file__)) + rel_path

        data = []
        with open(path) as f:
            for line in f:
                if 'number of wavelengths' not in line:
                    continue
                else:
                    line = f.next().split()
                    w = []
                    eps_1 = []
                    eps_2 = []
                    real_n = []
                    Im_n = []
                    while line:
                        try:
                            line = f.next().split()
                        except StopIteration:
                            line = []
                        if line:
                            w.append(float(line[0]))
                            Re_n.append(float(line['Re(n)-1']) + 1)
                            Im_n.append(float(line['Im(n)']))
                    data.append([
                        _np.array(w),
                        _np.array(Re_n),
                        _np.array(Im_n)])

             return data


class Crefin_Asil(SGPAHCrefin):

    def __init__(self):
        return


class Crefin_Gra(SGPAHCrefin):

    def __init__(self):
        return
