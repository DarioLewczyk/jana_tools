# Authorship: {{{
# Written by: Dario C. Lewczyk
# Date: 09-13-2024
#}}}
# imports: {{{ 
import os
from glob import glob
import pandas as pd
import numpy as np
from topas_tools.utils.topas_utils import Utils
from jana_tools.plotting.jana_plotting import JANA_Plot
import re
#}}}
# JANA_Tools: {{{ 
class JANA_Tools(Utils, JANA_Plot):
    # __init__: {{{ 
    def __init__(self, hklm_dir:str = None):
        Utils.__init__(self)
        JANA_Plot.__init__(self)
        # define internal variables: {{{
        self.jana_data = {} # This will store the relevant JANA data for you.
        #}}}
        # Get hklm data from directory: {{{
        if hklm_dir == None or not os.path.isdir(hklm_dir):
            print(f'There is no directory: {hklm_dir}. Navigate to the directory where your .prf files are located')
            hklm_dir = self.navigate_filesystem()
        os.chdir(hklm_dir)
        self.hklm_dir = hklm_dir
        #}}}
    #}}}
    # get_hklm_data: {{{
    def get_hklm_data(self,
            fileextension:str = 'prf', 
            modulated:bool = True, 
            num_cols:int = 17,
            lambda_angstrom:float = 1.540593,
            ):
        '''
        This function searches the directory given for 
        .prf files which are where JANA stores the hklm data
        along with other relevant information
        '''
        prf_files = glob(f'*.{fileextension}')
        self.jana_data['data_files'] = prf_files
        # collect and store the data: {{{
        for i, prf_fn in enumerate(prf_files):
            self.prf_file_parser(prf_fn, i,  modulated, num_cols, lambda_angstrom) 
        #}}} 
    #}}}
    # get_pattern_data: {{{
    def get_pattern_data(self, 
            fileextension:str = 'm90',  
            ):
        '''
        This function is used to collect patterns from JANA data files.
        '''
        files = glob(f'*.{fileextension}')
        # Loop through files: {{{
        for i, file in enumerate(files):

            # set the jana_data dictionary: {{{ 
            try:
                self.jana_data[i].update({
                    'pattern':{
                    'header': {},
                    'tth':[],
                    'q': [],
                    'yobs':[],
                    'error':[],
                    }
                })
            except:
                self.jana_data[i] = {
                    'pattern':{
                    'header': {},
                    'tth':[],
                    'q': [],
                    'yobs':[],
                    'error':[],
                    }
                } 

        #}}}
            with open(file, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    clean_line = self._clean_line(line)
                    last_entry = None # This is used to keep a record for updating the dictionary
                    # Loop through the entries in the cleaned line: {{{
                    for j, item in enumerate(clean_line):
                        string, value = self._is_string(item) # string is a bool and value is either an int, string, or float
                        if last_entry:
                            self.jana_data[i]['pattern']['header'][last_entry] = value # Set the paired value
                            last_entry = None
                        if j%2 == 0:
                            if string: 
                                self.jana_data[i]['pattern']['header'][value] = None # Set up to record the paired value
                                last_entry = value
                        if not string:
                            if j == 0:
                                self.jana_data[i]['pattern']['tth'].append(value)
                                # Try to convert tth to q: {{{
                                try:
                                    lam = self.jana_data[i]['pattern']['header']['lambda']
                                    q = self.convert_to_q(value, lam)
                                    self.jana_data[i]['pattern']['q'].append(q)
                                except:
                                    pass
                                #}}}
                            if j == 1:
                                self.jana_data[i]['pattern']['yobs'].append(value)
                            if j == 2:
                                self.jana_data[i]['pattern']['error'].append(value)
                    #}}}
                    


        #}}} 
        
    #}}}
    # prf_file_parser: {{{
    def prf_file_parser(self, 
            prf_fn:str = None, 
            idx:int = 0, 
            modulated:bool = True, 
            num_cols:int = 17,
            lambda_angstrom:float = 1.540593,
            ):
        '''
        This function takes a prf file (prf_fn) and will 
        pull out relevant information from it and drop it into a 
        dictionary and return that dictionary.

        modulated: tell the program to get hklm indices if true
        num_cols: Length of columns in the prf file
        lambda_angstrom: This is used to convert the 2theta into q for direct comparison
                        with other data
        '''
        # Initialize vars: {{{ 
        try:
            self.jana_data[idx].update({
                'hklm_data': {
                    'main': {'peaks': {}},
                    'satellite': {'peaks': {}},
                },
            })
        except:
            self.jana_data[idx] = {
                'hklm_data': {
                    'main': {'peaks': {}},
                    'satellite': {'peaks': {}},
                },
            }
        main_q = []
        main_tth = []
        main_ht = [] # holds strings for hover templates for main reflections
        
        satellite_q = []
        satellite_tth = []
        satellite_ht = [] # holds strings for hover templates for satellite reflections

        main_peaks = 0
        satellite_peaks = 0
        #}}}

        with open(prf_fn) as f:
            lines = f.readlines()
            for line in lines:
                # Split and clean the line: {{{
                clean_line = self._clean_line(line)
                #}}}
                if len(clean_line) == num_cols:
                    # Modulated case: {{{
                    if modulated:
                        h = int(clean_line[0])
                        k = int(clean_line[1])
                        l = int(clean_line[2])
                        m = int(clean_line[3])
                        hklm = f'{h} {k} {l} {m}'

                        fsq = float(clean_line[8])
                        fwhm = float(clean_line[9])
                        tth = float(clean_line[10])

                        d = lambda_angstrom / (2* np.sin(np.pi/180*tth)) # Get d spacing in angstrom
                        q = self.convert_to_q(tth = tth, lambda_angstrom=lambda_angstrom)

                        hklm_ht = '{}<br>hklm: ({})<br>d-spacing: {} {}<br>FSQ: {}<br>FWHM: {}' # format: type, hklm, d-spacing, fsq, fwhm
                        # Handle m=0: {{{
                        if m == 0: 
                            self.jana_data[idx]['hklm_data']['main']['peaks'][main_peaks] = {
                                'hklm': hklm,
                                'h': h,
                                'k': k,
                                'l': l,
                                'm': m,
                                'tth': tth,
                                'q': q,
                                'd-spacing': d,
                                'fsq': fsq,
                                'fwhm': fwhm,

                            }
                            intermediate_hklm_ht = hklm_ht.format('main', hklm, np.around(d,4), self._angstrom, fsq, fwhm)
                            final_hklm_ht = intermediate_hklm_ht+f'<br>tth: {np.around(tth,4)}<br>q: {np.around(q,4)}'
                            main_peaks += 1
                            main_q.append(q)
                            main_tth.append(tth)
                            main_ht.append(final_hklm_ht)
                        #}}} 
                        # m != 0: {{{
                        elif np.abs(m) > 0: 
                            self.jana_data[idx]['hklm_data']['satellite']['peaks'][satellite_peaks] = {
                                'hklm': hklm,
                                'tth': tth,
                                'h': h,
                                'k': k,
                                'l': l,
                                'm': m,
                                'q': q,
                                'd-spacing': d,
                                'fsq': fsq,
                                'fwhm': fwhm,

                            }
                            intermediate_hklm_ht = hklm_ht.format('satellite', hklm, np.around(d,4), self._angstrom, fsq, fwhm)
                            final_hklm_ht = intermediate_hklm_ht+f'<br>tth: {np.around(tth,4)}<br>q: {np.around(q,4)}'
                            satellite_peaks+= 1
                            satellite_q.append(q)
                            satellite_tth.append(tth)
                            satellite_ht.append(final_hklm_ht)
                        #}}}
                        
                    #}}} 
                    # Alternate case: {{{
                    else: 
                        print(f'You have elected to use the non-modulated case.\n'+
                        'This feature is not available yet.')
                    #}}}
            # Add tth and q arrays: {{{
            self.jana_data[idx]['hklm_data']['main']['tth'] = np.array(main_tth)
            self.jana_data[idx]['hklm_data']['main']['q'] = np.array(main_q)
            self.jana_data[idx]['hklm_data']['main']['hovertemplate'] = main_ht
            if modulated:
                self.jana_data[idx]['hklm_data']['satellite']['tth'] = np.array(satellite_tth)
                self.jana_data[idx]['hklm_data']['satellite']['q'] = np.array(satellite_q)
                self.jana_data[idx]['hklm_data']['satellite']['hovertemplate'] = satellite_ht
            #}}} 
    #}}}
    # categorize_composite_hklm: {{{ 
    def categorize_composite_hklm(self,index:int = 0,  modulation_axis:str = 'b'):
        '''
        Use the modulation axis to tell the program 
        which axis the modulation is along so that it can determine 
        the primary, secondary, common, and satellite indices for you.

        index: This is the index of the hklm dictionary data you want
        '''
        # Common axis setup: {{{
        common_h, common_k, common_l = (None, None, None)
        if modulation_axis == 'a':
            common_h = 0
        elif modulation_axis == 'b':
            common_k = 0
        elif modulation_axis == 'c':
            common_l = 0
        #}}}
        # define 4 dictionaries: {{{
        primary = {'tth':[], 'q':[]} # These are the hkl0 reflections
        secondary = {'tth':[], 'q':[]} # These are hklm axes where one of h,k,l is zero
        common = {'tth':[], 'q':[]} # These are m = 0 with a common h, k, or l = 0
        satellites = {'tth':[], 'q':[]} # h,k,l,m reflections

        primary_idx = 0
        secondary_idx = 0
        common_idx = 0
        satellite_idx = 0
        #}}}
        # loop through hklm dictionary: {{{
        for classification, peaks in self.jana_data[index]['hklm_data'].items():
            peaks = peaks['peaks'] # Get the actual peaks entry
            for i, peak in peaks.items():
                h = peak['h']
                k = peak['k']
                l = peak['l']
                m = peak['m']
                tth = peak['tth']
                q = peak['q']
                
                # Flags to categorize.
                add_primary, add_secondary, add_common, add_satellite = (False, False, False, False)

                # Common reflections: {{{
                if common_h == h and m == 0:
                    add_common = True
                if common_k == k and m == 0:
                    add_common = True
                if common_l == l and m == 0:
                    add_common = True
                #}}}
                # Primary reflections: {{{
                if m == 0:
                    add_primary = True
                #}}}
                # Secondary reflections: {{{
                if common_h == h:
                    add_secondary = True
                if common_k == k:
                    add_secondary = True
                if common_l == l:
                    add_secondary = True
                #}}}
                # Satellite reflections: {{{
                if m!=0:
                    if common_h != h and common_h != None:
                        add_satellite = True
                    elif common_k != k and common_k != None:
                        add_satellite = True
                    elif common_l != l and common_l != None:
                        add_satellite = True
                    
                #}}}
                # Update the dictionaries: {{{
                # primary: {{{
                if add_primary:
                    primary['tth'].append(tth)
                    primary['q'].append(q)
                    primary[primary_idx] = peak
                    primary_idx += 1
                #}}}
                # secondary: {{{
                if add_secondary:
                    secondary['tth'].append(tth)
                    secondary['q'].append(q)
                    secondary[secondary_idx] = peak
                    secondary_idx+= 1
                #}}}
                # common: {{{
                if add_common:
                    common['tth'].append(tth)
                    common['q'].append(q)
                    common[common_idx] = peak
                    common_idx+= 1
                #}}}
                # satellite: {{{
                if add_satellite:
                    satellites['tth'].append(tth)
                    satellites['q'].append(q)
                    satellites[satellite_idx] = peak
                    satellite_idx+= 1
                #}}}
                #}}}
        #}}} 
        # Get hovertemplates for plotting: {{{
        primary_ht, secondary_ht, common_ht, satellite_ht = self._get_composite_hovertemplates(primary, secondary, common, satellites)
        primary['ht'] = primary_ht
        secondary['ht'] = secondary_ht
        common['ht'] = common_ht
        satellites['ht'] = satellite_ht
        #}}}
        # Update the jana_data dictionary: {{{
        self.jana_data[index]['composite_hklm'] = {
                'primary': primary,
                'secondary': secondary,
                'common': common,
                'satellites': satellites,
        }
        #}}}
    #}}}
    # _get_composite_hovertemplates: {{{
    def _get_composite_hovertemplates(self, primary:dict = None, secondary:dict = None, common:dict = None, satellites:dict = None):
        '''
        This function will go through each of the dictionaries and generate the text 
        for informative hovering functionality with plots in plotly
        '''
        # initial definitions: {{{
        labels = ['primary', 'secondary', 'common', 'satellites']
        dicts = [primary, secondary, common, satellites]
        base_ht = '2theta: {tth}<br>q: {q}<br>label: {label}<br>hklm: ({hklm})<br>d: {d}<br>f^2: {fsq}<br>fwhm: {fwhm}<br>'
        primary_ht, secondary_ht, common_ht, satellite_ht = ([],[],[],[])
        #}}}
        # loop through the dictionaries: {{{ 
        for i, entry in enumerate(dicts):
            for j, reflection in entry.items():
                try:
                    label = labels[i]
                    hklm = reflection['hklm']
                    d = np.around(reflection['d-spacing'], 4)        
                    fsq = np.around(reflection['fsq'], 4)        
                    fwhm = np.around(reflection['fwhm'], 4)        
                    tth = np.around(reflection['tth'], 4)        
                    q = np.around(reflection['q'], 4)        
                    ht = base_ht.format(tth = tth, q = q, label=label, hklm = hklm, d = d, fsq = fsq, fwhm = fwhm)
                    # categorize the ht: {{{
                    if i == 0:
                        primary_ht.append(ht)
                    if i == 1:
                        secondary_ht.append(ht)   
                    if i == 2:
                        common_ht.append(ht)
                    if i == 3:
                        satellite_ht.append(ht)
                    #}}}
                except:
                    pass
        #}}} 
        return (primary_ht, secondary_ht, common_ht, satellite_ht)
    #}}}
    # _clean_line{{{ 
    def _clean_line(self, line):
        '''
        JANA's lines are a bit strange because they are delimited by spaces and usually also have 
        extra spaces for no apparent reason

        _clean_line() will allow you to pass a string and get a clean version of that string where only real text/values are retained. 
        '''
        splitline = line.split(' ') # split by spaces
        clean_line = [re.sub(r'\n', '', item) for item in splitline]
        clean_line = [item for item in clean_line if item.strip()] # Gets rid of any special or blank characters
        return clean_line
    #}}}
    # _is_string: {{{
    def _is_string(self, value:str = None):
        '''
        Sometimes it is useful to determine when you encounter a true string in a list 
        or if the string is really a float or int. 
        
        This will return (True, str) if the value is a string
        and (False, value) for other datatypes
        '''
        if isinstance(value, str):
            # integer case: {{{
            try: 
                integer = int(value)
                return(False, integer)
            except ValueError:
                pass
            #}}}
            # Float case: {{{
            try: 
                float_val = float(value)
                return (False, float_val)
            except ValueError:
                pass
            #}}}
            # String: {{{
            return (True, value)
            #}}}
        else:
            return False, value
    #}}}
#}}}
