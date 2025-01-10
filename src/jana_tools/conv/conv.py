# Authorship: {{{ 
# Dario C. Lewczyk
# 10-02-24
#}}}
# Imports: {{{
from topas_tools.utils.topas_utils import DataCollector
#}}}
# Converter: {{{
class Converter(DataCollector):
    # __init__: {{{
    def __init__(self):
        DataCollector.__init__(self,mode = 1)
    #}}} 
    # convert_topas_xy_to_dat: {{{
    def convert_topas_xy_to_dat(self, 
            extension:str = 'txt',
            lattice_prms:list = [1,1,1,90,90,90],
            ):
        '''
        This will convert all text files in a directory
        to the jana format and use the lattice params you give
        '''
        try:
            lattice_prms = tuple(lattice_prms)
            a, b, c, al, be, ga = lattice_prms # using tuple, unpack
        except:
            print('You must give a, b, c, alpha, beta, gamma')
        self.scrape_files(extension) # get the files
        for file in self.files:
            basename = file.strip(f'.{extension}') # Get the basename
            with open(file) as f:
                lines = f.readlines()
                with open(f'{basename}.dat', 'w') as f2:
                    for i, line in enumerate(lines):
                        if i == 0: 
                            f2.write(f'#cell {a} {b} {c} {al} {be} {ga}\n')
                        elif i == 1:
                            f2.write('#X\tI\n')
                        else:
                            splitline = line.split(',')
                            tth = splitline[0]
                            intensity = splitline[1] # contains \n
                            f2.write(f'{tth}\t{intensity}')

    #}}}
#}}}
