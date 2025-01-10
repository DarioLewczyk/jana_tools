# Authorship: {{{
# Written by: Dario C. Lewczyk
# Date: 09-13-2024
#}}}
# Imports: {{{
from topas_tools.plotting.plotting_utils import GenericPlotter
from topas_tools.utils.topas_utils import Utils, UsefulUnicode
import re
import os
import numpy as np
#}}}
# JANA_Plot: {{{
class JANA_Plot(GenericPlotter, UsefulUnicode):
    # __init__: {{{
    def __init__(self, hklm_dir:str = None): 
        GenericPlotter.__init__(self) 
        UsefulUnicode.__init__(self)
    #}}}
    # plot_pattern_with_hkl: {{{
    def plot_pattern_with_hkl(self,
            index:int = 0, 
            jana_data:dict = None, 
            plot_vs_q:bool = False, 
            composite:bool = False,
            hkl_offset:int = -1800, 
            hkl_names:list = ['main', 'satellites'],
            hkl_colors:list = ['blue', 'orange'],
            height:int = 650,
            width:int = 1000,
            marker_size:int = 8, 
            ):
        '''
        jana_data: This is a dictionary created by jana_tools.py
        plot_vs_q: do you want to be on a q scale or a 2theta scale? 
        composite: Do you want to plot hkls separated as in a composite structure or not?
        hkl_offset: the standard separation for hkl ticks
        hkl_names: a list of names for each of the hkls to be plotted.
        hkl_colors: a list of colors for the hkls
        height: height of the plot
        width: width of the plot
        marker_size: size of the hkl ticks
        '''
        # definitions: {{{
        pattern = jana_data[index]['pattern'] 
        hklm_data = jana_data[index]['hklm_data']
        yobs = pattern['yobs']
        # axis selection for pattern: {{{
        if plot_vs_q:
            pattern_x = pattern['q']
            hovertemplate = 'q: %{x}<br>Intensity: %{y}'
            xaxis_title = 'q (Å^-1)'
        else:
            pattern_x = pattern['tth']
            hovertemplate = '2theta: %{x}<br>Intensity: %{y}'
            xaxis_title = f'2{self._theta}{self._degree_symbol}'
        #}}}
        # dictionaries for hkls: {{{
        if composite:
            composite_dicts = jana_data[index]['composite_hklm']
            dictionaries = [
                composite_dicts['primary'], 
                composite_dicts['secondary'],
                composite_dicts['common'],
                composite_dicts['satellites'],
            ]
        else:
            dictionaries = [hklm_data['main'], hklm_data['satellite']]
        #}}}
        #}}}
        # Plot the pattern: {{{
        self.plot_data(
            pattern_x,
            yobs,
            name = 'Observed',
            color = 'black',
            mode = 'lines',
            xaxis_title=xaxis_title,
            yaxis_title='Intensity',
            width=width,
            height = height,
            hovertemplate = hovertemplate,
        )
        #}}}
        # plot the dictionaries of hkls: {{{
        hkl_plot_count = 1 # Allows the position of the hkls to be adjusted.
        for i, hklm_dict in enumerate(dictionaries):
            try:
                ht = hklm_dict['ht'] # for the composite dictionary
            except:
                ht = hklm_dict['hovertemplate'] # for the basic dictionary
            name = hkl_names[i]  
            color = hkl_colors[i]
            if plot_vs_q:
                hklm_x = hklm_dict['q']
            else:
                hklm_x = hklm_dict['tth']
            y = np.ones(len(hklm_x)) * (hkl_offset*hkl_plot_count)
            hkl_plot_count += 1
            self.add_data_to_plot(
                hklm_x,
                y,
                name = name, 
                symbol = 'line-ns',
                color = color,
                hovertemplate = ht,
                marker_size=marker_size
            )
        #}}}
        self.show_figure()
    #}}}

#}}}
