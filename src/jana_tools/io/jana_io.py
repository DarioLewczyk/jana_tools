# imports: {{{ 
import os
import pandas as pd
#}}}
# export_dataframe: {{{
def export_dataframe(df:pd.DataFrame = None, filename:str = None, filepath:str = None, sheet_name:str = 'Sheet1', overwrite:bool = False):
    '''
    filename: the filename. will add .xlsx inside of the code
    filepath: the path where you want the file to end up. Otherwise, will save in current directory.
    sheet_name: the name of the sheet you want to save the data to if applicable
    '''
    home = os.getcwd()
    filename = f'{filename}.xlsx'
    if filepath:
        if os.path.exists(filepath):
            os.chdir(filepath)
        else:
            os.mkdir(filepath)
            os.chdir(filepath) 
    if filename in os.listdir() and not overwrite:  
        with pd.ExcelWriter(filename, engine='openpyxl', mode = 'a') as writer:  
            df.to_excel(writer, sheet_name=sheet_name)
            
    else:
        with pd.ExcelWriter(filename) as writer:
            df.to_excel(writer, sheet_name=sheet_name)
    print(f'Your file was saved to: {os.path.join(filepath,filename)} in sheet: {sheet_name}')
    os.chdir(home)
#}}}
