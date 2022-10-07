#Load database
print('Loading Database')
from core_tools.data.SQL.connect import SQL_conn_info_local, set_up_remote_storage, sample_info, set_up_local_storage
set_up_remote_storage('131.180.205.81', 5432, 'xld_measurement_pc', 'XLDspin001', 'sixdots', "6dot", "XLD", "6D2S - SQ21-XX-X-XX-X")


from core_tools.data.gui.qml.data_browser import data_browser
from core_tools.data.gui.plot_mgr import data_plotter
from core_tools.data.ds.data_set import load_by_uuid

data_browser() 



#Libs
import matplotlib.pyplot as plt
import numpy as np

#Extract Data
import pandas as pd
import csv
import os 
def data_dump(data, folder, name):
    '''
    Dumps the data into a .csv file

    Parameters
    ----------
    data : dictionary
        Dictionary with keys for column names. example: dic_data = {'Temp':temps,'T1':np.empty(len(temps)),'std':np.empty(len(temps))}.
    folder : str
        Folder naming the kind of experiment
    name : str
        Name of the file

    Returns
    -------
    None.

    '''
    dire = os.path.join(r'M:\tnw\ns\qt\spin-qubits\data\Oriol Pietx\Data', folder)
    if not os.path.exists(dire):
        os.makedirs(dire)
        
        
    df = pd.DataFrame.from_dict(data, orient='index').transpose()
    df.to_csv(os.path.join(dire,name + ".csv"), index=False, mode='w+',header = list(data.keys()))
    
