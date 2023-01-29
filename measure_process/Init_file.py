from core_tools.data.SQL.connect import SQL_conn_info_local, set_up_remote_storage, sample_info, set_up_local_storage
set_up_remote_storage('131.180.205.81', 5432, 'xld_measurement_pc', 'XLDspin001', 'sixdots', "6dot", "XLD", "any")


from core_tools.data.gui.qml.data_browser import data_browser
from core_tools.data.gui.plot_mgr import data_plotter
from core_tools.data.ds.data_set import load_by_uuid

from core_tools.data.SQL.queries.dataset_gui_queries import query_for_measurement_results




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
    dire = os.path.join(r'M:\tnw\ns\qt\spin-qubits\projects\Hot qubit SiGe\Data', folder)
    if not os.path.exists(dire):
        os.makedirs(dire)
        
        
    df = pd.DataFrame.from_dict(data, orient='index').transpose()
    df.to_csv(os.path.join(dire,name + ".csv"), index=False, mode='w+',header = list(data.keys()))

def query_database(init_time = '2022-11-16 11:57:50', final_time='2022-11-16 23:59', name='', echo=True):
    '''
    queries the remote DB to extract the info about the measurements satisfyign the query.
    for now we query time-range but other things can be added    


    Parameters
    ----------
    init_time : str, optional
        initial time, see format for default. The default is '2022-11-16 11:57:50'.
    final_time : str, optional
        final query time. The default is '2022-11-16 23:59'.

    Returns
    -------
    res : dictionary
        contains uuid, name, date of all measurements satisfying the query 

    '''
    
    res = query_for_measurement_results.search_query(
            start_time=init_time, end_time=final_time,
            name = name,
            remote=True)
    
    uuids = np.array([], dtype='int')
    names = np.array([])
    dates = np.array([])
    for ii,exp in enumerate(res):
        uuids = np.append(uuids, int(exp.uuid))
        names = np.append(names,exp.name)
        dates = np.append(dates,exp.start_time)
    
    
    p = dates.argsort()
    dates = dates[p]
    names = names[p]
    uuids = uuids[p]
    
    datasets = []
    for ii in range(len(dates)):
        datasets.append({'uuid':uuids[ii],'name':names[ii],'date':dates[ii]})
    
    
    if echo:
        print("{:<3}\t{:<20}\t{:<20}".format('num','UUID','name'))
        for ii, ds in enumerate(datasets):
            print("{:<3}\t{:<20}\t{:<20}".format(ii,ds['uuid'],ds['name']))
        
    return datasets, names, uuids
 
 


if __name__ == "__main__":
    _db = data_browser() 
