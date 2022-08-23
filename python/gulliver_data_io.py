"""
    author: Suhas Vittal
    date:   23 August 2022
"""

import numpy as np
import pandas as pd

def read_timing_analysis_file(file_path, max_time=np.inf):
    reader = open(file_path, 'r')
    # Pandas data columns 
    syndrome_column = []
    hamming_weight_column = []
    time_taken_column = []

    line = reader.readline()
    while line != '':
        data_tuple = line.split(' ')
        line = reader.readline()
        # Extract data
        syndrome = data_tuple[0]
        time_taken = int(data_tuple[1])
        hamming_weight = 0
        # Skip this line if time taken is too large
        if time_taken > max_time:
            continue
        for i in range(len(syndrome)):
            if syndrome[i] != '0':
                hamming_weight += 1
        syndrome_column.append(syndrome)
        hamming_weight_column.append(hamming_weight)
        time_taken_column.append(time_taken)
    reader.close()
    # Sort the data.
    zipped_data = zip(syndrome_column, hamming_weight_column, time_taken_column)
    zipped_data = list(zipped_data)
    zipped_data.sort(key=lambda p: p[2])
    # Unzip the data
    syndrome_column = [x for (x,_,_) in zipped_data]
    hamming_weight_column = [x for (_,x,_) in zipped_data]
    time_taken_column = [x for (_,_,x) in zipped_data]
    index_column = list(range(len(syndrome_column)))

    df = pd.DataFrame({
        'Syndrome Index': index_column,
        'Syndrome': syndrome_column,
        'Hamming Weight': hamming_weight_column,
        'Time Taken': time_taken_column
    })
    return df

