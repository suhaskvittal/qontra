import numpy as np
import scipy as sp

import os

from sys import argv

def read_data_file(data_file, prefix):
    e_phys, e_log = [], []
    with open(data_file, 'r') as reader:
        lines = reader.readlines()
        for line in lines[1:]: # Skip first line as it is a header.
            line_data = line.split(',')
            x, y = float(line_data[0]), float(line_data[-1])
            if y > 0:
                e_phys.append(x)
                e_log.append(y)
            if len(e_phys) == 6:
                break
    if len(e_phys) < 5:
        return
    e_phys, e_log = np.array(e_phys), np.array(e_log)
    log_e_phys, log_e_log = np.log(e_phys), np.log(e_log)
    # Perform curve fit between the logs.
    res = sp.stats.linregress(log_e_phys, log_e_log)
    a, b = res.slope, res.intercept
    b = np.e**b
    print(f'{prefix}:\tLER = %.3f * p**%.3f' % (b, a))
    # Pseudo-threshold
    pthres = (1/b)**(1/(a-1))
    print('pseudothreshold:\t%.3e' % pthres)

code = argv[1]

if code == 'rsc':
    for d in [3, 5, 7, 9]:
        data_file = f'data/protean/rsc/output/d{d}.csv'
        if not os.path.exists(data_file):
            continue
        # Read the data from the file.
        read_data_file(data_file, f'rsc d={d}')
else:
    for version in ['v4.3', 'v4.4']:
        folder = f'data/protean/{code}/{version}'
        for m in ['x', 'z']:
            data_file = f'{folder}/output/basic_memory_{m}.csv'
            if not os.path.exists(data_file):
                continue
            # Read the data from the file.
            read_data_file(data_file, f'{version}, {m}')
