import os

from sys import argv

code = argv[1]
decoder = argv[2]
pmin = float(argv[3])
pmax = float(argv[4])
shots = int(argv[5])
n_proc = int(argv[6])

fix_error_opt = '-fix-error' if '-fix-error' in argv else ''

if decoder != 'mwpm' and decoder != 'restriction':
    print('Decoder should be one of \"mwpm\" or \"restriction\".')
    exit()

for version in ['v3.2']:
    folder = f'../data/protean/{code}/{version}'
    for mem_flag in ['', '-mx']:
        cmd = fr'''
                cd Release
                mkdir -p {folder}/output
                mpirun -np {n_proc} ./pr_base_memory {folder} --decoder {decoder} --s {shots}\
                        --pmin {pmin} --pmax {pmax} {fix_error_opt} {mem_flag}
                '''
        print('----------------------------')
        print(cmd)
        os.system(cmd)
