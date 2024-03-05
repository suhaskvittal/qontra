import os

from sys import argv

code = argv[1]
decoder = argv[2]
shots = int(argv[3])

n_proc = int(argv[4])

fix_error_opt = '-fix-error' if '-fix-error' in argv else ''

if decoder != 'mwpm' and decoder != 'restriction':
    print('Decoder should be one of \"mwpm\" or \"restriction\".')
    exit()

for version in ['v1', 'v3.2', 'v4.2']:
    folder = f'../data/protean/{code}/{version}'
    cmd = fr'''
            cd Release
            mkdir -p {folder}/output
            mpirun -np {n_proc} ./pr_base_memory {folder} --decoder {decoder} --s {shots} {fix_error_opt}
            '''
    print('----------------------------')
    print(cmd)
    os.system(cmd)
