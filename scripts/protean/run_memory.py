import os

from sys import argv

code = argv[1]
decoder = argv[2]

n_proc = int(argv[3])

folder = '../data/protean/%s' % code

if decoder != 'mwpm' and decoder != 'restriction':
    print('Decoder should be one of \"mwpm\" or \"restriction\".')
    exit()

shots = 100000

cmd_1 = r'''
        cd Release
        mkdir -p %s/output
        mpirun -np %d ./pr_base_memory %s --decoder %s --s %d
        ''' % (folder, n_proc, folder, decoder, shots)
cmd_2 = r'''
        cd Release
        mkdir -p %s/nn/data
        mpirun -np %d ./pr_nn_data %s --decoder %s --s %d --p 5e-4
        ''' % (folder, n_proc, folder, decoder, 10*shots)
cmd_3 = r'''
        cd Release
        ./pr_nn_train %s 1000
        ''' % (folder)
cmd_4 = r'''
        cd Release
        mkdir -p %s/output
        mpirun -np %d ./pr_nn_memory %s --decoder %s --s %d
        ''' % (folder, n_proc, folder, decoder, shots)

for cmd in [cmd_1]:
    print('----------------------------')
    print(cmd)
    os.system(cmd)
