import os

from sys import argv

proc = argv[1]

code = 'hysc/4_5/60_8_6_4'

if not os.path.exists(f'data/protean/{code}'):
    os.system(f'python scripts/protean/make_arch.py {code} 4')
if not os.path.exists(f'data/protean/{code}/v1/output/base_memory.csv'):
    os.system(
        f'''
            mkdir data/protean/{code}/v1/output
            mkdir data/protean/{code}/v2/output
            cd Release
            mpirun -np {proc} ./pr_base_memory ../data/protean/{code}/v1 --decoder mwpm\
                    --s 100000 --pmin 5e-4 --pmax 7e-4 -fix-error
            mpirun -np {proc} ./pr_base_memory ../data/protean/{code}/v1 --decoder mwpm\
                    --s 10000 --pmin 8e-4 --pmax 3e-3 -fix-error
            mpirun -np {proc} ./pr_base_memory ../data/protean/{code}/v2 --decoder mwpm\
                    --s 1000000 --pmin 5e-4 --pmax 7e-4 -fix-error
            mpirun -np {proc} ./pr_base_memory ../data/protean/{code}/v2 --decoder mwpm\
                    --s 100000 --pmin 8e-4 --pmax 1e-3 -fix-error
            mpirun -np {proc} ./pr_base_memory ../data/protean/{code}/v2 --decoder mwpm\
                    --s 10000 --pmin 2e-3 --pmax 3e-3 -fix-error
        ''')
