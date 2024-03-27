import os
from sys import argv
nproc = argv[1]

#for r_s in ['4_5', '5_5']:
#    os.system(f'python scripts/protean/evals/arch_evals.py hysc {r_s} 1 {nproc}')
for r_s in ['3_8', '3_10', '4_6']:
    os.system(f'python scripts/protean/evals/arch_evals.py hycc {r_s} 1 {nproc}')
