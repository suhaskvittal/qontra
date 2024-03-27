import os

#for r_s in ['3_8', '4_5', '4_6', '5_5', '5_6']:
#    os.system(f'python scripts/protean/evals/arch_evals.py hysc {r_s} 0')
for r_s in ['4_6']:
    os.system(f'python scripts/protean/evals/arch_evals.py hycc {r_s} 0')
