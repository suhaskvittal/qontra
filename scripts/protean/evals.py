
from sys import argv
import os

CODES = {
    'hexcc': { '': [ '7_1_3_3', '19_1_5_5', '37_1_7_7', '61_1_9_9' ] },
    'hysc': { 
        '4_5': [ '60_8_6_4', '160_18_8_6', '360_38_8_8', '660_68_10_8', '1800_182_10_10', '1920_194_12_10' ],
        '4_6': [ '36_8_4_4', '96_18_6_4', '324_56_6_6', '336_58_8_6', '720_122_8_8', '864_146_10_8', '2448_410_12_8' ],
        '5_5': [ '30_8_3_3', '40_10_4_4', '80_18_5_5', '150_32_6_6', '900_182_8_8' ],
        '5_6': [ '60_18_4_3', '120_34_6_5', '600_162_6_6', '2520_674_8_6' ]
    },
    'hycc': {
        '4_6': [ '24_8_4_4', '120_24_6_6', '216_40_8_8', '1320_224_10_10', '1440_244_12_12'],
        '4_8': [ '32_12_4_4', '400_104_8_8', '2688_676_12_12' ],
        '4_10': [ '40_16_4_4', '1000_304_8_8' ],
        '5_8': [ '320_116_4_4', '1920_676_8_8' ]
    }
}

def code_params(code):
    data = code.split('_')
    return tuple([int(x) for x in data])

def make_arch(family, subfamily):
    if len(argv) > 4:
        min_size = int(argv[4])
        max_size = int(argv[5])
    else:
        min_size = 1
        max_size = 10000000

    code_list = CODES[family][subfamily]
    arch_folder_prefix = f'../data/protean/{family}/{subfamily}'

    jid_blk = 'jid.ral' if 'cc' in family else ''
    version_list = [
        ('3.4', f'{jid_blk}.fla.ral.(prx.ral)+.rlb.rcr', 4, '-fno-opt-flags -fflag-jid'),  # Naive connectivity 4 implementation
        ('4.4', f'{jid_blk}.fla.ral({jid_blk}.prx.ral)+rlb.rcr', 4, '-fflag-jid')   # Optimized connectivity 4
    ]

    for code in code_list:
        n, _, dx, dz = code_params(code)
        if n < min_size or n > max_size:
            continue
        rounds = min(dx, dz)
        arch_folder = f'data/protean/{family}/{subfamily}/{code}'
        if os.path.exists(arch_folder):
            continue
        tanner_file = f'../data/tanner/{family}/{subfamily}/{code}.txt'
        # Define remaining flags:
        f_color_checks = '-color-checks' if 'cc' in family else ''
        for (vno, passes, conn, flags) in version_list:
            f_skip_schedule = '-skip-schedule' if n > 500 or vno == '3.4' else ''

            output_folder = f'../{arch_folder}/v{vno}'
            print('----------------------------')
            protean_cmd = fr'''
                        cd Release
                        mkdir -p {output_folder}
                        ./protean {tanner_file} {output_folder} \
                                --passes "{passes}"\
                                --s-rounds {rounds}\
                                --max-conn {conn} \
                                {flags} {f_color_checks} {f_skip_schedule} -skip-render
                      '''
            print(protean_cmd)
            os.system(protean_cmd)

def run_experiment_1(family, subfamily):
    code_list = CODES[family][subfamily]
    decoder = 'restriction2' if family == 'hycc' else 'mwpm'
    nproc = int(argv[4])
    
    call_list = [
        # pmin, pmax, shots
        (1e-4, 4e-4, 1_000_000),
        (5e-4, 8e-4, 100_000),
        (9e-4, 2e-3, 10_000)
    ]

    for code in code_list:
        n, k, _, _ = code_params(code)
        if n > 300 or k > 64:
            continue
        for vno in ['4.4']:
            folder = f'../data/protean/{family}/{subfamily}/{code}/v{vno}'
            for mflag in ['', '-mx']:
                for (pmin, pmax, shots) in call_list:
                    cmd = fr'''
                            cd Release
                            mkdir -p {folder}/output
                            mpirun -np {nproc} ./pr_base_memory {folder} --decoder {decoder} --s {shots} --pmin {pmin} --pmax {pmax} -fix-error {mflag}
                            '''
                    print('----------------------------')
                    print(cmd)
                    os.system(cmd)

def run_experiment_2(family, subfamily):
    if len(argv) > 4:
        min_size = int(argv[4])
        max_size = int(argv[5])
    else:
        min_size = 1
        max_size = 1000
    max_size = min(max_size, 1000)

    code_list = CODES[family][subfamily]
    arch_folder_prefix = f'../data/protean/{family}/{subfamily}'

    jid_blk = 'jid.ral' if 'cc' in family else ''
    version_list = [
        ('1', f'ral.rlb.rcr', 1000, '') # Basic implementation for testing scheduling
    ]

    for code in code_list:
        n, _, dx, dz = code_params(code)
        if n < min_size or n > max_size:
            continue
        rounds = min(dx, dz)
        arch_folder = f'data/protean/{family}/{subfamily}/{code}'
        if os.path.exists(f'arch_folder/v1'):
            continue
        tanner_file = f'../data/tanner/{family}/{subfamily}/{code}.txt'
        # Define remaining flags:
        f_color_checks = '-color-checks' if 'cc' in family else ''
        f_skip_schedule = ''
        for (vno, passes, conn, flags) in version_list:
            output_folder = f'../{arch_folder}/v{vno}'
            print('----------------------------')
            protean_cmd = fr'''
                        cd Release
                        mkdir -p {output_folder}
                        ./protean {tanner_file} {output_folder} \
                                --passes "{passes}"\
                                --s-rounds {rounds}\
                                --max-conn {conn} \
                                {flags} {f_color_checks} {f_skip_schedule} -skip-render
                      '''
            print(protean_cmd)
            os.system(protean_cmd)

def run_experiment_3(fam):
    if fam == 'hysc':
        code = 'hysc/5_5/30_8_3_3'
    else:
        code = 'hycc/4_6/24_8_4_4'
    nproc = int(argv[4])

    if fam == 'hysc':
        folders = [f'../data/protean/{code}/v1', f'../data/protean/{code}/v4.4']
        decs = ['pym', 'mwpm']
    else:
        folders = [f'../data/protean/{code}/v4.4', f'../data/protean/{code}/v4.4']
        decs = ['chamberland', 'restriction']
    for (i,f) in enumerate(folders):
        dc = decs[i]
        for mflag in ['', '-mx']:
            cmd = fr'''
                    cd Release
                    mkdir -p {f}/output
                    mpirun -np {nproc} ./pr_base_memory {f} --decoder {dc} --e 100 --pmin 5e-4 --pmax 1e-3 -fix-error {mflag}
                    '''
            print('----------------------------')
            print(cmd)
            os.system(cmd)


family = argv[1]
subfamily = argv[2]
exno = int(argv[3])

if subfamily == 'DNE':
    subfamily = ''

if exno == 0:
    # Make architectures.
    make_arch(family, subfamily)
elif exno == 1:
    run_experiment_1(family, subfamily)
elif exno == 2:
    run_experiment_2(family, subfamily)
elif exno == 3:
    run_experiment_3(family)