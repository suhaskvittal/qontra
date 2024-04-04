
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
        '4_6': [ '24_8_4_4', '96_20_6_6', '192_36_8_8'],
        '4_8': [ '336_88_8_8', '3456_868_12_12' ],
        '4_10': [ '40_16_4_4', '480_148_6_6', '640_196_8_8' ],
        '5_6': [ '120_36_4_4', '600_164_8_8' ]
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
        f_skip_schedule = '-skip-schedule' if n > 1000 else ''
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

def run_experiment_1(family, subfamily):
    code_list = CODES[family][subfamily]
    decoder = 'restriction' if family == 'hycc' else 'mwpm'
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

