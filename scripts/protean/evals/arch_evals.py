CODES = {
    'hexcc': {
        '': [
            '7_1_3_3',
            '19_1_5_5',
            '37_1_7_7',
            '61_1_9_9'
        ]
    },
    'hysc': {
        '3_8': [
            '48_6_6_3',
            '96_10_10_4',
            '216_20_12_5',
            '504_44_14_6',
            '2448_206_16_8',
        ],
        '4_5': [
            '60_8_6_4',
            '160_18_8_6',
            '960_98_10_8',
            '1920_194_12_10'
        ],
        '4_6': [
            '36_8_4_4',
            '336_58_8_6',
            '864_146_10_8'
        ],
        '5_5': [
            '30_8_3_3',
            '40_10_4_4',
            '80_18_5_5',
            '150_32_6_6',
            '900_182_8_8'
        ],
        '5_6': [
            '60_18_4_3',
            '120_34_6_5',
            '2520_674_8_6'
        ]
    },
    'hycc': {
        '3_8': [
            '48_8_4_4',
            '96_12_6_6',
            '192_20_8_8',
            '432_40_10_10',
            '1296_112_12_12'
        ],
        '3_10': [
            '120_20_4_4',
            '150_24_6_6',
            '360_52_8_8'
        ],
        '4_6': [
            '24_8_4_4',
            '96_20_6_6',
            '192_36_8_8',
            '240_44_8_8'
#           '1536_260_12_12'
        ],
        '4_10': [
            '40_16_4_4',
            '480_148_6_6',
            '640_196_8_8'
        ],
    }
}

from sys import argv
import os

def code_params(code):
    data = code.split('_')
    return tuple([int(x) for x in data])

def run_experiment_0(family, subfamily):
    code_list = CODES[family][subfamily]

    jid_blk = 'jid.ral' if 'cc' in family else ''
    other_flags = '-color-checks' if 'cc' in family else ''

    arch_folder_prefix = f'../data/protean/{family}/{subfamily}'

    version_list = [
#       ('1', f'{jid_blk}.rlb.rcr', 4, ''),          # Baseline: no modifications
#       ('2', f'{jid_blk}.fla.ral.rlb.rcr', 4, ''),  # Only flags
#       ('3.3', f'{jid_blk}.fla.ral.(prx.ral)+.rlb.rcr', 3, '-fno-opt-flags -fflag-jid'),  # Naive connectivity 3 implementation
#       ('3.4', f'{jid_blk}.fla.ral.(prx.ral)+.rlb.rcr', 4, '-fno-opt-flags -fflag-jid'),  # Naive connectivity 4 implementation
#       ('4.3', f'{jid_blk}.fla.ral({jid_blk}.prx.ral)+rlb.rcr', 3, '-fflag-jid'),  # Optimized connectivity 3
        ('4.4', f'{jid_blk}.fla.ral({jid_blk}.prx.ral)+rlb.rcr', 4, '-fflag-jid')   # Optimized connectivity 4
    ]

    for code in code_list:
        _, _, dx, dz = code_params(code)
        rounds = 1#min(dx, dz)
        arch_folder = f'../data/protean/{family}/{subfamily}/{code}'
        tanner_file = f'../data/tanner/{family}/{subfamily}/{code}.txt'
        for (vno, passes, conn, flags) in version_list:
            output_folder = f'{arch_folder}/v{vno}'
            render_folder = f'{output_folder}/render'
            print('----------------------------')
            protean_cmd = fr'''
                        cd Release
                        mkdir -p {output_folder}
                        mkdir -p {render_folder}
                        ./protean {tanner_file} {output_folder} --passes "{passes}" --s-rounds {rounds} --max-conn {conn} {other_flags}
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
    run_experiment_0(family, subfamily)
elif exno == 1:
    run_experiment_1(family, subfamily)

