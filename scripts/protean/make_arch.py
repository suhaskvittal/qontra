"""
    author: Suhas Vittal
    date:   16 January 2024
"""

import os

from sys import argv

code = argv[1]
rounds = argv[2]

tanner_graph_file = '../data/tanner/%s.txt' % code
arch_folder_prefix = '../data/protean/%s' % code

# Version list syntax:
#   (version-string, pass-string, max-connectivity, options)
jid_blk = 'jid.ral.' if '-use-jid' in argv else ''

version_list = [
    ('1', 'ral.rlb.rcr', 4, ''),
#   ('2.1', 'fla.ral.rlb.rcr', 4, ''),
#   ('2.2', f'{jid_blk}fla.ral.{jid_blk}rlb.rcr', 4, '-flag-jid'),
#   ('3.1', 'fla.ral(con.ral.prx.ral)+rlb.rcr', 4, ''),
    ('3.2', f'{jid_blk}fla.ral.{jid_blk}(con.ral.prx.ral)+rlb.rcr', 4, '-flag-jid'),
#   ('4.1', 'fla.ral(con.ral.prx.ral)+rlb.rcr', 3, ''),
    ('4.2', f'{jid_blk}fla.ral.{jid_blk}(con.ral.prx.ral)+rlb.rcr', 3, '-flag-jid'),
]
color_opt = '-color-checks' if '-color-checks' in argv else ''

for (version, pass_string, max_conn, other_opt) in version_list:
    output_folder = f'{arch_folder_prefix}/v{version}'
    render_folder = f'{output_folder}/render'
    print('----------------------------')
    
    protean_cmd = fr'''
                cd Release
                mkdir -p {output_folder}
                mkdir -p {render_folder}
                ./protean {tanner_graph_file} {output_folder}\
                        --passes "{pass_string}"\
                        --s-rounds {rounds}\
                        --max-conn {max_conn}\
                        {color_opt} {other_opt}
              '''
    print(protean_cmd)
    os.system(protean_cmd)
