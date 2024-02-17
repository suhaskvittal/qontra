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
#   (version-string, pass-string, max-connectivity)
if '--use-jid' in argv:
    version_list = [
        ('1', 'ral.rlb.rcr', 4),
        ('2.1', 'fla.ral.rlb.rcr', 4),
        ('2.2', 'jid.ral.fla.ral.jid.ral.rlb.rcr', 4),
        ('3.1', 'fla.ral(con.ral,prx.ral)+rlb.rcr', 4),
        ('3.2', 'jid.ral.fla.ral.jid.ral(con.ral.prx.ral)+rlb.rcr', 4),
        ('4.1', 'fla.ral(con.prx.ral)+rlb.rcr', 3),
        ('4.2', 'jid.ral.fla.ral.jid.ral(con.ral.prx.ral)+rlb.rcr', 3),
    ]
else:
    version_list = [
        ('1', 'ral.rlb.rcr', 4),
        ('2', 'fla.ral.rlb.rcr', 4),
        ('3', 'fla.ral(con.ral.prx.ral)+rlb.rcr', 4),
        ('4', 'fla.ral(con.ral.prx.ral)+rlb.rcr', 3),
    ]

for (version, pass_string, max_conn) in version_list:
    output_folder = '%s/v%s' % (arch_folder_prefix, version)
    render_folder = '%s/render' % output_folder
    print('----------------------------')
    
    protean_cmd = r'''
                cd Release
                mkdir -p %s
                mkdir -p %s
                ./protean --tanner %s --out %s --passes "%s" --s-rounds %s --render %s --max-conn %d -v
              ''' % (output_folder,
                     render_folder,
                     tanner_graph_file,
                     output_folder,
                     pass_string,
                     rounds,
                     render_folder,
                     max_conn)
    print(protean_cmd)
    os.system(protean_cmd)
