""" 
    author: Suhas Vittal
    date:   28 August 2023
"""

from collections import defaultdict

# Each CSS generating function should return:
#   (1) The data qubits
#   (2) The parity qubits
#   (3) The flag qubits
#   (4) The check dictionary
#   (5) Z check list
#   (6) X check list
#   (7) Z observable
#   (8) X observable
# Plus supplemental information in a dict (which can be empty)
# keys: 
#   (a) "color" : dictionary of checks to their color.

def concat(array):
    return ', '.join([str(x) for x in array])

def rsc(d):
    cid = d**2

    checks = {}
    zlist = []
    xlist = []
    
    # Do boundary checks first
    for r in range(0, d-1, 2):  # Z checks
        checks[cid] = [r+1, r, None, None]
        checks[cid+1] = [None, None, d*(d-1) + r+2, d*(d-1) + r+1]
        zlist.append(cid)
        zlist.append(cid+1)
        cid += 2
    for c in range(0, d-1, 2):  # X checks
        checks[cid] = [None, None, (c+1)*d + (d-1), c*d + (d-1)]
        checks[cid+1] = [(c+2)*d, (c+1)*d, None, None]
        xlist.append(cid)
        xlist.append(cid+1)
        cid += 2
    # Do the internal checks now
    for r in range(0, d-1):
        for c in range(0, d-1):
            x, y = c*d + r, (c+1)*d + r
            z, w = c*d + (r+1), (c+1)*d + (r+1)
            #
            #   x   y
            #    \ /
            #    / \
            #   z   w
            #
            if (r + c) & 0x1:
                # This is a Z check
                order = [w, y, z, x]
                zlist.append(cid)
            else:
                order = [w, z, y, x]
                xlist.append(cid)
            checks[cid] = order
            cid += 1
    data_qubits = list(range(d*d))
    parity_qubits = list(range(d*d, 2*d*d-1))
    z_obs = [(k+1)*d-1 for k in range(d)]
    x_obs = [k for k in range(d)]
    return data_qubits, parity_qubits, [], checks, zlist, xlist, [z_obs], [x_obs], {}

def hycc_d4(gr, gf, both_at_once=True, flags=True):
    checks = {}
    """
    Checks
        RED:
            0, 1, 2, 3, 4, 5
            6, 7, 8, 9, 10, 11
            12, 13, 14, 15, 16, 17
            18, 19, 20, 21, 22, 23
        GREEN:
            4, 5, 10, 11, 12, 13, 18, 19 (redundant)
            7, 9, 20, 22, 1, 3, 14, 16
            6, 8, 0, 2, 21, 23, 15, 17
        BLUE:
            9, 11, 18, 20, 3, 5, 12, 14
            2, 4, 8, 10, 19, 21, 13, 15
            0, 1, 6, 7, 16, 17, 22, 23 (redundant)
    """
    red, blue, green = [], [], []

    if both_at_once:
        oct_z_order = [7, 2, 4, 5, 6, 0, 3, 1]
        oct_x_order = [2, 1, 3, 0, 5, 7, 6, 4]
        hex_z_order = [3, 1, 0, 5, 4, 2]
        hex_x_order = [4, 0, 5, 2, 3, 1]
        def oct_order(rot, flip, *qubits):
            qubits = [qubits[(i+rot) % 8] for i in range(8)]
            if flip:
                qubits = qubits[::-1]
            z_sch = [None if i is None else qubits[i] for i in oct_z_order]
            x_sch = [None if i is None else qubits[i] for i in oct_x_order]
            nones = [None] * len(hex_z_order)
            z_sch.extend(nones)
            x_sch.extend(nones)
            return z_sch, x_sch

        def hex_order(rot, flip, *qubits):
            qubits = [qubits[(i+rot) % 6] for i in range(6)]
            if flip:
                qubits = qubits[::-1]
            z_sch = [None if i is None else qubits[i] for i in hex_z_order]
            x_sch = [None if i is None else qubits[i] for i in hex_x_order]
            nones = [None] * len(oct_z_order)
            z_sch = [*nones, *z_sch]
            x_sch = [*nones, *x_sch]
            return z_sch, x_sch

        red.extend(hex_order(0, 1, 0, 1, 3, 5, 4, 2))
        red.extend(hex_order(0, 1, 23, 22, 20, 18, 19, 21))
        red.extend(hex_order(0, 1, 7, 6, 8, 10, 11, 9,))
        red.extend(hex_order(0, 1, 16, 17, 15, 13, 12, 14))

        green.extend(oct_order(gr, gf, 14, 3, 1, 7, 9, 20, 22, 16))
        green.extend(oct_order(gr, gf, 2, 8, 6, 23, 21, 15, 17, 0))
        
        blue.extend(oct_order(0, False, 11, 18, 20, 14, 12, 5, 3, 9))
        blue.extend(oct_order(0, False, 19, 13, 15, 2, 4, 10, 8, 21))
    else:
        green_pl = [
            [2, 0, 23, 21, 17, 8, 15, 6],
            [3, 1, 22, 20, 7, 14, 9, 16],
        ]
        red_pl = [
            [0, None, 1, None, 2, 3, 4, 5],
            [None, 10, 11, 8, None, 9, 6, 7],
            [16, 17, 14, 15, 12, None, 13, None],
            [18, 19, 20, None, None, 21, 22, 23]
        ]
        blue_pl = [
            [5, 12, 3, 14, 9, 20, 11, 18],
            [10, 4, 8, 2, 21, 15, 19, 13]
        ]
        def make_checks(order):
            z_order = []
            x_order = []
            for i in range(16):
                if i < 8:
                    z_order.append(order[i])
                    x_order.append(None)
                else:
                    z_order.append(None)
                    x_order.append(order[i-8])
            return z_order, x_order
        for order in red_pl:
            red.extend(make_checks(order))
        for order in green_pl:
            green.extend(make_checks(order))
        for order in blue_pl:
            blue.extend(make_checks(order))

    cid = 24

    data_qubits = list(range(24))
    parity_qubits = []
    flag_qubits = []

    xlist = []
    zlist = []

    flag_dict = {}

    color_map = {}
    all_checks = [(red, 0), (green, 1), (blue, 2)]
    for (check_list, c) in all_checks:
        for (k, arr) in enumerate(check_list):
            if k & 0x1:
                xlist.append(cid)
            else:
                zlist.append(cid)
            parity_qubits.append(cid)
            checks[cid] = arr
            color_map[cid] = c
            # Figure out what flags we need to introduce.
            flag_dict[cid] = {}
            fid = cid+1
            all_fids = []

            if flags:
                w = len([x for x in arr if x is not None])
                for i in range(0, w, 2):
                    q1 = arr[i]
                    q2 = arr[i+1]
                    flag_dict[cid][q1] = fid
                    flag_dict[cid][q2] = fid
                    flag_qubits.append(fid)
                    all_fids.append(fid)
                    color_map[fid] = c
                    fid += 1
            flag_dict[cid]['all'] = all_fids
            cid = fid
    x_obs = [
                [20, 22, 7, 9],
                [15, 17, 21, 23],
                [0, 2, 6, 8],
                [1, 3, 14, 16],
                [2, 4, 8, 10],
                [3, 5, 12, 14],
                [9, 11, 18, 20],
                [13, 15, 19, 21]
            ]
    z_obs = x_obs

    return data_qubits, parity_qubits, flag_qubits, checks, zlist, xlist, z_obs, x_obs,\
            {'color': color_map, 'flag_dict': flag_dict}

def psc_d3():
    data_qubits = list(range(8))
    parity_qubits = list(range(8, 14))

    z_check_list = [
        [0, 1, None, None, None, None],
        [None, None, None, None, 7, 6],
        [6, 5, 4, 3, 2, 1]
    ]
    x_check_list = [
        [7, 4, 6, 1, 3, 0],
        [5, 3, None, None, None, None],
        [None, None, None, None, 4, 2]
    ]

    flag_qubits = []
    flag_dict = {}

    checks = {}
    zlist, xlist = [], []
    cid = 8
    fid = 14
    for ch in z_check_list:
        zlist.append(cid)
        checks[cid] = ch
        # Add flags if necessary
        all_fids = []
        flag_dict[cid] = {}
        if len(ch) == 6:
            for i in range(0, len(ch), 2):
                q1, q2 = ch[i], ch[i+1]
                flag_dict[cid][q1] = fid
                flag_dict[cid][q2] = fid
                flag_qubits.append(fid)
                all_fids.append(fid)
                fid += 1
        flag_dict[cid]['all'] = all_fids
        cid += 1

    for ch in x_check_list:
        xlist.append(cid)
        checks[cid] = ch
        # Add flags if necessary
        all_fids = []
        flag_dict[cid] = {}
        if len(ch) == 6:
            for i in range(0, len(ch), 2):
                q1, q2 = ch[i], ch[i+1]
                flag_dict[cid][q1] = fid
                flag_dict[cid][q2] = fid
                flag_qubits.append(fid)
                all_fids.append(fid)
                fid += 1
        flag_dict[cid]['all'] = all_fids
        cid += 1

    z_obs = [
                [0, 3, 5],
                [2, 4, 7]
            ]
    x_obs = [
                [0, 1, 2],
                [5, 6, 7]
            ]

    return data_qubits, parity_qubits, flag_qubits, checks, zlist, xlist, z_obs, x_obs, {}

def check_for_scheduling_conflict(code_data):
    data_qubits, parity_qubits, flag_qubits, checks, zlist, xlist, _, _, _ = code_data
    maxlen = max(len(checks[ch]) for ch in checks)
    conflict_list = []
    for i in range(maxlen):
        visited = set()
        for (ch, arr) in checks.items():
            j = arr[i]
            if j is None:
                continue
            if j in visited:
                conflict_list.append((j, i))
            visited.add(j)
#   print('conflicts:')
#   for (q, t) in conflict_list:
#       print('\tqubit %d at time %d' % (q, t))
    return len(conflict_list) > 0
    
def dump_info_about_code(code_data, only_checks=None):
    data_qubits, parity_qubits, flag_qubits, checks, zlist, xlist, _, _, _ = code_data
    print('checks:')
    for (ch, arr) in checks.items():
        if only_checks is None or ch in only_checks:
            print('\t%d:\t%s' % (ch, '\t'.join(['_' if x is None else str(x) for x in arr])))
    check_for_scheduling_conflict(code_data)

def write_syndrome_extraction_ops(writer, code_data, only_do_checks=None):
    data_qubits, parity_qubits, flag_qubits, checks, zlist, xlist, _, _, suppl = code_data
    if only_do_checks is not None:
        checks = {ch : x for (ch, x) in checks.items() if ch in only_do_checks}
        zlist = [ch for ch in zlist if ch in only_do_checks]
        xlist = [ch for ch in xlist if ch in only_do_checks]
    flag_dict = None if 'flag_dict' not in suppl else suppl['flag_dict']
    # Does all but measurement+reset and events.
    writer.write('@annotation inject_timing_error\n')
    if len(xlist) > 0:
        writer.write('h %s;\n' % concat(xlist))

    if flag_dict is not None:
        zflags = []
        for ch in zlist:
            zflags.extend(flag_dict[ch]['all'])
        if len(zflags) > 0:
            writer.write('h %s;\n' % concat(zflags))
        # Entangle flags with their corresponding parity qubits.
        max_flags = max(len(flag_dict[ch]['all']) for ch in parity_qubits)
        for i in range(max_flags):
            cnots = []
            for ch in parity_qubits:
                if i >= len(flag_dict[ch]['all']):
                    continue
                f = flag_dict[ch]['all'][i]
                if ch in xlist:
                    cnots.append(ch)
                    cnots.append(f)
                else:
                    cnots.append(f)
                    cnots.append(ch)
            if len(cnots) > 0:
                writer.write('cx %s;\n' % concat(cnots))
        
    max_weight = max(len(arr) for (_, arr) in checks.items())
    for i in range(max_weight):
        cnots = []
        for ch in xlist: 
            x = checks[ch][i]
            if x is None:
                continue
            # Check if x has a corresponding flag.
            if flag_dict is not None and x in flag_dict[ch]:
                f = flag_dict[ch][x]
                cnots.append(f)
                cnots.append(x)
            else:
                cnots.append(ch)
                cnots.append(x)
        for ch in zlist:
            x = checks[ch][i]
            if x is None:
                continue
            if flag_dict is not None and x in flag_dict[ch]:
                f = flag_dict[ch][x]
                cnots.append(x)
                cnots.append(f)
            else:
                cnots.append(x)
                cnots.append(ch)
        if len(cnots) > 0:
            writer.write('cx %s;\n' % concat(cnots)) 
    if flag_dict is not None:
        for i in range(max_flags):
            cnots = []
            for ch in parity_qubits:
                if i >= len(flag_dict[ch]['all']):
                    continue
                f = flag_dict[ch]['all'][i]
                if ch in xlist:
                    cnots.append(ch)
                    cnots.append(f)
                else:
                    cnots.append(f)
                    cnots.append(ch)
            if len(cnots) > 0:
                writer.write('cx %s;\n' % concat(cnots))
        if len(zflags) > 0:
            writer.write('h %s;\n' % concat(zflags))
    if len(xlist) > 0:
        writer.write('h %s;\n' % concat(xlist))

def write_asm_file_for_check_sch(code_data, output_file, checks=None, memory='z'):
    data_qubits, _, _, _, _, _, _, _, _ = code_data
    writer = open(output_file, 'w')
    write_syndrome_extraction_ops(writer, code_data, only_do_checks=checks)
    if memory == 'z':
        writer.write('cmpx %s;\n' % concat(data_qubits))
    else:
        writer.write('cmpz %s;\n' % concat(data_qubits))
    writer.close()

def write_asm_file_for_css(code_data, rounds, output_file, memory='z'):
    data_qubits, parity_qubits, flag_qubits, checks, zlist, xlist, z_obs_list, x_obs_list, suppl = code_data

    writer = open(output_file, 'w')
    # Initialization.
    writer.write('\n#\n# PROLOGUE\n#\n\n')

    writer.write('reset %s;\n' % concat(data_qubits))
    writer.write('@annotation no_tick\n')
    writer.write('reset %s;\n' % concat(parity_qubits))
    if flag_qubits:
        writer.write('@annotation no_tick\n')
        writer.write('reset %s;\n' % concat(flag_qubits))
    if memory == 'x':
        writer.write('h %s;\n' % concat(data_qubits))
    # Execute syndrome extraction.
    ectr = 0
    mctr = 0

    zflags, xflags = [], []
    flag_dict = None if 'flag_dict' not in suppl else suppl['flag_dict']
    if flag_dict is not None:
        for ch in xlist:
            xflags.extend(flag_dict[ch]['all'])
        for ch in zlist:
            zflags.extend(flag_dict[ch]['all'])
    all_checks = [*zlist, *xlist, *zflags, *xflags]

    writer.write('\n#\n# SYNDROME EXTRACTION\n#\n\n')
    for r in range(rounds):
        writer.write('# Round %d\n' % r)

        write_syndrome_extraction_ops(writer, code_data)
        
        writer.write('measure %s;\n' % concat(all_checks))
        writer.write('reset %s;\n' % concat(all_checks))
        
        m_per_round = len(all_checks)
        if memory == 'x':
            off = len(zlist)
            for (i, ch) in enumerate(xlist):
                if 'color' in suppl:
                    writer.write('@property color %d\n' % suppl['color'][ch])
                if r == 0:
                    writer.write('event %d, %d;\n' % (ectr, off+i))
                else:
                    writer.write('event %d, %d, %d;\n' % (ectr, mctr+off+i, mctr+off+i-m_per_round))
                ectr += 1
            # Also do zflags.
            off += len(xlist)
            for (i, fl) in enumerate(zflags):
                if 'color' in suppl:
                    writer.write('@property color %d\n' % suppl['color'][fl])
                if r == 0:
                    writer.write('event %d, %d;\n' % (ectr, off+i))
                else:
                    writer.write('event %d, %d;\n' % (ectr, mctr+off+i))
                ectr += 1
        else:
            for (i, ch) in enumerate(zlist):
                if 'color' in suppl:
                    writer.write('@property color %d\n' % suppl['color'][ch])
                if r == 0:
                    writer.write('event %d, %d;\n' % (ectr, i))
                else:
                    writer.write('event %d, %d, %d;\n' % (ectr, mctr+i, mctr+i-m_per_round))
                ectr += 1
            off = len(zlist) + len(xlist) + len(zflags)
            for (i, fl) in enumerate(xflags):
                if 'color' in suppl:
                    writer.write('@property color %d\n' % suppl['color'][fl])
                if r == 0:
                    writer.write('event %d, %d;\n' % (ectr, off+i))
                else:
                    writer.write('event %d, %d;\n' % (ectr, mctr+off+i))
                ectr += 1
        mctr += m_per_round
    # Epilogue
    writer.write('\n#\n# EPILOGUE\n#\n\n')
    if memory == 'x':
        writer.write('h %s;\n' % concat(data_qubits))
    writer.write('measure %s;\n' % concat(data_qubits))
    if memory == 'x':
        for (i, ch) in enumerate(xlist):
            pmeas = mctr - len(xlist) - len(zflags) - len(xflags) + i
            data_meas = [mctr + j for j in checks[ch] if j is not None]
            if 'color' in suppl:
                writer.write('@property color %d\n' % suppl['color'][ch])
            writer.write('event %d, %d, %s;\n' % (ectr, pmeas, concat(data_meas)))
            ectr += 1
    else:
        for (i, ch) in enumerate(zlist):
            pmeas = mctr - len(parity_qubits) - len(zflags) - len(xflags) + i
            data_meas = [mctr + j for j in checks[ch] if j is not None]
            if 'color' in suppl:
                writer.write('@property color %d\n' % suppl['color'][ch])
            writer.write('event %d, %d, %s;\n' % (ectr, pmeas, concat(data_meas)))
            ectr += 1
    obs_list = x_obs_list if memory == 'x' else z_obs_list
    for (i, obs) in enumerate(obs_list):
        data_meas = [mctr + j for j in obs if j is not None]
        writer.write('obs %d, %s;\n' % (i, concat(data_meas)))

if __name__ == '__main__':
    from sys import argv

    if len(argv) > 1 and argv[1] == '-s':
        for gr in range(8):
            for gf in [True, False]:
                code = hycc_d4(gr, gf)
                if not check_for_scheduling_conflict(code):
                    write_asm_file_for_css(code, 4, 'asm/qec/hycc/mzd4_%d_%d.asm' % (gr, gf))
    else:
#       code = hycc_d4(7, 0, True, False)
#       write_asm_file_for_css(code, 4, 'asm/qec/hycc/mzd4_both_at_once_no_flags.asm')
#       code = hycc_d4(7, 0)
#       write_asm_file_for_css(code, 4, 'asm/qec/hycc/mzd4_both_at_once.asm')
#       code = hycc_d4(0, 0, False)
#       write_asm_file_for_css(code, 4, 'asm/qec/hycc/mzd4_one_by_one.asm')
        for d in [3, 5, 7, 9]:
            write_asm_file_for_css(hexcc(d), d, 'asm/qec/hexcc/mzd%d.asm' % d)
