def __concat(arr, delimiter=', '):
    return delimiter.join(str(x) for x in arr)

def __w_op(writer, arr, opname):
    writer.write('%s %s;\n' % (opname, __concat(arr)))

def __w_property(writer, property_name, property_value):
    writer.write('@property %s %d\n' % (property_name, property_value))

def __w_annotation(writer, annotation_name):
    writer.write('@annotation %s\n' % annotation_name)

def w_newline(writer):
    writer.write('\n')

def w_big_comment(writer, text):
    writer.write('\n#\n# %s \n#\n\n' % text)

def w_comment(writer, text):
    writer.write('# %s\n' % text)

def w_reset(writer, qubits):
    __w_op(writer, qubits, 'reset')

def w_h(writer, qubits):
    __w_op(writer, qubits, 'h')

def w_cx(writer, qubits):
    __w_op(writer, qubits, 'cx')

def w_measure(writer, qubits, mctr_map):
    __w_op(writer, qubits, 'measure')
    for (i, q) in enumerate(qubits):
        mctr_map[q] = mctr_map['curr']+i
    mctr_map['curr'] += len(qubits)

def w_event(writer, event_id, *meas_id_array):
    __w_op(writer, [event_id, *meas_id_array], 'event')

def w_obs(writer, obs_id, *meas_id_array):
    __w_op(writer, [obs_id, *meas_id_array], 'obs')

def w_color(writer, color_id):
    __w_property(writer, 'color', color_id)

def write_color_code_asm(code_data, output_file, rounds, memory='z', flag_protocol=0):
    # Read all the code data.
    all_qubits = code_data['all_qubits']
    data_qubits = code_data['data_qubits']
    parity_qubits = code_data['parity_qubits']
    flag_qubits = code_data['flag_qubits']
    obs_list = code_data['obs_list'][memory]
    color_map = code_data['color_map']
    flag_map = code_data['flag_map']
    support_map = code_data['support_map']

    using_flags = len(flag_qubits) > 0
    # Convention: X stabilizer before Z.
    writer = open(output_file, 'w')
    #
    # Prologue: reset all qubits and initialize data in correct state.
    #
    w_big_comment(writer, 'PROLOGUE')
    w_reset(writer, all_qubits) 
    if memory == 'x':
        w_h(writer, data_qubits)
    #
    # Syndrome extraction round.
    #
    ectr, octr = 0, 0
    mctr_map = {'curr': 0}
    for r in range(rounds):
        prev_mctr_map = mctr_map.copy()

        w_big_comment(writer, 'ROUND %d' % r)
        for part in ['x', 'z']:
            __w_annotation(writer, 'inject_timing_error')
            # First, initialize all the parity qubits in |+>
            w_comment(writer, 'Parity qubit initialization (%s)' % part)
            if part == 'x':
                w_h(writer, parity_qubits)
            # Syndrome extraction operations are different whether we use flags or not.
            if using_flags:
                w_comment(writer, 'Flag qubit initialization (%s)' % part)
                # Now, initialize all flag qubits.
                if part == 'z':
                    w_h(writer, flag_qubits)
                fi = 0  # fi tracks the depth. We exit the below loop once all flags have been initialized.
                while True:
                    cx_list = []
                    for pq in parity_qubits:
                        flags = flag_map[pq]['all']
                        if fi >= len(flags):
                            continue
                        fq = flags[fi]
                        if part == 'x':
                            cx_list.extend([pq, fq])
                        else:
                            cx_list.extend([fq, pq])
                    if len(cx_list) == 0:
                        break
                    w_cx(writer, cx_list)
                    fi += 1
            # Now, we perform the CX operations.
            di = 0 # Like with the flag qubits, we do this iteratively until all CX gates are done.
            w_comment(writer, 'Stabilizer CX (%s)' % part)
            while True:
                cx_list = []
                for pq in parity_qubits:
                    if di >= len(support_map[pq]):
                        continue
                    dq = support_map[pq][di]
                    if dq is None:
                        continue
                    if part == 'x':
                        if using_flags:
                            fq = flag_map[pq][dq]
                            cx_list.extend([fq, dq])
                        else:
                            cx_list.extend([pq, dq])
                    else:
                        if using_flags:
                            fq = flag_map[pq][dq]
                            cx_list.extend([dq, fq])
                        else:
                            cx_list.extend([dq, pq])
                if len(cx_list) == 0:
                    break
                w_cx(writer, cx_list)
                di += 1
            # Now, we need to teardown the flags.
            if using_flags:
                w_comment(writer, "Flag teardown (%s)" % part)
                # We just repeat the initialization, but in reverse order.
                fi = 0
                while True:
                    cx_list = []
                    for pq in parity_qubits:
                        flags = flag_map[pq]['all']
                        if fi >= len(flags):
                            continue
                        fq = flags[fi]
                        if part == 'x':
                            cx_list.extend([pq, fq])
                        else:
                            cx_list.extend([fq, pq])
                    if len(cx_list) == 0:
                        break
                    w_cx(writer, cx_list)
                    fi += 1
                if part == 'z':
                    w_h(writer, flag_qubits)
                w_measure(writer, flag_qubits, mctr_map)
                w_reset(writer, flag_qubits)
                # Depending on the flag protocol, do something with the flags.
                w_comment(writer, 'Flag protocol (%s)' % part)
                if flag_protocol == 0:
                    # Make detection events for the flags.
                    #
                    # Note that we only care about flags during the "opposite" stabilizer measurements.
                    if memory != part:
                        for fq in flag_qubits:
                            w_color(writer, color_map[fq])
                            w_event(writer, ectr, mctr_map[fq])
                            ectr += 1
                elif flag_protocol == 1:
                    # Conditionally flip the data qubits.
                    clop_list = []
                    for fq in flag_qubits:
                        _, dq = support_map[fq]
                        if dq is None:
                            continue
                        clop_list.extend([mctr_map[fq], dq])
                    if part == 'x':
                        __w_op(writer, clop_list, 'clx')
                    else:
                        __w_op(writer, clop_list, 'clz')
            # Finally, measure the checks.
            w_comment(writer, 'Stabilizer measurement (%s)' % part)
            if part == 'x':
                w_h(writer, parity_qubits)
            w_measure(writer, parity_qubits, mctr_map)
            w_reset(writer, parity_qubits)
            if memory == part:
                # Create detection events.
                for pq in parity_qubits:
                    w_color(writer, color_map[pq])
                    if r == 0:
                        w_event(writer, ectr, mctr_map[pq])
                    else:
                        w_event(writer, ectr, mctr_map[pq], prev_mctr_map[pq])
                    ectr += 1
    #
    # Epilogue
    #
    w_big_comment(writer, 'EPILOGUE')
    # Measure all data qubits
    if memory == 'x':
        w_h(writer, data_qubits)
    w_measure(writer, data_qubits, mctr_map)
    # Now, we create the final detection events and the logical observables.
    for pq in parity_qubits:
        meas_list = [mctr_map[pq]]
        for dq in support_map[pq]:
            if dq is None:
                continue
            meas_list.append(mctr_map[dq])
        w_color(writer, color_map[pq])
        w_event(writer, ectr, *meas_list)
        ectr += 1
    for obs in obs_list:
        meas_list = [mctr_map[dq] for dq in obs]
        w_obs(writer, octr, *meas_list)
        octr += 1
    writer.close()

def hexcc(d, using_flags=True):
    offset = 2
    side_len = (3*d - 1) // 2

    loc = {}
    check_locs = []
    x_obs = []
    z_obs = []
    n = 0
    # As an initial pass, we will go through the structure
    # and mark the location of each qubit.
    for r in range(side_len):
        row_offset = offset
        for c in range(r+1):
            if row_offset == 0:
                check_locs.append((r, c))
            else:
                loc[(r, c)] = n
                if c == 0:  # This is the left edge of the triangle.
                    x_obs.append(n)
                    z_obs.append(n)
                n += 1
            row_offset = (row_offset+1) % 3
        offset = (offset+1) % 3
    # Now, compile the data qubits, parity qubits, and flag qubits.
    data_qubits = list(range(n))
    parity_qubits = []
    flag_qubits = []
    # There will be one parity qubit per plaquette.
    flag_map = {}
    support_map = {}
    color_map = {}

    get_loc = lambda _r, _c: loc[(_r, _c)] if (_r, _c) in loc else None

    for (i, j) in check_locs:
        pq = n
        # Compute support.
        #
        #     a   b           a   b
        #   f   P   c    -->  f   P   c
        #     e   d               e   d
        #
        a = get_loc(i-1, j-1)
        b = get_loc(i-1, j)
        c = get_loc(i, j+1)
        d = get_loc(i+1, j+1)
        e = get_loc(i+1, j)
        f = get_loc(i, j-1)
        n += 1
        # Populate data structures.
        # 
        # Flags:
        #   (b, c) have a flag, (e, d) have a flag, and (a, f) have a flag.
        # CX Order: b, c, e, d, a, f
        color_map[pq] = i % 3   # r = 0, g = 1, b = 2
        support_map[pq] = [b, c, e, d, a, f]
        parity_qubits.append(pq)

        if using_flags:
            flag_pair_array = [[b, c], [e, d], [a, f]]
            flag_map[pq] = { 'all': [] }
            for (q1, q2) in flag_pair_array:
                fq = n
                if q1 is None or q2 is None:
                    continue
                color_map[fq] = i % 3
                flag_map[pq]['all'].append(fq)
                flag_map[pq][q1] = fq
                flag_map[pq][q2] = fq
                support_map[fq] = [q1, q2]
                flag_qubits.append(fq)
                n += 1
    all_qubits = [*data_qubits, *parity_qubits, *flag_qubits]
    code_data = {
        'all_qubits': all_qubits,
        'data_qubits': data_qubits,
        'parity_qubits': parity_qubits,
        'flag_qubits': flag_qubits,
        'obs_list': {
            'x': [x_obs],
            'z': [z_obs]
        },
        'color_map': color_map,
        'flag_map': flag_map,
        'support_map': support_map
    }
    return code_data

if __name__ == '__main__':
    from sys import argv
    output_file = argv[1]
    d = int(argv[2])
    r = int(argv[3])
    write_color_code_asm(hexcc(d, True), output_file, r, memory='z', flag_protocol=1)

