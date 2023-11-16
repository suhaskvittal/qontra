from asm_gen_base import *

def hexcc(d):
    checks = defaultdict(list)

    offset = 2
    side_len = (3*d - 1) // 2

    loc = {}
    check_locs = []
    x_obs = []
    z_obs = []
    k = 0
    # As an initial pass, we will go through the structure
    # and mark the location of each qubit.
    for r in range(side_len):
        row_offset = offset
        for c in range(r+1):
            if row_offset == 0:
                check_locs.append((r, c))
            else:
                loc[(r, c)] = k 
                if c == 0:  # This is the left edge of the triangle.
                    x_obs.append(k)
                    z_obs.append(k)
                k += 1
            row_offset = (row_offset+1) % 3
        offset = (offset+1) % 3
    # Now, go through and create the remaining structures.
    data_qubits = list(range(k))
    parity_qubits = []
    flag_qubits = []
    checks = {}
    zlist, xlist = [], []

    get_loc = lambda _r, _c: loc[(_r, _c)] if (_r, _c) in loc else None

    print(check_locs)

    red_ctr, blue_ctr, green_ctr = 0, 0, 0
    color_map = {}
    print(check_locs)
    for (i, j) in check_locs:
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
        if i % 3 == 1:
            chz, chx = k + 1 + 3*green_ctr, k + 1 + 3*(green_ctr+1)
            green_ctr += 2
        elif i % 3 == 2:
            chz, chx = k + 2 + 3*blue_ctr, k + 2 + 3*(blue_ctr+1)
            blue_ctr += 2
        else:
            chz, chx = k + 3*red_ctr, k + 3*(red_ctr+1)
            red_ctr += 2
        color_map[chz] = i % 3  # r = 0, g = 1, b = 2
        color_map[chx] = i % 3

        checks[chz] = [b, c, d, a, f, e, None]
        checks[chx] = [None, b, a, f, c, d, e]
        parity_qubits.extend([chz, chx])
        zlist.append(chz)
        xlist.append(chx)
    return data_qubits, parity_qubits, flag_qubits, checks, zlist, xlist, [z_obs], [x_obs], {'color': color_map}

if __name__ == '__main__':
    from sys import argv
    output_file = argv[1]
    d = int(argv[2])
    r = int(argv[3])
    write_asm_file_for_css(hexcc(d), r, output_file)
