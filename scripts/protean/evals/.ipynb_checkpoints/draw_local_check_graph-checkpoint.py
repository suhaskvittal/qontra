import graphviz as gv
import networkx as nx

FOLDER = 'data/protean/hycc_d8_ex'

def split_csv_line(line: str, delim=',') -> list[str]:
    dat = []
    in_quote = False
    curr = []
    for c in line:
        if c == '\"':
            in_quote = not in_quote
        elif c == ',':
            if in_quote:
                curr.append(c)
            else:
                dat.append(''.join(curr))
        else:
            curr.append(c)
    return curr

def get_detector_to_check_map() -> dict[int, str]:
    qes_file = f'{FOLDER}/memory/zrm1.qes'
    with open(qes_file, 'r') as reader:
        lines = reader.readlines()
    # Go through lines and map detectors to checks
    d2c = {}    # detector to check
    check = None
    color = None
    for line in lines:
        if check is None:
            if '@annotation' not in line:
                continue
            dat = line.split(' ')
            name = dat[-1]
            if name[0] != 'x' and name[0] != 'z':
                continue
            check = name.strip()
        elif color is none:
            if '@property color' not in line:
                continue
            dat = line.split(' ')
            color = int(dat[-1]) 
        else:
            if 'event' not in line:
                continue
            dat = line.split(' ')
            args = dat[-1].split(',')
            eventno = int(args[0])
            d2c[eventno] = (check, color)
            check = None
    return d2c

def get_support_map(checks: list[str]) -> dict[str, list[int]]:
    tg_file = f'{FOLDER}/tanner_graph.txt'
    with open(tg_file, 'r') as reader:
        lines = reader.readlines()[1:] # Ignore header.
    support_map = {}
    for line in lines:
        dat = split_csv_line(line)
        check = dat[0]
        if check not in checks:
            continue
        support_map[check] = set(dat[2].split(','))
    return support_map

def make_interaction_graph(checks: list[str]) -> nx.Graph:
    support_map = get_support_map(checks)

    gr = nx.Graph()
    for (check, color) in checks:
        gr.add_node(check, color=color)
    # Add edges. Label them by common data qubits.
    for (i, (ch1, _)) in enumerate(checks):
        supp1 = support_map[ch1]
        for (j, (ch2, _)) in enumerate(checks):
            if i >= j:
                continue
            supp2 = support_map[ch2]
            comm = supp1 & supp2
            if len(comm) > 0:
                gr.add_edge(ch1, ch2, data=comm)
    return gr

def draw_interaction_graph(gr: nx.Graph, output_file: str):
    g = gv.Graph(filename=output_file)


from sys import argv

detectors = [int(d) for d in argv[1:]]

d2c = get_detector_to_check_map()
checks = [d2c[d] for d in detectors]

print(f'Checks: {checks}')
