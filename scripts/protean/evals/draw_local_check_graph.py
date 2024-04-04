


from sys import argv

detectors = [int(d) for d in argv[1:]]

d2c = get_detector_to_check_map()
checks = [d2c[d] for d in detectors]

print(f'Checks: {checks}')
