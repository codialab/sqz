#!/usr/bin/env python
from sys import stdout, stderr, exit
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF

import logging

#
# third party libraries
#

# logging
LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

def length_pline(line):
    splits = line.split('\t')
    haptype = splits[1]
    n = splits[2].count(',') + 1
    return f'{haptype}\t{n}'

def length_wline(line):
    splits = line.split('\t')
    haptype = '#'.join(splits[1:4])
    n = splits[6].count('>') + splits[6].count('<')
    return f'{haptype}\t{n}'

def parse_gfa(data, line_type, out):

    for line in data:
        if line.startswith(f'{line_type}\t'):
            if line_type == 'P':
                print(length_pline(line), file=out)
            elif line_type == 'W':
                print(length_wline(line), file=out)
            elif line_type == 'Z':
                print(length_wline(line), file=out)
            else:
                raise NotImplementedError
            
if __name__ == '__main__':

    description = '''count number of nodes in P/W/Z lines'''
    parser = ArgumentParser(formatter_class=ADHF, description=description)
    parser.add_argument('gfa_file', type=open, help='GFA file')
    parser.add_argument('-P', action='store_true', help='parse P-lines')
    parser.add_argument('-Z', action='store_true', help='parse Z-lines')
    parser.add_argument('-W', action='store_true', help='parse W-lines')
    args = parser.parse_args()

    if (args.P and args.Z) or (args.P and args.W) or (args.Z and args.W):
        LOG.error('I refuse to parse P and Z lines simultaneously')
        exit(1)

    if not (args.P or args.Z or args.W):
        LOG.error('You must specify either -P or -Z')
        exit(1)

    line_type = args.P and 'P' or args.Z and 'Z' or 'W'

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(ch)

    parse_gfa(args.gfa_file, line_type, stdout) 
