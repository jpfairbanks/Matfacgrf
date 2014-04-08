from __future__ import print_function
import argparse
import sys
from random import random as rand

def getargs():
    """getargs: for handling command line parameters"""
    parser = argparse.ArgumentParser(
        description="Randomly select lines from an mtx file.\
            The line counts are not adjusted.")
    parser.add_argument('-p',
                       help='the probabilty of selecting each line',
                       type=float,
                       metavar='P',
                       required=True,
                       dest='p')
    parser.add_argument('file',
                        help="the file containing edges; - uses stdin.",
                        metavar='graph.mtx',
                        )
    args = parser.parse_args()
    return args

def writeln(fp, line):
    fp.write(line)

if __name__ == '__main__':
    args = getargs()
    p = args.p
    sizeinfo = False
    outcount = 0
    errcount = 0
    #open the file for processing.
    if args.file == "-":
        fp = sys.stdin
    else:
        fp = open(args.file, 'r')
    for line in fp:
        #cludge processing MTX files.
        if line.startswith('%'): sys.stderr.write(line)
        elif not sizeinfo:
            sizeinfo = True
            words = line.split(" ")
            sys.stderr.write("%{0}".format(line))
            sys.stdout.write("%{0}".format(line))
            m, n, numedges = words[0], words[1], words[2]
        else:
            #processing the data lines.
            r = rand()
            if r < p:
                outcount += 1
                writeln(sys.stdout, line)
            else: #r >=p:
                errcount += 1
                writeln(sys.stderr, line)
    print("%linecount:{0}".format(errcount), file=sys.stderr)
    print("%linecount:{0}".format(outcount), file=sys.stdout)
