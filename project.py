import sys
from sys import argv
from FASTAReader import *

if __name__ == '__main__':
    fname = argv[1]
    seqs = FASTAFile(fname)
    for seq in seqs:
        print seq.name
        for s in seq:
            sys.stdout.write(s)
        print '\n'
