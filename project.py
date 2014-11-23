#!/usr/bin/env python
from sys import argv
from FASTAReader import *
from HMM import ProfileHMMHelper

if __name__ == '__main__':
    fname = argv[1]
    seqs = FASTAFile(fname)
    helper = ProfileHMMHelper(seqs)
    print helper.a
    print helper.b