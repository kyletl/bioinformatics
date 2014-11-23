#!/usr/bin/env python
from sys import argv
from FASTAReader import *
from HMM import AlignedProfileHMM
from HMM import ProfileHMM

# possible commands: init, train
if __name__ == '__main__':
	fname = argv[2]
	command = argv[1]
	seqs = FASTAFile(fname)

	if command == "init":
		helper = AlignedProfileHMMInit(seqs)
		print helper.a
		print helper.b

	elif command == "train":
		helper = AlignedProfileHMMInit(seqs)
		phmm = ProfileHMM(helper)
		print phmm.a
		print phmm.b
