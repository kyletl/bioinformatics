#!/usr/bin/env python
import sys
from sys import argv
from FASTAReader import *
from HMM import AlignedProfileHMMInit
from HMM import ProfileHMM


VOCAB = set("ARNDCQEGHILKMFPSTWYV")

def formatHelper(helper):

	for v in VOCAB:
		sys.stdout.write('    '+v+'    ')
	print '\n   m->m    m->i    m->d    i->m    i->i    d->m    d->d'

	for i in range(len(helper.a)-1):
		sys.stdout.write(str(i))
		for v in VOCAB:
			sys.stdout.write(' '+str(round(helper.b[i]['I'][v], 4))+' ')
		print str(i)+'\n        '
		if 'M' in helper.b[i].keys():
			for v in VOCAB:
				sys.stdout.write(' '+str(round(helper.b[i]['M'][v], 4)))
			print '\n'
		if 'M' in helper.a[i].keys():
			sys.stdout.write('     '+str(round(helper.a[i]['M']['M'], 4)))
			sys.stdout.write(' '+str(round(helper.a[i]['M']['I'], 4)))
			sys.stdout.write(' '+str(round(helper.a[i]['M']['D'], 4)))
			sys.stdout.write(' '+str(round(helper.a[i]['I']['M'], 4)))
			sys.stdout.write(' '+str(round(helper.a[i]['I']['I'], 4)))
			sys.stdout.write(' '+str(round(helper.a[i]['D']['M'], 4)))
			sys.stdout.write(' '+str(round(helper.a[i]['D']['D'], 4)))
			print '\n'
# possible commands: init, train
if __name__ == '__main__':
	fname = argv[2]
	command = argv[1]
	seqs = FASTAFile(fname)

	if command == "init":
		helper = AlignedProfileHMMInit(seqs)
		for h in helper.a:
			print h

	elif command == "train":
		helper = AlignedProfileHMMInit(seqs)
		phmm = ProfileHMM(helper)
		# print phmm.a
		for a in phmm.a:
			print a
		for b in phmm.b:
			print b
