#!/usr/bin/env python
"""
Taylor Cathcart, Cameron Setian, Kyle Tessier-Lavigne
CS 75 Final Project
11/20/14

Implementation of Hidden Markov Model.
"""
 
# imports
import math
import numpy as np
import sys
import copy

# constants
GAP = '-'
B = 'B'
D = 'D'
I = 'I'
M = 'M'
N = 'N'
PSEUDOCOUNT = 10
INS_DEL_SCALE = 0.2
BEGINNING = {
	I: {
		I: PSEUDOCOUNT,
		M: PSEUDOCOUNT,
		D: PSEUDOCOUNT * INS_DEL_SCALE,
	},
	B: {
		I: PSEUDOCOUNT,
		M: PSEUDOCOUNT,
		D: PSEUDOCOUNT,
	}
}
TRANSITIONS = {
	D: {
		D: PSEUDOCOUNT,
		M: PSEUDOCOUNT,
		I: PSEUDOCOUNT * INS_DEL_SCALE,
	},
	I: {
		I: PSEUDOCOUNT,
		M: PSEUDOCOUNT,
		D: PSEUDOCOUNT * INS_DEL_SCALE,
	},
	M: {
		I: PSEUDOCOUNT,
		M: PSEUDOCOUNT,
		D: PSEUDOCOUNT,
	}
}
END = {
    D: {
        N: PSEUDOCOUNT,
        I: PSEUDOCOUNT,
    },
    I: {
        I: PSEUDOCOUNT,
        N: PSEUDOCOUNT
    },
    M: {
        I: PSEUDOCOUNT,
        N: PSEUDOCOUNT
    },
}
VOCAB = set("ARNDCQEGHILKMFPSTWYV")
STATES = set([B, I, D, M, N])
BODY_STATES = set([I, D, M])
EMITS = set([I, M])


"""A class to define a profile-Hidden Markov Model.

Input:
    states - set of possible states
    transitions - sources: their possible destination set
    vocab - set of all possible emissions
    start (optional) - start state, must be in states
"""

class ProfileHMM:

	def __init__(self, hmm):
		assert len(hmm.a) == len(hmm.b)
		self.length = len(hmm.a) - 1
		self.vocab = VOCAB

		# a[t][i][j] 
		# i.e. chance of transition from state i at node t to state j
		self.a = copy.deepcopy(hmm.a)

		# b[t][i][e] = chances of emission e in state i at node t
		self.b = copy.deepcopy(hmm.b)

	def probability(self, seq):
		sequence = self.clean_sequence(seq)
		# alpha = self.check_forward(sequence)
		beta = self.check_backward(sequence)
		# print alpha[len(sequence)][(self.length + 1, N)] 
		# print beta[0][(0, B)]
		return beta[0][(0, B)]

	def train_once(self, sequence):
		pass

	def check_forward(self, sequence):
		alpha = [{} for t in range(len(sequence) + 1)] # alpha[t][state_tuple] = prob of being in particular state at time t

		# init data structure/base case
		alpha[0][(0, B)] = 1
		alpha[0][(0, I)] = 0
		for state in BODY_STATES:
			for k in range(1, self.length + 1):
				alpha[0][(k, state)] = 0
		alpha[0][(self.length + 1, N)] = 0

		# record all possible states, represented as tuples (node_number, state_type)
		all_states = self.get_all_states()
		all_states.reverse()

		# for increasing time values, find probability of particular state
		for t in range(1, len(sequence)+1):
			for (d, dest_state) in all_states:
				# get raw probability based on source states
				prob = 0
				for (s, src_state) in self.source_states(d, dest_state):
					prob += self.a[s][src_state][dest_state] * alpha[t-1][(s, src_state)]
				# scale probability by emission probability if needed
				if dest_state in EMITS:
					prob *= self.b[d][dest_state][sequence[t-1]] # adjust for sequence 0-index

				alpha[t][(d, dest_state)] = prob

		return alpha

	def check_backward(self, sequence):
		beta = [{} for i in range(len(sequence)+1)]

		# init data structures/base case
		all_states = self.get_all_states()
		for t, state in all_states:
			beta[len(sequence)][(t, state)] = 0
		beta[len(sequence)][(self.length+1, N)] = 1
		beta[len(sequence)][(self.length, I)] = self.a[self.length][I][N]
		beta[len(sequence)][(self.length, M)] = self.a[self.length][M][N]
		beta[len(sequence)][(self.length, D)] = self.a[self.length][D][N]

		# for decreasing time values, find probability of state given suffix
		for t in range(len(sequence)-1, -1, -1):
			for (s, src_state) in all_states:
				prob = 0
				# take a sigma over all destination states from this one
				for (d, dest_state) in self.dest_states(s, src_state):
					sigma = self.a[s][src_state][dest_state]
					# scale by emission probability for non-silent states
					if dest_state in EMITS:
						sigma *= beta[t+1][(d, dest_state)] * self.b[d][dest_state][sequence[t]]
					else:
						sigma *= beta[t][(d, dest_state)]
					prob += sigma

				beta[t][(s, src_state)] = prob

		return beta

	def source_states(self, pos, state):
		"""Return the set of all source states, in the form of coordinates, 
		to the specified state.
		"""
		assert state in STATES
		assert pos <= self.length or (pos == self.length + 1 and state is N)
		sources = set()

		if (state is not N or pos == self.length + 1) and state is not B:
			if pos == 0:
				if state is I:
					bgn = (0, B)
					sources.add(bgn)
					return sources
				else:
					return sources

			elif pos == 1 and state is not I:
				sources.add((0, B))
				sources.add((0, I))
				return sources

			elif pos > self.length + 1 or pos < 0:
				return source
			
			elif state is not I:
				src_pos = pos - 1

			else:
				src_pos = pos
			
			d_src = (src_pos, D)
			m_src = (src_pos, M)
			i_src = (src_pos, I)
			sources.add(d_src)
			sources.add(m_src)
			sources.add(i_src)

			return sources

		return sources

	def dest_states(self, pos, state):
		"""Inverse of source_states function."""
		assert state in STATES
		assert pos <= self.length or (pos == self.length + 1 and state is N)
		sources = set()

		if (state is not B or pos == 0) and state is not N:
			if pos == self.length:
				sources.add((self.length, I))					
				sources.add((self.length + 1, N))
				return sources

			if pos == 0:
				sources.add((0, I))
				sources.add((1, M))
				sources.add((1, D))
				return sources

			elif pos > self.length or pos < 0:
				return sources

			dest_pos = pos + 1
			d_dest = (dest_pos, D)
			m_dest = (dest_pos, M)
			i_dest = (pos, I)
			sources.add(d_dest)
			sources.add(m_dest)
			sources.add(i_dest)

			return sources

		return sources

	def get_all_states(self):
		"""Return list of states in reverse topological order."""
		ret = []
		ret.append((0, B))
		ret.append((0, I))
		for k in range(1, self.length + 1):
			ret.append((k, M))
			ret.append((k, D))
			ret.append((k, I))
		ret.append((self.length + 1, N))
		ret.reverse()
		return ret

	def clean_sequence(self, sequence):
		# return sequence
		new = ""
		for char in sequence:
			if char is not '-':
				new += char
		return new


class AlignedProfileHMMInit:

	def __init__(self, profiles):
		vocab = VOCAB
		# initialize model parameters data structures
		empty_count = {emit: PSEUDOCOUNT for emit in vocab}
		# empty_count[GAP] = 0
		no_count = {emit: 0 for emit in vocab}
		# no_count[GAP] = 0
		full_width = int(np.mean([len(p) for p in profiles]))

		# pi[i] = chances of starting with state i

		# init temporary data structures
		col_status = [None] * full_width
		col_count = [None] * full_width
		threshold = len(profiles) / 2

		# analyze input sequence to determine each column's state
		for t in range(full_width):
			count = dict(no_count)
			gaps = 0
			# check every profile
			for profile in profiles:
				emit = profile[t]
				if emit is GAP:
					gaps += 1
				else:
					count[emit] += 1
			# if more than half gaps, use an insert state
			if gaps > threshold:
				col_status[t] = I
			# otherwise it's a match state
			else:
				col_status[t] = M
			col_count[t] = count

		### EMISSION MATRIX ###

		# b[t][i][e] = chances of emission e in state i at node t
		b = [{I: dict(empty_count)}]

		# iterate through sequence to construct emission matrix using counts
		node = 0
		col = 0
		while col < full_width:
			# for insert state, add count to current node insert state
			if col_status[col] is I:
				AlignedProfileHMMInit._add_count(b[node][I], col_count[col])
				col += 1
			# otherwise, create a new node and add to its match state
			else:
				b.append({I: dict(empty_count), M: dict(empty_count)})
				node += 1
				AlignedProfileHMMInit._add_count(b[node][M], col_count[col])
				col += 1

		model_width = len(b)

		### TRANSITION MATRIX ###

		# a[t][i][j] 
		# i.e. chance of transition from state i at node t to state j
		a = [copy.deepcopy(BEGINNING)]
		a.extend([copy.deepcopy(TRANSITIONS) for i in range(model_width - 2)])
		a.append(copy.deepcopy(END))

		# iterate through sequence to construct transition matrix
		# scan through each sequence
		for profile in profiles:
			node = 0
			col = 0
			last_state = B
			while col < full_width:
				letter = profile[col]
				state = col_status[col]
				# if we gap on a match, record delete
				if state is M:
					if letter is GAP:
						a[node][last_state][D] += 1
						last_state = D
					else:
						a[node][last_state][M] += 1
						last_state = M
					node += 1
				elif state is I:
					if letter is not GAP:
						a[node][last_state][I] += 1
						last_state = I
				col += 1

		### NORMALIZATION ###

		for node in range(model_width):
			# normalize transition matrix
			for source in a[node]:
				norm = float(sum(a[node][source].values()))
				for dest in a[node][source]:
					a[node][source][dest] = math.log(a[node][source][dest] / norm)

			# normalize emission matrix
			for state in b[node]:
				norm = float(sum(b[node][state].values()))
				for emit in b[node][state]:
					b[node][state][emit] = math.log(b[node][state][emit] / norm)

		# hmm = ProfileHMM(width, vocab)
		self.a = a
		self.b = b
		# self.hmm = ProfileHMM()

	@staticmethod
	def _add_count(ledger, addition):
		for bucket in ledger:
			ledger[bucket] += addition[bucket]
		return ledger
