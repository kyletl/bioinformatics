#!/usr/bin/env python
"""
Taylor Cathcart, Cameron Setian, Kyle Tessier-Lavigne
CS 75 Final Project
11/20/14

Implementation of Hidden Markov Model.
"""
 
# imports
import numpy as np
import sys

# constants
GAP = '-'
B = 'B'
D = 'D'
I = 'I'
M = 'M'
N = 'N'
PSEUDOCOUNT = 1
BEGINNING = {
	I: {
		I: PSEUDOCOUNT,
		M: PSEUDOCOUNT,
		D: PSEUDOCOUNT,
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
		I: PSEUDOCOUNT,
	},
	I: {
		I: PSEUDOCOUNT,
		M: PSEUDOCOUNT,
		D: PSEUDOCOUNT,
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


"""A class to define a profile-Hidden Markov Model.

Input:
	states - set of possible states
	transitions - sources: their possible destination set
	vocab - set of all possible emissions
	start (optional) - start state, must be in states
"""

class ProfileHMM:

	def __init__(self, length, vocab):
		self.length = length
		self.vocab = vocab

		self.vocab_size = len(vocab)

		# states = set('DIM')
		# self.states = [set('I')]
		# self.states.extend([set('DIM')] * length) ## assume states is stored as states[column][set(states)]
		# self.transitions = transitions
		# self.num_states = len(states)

		# # pi[i] = chances of starting with state i
		# self.pi = {state: 1.0 / self.num_states for state in self.states}

		# # a[i][j] = "time-independent stochastic transition matrix"
		# # i.e. chance of transition from state i to state j
		# self.a = {
		# 	in_state: {
		# 		out_state: 1.0 / len(self.transitions[in_state]) 
		# 		for out_state in self.transitions[in_state]
		# 	} for in_state in self.transitions
		# }

		# # b[i][j] = chances of emission j in state i
		# if emits is not None:
		# 	self.b = {
		# 		emit_state: {
		# 			emit: 1.0 / self.vocab_size for emit in self.vocab
		# 		} for emit_state in emits
		# 	}
		# else:
		# 	self.b = {
		# 		state: {
		# 			emit: 1.0 / self.vocab_size for emit in self.vocab
		# 		} for state in self.states
		# 	}

	def converge(self):
		pass

	"""Function to scan forward through a sequence and find probabilities."""
	def scan_forward(self, sequence):
		width = len(sequence)
		result = {state: [0] * width for state in self.states}

		# init time 0
		for state in self.states:
			if state in self.emits:
				result[state][0] = self.pi[state] * self.b[state][sequence[0]]
			else:
				result[state][0] = self.pi[state]

		# induction from base case
		for t in range(1, width):
			# calculate prob of each destination state
			for dest in self.states:
				# for every source state in previous time
				for source in self.states:
					# add prob of transitioning into dest state
					result[dest][t] += result[source][t-1] * self.a[source][dest]
				# scale by emission probability if needed
				if dest in self.emits:
					result[dest][t] *= self.b[dest][sequence[t]]

		# result holds prob of each state/time given prefix
		return result

	"""Function to scan backward through a sequence and find probabilities."""
	def scan_backward(self, sequence):
		width = len(sequence)
		result = {state: [0] * width for state in self.states}

		# init final time
		for state in self.states[-1]:
			result[state][width-1] = 1

		result[I][width-1] = 1

		# induction from base case
		# for t in range(width-2, -1, -1):
			# check delete state

			# check match state

			# check insert state

			# calculate prob of each source state
			# for source in self.states:
			# 	# for every destination state
			# 	for dest in self.states:
			# 		# add prob of transitioning from source state
			# 		prob = result[dest][t+1] * self.a[source][dest]
			# 		if dest in self.emits:
			# 			prob *= self.b[dest][sequence[t+1]]
			# 		result[source][t] += prob

		# result holds prob of each state/time given suffix
		return result

class ProfileHMMHelper:

	def __init__(self, profiles):
		vocab = VOCAB
		# initialize model parameters data structures
		empty_count = {emit: PSEUDOCOUNT for emit in vocab}
		no_count = {emit: 0 for emit in vocab}
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

		# b[t][i][j] = chances of emission j in state i at node t
		b = [{I: dict(empty_count)}]

		# iterate through sequence to construct emission matrix using counts
		node = 0
		col = 0
		while col < full_width:
			# for insert state, add count to current node insert state
			if col_status[col] is I:
				ProfileHMMHelper._add_count(b[node][I], col_count[col])
				col += 1
			# otherwise, create a new node and add to its match state
			else:
				b.append({I: dict(empty_count), M: dict(empty_count)})
				node += 1
				ProfileHMMHelper._add_count(b[node][M], col_count[col])
				col += 1

		### TRANSITION MATRIX ###

		# a[t][i][j] 
		# i.e. chance of transition from state i at node t to state j
		a = [BEGINNING]
		a.extend([TRANSITIONS] * (len(b) - 2))
		a.append(END)

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

		# normI = sum(b[t][I].values())
		# normM = sum(b[t][M].values())
		# for emit in vocab:
		# 	b[t][I][emit] = b[t][I][emit] / normI
		# 	b[t][M][emit] = b[t][M][emit] / normM

		# hmm = ProfileHMM(width, vocab)

		self.a = a
		self.b = b


	@staticmethod
	def _add_count(ledger, addition):
		for bucket in ledger:
			ledger[bucket] += addition[bucket]
		return ledger
