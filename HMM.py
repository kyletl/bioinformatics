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
import copy

# constants
GAP = '-'
B = 'B'
D = 'D'
I = 'I'
M = 'M'
N = 'N'
PSEUDOCOUNT = 100
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


"""A class to define a profile-Hidden Markov Model.

Input:
	states - set of possible states
	transitions - sources: their possible destination set
	vocab - set of all possible emissions
	start (optional) - start state, must be in states
"""

class ProfileHMM:

	def __init__(self, hmm):
		# self.length = length
		self.vocab = VOCAB

		# a[t][i][j] 
		# i.e. chance of transition from state i at node t to state j
		self.a = copy.deepcopy(hmm.a)

		# b[t][i][e] = chances of emission e in state i at node t
		self.b = copy.deepcopy(hmm.b)

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

			elif pos > self.length + 1:
				return sources

			elif pos < 0:
				return sources
			
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
		assert pos <= self.length
		sources = set()

		if (state is not B or pos == 0) and state is not N:
			if pos == self.length:
					end = (self.length + 1, N)
					sources.add(end)
					return sources

			elif pos > self.length:
				return sources

			elif pos < 0:
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

	def check_forward(self, sequence):
		touched = set()
		alpha = [{}]

		# init data structure for all states
		alpha[0][(0, B)] = 1
		alpha[0][(0, I)] = 0
		for state in BODY_STATES:
			for k in range(1, self.length + 1)
				alpha[0][(k, state)] = 0
		alpha[0][(self.length + 1, N)] = 0

		# record all possible states, represented as tuples (node_number, state_type)
		all_states = set(alpha[0].values())

		for t in range(1, len(sequence)+1)

# 		self.vocab_size = len(vocab)

# 		# states = set('DIM')
# 		# self.states = [set('I')]
# 		# self.states.extend([set('DIM')] * length) ## assume states is stored as states[column][set(states)]
# 		# self.transitions = transitions
# 		# self.num_states = len(states)

# 		# # pi[i] = chances of starting with state i
# 		# self.pi = {state: 1.0 / self.num_states for state in self.states}

# 		# # a[i][j] = "time-independent stochastic transition matrix"
# 		# # i.e. chance of transition from state i to state j
# 		# self.a = {
# 		# 	in_state: {
# 		# 		out_state: 1.0 / len(self.transitions[in_state]) 
# 		# 		for out_state in self.transitions[in_state]
# 		# 	} for in_state in self.transitions
# 		# }

# 		# # b[i][j] = chances of emission j in state i
# 		# if emits is not None:
# 		# 	self.b = {
# 		# 		emit_state: {
# 		# 			emit: 1.0 / self.vocab_size for emit in self.vocab
# 		# 		} for emit_state in emits
# 		# 	}
# 		# else:
# 		# 	self.b = {
# 		# 		state: {
# 		# 			emit: 1.0 / self.vocab_size for emit in self.vocab
# 		# 		} for state in self.states
# 		# 	}

# 	def converge(self):
# 		pass

	"""Function to scan forward through a sequence and find probabilities."""
	def eval_forward(self, sequence):
		pass

	# 	width = len(sequence)
	# 	result = {state: [0] * width for state in self.states}

	# 	# init time 0
	# 	for state in self.states:
	# 		if state in self.emits:
	# 			result[state][0] = self.pi[state] * self.b[state][sequence[0]]
	# 		else:
	# 			result[state][0] = self.pi[state]

	# 	# induction from base case
	# 	for t in range(1, width):
	# 		# calculate prob of each destination state
	# 		for dest in self.states:
	# 			# for every source state in previous time
	# 			for source in self.states:
	# 				# add prob of transitioning into dest state
	# 				result[dest][t] += result[source][t-1] * self.a[source][dest]
	# 			# scale by emission probability if needed
	# 			if dest in self.emits:
	# 				result[dest][t] *= self.b[dest][sequence[t]]

	# 	# result holds prob of each state/time given prefix
	# 	return result

	"""Function to scan backward through a sequence and find probabilities."""
	def eval_backward(self, sequence):
		pass
		# width = len(sequence)
		# result = {state: [0] * width for state in self.states}

		# # init final time
		# for state in self.states[-1]:
		# 	result[state][width-1] = 1

		# result[I][width-1] = 1

		# # induction from base case
		# for t in range(width-2, -1, -1):
		# 	# check delete state

		# 	# check match state

		# 	# check insert state

		# 	# calculate prob of each source state
		# 	for source in self.states:
		# 		# for every destination state
		# 		for dest in self.states:
		# 			# add prob of transitioning from source state
		# 			prob = result[dest][t+1] * self.a[source][dest]
		# 			if dest in self.emits:
		# 				prob *= self.b[dest][sequence[t+1]]
		# 			result[source][t] += prob

		# # result holds prob of each state/time given suffix
		# return result

class AlignedProfileHMMInit:

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
		a = [BEGINNING]
		a.extend([TRANSITIONS] * (model_width - 2))
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

		for node in range(model_width):
			# normalize transition matrix
			for source in a[node]:
				norm = float(sum(a[node][source].values()))
				for dest in a[node][source]:
					a[node][source][dest] = a[node][source][dest] / norm

			# normalize emission matrix
			for state in b[node]:
				norm = float(sum(b[node][state].values()))
				for emit in b[node][state]:
					b[node][state][emit] = b[node][state][emit] / norm

		# hmm = ProfileHMM(width, vocab)

		self.a = a
		self.b = b
		# self.hmm = ProfileHMM()

	@staticmethod
	def _add_count(ledger, addition):
		for bucket in ledger:
			ledger[bucket] += addition[bucket]
		return ledger
