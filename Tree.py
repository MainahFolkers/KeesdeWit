#!/usr/bin/env python
from Protein_class_MF import *
from visualize_MF import *
import matplotlib.pyplot as plt
from copy import deepcopy
from collections import deque

sprotein = Protein("HHH")
maxdepth = sprotein.n 
path = []
all_direcs = ['r', 'd', 'u', 'l']

def depth_path(protein, depth, maxdepth):
	avail_direcs = all_direcs
	if depth == maxdepth:
		print(path)
		new = Protein.fold(sprotein)
		while not new:
			new = Protein.fold(sprotein)

		best = deepcopy(new)

		best.score = 0

		scores = []

		# if fold is valid after mutation continue to next iteration
		if new:
			new.make_grid()
			new.calc_score()
			scores.append(best.score)

			if new.score < best.score:
				print(new.score, "<", best.score)
				best = deepcopy(new)

	if depth < maxdepth:
		avail_direcs = deepcopy(all_direcs)
		if depth > 0:
			if path[depth-1] == 'r':
				avail_direcs.remove('l')
			elif path[depth-1] == 'd':
				avail_direcs.remove('u')
			elif path[depth-1] == 'l':
				avail_direcs.remove('r')
			elif path[depth-1] == 'u':
				avail_direcs.remove('d')
		for direc in avail_direcs:
			path.append(direc)
			depth_path(path, depth+1, maxdepth)
			path.pop()
	return path
depth_path(path, 0, maxdepth)

        