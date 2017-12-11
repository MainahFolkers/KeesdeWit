#!/usr/bin/env python
from Protein_class_MF import *
from visualize_MF import *
import matplotlib.pyplot as plt
from copy import deepcopy
from collections import deque

sprotein = Protein("HHPHHHPH")
maxdepth = (sprotein.n-2)
path = ['r']
all_direcs = ['r', 'd', 'u', 'l']
best_score = 1

def depth_path(protein, depth, maxdepth):
	if depth == maxdepth:

		x, y = 0, 0
		sprotein.coordinates = [[x, y]]

		aa = 0

		# iterate over amino acids
		while aa < sprotein.n - 1:
			[nx, ny] = sprotein.bend(path[aa], x, y)
			if [nx, ny] not in sprotein.coordinates:
				# update x and y
				[x, y] = [nx, ny]
				# new coordinates are safed
				sprotein.coordinates.append([x, y])
				# continue with next amino acid
				aa += 1
			elif [nx, ny] in sprotein.coordinates:
				return None

		grid = sprotein.make_grid()
		sprotein.calc_score()
		new_score = sprotein.score
		# if new_score < best_score:
		best_score = deepcopy(new_score)
		print(best_score)
		print(path)
		print(sprotein.coordinates)

	if depth < maxdepth:
		avail_direcs = deepcopy(all_direcs)
		if depth == 0:
			avail_direcs.remove('d')
		if depth >= 0:
			if path[depth] == 'r':
				avail_direcs.remove('l')
			elif path[depth] == 'd':
				avail_direcs.remove('u')
			elif path[depth] == 'l':
				avail_direcs.remove('r')
			elif path[depth] == 'u':
				avail_direcs.remove('d')
		for direc in avail_direcs:
			path.append(direc)
			depth_path(path, depth+1, maxdepth)
			path.pop()
	return path
depth_path(path, 0, maxdepth)

        