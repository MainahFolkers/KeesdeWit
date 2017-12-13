#!/usr/bin/env python
from Protein_class_MF import *
from visualize_MF import *
import matplotlib.pyplot as plt
from copy import deepcopy
from collections import deque

sprotein = Protein("HHPHHHPHHHHHHHHHHHHH")
maxdepth = (sprotein.n-2)
path = ['r']
all_direcs = ['r', 'd', 'u', 'l']
best_score = 0

def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
             for i in range(wanted_parts) ]

if maxdepth > 10:
	list_length = split_list(sprotein.chain, wanted_parts=(int(sprotein.n/10)))
	for i in list_length:
		maxdepth = (len(i) - 2)
		part_protein = Protein(i)


def depth_path(protein, depth, maxdepth):
	#global best_score
	if depth == maxdepth:

		x, y = 0, 0
		sprotein.coordinates = [[x, y]]

		aa = 0

		# iterate over amino acids
		while aa < maxdepth + 1:
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

		grid = part_protein.make_grid()
		sprotein.calc_score()
		new_score = sprotein.score
		if new_score < best_score:
			best_score = new_score
			return best_score

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
for r in list_length:		
	depth_path(path, 0, maxdepth)
	print(best_score)
	print(path)
	print(sprotein.coordinates)

        