#!/usr/bin/env python
from Protein_class import *
from visualize_MF import *
import matplotlib.pyplot as plt
from copy import deepcopy
from collections import deque
import math

sprotein = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")
maxdepth = (sprotein.n-2)
path = ['r']
all_direcs = ['r', 'd', 'u', 'l']
best_score = 0
best_score_list = []
best_direc_list = []


# splits longer proteins in seperate lists
def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
             for i in range(wanted_parts) ]

# actual depth-first search
def depth_path(protein, depth, maxdepth):
	global best_score

	# if max depth is reached a chosen directions are turned into coordinates
	if depth == maxdepth:

		x, y = 0, 0
		sprotein.coordinates = [[x, y]]

		aa = 0
 
		# iterate over amino acids --> TO DO: replace next block with calling check_val() method of protein
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

		# creates a grid and calculates scores from that
		grid = sprotein.make_grid()
		sprotein.calc_score()

		# evaluates whether current score is better than the best score found
		# if new score is better, it becomes new best score and its' directions are saved
		new_score = sprotein.score
		new_coords = sprotein.coordinates
		if new_score < best_score:
			best_score = new_score
			best_directions = path
			print(best_directions)
			# return best_score

	# if current depth is smaller than maxdepth, new directions are chosen constructively
	elif depth < maxdepth:
		avail_direcs = deepcopy(all_direcs)
		# because the first direction is already set from the start (see line 11) first direction to choose is
		# actually the second bond
		if depth == 0:
			avail_direcs.remove('d')
		# from third bond forward, directions are chosen constructively
		if depth >= 0:
			if path[depth] == 'r':
				avail_direcs.remove('l')
			elif path[depth] == 'd':
				avail_direcs.remove('u')
			elif path[depth] == 'l':
				avail_direcs.remove('r')
			elif path[depth] == 'u':
				avail_direcs.remove('d')
		# the actual recursive part of the depth-first search
		for direc in avail_direcs:
			path.append(direc)
			depth_path(path, depth+1, maxdepth)
			path.pop()
 
# if a protein is longer than 10 amino acids (so a maxdepth of 8) the protein is split
# into lists with a max length of 10 amino acids
if maxdepth > 9:
	# divides the length of the protein and rounds it up to the nearest integer
	num_of_list = math.ceil(sprotein.n/10)
	# splits the protein into the number of wanted lists (num_of_list)
	list_length = split_list(sprotein.chain, wanted_parts= num_of_list)
	print(list_length)
	# runs the depth first search for the amount of wanted lists, visualizes the folds,
	# prints the best score, coordinates and directions 
	for i in list_length:
		print(i)
		maxdepth = (len(i) - 2)
		sprotein = Protein(i)
		print(best_score)
		best_score = 0 
		depth_path(path, 0, maxdepth)
		best_score_list.append(best_score)
		print(best_score_list)
		#visualize(sprotein)
else:
	# runs the depth first search for the protein visualizes the folds,
	# prints the best score, coordinates and directions
	depth_path(path, 0, maxdepth)
	#visualize(sprotein)
	print(best_score)
	print(sprotein.coordinates)

        