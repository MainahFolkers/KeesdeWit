#!/usr/bin/env python
from Protein_class import *
from visualize_MF import *
import matplotlib.pyplot as plt
from copy import deepcopy
from collections import deque
import math

protein = "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"
sprotein = Protein(protein)
maxdepth = (sprotein.n-2)
path = ['r']
all_direcs = ['r', 'd', 'u', 'l']
best_score = 0
best_score_list = []
best_directions = deque()
best_coords = deque()
total_best_direc = []


# splits longer proteins in seperate lists
def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
             for i in range(wanted_parts) ]

# actual depth-first search
def depth_path(protein, depth, maxdepth, bestscore):
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
		# optimal_score = H_count
		if new_score < best_score :
			best_score = new_score
			best_directions.appendleft(deepcopy(path))
			best_coords.appendleft(deepcopy(sprotein.coordinates))
			# print(best_directions)
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
			# H_count = 0
			# for i in range(0,depth+1):
			# #for H in sprotein.chain:
			# 	if sprotein.chain[i] == "H":
			# 		H_count += 1
			# print(H_count)
			depth_path(path, depth+1, maxdepth, 0)
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
	num = 0
	for i in list_length:
		maxdepth = (len(i) - 2)
		sprotein = Protein(i)
		best_score = 0 

		depth_path(path, 0, maxdepth, best_score)
		best_score_list.append(best_score)

		print(best_score_list)

		sprotein = Protein(i)

		depth_path(sprotein, 0, maxdepth, best_score_list[num])

		best_directions.rotate(-1)

		while (num + 1) < len(best_directions):
			best_directions.popleft()

		d = 0
		best_direc = best_directions[num]
		while d < len(best_direc):
			total_best_direc.append(best_direc[d])
			# if total_best_direc[d] == 'r' and total_best_direc[d-1] == 'l':
				# while d < len(best_direc):
				# 	if total_best_direc[d] == 'r':
				# 		total_best_direc[d] = 'd'
			d += 1
		total_best_direc.append('d')
		#visualize(sprotein)
		num += 1
		sprotein = Protein(protein)
		sprotein.directions = total_best_direc
		print("DIT IS DE BESTE VOUWING: "+str(sprotein.directions))
		print(len(sprotein.directions))
		print(sprotein.chain)

else:

	# runs the depth first search for the protein visualizes the folds,
	# prints the best score, coordinates and directions
	depth_path(path, 0, maxdepth, best_score)
	best_score_list.append(best_score)
	print(best_score_list)
	depth_path(path, 0, maxdepth, best_score_list)
	print(best_coords)
	print(best_score_list)
	print(best_directions)
	visualize(sprotein)
