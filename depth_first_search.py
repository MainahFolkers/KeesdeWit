#!/usr/bin/env python
from Protein_class import *
import matplotlib.pyplot as plt
from copy import deepcopy
from collections import deque
import math

def depth_first_search(protein):
	maxdepth = (protein.n-2)
	all_direcs = ['r', 'd', 'u', 'l']
	scores = []
	new = protein
	best = deepcopy(new)
	best.score = 0
	print(best)

	# # splits longer proteins in seperate lists
	# def split_list(alist, wanted_parts=1):
	#     length = len(alist)
	#     return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
	#              for i in range(wanted_parts) ]

	# actual depth-first search
	def depth_path(protein, depth, maxdepth, bestscore):
		# if current depth is smaller than maxdepth, new directions are chosen constructively
		if depth < maxdepth:
			avail_direcs = deepcopy(all_direcs)

			# Because we want to eliminate as many mirror immages as possible, first direction
			# is always set to right
			if depth == 0:
				protein.directions.append('r')
				avail_direcs.remove('d')
				avail_direcs.remove('l')
			# from third bond forward, directions are chosen constructively
			elif depth > 0:
				if protein.directions[depth] == 'r':
					avail_direcs.remove('l')
				elif protein.directions[depth] == 'd':
					avail_direcs.remove('u')
				elif protein.directions[depth] == 'l':
					avail_direcs.remove('r')
				elif protein.directions[depth] == 'u':
					avail_direcs.remove('d')
			# the actual recursive part of the depth-first search
			for direc in avail_direcs:
				protein.directions.append(direc)
				depth_path(protein, depth+1, maxdepth, 0)
				protein.directions.pop()

		# if max depth is reached a chosen directions are turned into coordinates
		elif depth == maxdepth:

			H_count = 0
			for i in range(0,depth+1):
					if protein.chain[i] == "H":
						H_count += 1

			if H_count%2 == 0:
				max_H_score = H_count/2
			else:
				max_H_score = (H_count/2) - 0.5

			C_count = 0
			for i in range(0,depth+1):
					if protein.chain[i] == "C":
						C_count += 1

			if C_count%2 == 0:
				max_C_score = ((C_count/2) * 5) + max_H_score
			else:
				max_C_score = (((C_count/2) * 5) - 2.5) + max_H_score

			x, y = 0, 0
			protein.coordinates = [[x, y]]

			aa = 0
	 
			# iterate over amino acids --> TO DO: replace next block with calling check_val() method of protein
			while aa < maxdepth+1:
				[nx, ny] = protein.bend(protein.directions[aa], x, y)
				if [nx, ny] not in protein.coordinates:
					# update x and y
					[x, y] = [nx, ny]
					# new coordinates are safed
					protein.coordinates.append([x, y])
					# continue with next amino acid
					aa += 1
				elif [nx, ny] in protein.coordinates:
					return None

			new = protein
			# creates a grid and calculates scores from that
			protein.make_grid()
			protein.calc_score()

			scores.append(best.score)

			# evaluates whether current score is better than the best score found
			# if new score is better, it becomes new best score and its' directions are saved
			# optimal_score = H_count
			if new.score < best.score:
				if C_count == 0:
					print("AANTAL H'S IS: "+str(H_count))
					print(max_H_score)
					if abs(new.score) > max_H_score:
						print("SCORE IS GROTER DAN HET AANTAL H'S GEDEELD DOOR 2 NAMELIJK: "+str(new.score))
						best = deepcopy(new)
					elif abs(new.score) == max_H_score:
						print("SCORE IS GELIJK AAN HET AANTAL H'S GEDEELD DOOR 2 NAMELIJK: "+str(new.score))
					elif abs(new.score) < max_H_score:
						print("SCORE IS LAGER DAN HET AANTAL H'S GEDEELD DOOR 2 NAMELIJK: "+str(new.score))
				else:
					print("AANTAL H'S IS: "+str(H_count))
					print(max_H_score)
					print("AANTAL C'S IS: "+str(C_count))
					print(max_C_score)
					if abs(new.score) > max_C_score:
						print("SCORE IS GROTER DAN HET AANTAL C'S GEDEELD DOOR 2 NAMELIJK: "+str(new.score))
						best = deepcopy(new)
					elif abs(new.score) == max_C_score:
						print("SCORE IS GELIJK AAN HET AANTAL C'S GEDEELD DOOR 2 NAMELIJK: "+str(new.score))
					elif abs(new.score) < max_C_score:
						if abs(new.score) >= (C_count + H_count)/2:
							print("SCORE IS GELIJK AAN HET AANTAL GELADEN MOLECULEN NAMELIJK: "+str(new.score))
				# best_score = new.score
				# best_directions.appendleft(deepcopy(path))
				# best_coords.appendleft(deepcopy(protein.coordinates))

	# if a protein is longer than 10 amino acids (so a maxdepth of 8) the protein is split
	# into lists with a max length of 10 amino acids
	# if maxdepth > 9:

	# 	# divides the length of the protein and rounds it up to the nearest integer
	# 	num_of_list = math.ceil(protein.n/10)



	# else:


	# runs the depth first search for the protein visualizes the folds,
	# prints the best score, coordinates and directions
	depth_path(protein, 0, maxdepth, best.score)
	print(best.score)
	print(best.directions)

	return best