
from Protein_class import *
import matplotlib.pyplot as plt
from copy import deepcopy
from collections import deque
from itertools import islice
import math

# splits longer proteins in seperate lists
def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
    	for i in range(wanted_parts) ]

# actual depth-first search
def depth_fold(directions, depth, maxdepth, folds):
	all_direcs = ['r', 'd', 'u', 'l']
	# if current depth is smaller than maxdepth, new directions are chosen constructively
	if depth < maxdepth:
		# if last_direc == None:
		avail_direcs = deepcopy(all_direcs)
		# Because we want to eliminate as many mirror immages as possible, first direction
		# is always set to right
		if depth == 0:
			directions.append('r')
			avail_direcs.remove('d')
			avail_direcs.remove('l')
		# from third bond forward, directions are chosen constructively
		elif depth > 0:
			if directions[depth] == 'r':
				avail_direcs.remove('l')
			elif directions[depth] == 'd':
				avail_direcs.remove('u')
			elif directions[depth] == 'l':
				avail_direcs.remove('r')
			elif directions[depth] == 'u':
				avail_direcs.remove('d')
		# the actual recursive part of the depth-first search
		for direc in avail_direcs:
			directions.append(direc)
			depth_fold(directions, depth+1, maxdepth, folds)
			directions.pop()

		# elif last_direc != None:
		# 	avail_direcs = deepcopy(all_direcs)
		# 	if depth == 0:
		# 		if last_direc == 'r':
		# 			avail_direcs.remove('l')
		# 		elif last_direc == 'd':
		# 			avail_direcs.remove('u')
		# 		elif last_direc == 'l':
		# 			avail_direcs.remove('r')
		# 		elif last_direc == 'u':
		# 			avail_direcs.remove('d')
		# 	elif depth > 0:
		# 		if protein.directions[depth] == 'r':
		# 			avail_direcs.remove('l')
		# 		elif protein.directions[depth] == 'd':
		# 			avail_direcs.remove('u')
		# 		elif protein.directions[depth] == 'l':
		# 			avail_direcs.remove('r')
		# 		elif protein.directions[depth] == 'u':
		# 			avail_direcs.remove('d')
		# 	# the actual recursive part of the depth-first search
		# 	for direc in avail_direcs:
		# 		protein.directions.append(direc)
		# 		depth_fold(protein, depth+1, maxdepth)
		# 		protein.directions.pop()

	# if max depth is reached a chosen directions are turned into coordinates
	else:
		folds.append(deepcopy(directions))
	return folds

def depth_score(folds, protein, best):
	print("WE ZIJN IN DEPTH SCORE YAAAAAY")
	print(len(folds))
	for fold in folds:
		protein.directions = deepcopy(fold)
		# iterate over amino acids --> TO DO: replace next block with calling check_val() method of protein

		if not protein.mut_fold():
			print("ONGELDIGE VOUWING BOOOOEEE")
			continue
		elif protein.mut_fold():
			# creates a grid and calculates scores from that
			protein.make_grid(protein.n)
			protein.calc_score()

		# evaluates whether new score is better than the best score found
		# if new score is better, it becomes new best score and its' directions are saved
		if protein.score < best.score:
			best = deepcopy(protein)
	return best

def depth_first_search(protein):
	maxdepth = (protein.n-2)
	new = deepcopy(protein)
	best = deepcopy(protein)
	best.score = 1
	last_direc = None
	directions = []
	folds = []

	# division = 2
	# if maxdepth > 8:
	# 	maxdepth = ceil((new.n / division) -1)
	# 	division += 1
	# 	wanted_parts = int(new.n / (maxdepth + 1))

	# 	frag_chain = new.chain[:maxdepth + 2]
	# 	fragment = Protein(frag_chain)

	# 	for num in range(wanted_parts + 1):
	# 		depth_fold(0, maxdepth)
	# 		best = depth_score(folds, fragment, best)
	# 		folds = []
	# 		last_direc = new.directions[maxdepth]

	# else:
	folds = depth_fold(directions, 0, maxdepth, folds)
	best = depth_score(folds, protein, best)

	print("depth first search: Best score = ", best.score)
	return best

