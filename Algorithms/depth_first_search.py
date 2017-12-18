
from Protein_class_2D import *
from copy import deepcopy

MAX_LEN = 10

def depth_fold(directions, depth, maxdepth, folds, last_direc):
	# directions a bond can move to (right, down, left, up)
	all_direcs = ['r', 'd', 'u', 'l']

	if depth < maxdepth:
		# last direc is last direction of previous fragment
		if last_direc == None:
			# reset available directions to choose from to all directions
			avail_direcs = deepcopy(all_direcs)
			if depth == 0:
				# first bond is to right
				directions.append('r')
				# the chain cannot go back where it came from
				avail_direcs.remove('l')
				# second bond can only bend up or to right
				avail_direcs.remove('d')
			# from third bond onward, directions are chosen constructively
			elif depth > 0:
				# the chain cannot go back where it came from
				if directions[depth] == 'r':
					avail_direcs.remove('l')
				elif directions[depth] == 'd':
					avail_direcs.remove('u')
				elif directions[depth] == 'l':
					avail_direcs.remove('r')
				elif directions[depth] == 'u':
					avail_direcs.remove('d')
			# iterate over available directions
			for direc in avail_direcs:
				# add next direction
				directions.append(direc)
				# continue to next bond
				depth_fold(directions, depth+1, maxdepth, folds, last_direc)
				# delete last direction
				directions.pop()

		elif last_direc != None:
			# reset available directions to choose from to all directions
			avail_direcs = deepcopy(all_direcs)
			if depth == 0:
				# first bond is to right
				directions.append('r')
				# the chain cannot go back where it came from
				avail_direcs.remove('l')
				# second bond can only bend up or to right
				avail_direcs.remove('d')
			# from third bond onward, directions are chosen constructively
			elif depth > 0:
				# the chain cannot go back where it came from
				if directions[depth] == 'r':
					avail_direcs.remove('l')
				elif directions[depth] == 'd':
					avail_direcs.remove('u')
				elif directions[depth] == 'l':
					avail_direcs.remove('r')
				elif directions[depth] == 'u':
					avail_direcs.remove('d')
			# iterate over available directions
			for direc in avail_direcs:
				# add next direction
				directions.append(direc)
				# continue to next bond
				depth_fold(directions, depth+1, maxdepth, folds, last_direc)
				# delete last direction
				directions.pop()

	# if max depth is reached
	else:
		folds.append(deepcopy(directions))

	return folds

def depth_score(folds, protein, best):

	fold_number = 0
	# iterate over all possible folds
	for fold in folds:
		# set directions in protein
		protein.directions = deepcopy(fold)

		# if the fold is invalid
		if not protein.mut_fold():
			# remove current fold
			folds.pop(fold_number)
			continue

		# if fold is valid
		elif protein.mut_fold():
			# place protein on grid
			protein.make_grid(protein.n)
			# calculate score
			protein.calc_score()

		# if score improved
		if protein.score < best.score:
			# save current best fold
			best = deepcopy(protein)
		# if current score is worse
		if protein.score > best.score:
			#remove current fold
			folds.pop(fold_number)

	fold_number += 1

	return best

def depth_first_search(protein):

	# maximum depth is amount of bonds - 1
	maxdepth = protein.n - 2

	best = deepcopy(protein)

	last_direc = None

	directions = []
	folds = []
	total_folds = []

	# if protein is longer than 10 amino acids
	if protein.n > MAX_LEN:

		fragments = [] 
		complete_fold = []

		# iterate over amount of fragments
		for i in range(protein.n/MAX_LEN):

			# selects part of the protein chain
			frag_chain = protein.chain[int(protein.n / MAX_LEN * i) : int(protein.n / MAX_LEN * (i + 1))]
			# turns fragment into protein
			fragment = Protein(frag_chain)
			# maximum depth is amount of bonds - 1
			maxdepth = fragment.n - 1
			# add fragments to one list
			fragments.append(fragment)

		# iterate over all fragments
		for fragment in fragments:
			best.score = 1

			# fill folds with all possible direction arrays starting at depth = 0
			folds = depth_fold(directions, 0, maxdepth, folds, last_direc)
			# calculate best score for fragment
			best = depth_score(folds, fragment, best)
			# retrieve all folds with maximum scores
			best = depth_score(folds, fragment, best)
			# save best folds
			total_folds.append(folds)
			# empty list for next fragment
			folds = []
			# remember last direction of fragment
			last_direc = best.directions[maxdepth]


	else:

		best.score = 1
		# fill folds with all possible direction arrays starting at depth = 0
		folds = depth_fold(directions, 0, maxdepth, folds, last_direc)
		# calculate best score for protein
		best = depth_score(folds, protein, best)

	print("depth first search: Best score = ", best.score)

	return best