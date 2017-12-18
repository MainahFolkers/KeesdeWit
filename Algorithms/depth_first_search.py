from Protein_class_2D import *
from copy import deepcopy

def depth_fold(directions, depth, max_depth, folds):
	# directions a bond can move to (right, down, left, up)
	all_direcs = ['r', 'd', 'u', 'l']

	if depth < max_depth:
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
			depth_fold(directions, depth+1, max_depth, folds)
			# delete last direction
			directions.pop()
	# if max depth is reached
	elif depth is max_depth:
		folds.append(deepcopy(directions))

	return folds

def depth_first_search(protein):
	algo = 'df'

	scores = []
	# maximum depth is amount of bonds - 1
	max_depth = protein.n - 2

	new = deepcopy(protein)
	best = deepcopy(protein)

	# set impossible best score (minimum score = 0) so first fold always improves
	best.score = 1

	directions = []
	folds = []

	# fill folds with all possible direction arrays starting at depth = 0
	folds = depth_fold(directions, 0, max_depth, folds)

	# iterate over directions arrays
	for fold in folds:
		# set directions in protein
		protein.directions = deepcopy(fold)

		# if the fold is valid
		if protein.mut_fold():
			# place protein on grid
			protein.make_grid()
			# caluclate score
			protein.calc_score()
			# save best score per iteration
			scores.append(best.score)

		# if score improved
		if protein.score < best.score:
			# save current best fold
			best = deepcopy(protein)

	print("Depth first search: Best score = ", best.score)

	# open output file
	with open("ALGOS_" + protein.chain + ".txt", 'a+') as ofile:
		# write algorithm name as label
		ofile.write(algo + ',')
		# write comma seperated scores
		for score in scores:
			ofile.write(str(score) + ',')
		for i in range(2000 - len(scores)):
			ofile.write('-3,')
		# write new line
		ofile.write('\n')
	# close output file
	ofile.close()

	return best
