
from Protein_class import *
import matplotlib.pyplot as plt
from copy import deepcopy
from collections import deque
from itertools import islice
import math

def depth_first_search(protein):
	maxdepth = (protein.n-2)
	all_direcs = ['r', 'd', 'u', 'l']
	folds = []
	new = deepcopy(protein)
	last_direc = None
	# print("MAXDEPTH IS NU: "+str(maxdepth))

	# splits longer proteins in seperate lists
	def split_list(alist, wanted_parts=1):
	    length = len(alist)
	    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
	             for i in range(wanted_parts) ]

	# actual depth-first search
	def depth_fold(protein, depth, maxdepth, bestscore):
		# if current depth is smaller than maxdepth, new directions are chosen constructively
		if depth < maxdepth:
			if last_direc == None:
				avail_direcs = deepcopy(all_direcs)
				# Because we want to eliminate as many mirror immages as possible, first direction
				# is always set to right
				if depth == 0:
					new.directions.append('r')
					avail_direcs.remove('d')
					avail_direcs.remove('l')
				# from third bond forward, directions are chosen constructively
				elif depth > 0:
					if new.directions[depth] == 'r':
						avail_direcs.remove('l')
					elif new.directions[depth] == 'd':
						avail_direcs.remove('u')
					elif new.directions[depth] == 'l':
						avail_direcs.remove('r')
					elif new.directions[depth] == 'u':
						avail_direcs.remove('d')
				# the actual recursive part of the depth-first search
				for direc in avail_direcs:
					new.directions.append(direc)
					depth_fold(new, depth+1, maxdepth, new.score)
					new.directions.pop()

			elif last_direc != None:
				avail_direcs = deepcopy(all_direcs)
				if depth == 0:
					if last_direc == 'r':
						avail_direcs.remove('l')
					elif last_direc == 'd':
						avail_direcs.remove('u')
					elif last_direc == 'l':
						avail_direcs.remove('r')
					elif last_direc == 'u':
						avail_direcs.remove('d')
				elif depth > 0:
					if new.directions[depth] == 'r':
						avail_direcs.remove('l')
					elif new.directions[depth] == 'd':
						avail_direcs.remove('u')
					elif new.directions[depth] == 'l':
						avail_direcs.remove('r')
					elif new.directions[depth] == 'u':
						avail_direcs.remove('d')
				# the actual recursive part of the depth-first search
				for direc in avail_direcs:
					new.directions.append(direc)
					depth_fold(new, depth+1, maxdepth, new.score)
					new.directions.pop()						

		# if max depth is reached a chosen directions are turned into coordinates
		elif depth == maxdepth:

			# iterate over amino acids --> TO DO: replace next block with calling check_val() method of protein
			new.x, new.y = 0, 0
			new.coordinates = [[new.x, new.y]]

			# start at second amino acid
			aa = 1
			while aa < maxdepth+2:
				# bond index
				b = aa - 1
				# check whether the direction yields a valid coordinate
				if new.check_val(aa, b):
					aa += 1
				# if mutation yiels invalid coordinate
				else:
					return None
			folds.append(deepcopy(new.directions))

	def depth_score(folds, protein):
		H_count = 0

		for i in range(0, maxdepth+2):
				if new.chain[i] == "H":
					H_count += 1

		if H_count%2 == 0:
			max_H_score = (H_count/2) * -1
		else:
			max_H_score = ((H_count/2) - 0.5) * -1

		C_count = 0
		for i in range(0, maxdepth+2):
				if new.chain[i] == "C":
					C_count += 1

		if C_count%2 == 0:
			max_C_score = (((C_count/2) * 5) * -1) + max_H_score
		else:
			max_C_score = ((((C_count/2) * 5) - 2.5) * -1) + max_H_score

		num = 0
		for fold in folds:
			num += 1
			# print(len(folds))
			# print(num)
			new.directions = deepcopy(fold)
			# iterate over amino acids --> TO DO: replace next block with calling check_val() method of protein
			new.x, new.y = 0, 0
			new.coordinates = [[new.x, new.y]]

			# start at second amino acid
			aa = 1
			while aa < maxdepth+2:
				# bond index
				b = aa - 1
				# check whether the direction yields a valid coordinate
				if new.check_val(aa, b):
					aa += 1
			# print(new.directions)
			# print(new.coordinates)
			

			# creates a grid and calculates scores from that
			new.make_grid()
			new.calc_score()

			# scores.append(best.score)

			# evaluates whether new score is better than the best score found
			# if new score is better, it becomes new best score and its' directions are saved
			# optimal_score = H_count
			# print(new.directions)
			# print(new.chain)

			if new.score < protein.score:
				protein = deepcopy(new)
				# print("PROTEIN SCORE IS NU: "+str(protein.score))
				if C_count == 0:
					print("AANTAL H'S IS: "+str(H_count))
					# print("SCORE IS: "+str(max_H_score))
					if new.score < max_H_score:
						print("SCORE IS BETER DAN HET AANTAL H'S GEDEELD DOOR 2 NAMELIJK: "+str(new.score))
						print(protein.score)
						print(protein.directions)
						return protein
					elif new.score == max_H_score:
						print("SCORE IS GELIJK AAN HET AANTAL H'S GEDEELD DOOR 2 NAMELIJK: "+str(new.score))
						# print(protein.score)
						print(protein.directions)
						return protein
					elif new.score > max_H_score:
						print("SCORE IS SLECHTER DAN HET AANTAL H'S GEDEELD DOOR 2 NAMELIJK: "+str(new.score))
				else:
					print("AANTAL H'S IS: "+str(H_count))
					print("SCORE IS: "+str(max_H_score))
					print("AANTAL C'S IS: "+str(C_count))
					print("SCORE IS: "+str(max_C_score))
					if new.score < max_C_score:
						print("SCORE IS BETER DAN HET AANTAL C'S GEDEELD DOOR 2 NAMELIJK: "+str(new.score))
					elif new.score == max_C_score:
						print("SCORE IS HET AANTAL C'S GEDEELD DOOR 2 NAMELIJK: "+str(new.score))
					elif new.score > max_C_score:
						if new.score <= (C_count + H_count)/2:
							print("SCORE IS GELIJK AAN HET AANTAL GELADEN MOLECULEN NAMELIJK: "+str(new.score))

	division = 2
	while maxdepth > 8:
		maxdepth = int((new.n/division) -1)	
		division += 1
		num_of_list = int(new.n/(maxdepth+1))
		wanted_parts = int(new.n/10)

print (split_list(A, wanted_parts))
	# print("DIT IS NUM OF LIST: "+str(num_of_list))
	# while num < range(num_of_list):
	# 	print(new.chain[2])
	for num in range(num_of_list):
		new.n = (maxdepth + 2) * num
		depth_fold(protein, 0, maxdepth, protein.score)
		depth_score(folds, protein)
		print(len(folds))
		if num < num_of_list-1:
			del folds[:]
		last_direc = new.directions[maxdepth]
	depth_score(folds,protein)
	
	
	


	# if a protein is longer than 10 amino acids (so a maxdepth of 8) the protein is split
	# into lists with a max length of 10 amino acids
	# if maxdepth > 9:

	# 	# divides the length of the protein and rounds it up to the nearest integer
	# 	num_of_list = math.ceil(protein.n/10)



	# else:


	# runs the depth first search for the protein visualizes the folds,
	# prints the best score, coordinates and directions
