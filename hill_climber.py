from Protein_class_MF import *
from fold_MF import *
# import matplotlib.pyplot as plt
# import numpy as np


# ------------------------------------------------------------------------------------
# EXTREMELY IMPORTANT
# 
# File "Protein_class_MF.py" now returns a tuple with (Hbonds, score)
# 
# for plotting in file "main.MF.py": store result of calc_score() in variable "Hbondsvar"
# accessing the actual list of Hbonds from the tuple that calc_score returns:
# --> Hbondsvar[0]
# 
# ------------------------------------------------------------------------------------



def hill_climb(protein):

	# ---------------------------------------------------------
	# init beginning
	# ---------------------------------------------------------

	# defining start points

	# init, first folding is the only folding --> hence also the best one
	bestFold.coordinates = fold(protein)

	# calculate the score of the first folding
	# calc_score yields tuple -> score is 2nd thing in the tuple
	bestScoreVar = bestFold.calc_score()
	bestScore = bestScoreVar[1]
	

	# ---------------------------------------------------------
	# climb that hill!    
	# --------------------------------------------------------- 

	# for now, just some iterations. Later: climbing until no improvements for x iterations
	for i in range(50):

		# calculate current score
		newFold = fold(protein)
		newScoreVar = newFold.calc_score()
		print(newScoreVar[1])

		# check if improvement
		if newScoreVar[1] < bestScore:
			print("JOEPIIIEEE we zitten in de if")
			# keep that folding & score
			bestScore = newScoreVar[1]
			bestFold = newFold
		print(newFold.coordinates)

	# return best score
	return bestScore, bestFold

# result = hill_climb(bestFold)
# print("score result = ")
# print(result)







#     fold(protein)
#     protein.calc_score()
#     cur_score = protein.score

#     fold(protein)
#     protein.calc_score()
#     new_score = protein.score

#     if cur_score > new_score:
#         cur_score = new_score

# # x = choice(for i in range(protein.n))
# # directions = [r, d, u, l, u, d, r]
# # direction[x] = andere richting

# protein = Protein("HHPHHHPH")
# hill_climb(protein)
