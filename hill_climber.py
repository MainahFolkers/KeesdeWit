from Protein_class_MF import *
from fold_MF import *
import matplotlib.pyplot as plt
import numpy as np


# --------------------------------------------------------------------------------------------
# EXTREMELY IMPORTANT
# 
# File "Protein_class_MF.py": function "calc_score" now returns a tuple with (Hbonds, score)
# --------------------------------------------------------------------------------------------

# this is what we're gonna fold optimally
protein = Protein("HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH")

# hillclimber
def hill_climb(protein):

	# ---------------------------------------------------------
	# init beginning
	# ---------------------------------------------------------
	
	# so that first climb is always an improvement
	bestScore = 1000000

	# tell the climber that newFold and bestFold are instances of class Protein
	# (hier zijn credits voor Bas in place -> we misten een instance van class protein.... )
	newFold = protein
	bestFold = protein

	# ---------------------------------------------------------
	# climb!    
	# --------------------------------------------------------- 

	# for now, just some iterations. Later: climbing until no improvements for x iterations
	for i in range(10000):

		# calculate current score
		newFold.coordinates = fold(protein)

		xs = []
		ys = []
		# iterate over amino acids
		for i in range(newFold.n):
		    # seperate coordinates in xs and ys
		    xs.append(newFold.coordinates[i][0])
		    ys.append(newFold.coordinates[i][1])

		xmin = min(xs)
		ymin = min(ys)
		# calculate ranges voor x and y
		xran = max(xs) - xmin
		yran = max(ys) - ymin

		# make grid with minimal ranges
		newFold.make_grid(xran, yran)

		# iterate over amino acids
		for i in range(newFold.n):
		    # transform coordinates onto minimal grid
		    xs[i] = xs[i] - xmin
		    ys[i] = ys[i] - ymin

		    # place amino acids onto grid
		    newFold.grid[xs[i]][ys[i]].aa = newFold.chain[i]
		    newFold.grid[xs[i]][ys[i]].i = i

		newScoreVar = bestFold.calc_score()
		# print(newScoreVar[1])

		# check if improvement
		if newScoreVar[1] < bestScore:
			# print("JOEPIIIEEE we zitten in de if")
			# keep that folding & score
			bestScore = newScoreVar[1]
			bestFold = newFold
			print(str(bestScore), " verbetering! :)" )
		# print(bestFold.coordinates)

	# return best score
	return bestScore, bestFold

result = hill_climb(protein)
# print("score result = ")
print(result[0])

# # determine color for point on plot
# if result.chain[i] == 'H':
#     col = 'red'
# else:
#     col = 'blue'

# plt.scatter(xs[i], ys[i], s=120, zorder=2, color=col)
# plt.annotate(i, xy=(xs[i], ys[i]), xytext=(xs[i] + 0.05, ys[i] + 0.05), fontsize=20)

# # plot black line behind / between points
# plt.plot(xs, ys, lw=3, zorder=1, color='black')

# # show grid on plot
# plt.grid(b=True)

# # set tick spacing to 1
# plt.xticks(np.arange(min(xs), xran + 1, 1))
# plt.yticks(np.arange(min(ys), yran  + 1, 1))

# # plot dotted line to visualize H-bond
# for i in range(abs(protein.score)):
#     plt.plot(Hbondsvar[0][i], Hbondsvar[0][i + 1], lw=3, zorder=3, color='black', linestyle='--')

# # show plot
# plt.show()







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
