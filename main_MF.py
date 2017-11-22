from Protein_class_MF import *
from fold_MF import *
from hill_climber import *
import matplotlib.pyplot as plt
import numpy as np

# protein to optimize
protein = Protein("HHPHHHPH")

# optimize folding with hillclimber algorithm
# resultsTuple = hill_climb(protein)

# fold protein until valid folding
protein.coordinates = fold(protein)

xs = []
ys = []
# iterate over amino acids
for i in range(protein.n):
    # seperate coordinates in xs and ys
    xs.append(protein.coordinates[i][0])
    ys.append(protein.coordinates[i][1])

xmin = min(xs)
ymin = min(ys)
# calculate ranges voor x and y
xran = max(xs) - xmin
yran = max(ys) - ymin

# make grid with minimal ranges
protein.make_grid(xran, yran)

plt.figure()

# iterate over amino acids
for i in range(protein.n):
    # transform coordinates onto minimal grid
    xs[i] = xs[i] - xmin
    ys[i] = ys[i] - ymin

    # place amino acids onto grid
    protein.grid[xs[i]][ys[i]].aa = protein.chain[i]
    protein.grid[xs[i]][ys[i]].i = i

    # determine color for point on plot
    if protein.chain[i] == 'H':
        col = 'red'
    else:
        col = 'blue'

    plt.scatter(xs[i], ys[i], s=120, zorder=2, color=col)
    plt.annotate(i, xy=(xs[i], ys[i]), xytext=(xs[i] + 0.05, ys[i] + 0.05), fontsize=20)

# plot black line behind / between points
plt.plot(xs, ys, lw=3, zorder=1, color='black')

# show grid on plot
plt.grid(b=True)

# set tick spacing to 1
plt.xticks(np.arange(min(xs), xran + 1, 1))
plt.yticks(np.arange(min(ys), yran  + 1, 1))

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
Hbondsvar = protein.calc_score()
print(Hbondsvar[1])
print("Score = ", protein.score)
print protein.directions

# plot dotted line to visualize H-bond
for i in range(abs(protein.score)):
    plt.plot(Hbondsvar[0][i], Hbondsvar[0][i + 1], lw=3, zorder=3, color='black', linestyle='--')

# show plot
plt.show()
