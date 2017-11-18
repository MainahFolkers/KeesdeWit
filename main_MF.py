from Protein_class_MF import *
from fold_MF import *
import matplotlib.pyplot as plt
import numpy as np

protein = Protein("HHPHHHPH")

# fold protein until valid folding
protein.coordinates = fold(protein)
while not protein.coordinates:
    protein.coordinates = fold(protein)

xs = []
ys = []
# iterate over amino acids
for i in range(protein.n):
    # seperate coordinate sin xs and ys
    xs.append(protein.coordinates[i][0])
    ys.append(protein.coordinates[i][1])

xmin = min(xs)
ymin = min(ys)
# calculate ranges voor x and y
xran = max(xs) - xmin + 1
yran = max(ys) - ymin + 1

# make grid with minimal ranges
protein.make_grid(xran, yran)

# iterate over amino acids
for i in range(protein.n):
    # transform coordinates onto minimal grid
    xs[i] = xs[i] - xmin
    ys[i] = ys[i] - ymin

    x = xs[i]
    y = ys[i]

    # place amino acids onto grid
    protein.grid[x][y].aa = protein.chain[i]
    protein.grid[x][y].i = i

    # determine color for point on plot
    if protein.chain[i] == 'H':
        col = 'red'
    else:
        col = 'blue'

    plt.scatter(x, y, color=col)
    plt.annotate(i, xy=(x, y))
# improvement: draw line behind points
plt.plot(xs, ys, color='black')
# show grid on plot
plt.grid(b=True)

xmax = max(xs)
ymax = max(ys)

# set tick spacing to 1
plt.xticks(np.arange(min(xs), xmax+1, 1))
plt.yticks(np.arange(min(ys), ymax+1, 1))

# calculte and print score
# does not calculate score correctly!
protein.calc_score(xmax, ymax)
print "Score = ", protein.score

# show plot
plt.show()
