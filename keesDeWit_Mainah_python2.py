# state space
def complexity(n):
    return 3**(n - 2)

# protein with subclass amino acid
class Protein:
    class AA:
        aa = None
        i = None
        # up down left right
        direc = None
        
    def __init__(self, chain):
        # chain is string of amino acid sequence
        self.chain = chain
        # n is amount of amino acids in protein chain
        self.n = len(self.chain)
        # create a grid of 2n x 2n (can later expand to 3D)
        self.grid = [[self.AA() for y in range(2*self.n)] for x in range(2*self.n)]
        # stability score is calculated with score function
        self.score = self.score()
    
    # still needs to be implemented in Score_function.ipynb
    def score(self):
        score = 0
        return score

# demo protein
protein = Protein("HHPHHHPH")
print protein.chain
print protein.n
print protein.score
print protein.AA.aa
print protein.AA.i

# fold function
def fold(protein):
    # iterate over amino acids c = column
    for c in range(protein.n):
        # fold protein straight down
        protein.grid[0][c].aa = protein.chain[c]
        protein.grid[0][c].i = c

# demo fold function
fold(protein)
print protein.grid[0][0].aa

# plot folded protein
import matplotlib.pyplot as plt
# iterate over grid rows
for r in range(protein.n):
    # iterate over grid columns
    for c in range(protein.n):
        if protein.grid[r][c].i is not None:
            if protein.grid[r][c].aa == 'H':
                color = 'red'
            else:
                color = 'blue'
            # plot that shit 
            plt.plot(r, c, color=color)

# demo plot folded protein
plt.show()