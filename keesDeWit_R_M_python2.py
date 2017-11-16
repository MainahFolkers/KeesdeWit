# complexity
def complexity(n):
    return 3**(n - 1)


# class protein
class Protein:
    
    # maak een class protein en een class amino -> dan kun je inheriten 
    # en dan kun je makkelijk de eigenschappen aanroepen
        
    def __init__(self, chain):
        
        # chain is string of amino acid sequence
        self.chain = chain
        
        # n is amount of amino acids in protein chain
        self.n = len(self.chain)
        
        # create a grid of n x n (can later expand to 3D)
        # x = xcoordinate in coordinate system, y = y coordinate
        # for aa in range (self.n):
        self.grid = [[{'aa': None, 'x': None, 'y': None}] for n in range(self.n)] # for x in range(self.n)]
        
        # stability score is calculated with score function
        # self.score = self.score() 
#     def score(self):
#         score = 0
#         return score;


# demo protein
protein = Protein("HHPHHHPH")
print protein.chain
print protein.n
print protein.grid
# print protein.score


# folding function
import random

# function that randomly selects the direction of a chain: up, down, left, right
def directionChooser(chain):
    direction = random.randint(0,3)
    if direction == 0:
        return 'u'
    elif direction == 1:
        return 'r'
    elif direction == 2:
        return 'd'
    else: 
        return 'l'

# folds the protein. calls functions 'directionChooser', coordMaker'.
def fold(protein):
    
    # array with 'path' of protein: first chain always goes up
    foldingArray = ['u']
    
    # select direction of second chain. only up/right because mirroring
    direction = random.randint(0,1)
    if direction == 0:
        foldingArray.append('u')
    else:
        foldingArray.append('r')
        
    # third amino acid and up: random direction
    # fold every chain (chain between amino acids: #chains = (#aa's - 1)
    # for c in range (2, protein.n - 2):
    #     foldingArray.append(directionChooser(c))
    
    # def coordMaker(foldingArray)

    for aa in range(protein.n):
#         print(amino)
        print(protein.grid[0]['aa'])# = Protein.chain[amino]
         
    #             while collision :#pseude\0
    #                 direction = directionChooser(c)
    #                 if direction = goed:
    #                     break en return right direction
    #                 and push in array
    
    # call directionChooser until no collisions with previous chain:
    # while true
    #     newDirection = directionChooser(c)
    #     if != collision
    #         break
                   
    # second chain either goes up or right (not left because of mirroring)

    # fold protein straight down
    # Protein.grid[0]['aa'] = Protein.chain[i]
#     currentX, currentY = 0, 0
    
#     for direction in folding:
#         if direction == 'u':
#             Protein.grid[0][c]['x'] = currentX
#             Protein.grid[0][c]['y'] = currentY + 1
#         elif direction== 'd':
#             Protein.grid[0][c]['x'] = currentX
#             Protein.grid[0][c]['y'] = currentY - 1
#         elif direction== 'r':
#             Protein.grid[0][c]['x'] = currentX + 1
#             Protein.grid[0][c]['y'] = currentY
#         else:
#             Protein.grid[0][c]['x'] = currentX - 1
#             Protein.grid[0][c]['y'] = currentY
        # pseudo: 0 = naar links relative aan eerdere chain
        # 1 = rechtdoor relative
        # 
                
        #Protein.grid[0][c]['i'] = c


# demo folding protein
fold(protein)
print protein.grid


# plot
import matplotlib.pyplot as plt

# iterate of grid rows
for r in range(protein.n):
    # iterate over grid columns
    for c in range(protein.n):
        if protein.grid[r][c]['i'] is not None:
            if protein.grid[r][c]['aa'] == 'H':
                color = 'red'
            else:
                color = 'blue'
                
            plt.plot(r, c, color=color)   


# plot that shit
plt.show()
