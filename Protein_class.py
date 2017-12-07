from random import *

class Protein:
    
    # amino acid properties
    class AA():

        # amino acid type
        aa = None

        # amino acid index in chain
        i = None

    # protein chain properties
    def __init__(self, chain):

        # sequence
        self.chain = chain

        # chain length (number of amino acids)
        self.n = len(self.chain)

        # folding and grid
        self.directions = ['' for i in range(self.n - 1)]
        self.coordinates = [[0, 0], [1, 0]]
        self.xs = []
        self.ys = []
        self.xran = 0
        self.yran = 0
        self.grid = []

        # stability score (calculated with score function)
        self.score = 0

    """
    fold() folds the protein randomly while adhering to constraints.
    fold() takes an argument of the class Protein, and an argument
    that specifies at what point folding should start.
    if no second argument is passed, folding starts from scratch.
    """
    def rand_fold(protein):
        
        # 1st AA coordinates are (x = 0, y = 0)
        x, y = 0, 0
        self.coordinates[0] = [x, y]

        # 1st direction always to the right
        self.directions[0] = 'r'

        # direction 'right' means y-coordinate doesn't change; and x = x + 1 
        x = x + 1

        # 2nd AA's coordinate placed in coordinate list (variable 'y' still 0)
        self.coordinates[1] = [x, y]    

        # 2nd direction can only be 'up' or 'right' (symmetry pruning)
        direcs = ['r', 'u']

        # choose randomly to go 'up' or 'right'
        direc = choice(direcs)

        # put 2nd direction in list
        self.directions[1] = direc

        # calculate coordinate of 3rd AA depending on the chosen direction
        self.coordinates[2] = self.bend(direc, x, y)

        # bend() returns [x, y] --> set current x and y to those values
        x = self.coordinates[2][0]
        y = self.coordinates[2][1]

        # continue folding from fourth amino acid
        aa_index = 3

        # all possible directions (right, down, left, up)
        all_direcs = ['r', 'd', 'l', 'u']

        # 2nd direction is most recently chosen before entering the loop
        chosen_direc = self.directions[1]

        # contunue folding until end of chain
        while aa_index < self.n:

            # available directions to choose from in current iteration
            avail_direcs = all_direcs

            # the chain cannot go back where it came from
            if chosen_direc == 'r':
                avail_direcs.remove('l')
            elif chosen_direc == 'd':
                avail_direcs.remove('u')
            elif chosen_direc == 'l':
                avail_direcs.remove('r')
            elif chosen_direc == 'u':
                avail_direcs.remove('d')

            # while there are still directions available to choose
            while avail_direcs:

                # random choice 1 of 3 available directions
                chosen_direc = choice(avail_direcs)

                # i - 1 is index of bond
                self.directions[aa_index - 1] = chosen_direc

                # sample without replacement, remove chosen direction from options
                avail_direcs.remove(chosen_direc)

                # calculate new coordinates
                [nx, ny] = self.bend(protein.directions[aa_index], x, y)

                # if the current new coordinate is not yet occupied by other amino acid
                if [nx, ny] not in self.coordinates:
                    # update old x and y
                    [x, y] = [nx, ny]
                    # new coordinates are safed
                    self.coordinates[aa_index] = [nx, ny]

                    # continue with next amino acid
                    aa_ndex += 1
                    break

                # if folding is invalid
                elif [nx, ny] in self.coordinates:
                    return None

            # if there are no directions left to choose from
            if not avail_direcs:
                return None

        # as usual, return class Protein :)        
        return self

    """
    bend() converts directions to coordinates
    'up' and 'down' affect the y coordinates, 
    'left' and 'right' affect the x coordinates.
    """
    def bend(self, direc, x, y):
       
        if direc == 'u':
            y = y + 1
        elif direc == 'r':
            x = x + 1
        elif direc == 'd':
            y = y - 1
        elif direc == 'l':
            x = x - 1
            
        return [x, y]

    """
    make_grid() works with the coordinates output from bend()
    """
    def make_grid(self):
        xs = []
        ys = []

        # iterate over amino acids
        for aa in range(self.n):

            # seperate coordinates in arrays for xs and for ys
            xs.append(self.coordinates[aa][0])
            ys.append(self.coordinates[aa][1])

        # define minima
        xmin = min(xs)
        ymin = min(ys)

        # calculate ranges voor x and y
        self.xran = max(xs) - xmin
        self.yran = max(ys) - ymin

        # instantiate dynamic grid (based on current x and y minima and maxima)
        self.grid = [[self.AA() for y in range(self.yran + 1)] for x in range(self.xran + 1)]

        # iterate over amino acids in chain
        for aa in range(self.n):

            # transform coordinates onto dynamic grid
            xs[aa] -= xmin
            ys[aa] -= ymin

            # place amino acids onto grid
            self.grid[xs[aa]][ys[aa]].aa = self.chain[aa]
            self.grid[xs[aa]][ys[aa]].i = aa

        # save coordinates in arrays for xs and ys
        self.xs = xs
        self.ys = ys

    """
    calc_score() calculates the score by traversing the grid and looking for adjacent Hs.
    """
    def calc_score(self):

        # lowerbound: worst score = 0
        score = 0

        # instantiate array of H-bonds
        Hbonds = []

        # iterate over rows in grid
        for r in range(self.xran + 1):

            # iterate over columns in grid
            for c in range(self.yran + 1):

                # check if there is a H in cell of grid
                if self.grid[r][c].aa == 'H':

                    # the cell index that the H is in
                    i = self.grid[r][c].i

                    # if cell index does not go outside of grid range
                    if r + 1 < self.xran + 1:

                        # ni is next index 
                        ni = self.grid[r + 1][c].i

                        # naa is next amino acid
                        naa = self.grid[r + 1][c].aa

                        # check if there is an H on the right that is not connected in the chain
                        if naa is 'H' and abs(ni-i) is not 1:

                            # save that H bond
                            Hbonds.append([i, ni])

                            # update score
                            score -= 1

                    # if index does not go outside of grid range
                    if c + 1 < self.yran + 1:

                        # ni is next index
                        ni = self.grid[r][c + 1].i

                        # naa is next amino acid
                        naa = self.grid[r][c + 1].aa

                        # check if there is an H above that is not connected in the chain
                        if naa is 'H' and abs(ni-i) is not 1:

                            # save the H bond
                            Hbonds.append([i, ni])

                            # update score
                            score -= 1

        # save ultimate score
        self.score = score

        # yield the places of the H-bonds (for visualization)
        return Hbonds
