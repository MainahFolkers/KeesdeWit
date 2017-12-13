from random import choice
from copy import deepcopy

class Protein:
    class AA():
        # aa is type of amino acid (H/P, can later expand with C)
        aa = None
        # i in index in amino acid chain
        i = None

    def __init__(self, chain):
        # chain is string of amino acid sequence
        self.chain = chain
        # n is amount of amino acids in protein chain
        self.n = len(self.chain)
        # directions and coordinates assigned in fold functions
        self.directions = ['' for i in range(self.n - 1)]
        self.coordinates = [[] for i in range(self.n)]
        # coordinates components x and y
        self.xs, self.ys = [], []
        # x and y ranges calculated in make_grid()
        self.xran, self.yran, self.grid = 0, 0, []
        # protein stability is calculated with score function
        self.score = 0

    """
    rand_fold() folds the protein randomly while adhering to constraints.
    """
    def rand_fold(self):

        # function only outputs valid folds -> while invalid, fold again
        while True:
            # reset coordinates and directions for every new fold
            self.coordinates = [[] for i in range(self.n)]
            self.directions = ['' for i in range(self.n - 1)]

            # first AA starts on (x = 0, y = 0)
            self.x, self.y = 0, 0
            # NOTE: the 'self'-part in x and y is necessary for check_val() because scope
            self.coordinates[0] = [self.x, self.y]

            # first bond is to right
            self.directions[0] = 'r'

            # second AA always to right (x = x + 1, y = 0) which is (x = 1, y = 0)
            self.x = self.x + 1
            self.coordinates[1] = [self.x, self.y]

            # third AA can only bend up or to right
            direcs = ['r', 'u']
            direc = choice(direcs)

            # second bond is saved
            self.directions[1] = direc

            # save coordinates of third AA
            self.coordinates[2] = self.bend(direc, self.x, self.y)
            self.x = self.coordinates[2][0]
            self.y = self.coordinates[2][1]

            # directions a bond can move to (right, down, left, up)
            all_direcs = ['r', 'd', 'l', 'u']

            # most recently chosen direction is now the second bond
            chosen_direc = self.directions[1]

            # first three AAs already have coordinates -> start at fourth AA
            amino_coord_start = 3

            # continue going through AAs in the entire chain
            for amino_coord in range(amino_coord_start, self.n):

                # reset available directions to choose from to all directions
                avail_direcs = deepcopy(all_direcs)

                # the chain cannot go back where it came from
                if chosen_direc == 'r':
                    avail_direcs.remove('l')
                elif chosen_direc == 'd':
                    avail_direcs.remove('u')
                elif chosen_direc == 'l':
                    avail_direcs.remove('r')
                elif chosen_direc == 'u':
                    avail_direcs.remove('d')

                # while there are still directions left to choose from
                while avail_direcs:

                    # random choice 1 of 3 directions
                    chosen_direc = choice(avail_direcs)

                    # 'amino_coord - 1' is index of bond
                    self.directions[amino_coord - 1] = chosen_direc

                    # sample without replacement, remove chosen direction from options
                    avail_direcs.remove(chosen_direc)

                    # if chosen direction yields a valid coordinate:
                    if self.check_val(amino_coord, amino_coord - 1):
                        # if end of chain is reached, output valid fold
                        if amino_coord == self.n - 1:
                            return self
                        # if not the last AA, onwards to calculate the next coordinate!
                        else:
                            break

                # if there are no directions left to choose from
                if not avail_direcs:
                    # stop current folding and start again at (0,0) to make new random fold)
                    break

    """
    check_val() checks if folding is valid.
    if the proposed coordinates (calculated from self.direction using self.bend()) are 
    not already occupied, the coordinate is assigned to current AA in self.coordinates.
    NOTE: be careful with input, as it differs between random sampler and hillclimber.
    """
    def check_val(self, amino_coord, amino_dir):
        # calculate next coordinates according to selected direction
        [nx, ny] = self.bend(self.directions[amino_dir], self.x, self.y)

        # if those coordinates on grid are not yet occupied by other amino acid
        if [nx, ny] not in self.coordinates:

            # update old x and y
            [self.x, self.y] = [nx, ny]

            # save the new coordinates
            self.coordinates[amino_coord] = [nx, ny]
            
            return True

    """
    mut_dir() randomly mutates (a) randomly selected position(s) in the chain.
    if no number of mutations is specified, one mutation is made.
    """
    def mut_dir(self, mutnum = 1):
        # mutate a number of positions in the protein
        for mutation in range(mutnum):

            # randomly choice position (not first bond because symmetry)
            mut = choice(range(1, self.n - 1))

            # original direction on the chosen position cannot be chosen again
            direcs = ['r', 'd', 'l', 'u']
            direcs.remove(self.directions[mut])

            # replace original direction with randomly chosen mutation
            self.directions[mut] = choice(direcs)

        return self

    """
    mut_fold() folds mutated protein, with mutation position randomly
    determined by function mut_dir().
    """
    def mut_fold(self):
        # start interpreting self.directions from coordinates (0,0)
        self.x, self.y = 0, 0
        self.coordinates = [[x, y]]

        # start at first bond
        amino_dir_start = 0

        # calculate AA coordinates for every bonds direction
        for amino_dir in range(amino_dir_start, self.n - 1):

            # check whether the direction yields a valid coordinate
            if self.check_val(amino_coord, amino_coord - 1):
                # upon reaching end of chain, output valid fold
                if amino_coord == self.n - 1:
                    return self
                # if not the last AA, onwards to calculate the next coordinate!
                else: 
                    return None

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
    make_grid() works with the coordinates
    """
    def make_grid(self):
        # empty arrays with x and y components of coordinates
        self.xs = []
        self.ys = []

        # start at the first amino acid in the chain (index = 0)
        start = 0

        # place whole chain on grid
        for aa in range(start, self.n):

            # seperate coordinates in arrays for xs and for ys
            self.xs.append(self.coordinates[aa][0])
            self.ys.append(self.coordinates[aa][1])

        # define minima
        xmin = min(self.xs)
        ymin = min(self.ys)

        # calculate ranges voor x and y
        self.xran = max(self.xs) - xmin
        self.yran = max(self.ys) - ymin

        # instantiate dynamic grid (based on current x and y minima and maxima)
        self.grid = [[self.AA() for y in range(self.yran + 1)] for x in range(self.xran + 1)]

        # iterate over amino acids in chain
        for aa in range(self.n):

            # transform coordinates onto minimal grid
            self.xs[aa] -= xmin
            self.ys[aa] -= ymin

            # place amino acids onto grid
            self.grid[self.xs[aa]][self.ys[aa]].aa = self.chain[aa]
            self.grid[self.xs[aa]][self.ys[aa]].i = aa

        # # save x and y components --> rianne did specify self.x/y already at top
        # self.xs = xs
        # self.ys = ys

    """
    calc_score() calculates the score by traversing the grid 
    looking for adjacent Hs that are not connected in the chain.
    """
    def calc_score(self):

        # lowerbound: worst score = 0
        score = 0
        # upperbound: best score = -H / 2; where H is amount of H's in chain

        # instantiate list of 1-point bonds (H-H / H-C)
        HHCbonds = []
        # isntantiate list of 5-point bonds (C-C)
        CCbonds = []

        # iterate over rows in grid
        for r in range(self.xran + 1):

            # iterate over columns in grid
            for c in range(self.yran + 1):

                # check if there is an AA of interest (H/C) in cell of grid
                if self.grid[r][c].aa == 'H' or self.grid[r][c].aa =='C':
                    # the cell index that the H is in
                    i = self.grid[r][c].i

                    # if cell index does not go outside of grid range
                    if r + 1 < self.xran + 1:
                        # ni is next index
                        ni = self.grid[r + 1][c].i
                        # naa is next amino acid
                        naa = self.grid[r + 1][c].aa

                        # 1-point bond scenario H-H or C-H: check to right and unconnected in chain
                        if (self.grid[r][c].aa == 'H' or self.grid[r][c].aa == 'C') \
                                and (naa is 'H' and abs(ni-i) is not 1):
                            # save that bond
                            HHCbonds.append([i, ni])
                            score -= 1

                        # 1-point bond scenario H-C: check to right and unconnected in chain
                        elif (self.grid[r][c].aa == 'H' and naa is 'C') and abs(ni-i) is not 1:
                            # save that bond
                            HHCbonds.append([i, ni])
                            score -= 1

                        # 5-point bond scenario C-C: check to right and unconnected in chain
                        elif (self.grid[r][c].aa == 'C' and naa is 'C') and abs(ni-i) is not 1:
                            # save that H-bond
                            CCbonds.append([i, ni])
                            score -= 5

                    # if index does not go outside of grid range
                    if c + 1 < self.yran + 1:
                        # ni is next index
                        ni = self.grid[r][c + 1].i
                        # naa is next amino acid
                        naa = self.grid[r][c + 1].aa

                        # 1-point bond scenario H-H or C-H: check to right and unconnected in chain
                        if (self.grid[r][c].aa == 'H' or self.grid[r][c].aa == 'C') \
                                and (naa is 'H' and abs(ni-i) is not 1):
                            # save that bond
                            HHCbonds.append([i, ni])
                            score -= 1

                        # 1-point bond scenario H-C: check to right and unconnected in chain
                        elif (self.grid[r][c].aa == 'H' and naa is 'C') and abs(ni-i) is not 1:
                            # save that bond
                            HHCbonds.append([i, ni])
                            score -= 1

                        # 5-point bond scenario C-C: check to right and unconnected in chain
                        if (self.grid[r][c].aa == 'C' and naa is 'C') and abs(ni-i) is not 1:
                            # save that H-bond
                            CCbonds.append([i, ni])
                            score -= 5

            # save ultimate score
            self.score = score

        # return the places of the bonds (for visualization)
        return [HHCbonds, CCbonds] 