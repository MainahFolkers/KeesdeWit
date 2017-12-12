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
        self.directions = ['' for i in range(self.n - 1)]
        # place first 2 AAs already known -> first bond always to right
        self.coordinates = [[] for i in range(self.n)]
        # so xs and ys also known
        self.xs = [0, 1]
        self.ys = [0, 0]
        # range x and y calculated in make_grid -> for dynamic grid
        self.xran = 0
        self.yran = 0
        self.grid = []
        # stability score is calculated with score function
        self.score = 0

    """
    rand_fold() folds the protein randomly while adhering to constraints.
    rand_fold() takes an argument of the class Protein, and an argument
    that specifies at what point folding should start.
    if no second argument is passed, folding starts from scratch.
    """
    def rand_fold(self):
        
        # first AA starts on (x = 0, y = 0)
        x, y = 0, 0
        self.coordinates[0] = [x, y]

        # first bond is to right
        self.directions[0] = 'r'

        # second AA always to right (x = x + 1, y = 0) which is (x = 1, y = 0)
        x = x + 1
        self.coordinates[1] = [x, y]

        # third AA can only bend up or to right
        direcs = ['r', 'u']
        direc = choice(direcs)

        # second bond is saved
        self.directions[1] = direc

        # save coordinates of third AA
        self.coordinates[2] = self.bend(direc, x, y)
        x = self.coordinates[2][0]
        y = self.coordinates[2][1]

        # directions a bond can move to (right, down, left, up)
        all_direcs = ['r', 'd', 'l', 'u']

        # most recently chosen direction is now the second bond
        chosen_direc = self.directions[1]

        # first three AAs already have coordinates -> start at fourth
        aa = 3

        # continue folding until end of chain
        while aa < self.n:
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
                # aa - 1 is index of bond
                self.directions[aa - 1] = chosen_direc
                # sample without replacement, remove chosen direction from options
                avail_direcs.remove(chosen_direc)

                ####### this line (function check_val: see line 127)
                self.check_val(x, y, aa - 1)

                ####### replaces the following lines:
                # # calculate new coordinates
                # [nx, ny] = self.bend(self.directions[aa - 1], x, y)
                # # if point on grid is not yet occupied by other amino acid
                # if [nx, ny] not in self.coordinates:
                #     # update x and y
                #     [x, y] = [nx, ny]
                #     # new coordinates are saved
                #     self.coordinates[aa] = [x, y]
                #     # continue with next amino acid
                #     aa += 1
                # # if folding is invalid
                # elif [nx, ny] in self.coordinates:
                #     return None

            # if there are no directions left to choose from
            if not avail_direcs:
                return None

        # return class protein
        return self

    """
    check_val() checks if folding is valid.
    if the coordinates (calculated with bend() using self.directions) 
    are not already occupied by another AA, the coordinate is assigned 
    to the current AA in self.coordinates.
    """
    def check_val(x, y, dir_to_check):

        # calculate next coordinates
        [nx, ny] = self.bend(self.directions[dir_to_check], x, y)

        # if point on grid is not yet occupied by other amino acid
        if [nx, ny] not in self.coordinates:
            # update old x and y
            [x, y] = [nx, ny]
            # new coordinates are saved
            self.coordinates[dir_to_check] = [x, y]

        # if folding is invalid
        elif [nx, ny] in self.coordinates:
            return None

    """
    mut_dir() randomly chooses a position in the chain to mutate
    thereby undoing the former folding.
    improvement: in different file, determine amount of mutations.
    """
    def mut_dir(self):

        # randomly choice position in protein to mutate
        # first bond cannot be chosen to avoid symmetry
        mut = choice(range(1, self.n - 1))

        # all direction options (right, down, left, up)
        direcs = ['r', 'd', 'l', 'u']
        # remove original direction on the chosen position from options
        direcs.remove(self.directions[mut])

        # random choice 1 of 3 directions
        direc = choice(direcs)
        # save mutated direction
        self.directions[mut] = direc

        # return class protein
        return self

    """
    mut_fold() folds mutated protein, with mutation position randomly
    determined by function mut_dir().

    --> --> --> --> maybe better to remove mut_fold() entirely and just 
    use protein.check_val(x, y, dir_to_check) in the hillclimber script?
    """
    def mut_fold(self):

        # first AA starts on (x = 0, y = 0)
        x, y = 0, 0
        self.coordinates = [[x, y]]

        # start at first bond
        b = 0

        # iterate over bonds
        while b < self.n - 1:

            ###### this line (see line 127 for function)
            self.check_val(x, y, b)

            ###### replaces the lines below:
            # # calculate next coordinates
            # [nx, ny] = self.bend(self.directions[b], x, y)
            # # if point on grid is not yet occupied by other amino acid
            # if [nx, ny] not in self.coordinates:
            #     # update old x and y
            #     [x, y] = [nx, ny]
            #     # new coordinates are safed
            #     self.coordinates[b] = [x, y]
            #     # continue with next bond
            #     b += 1
            # # if folding is invalid
            # elif [nx, ny] in self.coordinates:
            #     return None

        # return protein
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
    make_grid() works with the coordinates
    """
    def make_grid(self):

        # iterate over amino acids
        for aa in range(self.n):

            # seperate coordinates in arrays with xs and ys
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

    """
    calc_score() calculates the score by traversing the grid 
    and looking for adjacent Hs that are not connected in the chain.
    """
    def calc_score(self):

        # lowerbound: worst score = 0
        score = 0
        # upperbound: best score = -H / 2; where H is amount of H's in chain

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

                        # check if there is an H on right that is unconnected in chain
                        if naa is 'H' and abs(ni-i) is not 1:
                            # save that H-bond
                            Hbonds.append([i, ni])
                            # update score
                            score -= 1

                    # if index does not go outside of grid range
                    if c + 1 < self.yran + 1:
                        # ni is next index
                        ni = self.grid[r][c + 1].i
                        # naa is next amino acid
                        naa = self.grid[r][c + 1].aa

                        # check if there is an H above that is not connected in chain
                        if naa is 'H' and abs(ni-i) is not 1:
                            # save that H-bond
                            Hbonds.append([i, ni])
                            # update score
                            score -= 1

            # save ultimate score
            self.score = score

        # return the places of the H-bonds (for visualization)
        return Hbonds