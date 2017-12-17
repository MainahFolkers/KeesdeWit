from random import choice
from copy import deepcopy
import matplotlib.pyplot as plt

class Protein:
    class AA():
        # aa is type of amino acid (H/P, can later expand with C)
        aa_type = None
        # i is index in amino acid chain
        i = None

    def __init__(self, chain):
        # chain is string of amino acid sequence
        self.chain = chain
        # n is amount of amino acids in protein chain
        self.n = len(self.chain)
        # directions and coordinates assigned in fold functions
        self.directions = []
        self.coordinates = []
        # coordinates components x and y
        self.xs, self.ys = [], []
        # x and y ranges calculated in make_grid()
        self.xran, self.yran = 0, 0
        self.grid = []
        # protein stability is calculated with score function
        self.score = 0

    """
    rand_fold() folds the protein randomly while adhering to constraints.
    """
    def rand_fold(self):
        while True:
            # reset coordinates and directions for every new fold
            self.coordinates = []
            self.directions = ['' for i in range(self.n - 1)]

            # first amino acid starts on (x = 0, y = 0)
            self.x, self.y = 0, 0

            self.coordinates.append([self.x, self.y])

            # first bond is to right
            self.directions[0] = 'r'

            # second AA always to right (x = x + 1, y = 0) which is (x = 1, y = 0), to avoid symmetry
            self.x += 1
            self.coordinates.append([self.x, self.y])

            # third AA can only bend up or to right
            avail_direcs = ['r', 'u']
            chosen_direc = choice(avail_direcs)

            # second bond is saved
            self.directions[1] = chosen_direc

            # save coordinates of third amino acid
            self.coordinates.append(self.bend(chosen_direc, self.x, self.y))
            self.x = self.coordinates[2][0]
            self.y = self.coordinates[2][1]

            # directions a bond can move to (right, down, left, up)
            all_direcs = ['r', 'd', 'l', 'u']

            # first three amino acids already have coordinates -> start at fourth
            aa = 3

            # iterate over amino acids
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
                    b = aa - 1
                    self.directions[b] = chosen_direc

                    # sample without replacement, remove chosen direction from options
                    avail_direcs.remove(chosen_direc)

                    # if chosen direction yields valid coordinates
                    if self.check_val(aa, b):
                        if aa is self.n - 1:
                            return self
                        # continue to next amino acid
                        aa += 1
                        break

                # if there are no directions left to choose from
                if not avail_direcs:
                    # stop current folding
                    break

    """
    check_val() checks if folding is valid.
    If the proposed coordinates (calculated from self.directions using self.bend()) are
    not already occupied, the coordinate is assigned to current amino acid in self.coordinates.
    NOTE: be careful with input, as it differs between random sampler and hillclimber.
    """
    def check_val(self, aa, b):
        # calculate next coordinates according to selected direction
        [nx, ny] = self.bend(self.directions[b], self.x, self.y)

        # if those coordinates on grid are not yet occupied by other amino acid
        if [nx, ny] not in self.coordinates:

            # update old x and y
            [self.x, self.y] = [nx, ny]

            # save the new coordinates
            self.coordinates.append([nx, ny])
            return True

    """
    mut_dir() randomly mutates direction(s) in the chain.
    """
    def mut_dir(self, AOM):
        # mutate a number of positions in the protein
        for mutation in range(AOM):
            # randomly choice position (not first bond because of symmetry)
            mut = choice(range(1, self.n - 1))

            # directions a bond can move to (right, down, left, up)
            avail_direcs = ['r', 'd', 'l', 'u']
            # original direction cannot be chosen again
            avail_direcs.remove(self.directions[mut])

            # replace original direction with mutation
            self.directions[mut] = choice(avail_direcs)

    """
    mut_fold() folds mutated protein, with mutation position randomly
    determined by function mut_dir().
    """
    def mut_fold(self):
        # first amino acid starts on (x = 0, y = 0)
        self.x, self.y = 0, 0
        self.coordinates = [[self.x, self.y]]

        # start at second amino acid
        aa = 1

        # iterate over amino acids
        while aa < self.n:
            # bond index
            b = aa - 1
            # check whether the direction yields a valid coordinate
            if self.check_val(aa, b):
                aa += 1
            # if mutation yiels invalid coordinate
            else:
                return None

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
        # empty arrays with x and y components of coordinates
        self.xs = []
        self.ys = []

        # iterate over amino acids
        for aa in range(self.n):
            # seperate coordinates in arrays for xs and for ys
            self.xs.append(self.coordinates[aa][0])
            self.ys.append(self.coordinates[aa][1])

        # define minima
        xmin = min(self.xs)
        ymin = min(self.ys)

        # calculate ranges voor x and y
        self.xran = max(self.xs) - xmin
        self.yran = max(self.ys) - ymin

        # make minimal grid (based on current x and y ranges)
        self.grid = [[self.AA() for y in range(self.yran + 1)] for x in range(self.xran + 1)]

        # iterate over amino acids in chain
        for aa in range(self.n):

            # transform coordinates onto minimal grid
            self.xs[aa] -= xmin
            self.ys[aa] -= ymin

            # place amino acids onto grid
            self.grid[self.xs[aa]][self.ys[aa]].aa_type = self.chain[aa]
            self.grid[self.xs[aa]][self.ys[aa]].i = aa

    """
    calc_score() calculates the score by traversing the grid
    looking for adjacent Hs that are not connected in the chain.
    """
    def calc_score(self):
        # lowerbound: worst score = 0
        score = 0

        # iterate over rows in grid
        for r in range(self.xran + 1):
            # iterate over columns in grid
            for c in range(self.yran + 1):
                # check if there is an amino acid of interest (H/C)
                if self.grid[r][c].aa_type is 'H' or self.grid[r][c].aa_type is 'C':
                    # the cell index that the H/C is in
                    i = self.grid[r][c].i

                    # if row index does not go outside of grid range
                    if r + 1 < self.xran + 1:
                        # ni is next index
                        ni = self.grid[r + 1][c].i
                        # naa is next amino acid
                        naa = self.grid[r + 1][c].aa_type

                        # 1-point bond scenario H-H: check to right and unconnected in chain
                        if self.grid[r][c].aa_type is 'H' and naa is 'H':
                            if abs(ni-i) is not 1:
                                score -= 1

                        # 1-point bond scenario H-C: check to right and unconnected in chain
                        elif self.grid[r][c].aa_type is 'H' and naa is 'C':
                            if abs(ni-i) is not 1:
                                score -= 1

                        # 1-point bond scenario C-H: check to right and unconnected in chain
                        elif self.grid[r][c].aa_type is 'C' and naa is 'H':
                            if abs(ni-i) is not 1:
                                score -= 1

                        # 5-point bond scenario C-C: check to right and unconnected in chain
                        elif self.grid[r][c].aa_type is 'C' and naa is 'C':
                            if abs(ni-i) is not 1:
                                score -= 5

                    # if column index does not go outside of grid range
                    if c + 1 < self.yran + 1:
                        # ni is next index
                        ni = self.grid[r][c + 1].i
                        # naa is next amino acid
                        naa = self.grid[r][c + 1].aa_type

                        # 1-point bond scenario H-H: check above and unconnected in chain
                        if self.grid[r][c].aa_type is 'H' and naa is 'H':
                            if abs(ni-i) is not 1:
                                score -= 1

                        # 1-point bond scenario H-C: check above and unconnected in chain
                        if self.grid[r][c].aa_type is 'H' and naa is 'C':
                            if abs(ni-i) is not 1:
                                score -= 1

                        # 1-point bond scenario C-H: check above and unconnected in chain
                        elif self.grid[r][c].aa_type is 'C' and naa is 'H':
                            # ni could also be None
                            if abs(ni-i) is not 1:
                                score -= 1

                        # 5-point bond scenario C-C: check above and unconnected in chain
                        elif self.grid[r][c].aa_type is 'C' and naa is 'C':
                            if abs(ni-i) is not 1:
                                score -= 5

            # save ultimate score
            self.score = score

    def visualize(self, AOM, run):
        plt.figure()

        # iterate over amino acids
        for aa in range(self.n):
            # determine color for point on plot (later expand with C)
            if self.chain[aa] is 'H':
                col = 'red'
            elif self.chain[aa] is 'C':
                col = 'green'
            elif self.chain[aa] is 'P':
                col = 'blue'

            plt.scatter(self.xs[aa], self.ys[aa], s=120, zorder=2, color=col)
            plt.annotate(aa, xy=(self.xs[aa], self.ys[aa]), xytext=(self.xs[aa] + 0.05, self.ys[aa] + 0.05), fontsize=10)

        # plot black line behind / between points
        plt.plot(self.xs, self.ys, lw=3, zorder=1, color='black')

        # show grid on plot
        plt.grid(b=True)

        # set tick spacing to 1
        plt.xticks(range(min(self.xs), self.xran + 1))
        plt.yticks(range(min(self.ys), self.yran  + 1))

        plt.title("Score = " + str(self.score), fontsize=30)

        # save figure
        plt.savefig(self.chain + '_score=' + str(self.score) + '_AOM=' + str(AOM) + '_run=' + str(run) + '.jpg')
