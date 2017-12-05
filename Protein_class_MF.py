from random import choice

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
        self.coordinates = [[0, 0], [1, 0]]
        self.xs = []
        self.ys = []
        self.xran = 0
        self.yran = 0
        self.grid = []
        # stability score is calculated with score function
        self.score = 0

    def fold(self):
        x = 0
        y = 0
        self.coordinates = [[x, y]]
        # start at first bond
        b = 0
        # iterate over bonds
        while b < self.n - 1:

            # calculate next coordinates
            [nx, ny] = self.bend(self.directions[b], x, y)

            # if point on grid is not yet occupied by other amino acid
            if [nx, ny] not in self.coordinates:
                # update old x and y
                [x, y] = [nx, ny]
                # new coordinates are safed
                self.coordinates.append([x, y])

                # continue with next bond
                b += 1

            # if folding is invalid
            elif [nx, ny] in self.coordinates:
                return None

        return self

    def bend(self, direc, x, y):
        # u is up
        if direc == 'u':
            y = y + 1
        # r is right
        elif direc == 'r':
            x = x + 1
        # d is down
        elif direc == 'd':
            y = y - 1
        # l is left
        elif direc == 'l':
            x = x - 1
        return [x, y]

    def make_grid(self):
        xs = []
        ys = []
        # iterate over amino acids
        for aa in range(self.n):
            # seperate coordinates in xs and ys
            xs.append(self.coordinates[aa][0])
            ys.append(self.coordinates[aa][1])

        xmin = min(xs)
        ymin = min(ys)
        # calculate ranges voor x and y
        self.xran = max(xs) - xmin
        self.yran = max(ys) - ymin

        self.grid = [[self.AA() for y in range(self.yran + 1)] for x in range(self.xran + 1)]

        # iterate over amino acids
        for aa in range(self.n):
            # transform coordinates onto minimal grid
            xs[aa] -= xmin
            ys[aa] -= ymin

            # place amino acids onto grid
            self.grid[xs[aa]][ys[aa]].aa = self.chain[aa]
            self.grid[xs[aa]][ys[aa]].i = aa

        self.xs = xs
        self.ys = ys

    def calc_score(self):
        # lowerbound: worst score = 0
        score = 0

        Hbonds = []

        # iterate over rows in grid
        for r in range(self.xran + 1):
            # iterate over columns in grid
            for c in range(self.yran + 1):
                # check if there is a H in cell of grid
                if self.grid[r][c].aa == 'H':
                    i = self.grid[r][c].i
                    # if index does not go outside of grid range
                    if r + 1 < self.xran + 1:
                        # ni is next index
                        ni = self.grid[r + 1][c].i
                        # naa is next amino acid
                        naa = self.grid[r + 1][c].aa
                        # check if there is a H on the right
                        if naa is 'H' and abs(ni-i) is not 1:
                            Hbonds.append([i, ni])
                            score -= 1
                    # if index does not go outside of grid range
                    if c + 1 < self.yran + 1:
                        # ni is next index
                        ni = self.grid[r][c + 1].i
                        # naa is next amino acid
                        naa = self.grid[r][c + 1].aa
                        # check if there is a H above
                        if naa is 'H' and abs(ni-i) is not 1:
                            Hbonds.append([i, ni])
                            score -= 1
        self.score = score
        return Hbonds
