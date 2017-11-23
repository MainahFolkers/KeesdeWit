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
        self.directions = []
        self.coordinates = []
        self.xs = []
        self.ys = []
        self.xran = 0
        self.yran = 0
        self.grid = []
        # stability score is calculated with score function
        self.score = 0

    def make_grid(self):
        xs = []
        ys = []
        # iterate over amino acids
        for i in range(self.n):
            # seperate coordinates in xs and ys
            xs.append(self.coordinates[i][0])
            ys.append(self.coordinates[i][1])

        xmin = min(xs)
        ymin = min(ys)
        # calculate ranges voor x and y
        self.xran = max(xs) - xmin
        self.yran = max(ys) - ymin

        self.grid = [[self.AA() for y in range(self.yran + 1)] for x in range(self.xran + 1)]

        # iterate over amino acids
        for i in range(self.n):
            # transform coordinates onto minimal grid
            xs[i] = xs[i] - xmin
            ys[i] = ys[i] - ymin

            # place amino acids onto grid
            self.grid[xs[i]][ys[i]].aa = self.chain[i]
            self.grid[xs[i]][ys[i]].i = i

        self.xs = xs
        self.ys = ys

    def calc_score(self):
        # lowerbound: worst score = 0
        # over estimation upperbound: best score = -H/2 - 1
        score = 0
        Hbonds = []

        # iterate over rows in grid
        for r in range(self.xran):
            # iterate over columns in grid
            for c in range(self.yran):
                # check if there is a H in cell of grid
                if self.grid[r][c].aa == 'H':
                    i = self.grid[r][c].i
                    # ni is next index
                    ni = self.grid[r + 1][c].i
                    # naa is next amino acid
                    naa = self.grid[r + 1][c].aa
                    # check if there is a H on the right
                    if naa is 'H' and abs(ni-i) is not 1:
                        Hbonds.append([i, ni])
                        score = score - 1
                    # ni is next index
                    ni = self.grid[r][c + 1].i
                    # naa is next amino acid
                    naa = self.grid[r][c + 1].aa
                    # check if there is a H above
                    if naa is 'H' and abs(ni-i) is not 1:
                        Hbonds.append([i, ni])
                        score = score - 1

        self.score = score
        return Hbonds
