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
        self.coordinates = []
        self.xran = 0
        self.yran = 0
        self.grid = self.make_grid(self.xran, self.yran)
        # stability score is calculated with score function
        self.score = 0

    def make_grid(self, xran, yran):
        self.xran = xran
        self.yran = yran
        self.grid = [[self.AA() for y in range(yran + 1)] for x in range(xran + 1)]

    def calc_score(self):
        # lowerbound: worst score = 0
        # over estimation upperbound: best score = (-2H - 2)/2 = -H - 1
        score = 0
        Hbonds = []

        # iterate over rows in grid
        for r in range(1, self.xran + 1):
            # iterate over columns in grid
            for c in range(1, self.yran + 1):
                # check if there is a H in cell of grid
                if self.grid[r - 1][c - 1].aa == 'H':
                    i = self.grid[r - 1][c - 1].i
                    # ni is next index
                    ni = self.grid[r][c - 1].i
                    # naa is next amino acid
                    naa = self.grid[r][c - 1].aa
                    # check if there is a H on the right
                    if r < self.xran and naa is 'H' and abs(ni-i) is not 1:
                        Hbonds.append([r -1, r])
                        Hbonds.append([c - 1, c - 1])
                        score = score - 1
                    # ni is next index
                    ni = self.grid[r - 1][c].i
                    # naa is next amino acid
                    naa = self.grid[r - 1][c].aa
                    # check if there is a H above
                    if c < self.yran and naa is 'H' and abs(ni-i) is not 1:
                        Hbonds.append([r - 1, r - 1])
                        Hbonds.append([c - 1, c])
                        score = score - 1
        self.score = score
        return Hbonds
