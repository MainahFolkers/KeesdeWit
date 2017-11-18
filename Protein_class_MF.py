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
        self.grid = [[self.AA() for y in range(yran)] for x in range(xran)]


    def calc_score(self, xmax, ymax):
        score = 0
        # iterate over rows in grid
        for r in range(xmax):
            # iterate over columns in grid
            for c in range(ymax):
                # check if there is a H in cell of grid
                if self.grid[r][c].aa == 'H':
                    i = self.grid[r][c].i
                    ni = self.grid[r + 1][c].i
                    naa = self.grid[r + 1][c].aa
                    # check if there is a H on the right
                    #if r + 1 < xmax and naa is 'H' and abs(i-ni) is not 1:
                    if r + 1 < xmax and naa is 'H':
                        score = score - 1
                    ni = self.grid[r][c + 1].i
                    naa = self.grid[r][c + 1].aa
                    # check if there is a H above
                    #if c + 1 < ymax and naa is 'H' and abs(i-ni) is not 1:
                    if c + 1 < ymax and naa is 'H':
                        score = score - 1
        self.score = score
