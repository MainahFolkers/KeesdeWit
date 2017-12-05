from random import *

# should this be a function in the Protein class?
def fold(protein):
    protein.coordinates = [[] for i in range(protein.n)]
    protein.directions = ['' for i in range(protein.n - 1)]
    x = 0
    y = 0

    # first amino acid starts on (x=0, y=0)
    protein.coordinates[0] = [x, y]

    # second amino acid is always on right (x=x+1, y=0) which represents (x=1, y=0)
    x = x + 1
    protein.coordinates[1] = [x, y]
    # first bond is to right
    protein.directions[0] = 'r'

    # third amino acid can only bend up or to right
    direcs = ['r', 'u']
    direc = choice(direcs)
    protein.coordinates[2] = protein.bend(direc, x, y)
    x = protein.coordinates[2][0]
    y = protein.coordinates[2][1]
    # second bond is to right or up
    protein.directions[1] = direc

    # start at fourth amino acid
    aa = 3
    # iterate over amino acids
    while aa < protein.n:
        # improvement: remove one direction based on previous direction
        direcs = ['r', 'd', 'l', 'u']

        # while there are still directions left to choose
        while direcs:
            # random choice 1 of 4 directions
            direc = choice(direcs)
            # i - 1 is index of bond
            protein.directions[aa - 1] = direc
            # sample without replacement, remove chosen direction from options
            direcs.remove(direc)
            # calculate new coordinates
            [nx, ny] = protein.bend(direc, x, y)

            # if point on grid is not yet occupied by other amino acid
            if [nx, ny] not in protein.coordinates:
                # update old x and y
                [x, y] = [nx, ny]
                # new coordinates are safed
                protein.coordinates[aa] = [nx, ny]

                # continue with next amino acid
                aa += 1
                break

        # if there are no directions left to choose from
        if not direcs:
            return None
    return protein
