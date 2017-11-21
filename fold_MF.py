from random import *

# should this be a function in the Protein class?
def fold(protein):
    coordinates = [[] for i in range(protein.n)]
    protein.directions = ['' for i in range(protein.n - 1)]
    x = 0
    y = 0
    # first amino acid starts on (x=0, y=0)
    coordinates[0] = [x, y]
    # second amino acid is always on right (x=x+1, y=0) which represents (x=1, y=0)
    x = x + 1
    coordinates[1] = [x, y]
    # first bond is to right
    protein.directions[0] = 'r'

    # third amino acid can only bend up or to right
    direcs = ['r', 'u']
    direc = choice(direcs)
    coordinates[2] = bend(direc, x, y)
    x = coordinates[2][0]
    y = coordinates[2][1]
    # second bond is to right or up
    protein.directions[1] = direc

    # start at fourth amino acid
    i = 3
    # iterate over amino acids
    while i < protein.n:
        # improvement: remove one direction based on previous direction
        direcs = ['r', 'd', 'l', 'u']

        # while there are still directions left to choose
        while direcs:
            # random choice 1 of 4 directions
            direc = choice(direcs)
            # i - 1 is index of bond
            protein.directions[i - 1] = direc
            # sample without replacement, remove chosen direction from options
            direcs.remove(direc)
            # calculate new coordinates
            [nx, ny] = bend(direc, x, y)

            # if point on grid is not yet occupied by other amino acid
            if [nx, ny] not in coordinates:
                # new coordinates are safed
                [x, y] = [nx, ny]
                coordinates[i] = [x, y]

                # continue with next amino acid
                i = i + 1
                break

        # if there are no directions left to choose from
        if not direcs:
            return fold(protein)
    return coordinates

def bend(direc, x, y):
    # u is up
    if direc == 'u':
        y = y + 1
    # r is right
    elif direc == 'r':
        x = x + 1
    # d is down
    elif direc == 'd':
        y = y - 1
    # l is left, maybe use else instead of elif
    elif direc == 'l':
        x = x - 1
    return [x, y]
