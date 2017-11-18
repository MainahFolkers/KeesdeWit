from random import *

def fold(protein):
    coordinates = [[] for i in range(protein.n)]
    x = 0
    y = 0
    # first amino acid starts on (x=0, y=0)
    coordinates[0] = [x, y]
    # second amino acid is always on right (x=x+1, y=0) which represents (x=1, y=0)
    x = x + 1
    coordinates[1] = [x, y]

    # start at third amino acid
    i = 2
    # iterate over amino acids
    while i < protein.n:
        # improvement: remove one direction based on previous direction
        direcs = ['r', 'd', 'l', 'u']

        while direcs:
            # random choice 1 of 4 directions
            direc = choice(direcs)
            # sample without replacement, remove chosen direction from options
            direcs.remove(direc)

            # calculates new coordinates
            [nx, ny] = direction(direc, x, y)

            # if point on grid is not yet occupied by other amino acid
            if [nx, ny] not in coordinates:
                # new coordinates are safed
                [x, y] = [nx, ny]
                coordinates[i] = [x, y]

                # continue with next amino acid
                i = i + 1
                break

        if not direcs:
            return None
    return coordinates

def direction(direc, x, y):
    # u is up
    if direc == 'u':
        y = y - 1
    # r is right
    elif direc == 'r':
        x = x + 1
    # d is down
    elif direc == 'd':
        y = y + 1
    # l is left, maybe use else instead of elseif
    elif direc == 'l':
        x = x - 1
    return [x, y]
