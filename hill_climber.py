from Protein_class import *
from copy import deepcopy
import csv

# improvement: in different file, determine amount of mutations
def mut_dir(protein):
    # randomly choice position in protein to mutate
    # the first bond cannot be chosen to avoid symmetry
    mut = choice(range(1, protein.n - 1))
    # directions a bond can move to (right, down, left, up)
    direcs = ['r', 'd', 'l', 'u']
    # remove original direction from options
    direcs.remove(protein.directions[mut])
    # random choice 1 of 3 directions
    direc = choice(direcs)
    # save mutated direction
    protein.directions[mut] = direc

    return protein

def mut_fold(protein):
    # first amino acid starts at x=0 and y = 0
    x, y = 0, 0
    # save first coordinates
    protein.coordinates = [[x, y]]
    # start at first bond
    b = 0
    # iterate over bonds
    while b < protein.n - 1:

        # calculate next coordinates
        [nx, ny] = protein.bend(protein.directions[b], x, y)

        # if point on grid is not yet occupied by other amino acid
        if [nx, ny] not in protein.coordinates:
            # update old x and y
            [x, y] = [nx, ny]
            # new coordinates are safed
            protein.coordinates.append([x, y])

            # continue with next bond
            b += 1

    # if folding is invalid
    elif [nx, ny] in protein.coordinates:
            return None

    return protein

def hill_climb(protein, ITER):
    # open output file which stores best scores per iteration
    ofile  = open('hill_scores.csv', "w", newline="")
    writer = csv.writer(ofile, delimiter=",")

    # fold protein randomly
    new = protein.rand_fold()
    # until valid folding
    while not new:
        new = protein.rand_fold()

    best = deepcopy(new)
    # first folding will always be an improvement
    # because minimum score = 0
    best.score = 1

    scores = []

    i = 0
    # improvement: determine how many iterations
    while i < ITER:

        # print each 100th iteration
        if i % (ITER/100) == 0:
            print(i)

        new = deepcopy(best)
        # mutate direction
        new = mut_dir(new)

        # if fold is valid after mutation
        if mut_fold(new):
            # set protein on grid
            new.make_grid()
            # caluclate score
            new.calc_score()
            # add best score to scores array
            scores.append(best.score)

            # continue to next iteration
            i += 1

            if new.score < best.score:
                print("Climbing!", new.score, "<", best.score)
                # save new best folding
                best = deepcopy(new)

        # if mutated fold is invalid
    elif not mut_fold(new):
            # try to fold again
            continue

    # write scores to output file
    writer.writerow(scores)
    # close output file
    ofile.close()

    print("Hill climber: Best score = ", best.score)
    return best
