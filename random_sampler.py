from Protein_class import *
from copy import deepcopy
import csv

def rand_samp(protein, ITER):
    # open output file which stores best scores per iteration
    ofile  = open('scores.csv', "w", newline="")
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

    # improvement: determine how many iterations
    i = 0
    while i <ITER:

        # print each 100th iteration
        if i % (ITER/100) == 0:
            print(i)

        # fold protein randomly again
        new = protein.rand_fold()

        # if fold is valid
        if new:
            # set protein on grid
            new.make_grid()
            # caluclate score
            new.calc_score()
            # write i and best.score to output file
            writer.writerow([i, best.score])

            # continue to next iteration
            i += 1

            if new.score < best.score:
                print("Improving!", new.score, "<", best.score)
                # save new best folding
                best = deepcopy(new)

        # if fold is invalid
        elif not new:
            # try to fold again
            continue

    # close output file
    ofile.close()

    print("Random sampler: Best score = ", best.score)
    return best
