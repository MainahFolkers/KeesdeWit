from Protein_class import *
from copy import deepcopy
import csv

def hill_climb(protein, ITER):
    # open output file which stores best scores per iteration
    ofile  = open('hill_scores.csv', "w", newline="")
    writer = csv.writer(ofile, delimiter=",")

    # randomly fold the protein
    new = protein.rand_fold()

    # this first folding is the only folding, so the current best
    best = deepcopy(new)

    # best score = 1 so that even if first mutated fold = 0 so first fold is always improvement
    best.score = 1

    # array to see the score progress
    scores = []

    # improvement: determine how many iterations
    while i in range(ITER):

        # print each 100th iteration
        if i % (ITER/100) == 0:
            print(i)

        # new folding continues on current best fold
        new = deepcopy(best)

        # randomly mutate a direction
        new = new.mut_dir()
        # continue mutating and trying to fold until valid folding
        while not new.mut_fold:
            new = new.mut_dir()

        # # if fold is valid after mutation
        # if new.mut_fold():

        # when valid mutated fold has been made: place protein on grid
        new.make_grid()
        # caluclate score
        new.calc_score()
        # add best score to scores array
        scores.append(best.score)

        # # continue to next iteration
        # i += 1

        # if score improved, continue climbing with this fold
        if new.score < best.score:
            print("Climbing!", new.score, "<", best.score)
            # save new best fold
            best = deepcopy(new)

    # write scores to output file
    writer.writerow(scores)
    # close output file
    ofile.close()

    # print findings to terminal
    print("Hill climber: Best score = ", best.score)

    # output best protein
    return best