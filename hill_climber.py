from Protein_class import *
from copy import deepcopy
import csv

def hill_climb(protein, ITER):
    # # open output file which stores best scores per iteration
    # ofile  = open('hill_scores.csv', "w", newline="")
    # writer = csv.writer(ofile, delimiter=",")

    # fold protein randomly
    new = protein.rand_fold()

    # this first folding is the only folding, so the current best
    best = deepcopy(new)

    # set impossible best score (minimum score = 0) so first fold always improves
    best.score = 1
    

    # array to see the score progress
    scores = []

    # improvement: determine how many iterations
    for i in range(ITER):
        # # print each 100th iteration
        # if i % (ITER/100) == 0:
        #     print(i)

        # new folding continues on current best fold
        new = deepcopy(best)

        # randomly mutate a number of directions
        new = new.mut_dir(14)

        # continue mutating and trying to fold until valid folding
        while new.mut_fold == None:
            new = new.mut_dir()

        # when valid mutated fold has been made: place protein on grid
        new.make_grid()
        # caluclate score
        new.calc_score()
        # add best score to scores array
        scores.append(best.score)

        # if score improved, continue climbing with this fold
        if new.score < best.score:
            # print to terminal
            print("Climbing!", new.score, "<", best.score)
            # save new best fold to continue
            best = deepcopy(new)

    # # write scores to output file
    # writer.writerow(scores)
    # # close output file
    # ofile.close()

    # print findings to terminal
    print("Hill climber: Best score = ", best.score)

    # output best protein
    return best