from Protein_class import *
from copy import deepcopy

def hill_climb(protein, ITER, AOM):
    # algorithm is hill climber for plot legend label
    algo = 'hc'

    scores = []

    # fold protein randomly until valid fold
    new = protein.rand_fold()

    # first fold is the only folding, so the current best
    best = deepcopy(new)

    # set impossible best score (minimum score = 0) so first fold always improves
    best.score = 1

    i = 0
    while i < ITER:
        # new folding continues on current best fold
        new = deepcopy(best)

        # randomly mutate a number of bonds
        new.mut_dir(AOM)

        # if the mutated fold is valid
        if new.mut_fold():
            # place protein on grid
            new.make_grid()
            # caluclate score
            new.calc_score()

            # save best score per iteration
            scores.append(best.score)

            # if score improved
            if new.score < best.score:
                print("Climbing!", new.score, "<", best.score, "i =", i)
                # save current best fold
                best = deepcopy(new)

            # continue to next iteration
            i += 1

    print("Hill climber: Best score = ", best.score)

    # open output file
    with open("ALGOS_" + protein.chain + ".txt", 'a+') as ofile:
        # write algorithm name as label
        ofile.write(algo + ',')
        # write comma seperated scores
        for score in scores:
            ofile.write(str(score) + ',')
        # write new line
        ofile.write('\n')
    # close output file
    ofile.close()

    # output best protein
    return best
