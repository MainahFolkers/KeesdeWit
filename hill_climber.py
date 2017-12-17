from Protein_class import *
from copy import deepcopy

def hill_climb(protein, ITER, AOM):
    scores = []

    # fold protein randomly until valid folding
    new = protein.rand_fold()

    # this first folding is the only folding, so the current best
    best = deepcopy(new)

    # set impossible best score (minimum score = 0) so first fold always improves
    best.score = 1

    i = 0
    while i < ITER:
        # new folding continues on current best fold
        new = deepcopy(best)

        # randomly mutate a number of directions
        new.mut_dir(AOM)

        # if the folding is valid
        if new.mut_fold():
            # place protein on grid
            new.make_grid()
            # caluclate score
            new.calc_score()

            # save best score per iteration
            scores.append(best.score)

            # continue to next iteration
            i += 1

            # if score improved
            if new.score < best.score:
                print("Climbing!", new.score, "<", best.score, "i =", i)
                # save new best fold
                best = deepcopy(new)

    #print("Hill climber: Best score = ", best.score)

    with open(protein.chain +".txt", 'a+') as ofile:
        ofile.write(str(AOM) + ',')
        for score in scores:
            ofile.write(str(score) + ',')
        ofile.write('\n')
    ofile.close()

    # output best protein
    return best
