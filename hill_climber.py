from Protein_class import *
from copy import deepcopy

def hill_climb(protein, ITER, AOM = 1):
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

        # fold mutated protein
        new = new.mut_fold()

        # if the folding is valid
        if new:
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
                print("Climbing!", new.score, "<", best.score)
                # save new best fold
                best = deepcopy(new)

    print("Hill climber: Best score = ", best.score, "AOM = ", AOM)

    with open("AOM=" + str(AOM) + '_' + protein.chain +".txt", 'a+') as ofile:
        ofile.write("AOM=" + str(AOM) + ',')
        for score in scores:
            ofile.write(str(score) + ',')
        ofile.write('\n')
    ofile.close()

    # output best protein
    return best
