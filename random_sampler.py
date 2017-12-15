from Protein_class import *
from copy import deepcopy

def rand_samp(protein, ITER):
    scores = []

    # fold protein randomly until valid folding
    new = protein.rand_fold()

    # we have only one fold, so that is the current best
    best = deepcopy(new)

    # set impossible best score (minimum score = 0) so first fold always improves
    best.score = 1

    # improvement: determine how many iterations
    for i in range(ITER):
        # fold protein randomly again
        new = protein.rand_fold()

        # set protein on grid
        new.make_grid()
        # caluclate score
        new.calc_score()

        # save best score per iteration
        scores.append(best.score)

        # is score improved
        if new.score < best.score:
            print("Improving!", new.score, "<", best.score)
            # save new best folding
            best = deepcopy(new)

    print("Random sampler: Best score = ", best.score, "ITER = ", ITER)

    return best
