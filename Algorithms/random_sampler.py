from Protein_class import *
from copy import deepcopy

def rand_samp(protein, ITER):
    scores = []

    # fold protein randomly until valid fold
    new = protein.rand_fold()

    # we have only one fold, so that is the current best
    best = deepcopy(new)

    # set impossible best score (minimum score = 0) so first fold always improves
    best.score = 1

    for i in range(ITER):
        # fold protein randomly again
        new = protein.rand_fold()

        # place protein on grid
        new.make_grid()
        # caluclate score
        new.calc_score()

        # save best score per iteration
        scores.append(best.score)

        # is score improved
        if new.score < best.score:
            print("Improving!", new.score, "<", best.score)
            # save current best fold
            best = deepcopy(new)

    print("Random sampler: Best score = ", best.score)

    # open output file
    with open("ALGOS_"protein.chain +".txt", 'a+') as ofile:
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
