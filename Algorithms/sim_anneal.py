from Protein_class import *
from copy import deepcopy
from random import uniform

def sim_anneal(protein, ITER, AOM, COOL):
    algo = 'sa'

    scores = []

    # fold protein randomly until valid fold
    new = protein.rand_fold()

    # first fold is the only folding, so the current best
    best = deepcopy(new)
    all_time_best = deepcopy(new)

    # set impossible best score (minimum score = 0) so first fold always improves
    best.score = 1
    all_time_best.score = 1

    i = 1
    while i < ITER + 1:
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

            # temperature is acceptance probability: linear / hyperbolic cooling schedule
            # if the cooling schedule is linear
            if COOL is "linear":
                # temperature decreases linear
                temp = float(ITER - i)
            # if cooling schedule is hyperbolic
            elif COOL is "hyperbolic":
                # temperature decreases hyperbolic
                temp = ITER / i

            # if score improved
            if new.score < best.score:
                print("Improving!", new.score, "<", best.score)
                # save current best fold
                best = deepcopy(new)

            # if score deteriorated
            else:
                # acceptance chance for deteriorated fold
                chance = uniform(1, ITER)
                # if chance is below current temperature, accept deteriorated fold
                if chance < temp:
                    print("Fold accepted anyway!")
                    # save current best fold
                    best = deepcopy(new)

            # if best score improved
            if best.score < all_time_best.score:
                # save current all time best fold
                all_time_best = deepcopy(best)

            # continue to next iteration
            i += 1

    print("Simulated annealing: Best score = ", all_time_best.score)

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

    # output all time best protein
    return all_time_best
