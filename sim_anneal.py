from Protein_class import *
from copy import deepcopy
from random import uniform

def sim_anneal(protein, ITER, AOM, COOL):
    scores = []

    # fold protein randomly
    new = protein.rand_fold()

    # this first folding is the only folding, so the current best
    best = deepcopy(new)

    # set impossible best score (minimum score = 0) so first fold always improves
    best.score = 1
    best_score = 1

    # in temp: dividing by zero impossible, so start at 1
    i = 0
    while i < ITER:
        # mutate some directions in current  best fold
        new = deepcopy(best)
        new.mut_dir(AOM)

        # continue mutating and trying to fold until valid folding
        if new.mut_fold():
            # when valid mutated fold has been made: place protein on grid
            new.make_grid()
            # caluclate score
            new.calc_score()

            # save best score per iteration
            scores.append(best_score)

            # continue to next iteration
            i += 1

            # temperature is acceptance probability: linear / hyperbolic cooling schedule
            if COOL is "linear":
                temp = (ITER - i) / 100
            elif COOL is "hyperbolic":
                temp = 100 / i

            # if score improved, accept new fold as best fold
            if new.score < best.score:
                print("Improving!", new.score, "<", best.score)
                best = deepcopy(new)
                best_score = deepcopy(best.score)
            # if score deteriorated, accept new fold according to temperature
            else:
                # determine acceptance for this worse score
                chance = uniform(0, 100)
                # if change is below current temperature, accept worse fold
                if chance < temp:
                    print("Fold accepted anyway!")
                    best = deepcopy(new)

            if best_score < all_time_best.score:
                all_time_best = deepcopy(best)

    print("Simulated annealing: Best score = ", all_time_best.score)

    with open("COOL_" protein.chain + "AOM= " + str(AOM) +".txt", 'a+') as ofile:
        ofile.write(COOL + ',')
        for score in scores:
            ofile.write(str(score) + ',')
        ofile.write('\n')
    ofile.close()

    # output best protein
    return all_time_best
