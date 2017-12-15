from Protein_class import *
from copy import deepcopy
from random import randint

def sim_anneal(protein, ITER, AOM = 1):
    scores = []

    # fold protein randomly
    new = protein.rand_fold()

    # this first folding is the only folding, so the current best
    best = deepcopy(new)

    # set impossible best score (minimum score = 0) so first fold always improves
    best.score = 1

    i = 1

    # in temp: dividing by zero impossible, so start at 1
    while i < ITER + 1
        # new folding continues on current best fold
        new = deepcopy(best)

        # randomly mutate a number of directions
        new = new.mut_dir(AOM)

        # continue mutating and trying to fold until valid folding
        if new.mut_fold():
            # when valid mutated fold has been made: place protein on grid
            new.make_grid()
            # caluclate score
            new.calc_score()

            # save best score per iteration
            scores.append(best.score)

            # continue to next iteration
            i += 1

        # temperature is acceptance probability -> for now linearly declining
        temp = 100 / i

        # if score improved, accept new fold as best fold
        if new.score < best.score:
            print("Improving!", new.score, "<", best.score)
            best = deepcopy(new)
        # if score deteriorated, accept new fold according to temperature
        else:
            # determine acceptance for this worse score
            chance = randint(0, 100)
            # if change is below current temperature, accept worse fold
            if chance < temp:
                print("Fold accepted anyway!")
                best = deepcopy(new)

    print("Simulated annealing: Best score = ", best.score)

    # output best protein
    return best
