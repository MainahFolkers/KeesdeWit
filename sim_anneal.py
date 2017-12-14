from Protein_class import *
from copy import deepcopy
from random import randint
import csv

def sim_anneal(protein, ITER):

    # fold protein randomly
    new = protein.rand_fold()

    # this first folding is the only folding, so the current best
    best = deepcopy(new)

    # set impossible best score (minimum score = 0) so first fold always improves
    best.score = 1
    
    # array to see the score progress
    scores = []

    # in temp: dividing by zero impossible, so start at 1
    for i in range(1, ITER + 1):

        # mutate some directions in current best fold
        new = deepcopy(best)
        new = new.mut_dir(3)

        # calculate and save new score
        new.make_grid()
        new.calc_score()
        scores.append(best.score)

        # temperature is acceptance probability -> for now linearly declining
        temp = 100 / i

        # if score improved, accept new fold as best fold
        if new.score < best.score:
            print("score improved, new fold accepted as best fold!")
            print("new score is:")
            print(new.score)
            best = deepcopy(new)
        # if score deteriorated, accept new fold according to temperature
        else:
            # determine acceptance for this worse score
            chance = randint(0, 100)
            print("chance is: ")
            print(chance)
            # if change is below current temperature, accept worse fold
            if chance < temp:
                print("aangenomen yaaay")
                best = deepcopy(new)
            else:
                print("niet aangenomen aaahww")

    # print findings to terminal
    print("Simulated annealing: Best score = ", best.score)

    # output best protein
    return best