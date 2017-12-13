from Protein_class_MF import *
from visualize_MF import *
from fold_MF import *
import matplotlib.pyplot as plt
from copy import deepcopy

def mut_dir(protein):
    # randomly choice position in protein to mutate
    # can the first bond be chosen? no
    mut = choice(range(1, protein.n - 1))
    direcs = ['r', 'd', 'l', 'u']
    # remove original direction from options
    direcs.remove(protein.directions[mut])
    # random choice 1 of 3 directions
    direc = choice(direcs)
    protein.directions[mut] = direc

    return protein

def hill_climb(protein, ITER):
    # initial folding with fold_MF.py
    # improvement: do initial fold with class function
    new = fold(protein)
    while not new:
        new = fold(protein)

    best = deepcopy(new)
    # first folding will always be an improvement
    # because minimum score = 0
    best.score = 1

    scores = []

    i = 0
    # improvement: determine how many iterations
    while i < ITER:

        new = deepcopy(best)
        new = mut_dir(new)

        # if fold is valid after mutation continue to next iteration
        if new.fold():
            new.make_grid()
            new.calc_score()
            scores.append(best.score)
            print(round(i/ITER * 100) , '%')
            i += 1

            if new.score < best.score:
                print("Climbing!", new.score, "<", best.score)
                best = deepcopy(new)

        elif not new.fold():
            continue

    # plot best score per iteration
    plt.plot(range(i), scores)

    visualize(best)
    print("Best score = ", best.score)
    return best 
