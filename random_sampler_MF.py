from Protein_class_MF import *
from fold_MF import *
from visualize_MF import *
import matplotlib.pyplot as plt
from copy import deepcopy

def rand_samp(protein, ITER):
    new = fold(protein)
    while not new:
        new = fold(protein)

    best = deepcopy(new)
    # first folding will always be an improvement
    # because minimum score = 0
    best.score = 1

    scores = []

    # improvement: determine how many iterations
    i = 0
    while i <ITER:

        new = fold(protein)

        # if fold is valid after mutation continue to next iteration
        if new:
            new.make_grid()
            new.calc_score()
            scores.append(best.score)
            print(round(i/ITER * 100) , '%')
            i += 1

            if new.score < best.score:
                print("Improving!", new.score, "<", best.score)
                best = deepcopy(new)

        elif not new:
            continue

    # plot best score per iteration
    plt.plot(range(i), scores)

    visualize(best)
    print("Best score = ", best.score)
    return best
