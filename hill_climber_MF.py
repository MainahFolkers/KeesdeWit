from Protein_class_MF import *
from fold_MF import *
from plot_MF import *

def hill_climb(protein):

    best_score = 1

    # improvement: determine how many iterations
    for i in range(10):
        new = fold(protein)
        new.make_grid()
        new.calc_score()

        if best_score > new.score:
            print "Climbing!", best_score, ">", new.score
            best_score = new.score
            best = new

    visualize(best)
    print "Best score = ", best_score, best.score
    return best

# x = choice(for i in range(protein.n))
# directions = [r, d, u, l, u, d, r]
# direction[x] = andere richting

protein = Protein("HHPHHHPH")
protein = hill_climb(protein)
