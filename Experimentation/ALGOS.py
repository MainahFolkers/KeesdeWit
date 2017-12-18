from Protein_class import *
from random_sampler import *
from hill_climber import *
from sim_anneal import *
#from depth_first_search import *
from copy import deepcopy
import matplotlib.pyplot as plt
from time import time

start_time = time()

chains = [
    "HHPHHHPH",
    #"HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH",
    #"HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    ]

# amount of iterations
ITER = 10
# amount of runs
RUNS = 2
# amount of mutations
AOM = 5
# cooling schedule
COOL = "hyperbolic"

for chain in chains:
    best_score = 1

    # data series for y-axis in legend
    ys = []

    print("Protein chain =", chain)
    protein = Protein(chain)

    for algo in ['rs', 'hc', 'sa']: #, 'df']:
        # sample size for mean decline in best score per iteration
        for run in range(1, RUNS + 1):
            #print("Run =", run)

            if algo is 'rs':
                protein = rand_samp(protein, ITER)
                protein.visualize('NA', 'NA', run)
            elif algo is 'hc':
                protein = hill_climb(protein, ITER, AOM)
                protein.visualize(AOM, 'NA', run)
            elif algo is 'sa':
                protein = sim_anneal(protein, ITER, AOM, COOL)
                protein.visualize(AOM, COOL, run)
            #elif algo is 'df':
                #protein = depth_first_search(protein)
                #protein.visualize('NA', 'NA', run)

            if protein.score < best_score:
                best_score = deepcopy(protein.score)
                print("Current best score = ", best_score)

    with open("ALGOS_"protein.chain + ".txt", 'r') as ifile:
        data = ifile.readlines()
        for line in data:
            # remove comma
            scores = line.split(',')
            # remove \n
            scores = scores[:-1]
            # covert string to int
            scores[1:] = [ int(score) for score in scores[1:] ]
            ys.append(scores)
    ifile.close()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xlabel("Iterations")
    plt.ylabel("Best score")

    ys_mean = []
    # iterate over runs
    for i in range(len(ys) - RUNS + 1):
        # if new AOM starts
        if i % RUNS is 0:
            y_mean = []
            # append AOM legend label
            y_mean.append(ys[i][0])
            # iterate over scores of 1 run
            for j in range(1, len(ys[i])):
                score_sum = 0
                # iterate over ys with same AOM
                for k in range(i, i + RUNS):
                    score_sum += ys[k][j]
                score_mean = score_sum / RUNS
                y_mean.append(score_mean)
            ys_mean.append(y_mean)

    i = len(ys_mean)
    for y_mean in ys_mean:
        # x is amount of iterations
        x = range(len(y_mean) - 2)
        # first element of y is label
        ax.plot(x, y_mean[2:], label=y_mean[0], zorder=i)
        i -= 1

    plt.legend(title="Cooling schedule", bbox_to_anchor=(1.04,1), loc="upper left", borderaxespad=0)
    plt.savefig("ALGOS_"protein.chain + ".jpg", bbox_inches = 'tight')

minutes = (time() - start_time) / 60
print(minutes, " minutes total")
