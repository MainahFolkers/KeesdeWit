from Protein_class import *
from sim_anneal import *
from copy import deepcopy
import matplotlib.pyplot as plt
import time

start_time = time.time()

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
AOM = 1

start_chain = time.time()

for chain in chains:
    best_score = 1

    # data series for y-axis in legend
    ys = []

    print("Protein chain =", chain)
    protein = Protein(chain)

    for cool in ["linear", "hyperbolic"]:
        # sample size for mean decline in best score per iteration
        for run in range(1, RUNS + 1):
            #print("Run =", run)
            protein = sim_anneal(protein, ITER, AOM, cool)
            protein.visualize(AOM, cool, run)

            if protein.score < best_score:
                best_score = deepcopy(protein.score)
                print("Current best score = ", best_score)

    with open("COOL_"protein.chain +".txt", 'r') as ifile:
        data = ifile.readlines()
        for line in data:
            # remove comma
            scores = line.split(',')
            # remove \n
            scores = scores[:-1]
            # covert string to int
            scores = [ int(score) for score in scores[1:] ]
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
    plt.savefig("COOL_" + protein.chain + ".jpg", bbox_inches = 'tight')

    minutes = (time.time() - start_chain) / 60
    print(minutes, " minutes for ", protein.chain)
    start_chain = time.time()

minutes = (time.time() - start_time) / 60
print(minutes, " minutes total")
