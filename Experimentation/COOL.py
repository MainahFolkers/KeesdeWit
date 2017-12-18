from Protein_class import *
from sim_anneal import *
from copy import deepcopy
import matplotlib.pyplot as plt
from time import time

# save start time
start_time = time()

# tested protein chains
chains = [
    "HHPHHHPH",
    "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH",
    "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    ]

# amount of iterations
ITER = 2000
# amount of runs
RUNS = 5
# amount of mutations
AOM = 7

# iterate over chains
for chain in chains:
    # set impossible best score (minimum score = 0)
    best_score = 1

    # data series for y-axis in legend
    ys = []

    print("Protein chain =", chain)
    protein = Protein(chain)

    # iterate over cooling schedules
    for cool in ["linear", "hyperbolic"]:
        # iterate over runs
        for run in range(1, RUNS + 1):
            print("Run =", run)

            protein = sim_anneal(protein, ITER, AOM, cool)

            # if score improved
            if protein.score < best_score:
                # visualize protein
                protein.visualize(AOM, cool, run)
                best_score = deepcopy(protein.score)
                print("Current best score = ", best_score)

    # open input file
    with open("COOL_" + protein.chain + "_AOM=" + str(AOM) +".txt", 'r') as ifile:
        data = ifile.readlines()
        for line in data:
            # remove comma
            scores = line.split(',')
            # remove \n
            scores = scores[:-1]
            # covert score strings to integers
            scores[1:] = [ int(score) for score in scores[1:] ]
            ys.append(scores)
    # close input file
    ifile.close()

    # open figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # label axes
    plt.xlabel("Iterations")
    plt.ylabel("Best score")

    ys_mean = []

    # iterate over runs
    for i in range(len(ys) - RUNS + 1):
        # if new sample starts
        if i % RUNS is 0:
            y_mean = []
            # append legend label
            y_mean.append(ys[i][0])
            # iterate over scores of 1 run
            for j in range(1, len(ys[i])):
                score_sum = 0
                # iterate over ys in same sample
                for k in range(i, i + RUNS):
                    score_sum += ys[k][j]
                # calculate mean score per iteration
                score_mean = score_sum / RUNS
                y_mean.append(score_mean)
            ys_mean.append(y_mean)

    i = len(ys_mean)
    for y_mean in ys_mean:
        # x is amount of iterations
        x = range(len(y_mean) - 2)
        # plot line, first element of y is legend label
        ax.plot(x, y_mean[2:], label=y_mean[0], zorder=i)
        i -= 1

    # make plot legend
    plt.legend(title="Cooling schedule", bbox_to_anchor=(1.04,1), loc="upper left", borderaxespad=0)

    # save figure
    plt.savefig("COOL_" + protein.chain + ".jpg", bbox_inches = 'tight')

    # close figure
    plt.clf()

# amount of minutes program was running
minutes = (time() - start_time) / 60
print(minutes, " minutes total")
