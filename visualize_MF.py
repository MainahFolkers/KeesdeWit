import matplotlib.pyplot as plt
import numpy as np

def visualize(protein):
    plt.figure()

    # iterate over amino acids
    for i in range(protein.n):
        # determine color for point on plot (later expand with C)
        if protein.chain[i] == 'H':
            col = 'red'
        else:
            col = 'blue'

        plt.scatter(protein.xs[i], protein.ys[i], s=120, zorder=2, color=col)
        plt.annotate(i, xy=(protein.xs[i], protein.ys[i]), xytext=(protein.xs[i] + 0.05, protein.ys[i] + 0.05), fontsize=20)

    # plot black line behind / between points
    plt.plot(protein.xs, protein.ys, lw=3, zorder=1, color='black')

    # show grid on plot
    plt.grid(b=True)

    # set tick spacing to 1
    plt.xticks(range(min(protein.xs), protein.xran + 1))
    plt.yticks(range(min(protein.ys), protein.yran  + 1))

    plt.title("Score = " + str(protein.score), fontsize=40)

    # improvement: plot dotted lines for H-bonds

    # show plot
    plt.show()
