# Algorithms

To determine the optimal conformation of proteins (different in chain length and sequence), we have implemented several algorithms. 
* The random sampler folds the protein randomly, using the function protein.rand_fold(). The direction of the first bond is always to the right, to avoid symmetry. The second bond can only ever go up or down, also to prune some solutions. 
* The hillclimber folds a protein randomly on its first iteration. On following iterations, the direction of a number of bonds is mutated. Only if this yields a valid protein, the try is counted as an iteration. 
* The simulated annealing is an extension of the hillclimber. Here too, a protein is randomly folded and then mutated. Only, where a hillclimber can get stuck in a local minimum score, with simulated annealing two cooling schedules are provided in the algorithm, with a linear and a hyperbolic probability of accepting a deteriorated score. 
* The depth-first algorithm constructively folds a protein, thereby calculating every possibility. However, with chains > 13, the computer runs out of memory. Therefore, greedy choices are made to run the algorithm in reasonable time. Thereby the guarantee on the best score is given up.


## Getting Started

The algorithms are called in the files AOM.py, COOL.py, and ALGOS.py in the directory KeesdeWit/Experimentation. Run those files to call the algorithms and to visualize the results!

### Prerequisites

* Python 3.x. Installing: https://www.python.org/downloads/
* matplotlib.pyplot. Installing:
  * python -mpip install -U pip
  * python -mpip install -U matplotlib

## Running the tests

The algorithms are called in the files AOM.py, COOL.py, and ALGOS.py in the directory KeesdeWit/Experimentation. Run those files to call the algorithms and to visualize the results!

## Built With

* Sublime2 text editor https://www.sublimetext.com/
* Atom text editor https://atom.io/

## Versioning

Python 3.6.1

libraries:
* copy (deepcopy) https://docs.python.org/3/library/copy.html
* matplotlib.pyplot https://matplotlib.org/api/_as_gen/matplotlib.pyplot.plot.html
* time https://docs.python.org/3/library/time.html
* random (choice, uniform)  https://docs.python.org/3/library/random.html
* collections (deque) https://docs.python.org/3/library/collections.html?highlight=collections#module-collections
* math https://docs.python.org/3/library/math.html?highlight=math#module-math

## Authors

* Mainah Folkers, student number 10535845
* Milan Scholten, student number 10551891
* Rianne M. Schoon, student number 10742794

## Acknowledgments

* Bas Terwijn for the tech assist
* Daan van den Berg for feedback on our approach and the project presentation
* Bart for feedback on our practice presentations
* Zhang, J., Kou, S. C., & Liu, J. S. (2007). Biopolymer structure simulation and optimization via fragment regrowth Monte Carlo. The Journal of chemical physics, 126(22), 06B605.
* Custódio, F. L., Barbosa, H. J., & Dardenne, L. E. (2004). Investigation of the three-dimensional lattice HP protein folding model using a genetic algorithm. Genetics and Molecular Biology, 27(4), 611-615.
* Dill, K. A. et al. (1995). Principles of protein folding - A perspective from simple exact models. Protein Science, 4, 561-602.
