# Experimentation

To determine the optimal conformation of proteins (different in chain length and sequence), we have implemented several algorithms (see also KeesdeWit/Algorithms).
Several experiments can be performed:
* Firstly, the amount of mutations (AOM) can be determined for using the hillclimber and simulated annealing algorithms. This can be determined by running the file AOM.py. Choose the desired chain, number of runs and the number of iterations, and you're good to go! The best result is visualized and the outputted graph depicts the course of the determination. 
* Secondly, for the simulated annealing the optimal cooling schedule can be determined. Run COOL.py to see whether a linear or hyperbolic cooling schedule is optimal. 
* Lastly, the performance of all available algorithms can be determined with ALGOS.py. The outputted graph depicts the mean score achieved over a number of runs with a number of iterations. 

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
