# Classes

This project is about determining the optimal conformation of proteins (different in chain length and sequence). The proteins can fold with 90 degree angles along x / y, or x / y / z axes, for 2D or 3D, respectively.
The classes contains a class with characteristics of amino acids in the chain.
The classes contain several functions for folding the protein: rand_fold() for random folding, starting from zero. Here, the first bond is always to the right to avoid symnmetry, and the second bond can only go up or to the right. rand_fold() always yields a valid fold. 
mut_dir() mutates the direction of a bond. This function is called in the hillclimber and simulated annealing algorithms, after which mut_fold() is called. This function folds the mutated protein. 
Both rand_fold() and mut_fold() call check_val(), whick checks if a direction is valid so that the protein does not overlap itself.
After folding, the protein is placed on a grid. bend() converts bond directions to amino-acid coordinates on the grid. make_grid() then places the protein on a dynamic, protein-appropriate minimal grid.
To calculate the score of a protein, calc_score() walks over the grid and checks proximities of H-H, H-C, or C-C. 
In case of 2D, the protein can then be visualized on its grid with visualize(), with blue and red dots for P and H amino-acids, respectively.

## Getting Started

The classes can be called in all other files.

### Prerequisites

* Python 3.x. Installing: https://www.python.org/downloads/
* matplotlib.pyplot. Installing:
  * python -mpip install -U pip
  * python -mpip install -U matplotlib

## Running the tests

The classes are called in the files AOM.py, COOL.py, and ALGOS.py in the directory KeesdeWit/Experimentation. Run those files to call the algorithms and to visualize the results!

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
