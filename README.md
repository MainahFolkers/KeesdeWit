# KeesdeWit

Proteins are strings of amino acids. The spatial conformation of these chains is of utmost importance for a protein's functioning in biological processes.
It is known that hydrophobic amino acids (H) and cysteine amino acids (C) prefer proximity to one another. Such is not the case for polar amino acids (P). Two adjacent H's, two C's, or a C and an H: can all form a bond that increases the stabilty of the protein. 
The goal of this project is to fold the given proteins with maximal stability.

## Getting Started
There are three subrepositories. 'Algorithms' contains random_sampler.py, hill_climber.py, sim_anneal.py, depth_first_search.py; the names of which reveal the algorithm. 'Classes' contains Protein_class_2D.py and Protein_class_3D.py, for protein folding in 2D or 3D, respectively. 'Experimentation' contains AOM.py (to determine the optimal amount of mutations per iteration, relevant for the hillclimber and simulated annealing algorithms), COOL.py (to determine the optimal cooling schedule for the simualted annealing), ALGOS.py (to determine which algorithm yields the best score given the samen amount of iterations and the same AOM for several runs).
Running main.py gives examples of output that you get in the experimentation files. 

### Prerequisites
* Python 3.x. Installing: https://www.python.org/downloads/
* matplotlib.pyplot. Installing:
  * python -mpip install -U pip
  * python -mpip install -U matplotlib

## Running the tests
In main.py, specify the test you want to perform. Give input in the terminal to see the visualizations of folds!

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
