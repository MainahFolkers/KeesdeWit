# KeesdeWit

Proteins are strings of amino acids. The spatial conformation of these chains is of utmost importance for a protein's functioning in biological processes.
It is known that hydrophobic amino acids (H) and cysteine amino acids (C) prefer proximity to one another. Such is not the case for polar amino acids (P). Two adjacent H's, two C's, or a C and an H: can all form a bond that increases the stabilty of the protein. 
The goal of this project is to fold the given proteins with maximal stability.

## Getting Started
Download Protein_class.py (for 2D folding) along with main.py, random_sampler.py, hill_climber.py, sim_anneal.py and depth_first_search.py and run main.py.
If you wish to fold in 3D download Protein_class_3D.py and change every reference to Protein_class to Protein_class_3D.

### Prerequisites
Python 3


## Running the tests
By downloading and running AOM.py a test to see what the effect of varying the amount of mutations per iteration of the hill climber algoritm is can be performed.

## Built With

* Sublime2 text editor https://www.sublimetext.com/
* Atom text editor https://atom.io/

## Versioning

Python 3.6.1

libraries:
* copy (deepcopy)
* matplotlib.pyplot
* time
* random (choice, uniform)

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
