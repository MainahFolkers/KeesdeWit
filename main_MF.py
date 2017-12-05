from random_sampler_MF import *
from hill_climber_MF import *

# amount of iterations
ITER = 10000

# short protein declaration
sprotein = Protein("HHPHHHPH")
# random sampler algorithm
rand_samp(sprotein, ITER)
# hill climber algorithm
hill_climb(sprotein, ITER)

# long protein declaration
lprotein = Protein("HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH")
# random sampler algorithm
rand_samp(lprotein, ITER)
# hill climber algorithm
# improvement: make > 1 mutations per iteration (14 is best for 50 long protein)
hill_climb(lprotein, ITER)

# trie met array van directions, bij het maken van de trie gelijk de validiteit van de vouwing testen
# depth first, weinig geheugen, recursie is een optie: maxdepth is aantal bindingen

# visualisatie met tkinter -> erg ingewikkeld

# experimentatie: lijn grafiek van gemiddelde per algoritme
