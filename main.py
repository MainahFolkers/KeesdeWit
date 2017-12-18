from Protein_class_2D import *
from random_sampler import *
from hill_climber import *
from sim_anneal import *
# from depth_first_search import *

chain = input("Chain: (HHPHHHPH)")

protein = Protein(chain)

# ask user for amount of iterations
ITER = int(input("Amount of interations: "))

# random sampler algorithm
protein = rand_samp(protein, ITER)

vis = input("Visualize? (y/n) ")

if vis == 'y':
    # visualize best random sampled protein
    protein.visualize("NA", "NA", "NA")

# ask user for amount of iterations
AOM = int(input("Amount of mutations: "))

# hill climber algorithm
protein = hill_climb(protein, ITER, AOM)

vis = input("Visualize? (y/n) ")

if vis == 'y':
    # visualize best random sampled protein
    protein.visualize(AOM, "NA", "NA")
    
# ask user input for cooling schedule: linear / hyperbolic
cool = input("Cooling schedule: linear / hyperbolic: ")

# simmulated annealing algorithm
protein = sim_anneal(protein, ITER, AOM, cool)
vis = input("Visualize? (y/n) ")

if vis == 'y':
    # visualize best random sampled protein
    protein.visualize(AOM, cool, "NA")

# Depth first search algorithm 
protein = depth_first_search(protein)

vis = input("Visualize? (y/n) ")

if vis == 'y':
    # visualize best random sampled protein
    protein.visualize("NA", "NA", "NA")

