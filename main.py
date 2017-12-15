from random_sampler import *
from hill_climber import *
from sim_anneal import *
from visualize import *
import matplotlib.pyplot as plt

chain = input("Chain: (HHPHHHPH)")

protein = Protein(chain)

# ask user for amount of iterations
ITER = int(input("Amount of interations: "))

# random sampler algorithm
protein = rand_samp(protein, ITER)

vis = input("Visualize? (y/n) ")

if vis == 'y':
    # visualize best random sampled protein
    visualize(protein)

# ask user for amount of iterations
AOM = int(input("Amount of mutations: "))

# hill climber algorithm
protein = hill_climb(protein, ITER, AOM)

vis = input("Visualize? (y/n) ")

if vis == 'y':
    # visualize best hill climbed protein
    visualize(protein)

# simmulated annealing algorithm
protein = sim_anneal(protein, ITER, AOM)

vis = input("Visualize? (y/n) ")

if vis == 'y':
    # visualize best hill climbed protein
    visualize(protein)

# Depth first search algorithm 
protein = depth_path(path, 0, maxdepth, 0)

vis = input("Visualize? (y/n) ")

if vis == 'y':
    # visualize best hill climbed protein
    visualize(protein)