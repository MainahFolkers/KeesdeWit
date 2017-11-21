def hill_climb(protein):
    for i in range(50):
        fold(protein)
        protein.calc_score()
        cur_score = protein.score

        fold(protein)
        protein.calc_score()
        new_score = protein.score

        if cur_score > new_score:
            cur_score = new_score

# x = choice(for i in range(protein.n))
# directions = [r, d, u, l, u, d, r]
# direction[x] = andere richting

protein = Protein("HHPHHHPH")
hill_climb(protein)
