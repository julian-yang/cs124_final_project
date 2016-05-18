from greedy import generate_permutations


def get_solution_pairs(genotype):
    num_replacements = genotype.count('2')
    replacements = generate_permutations(num_replacements)
    solutions = [list(genotype) for i in range(len(replacements))]
    begin_find = 0
    for i in range(num_replacements):
        replace_index = genotype.find('2', begin_find)
        begin_find = replace_index + 1
        for j in range(len(replacements)):
            solutions[j][replace_index] = replacements[j][i]

    solutions[:] = [''.join(s) for s in solutions]
    l1 = solutions[:len(solutions)/2]
    l2 = solutions[len(solutions)/2:]
    l2.reverse()
    fills = zip(l1, l2)
    return fills


def determine_h1_h2_h3_h4(genotype, genotype_solutions, haplotype_resolution):
    (h1, h2) = genotype_solutions[genotype][0]
    (h3, h4) = genotype_solutions[genotype][1]
    if haplotype_resolution[h1] + haplotype_resolution[h2] != 3 \
            or haplotype_resolution[h3] + haplotype_resolution[h4] != 3:
        return h1, h2, h3, h4, False
    if haplotype_resolution[h1] == 2:
        temp = h1
        h1 = h2
        h2 = temp
    if haplotype_resolution[h3] == 2:
        temp = h3
        h3 = h4
        h4 = temp
    return h1, h2, h3, h4, True


def branch_and_bound(input, greedy_solution_count):
    # generate resolution coverage.
    # maps haplotypes to genotypes it covers
    unique_genotypes = set(input)
    # dict with key -> genotype; value -> list of tuple resolutions
    genotype_solutions = {}
    haplotype_resolution = {}

    genotypes_with_two_resolutions = []
    for genotype in unique_genotypes:
        genotype_solutions[genotype] = get_solution_pairs(genotype)
        for (h1, h2) in genotype_solutions[genotype]:
            if h1 in haplotype_resolution:
                haplotype_resolution[h1].add(genotype)
            else:
                haplotype_resolution[h1] = set([genotype])
            if h2 in haplotype_resolution:
                haplotype_resolution[h2].add(genotype)
            else:
                haplotype_resolution[h2] = set([genotype])
        if len(genotype_solutions[genotype]) == 2:
            genotypes_with_two_resolutions.append(genotype)

    # now we need to prune solutions
    # Prune Case 1: Remove all extra resolutions with only coverage 2
    for genotype, solutions in genotype_solutions.iteritems():
        # solutions is a list containing tuples of haplotype resolutions
        found_coverage2 = False
        for (h1, h2) in solutions:
            if haplotype_resolution[h1] + haplotype_resolution[h2] == 2:
                if not found_coverage2:
                    found_coverage2 = True
                else:
                    solutions.remove((h1, h2))

    # Prune Case 2: for all genotypes with 2 resolutions, such that:
    # Mi has (h1, h2) and (h4, h5) and Mj has (h2, h3) and (h5, h6)
    # if h1, h3, h4, h6 have coverage 1 and h2 and h4 have coverage 2:
    # keep only (h1, h2) and (h2, h3) [delete (h4, h5) and (h5, h6)
    for genotype1 in genotypes_with_two_resolutions:
        (h1, h2, h3, h4, valid) = \
            determine_h1_h2_h3_h4(genotype1, genotype_solutions,
                                  haplotype_resolution)
        if not valid:
            continue
        for genotype2 in genotypes_with_two_resolutions:
            if genotype1 == genotype2:
                continue
            (h5, h6, h7, h8, valid) = \
                determine_h1_h2_h3_h4(genotype2, genotype_solutions,
                                      haplotype_resolution)
            if not valid:
                continue
            if (h2 == h6 and h4 == h8) or (h2 == h8 and h4 == h6):
                # delete h3,h4 and h7,h8 (aka the latter 2 pairs of each
                # genotype)
                del genotype_solutions[genotype1][-1]
                del genotype_solutions[genotype2][-1]







