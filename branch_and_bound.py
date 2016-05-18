from greedy import generate_permutations
from genotype_generator import generate_genotype_input
from greedy import greedy_with_regret_solver
from greedy import convert_input
from solutions_checker import check_solution2


def get_solution_pairs(genotype):
    num_replacements = genotype.count('2')
    if num_replacements == 0:
        return [(genotype, genotype)]
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
    if len(haplotype_resolution[h1]) + len(haplotype_resolution[h2]) != 3 \
            or len(haplotype_resolution[h3]) + len(haplotype_resolution[h4]) != 3:
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


def branch_and_bound(input_genotypes, greedy_solution_count):
    # generate resolution coverage.
    # maps haplotypes to genotypes it covers
    unique_genotypes = list(set(input_genotypes))
    # dict with key -> genotype; value -> list of tuple resolutions
    genotype_resolutions = {}
    haplotype_resolution = {}

    genotypes_with_two_resolutions = []
    for genotype in unique_genotypes:
        genotype_resolutions[genotype] = get_solution_pairs(genotype)
        for (h1, h2) in genotype_resolutions[genotype]:
            if h1 in haplotype_resolution:
                haplotype_resolution[h1].add(genotype)
            else:
                haplotype_resolution[h1] = {genotype}
            if h2 in haplotype_resolution:
                haplotype_resolution[h2].add(genotype)
            else:
                haplotype_resolution[h2] = {genotype}

    # now we need to prune solutions
    # Prune Case 1: Remove all extra resolutions with only coverage 2
    for genotype, solutions in genotype_resolutions.iteritems():
        # solutions is a list containing tuples of haplotype resolutions
        found_coverage2 = False
        for (h1, h2) in solutions:
            if len(haplotype_resolution[h1]) + len(haplotype_resolution[h2]) == 2:
                if not found_coverage2:
                    found_coverage2 = True
                else:
                    solutions.remove((h1, h2))
        if len(solutions) == 2:
            genotypes_with_two_resolutions.append(genotype)

    # Prune Case 2: for all genotypes with 2 resolutions, such that:
    # Mi has (h1, h2) and (h4, h5) and Mj has (h2, h3) and (h5, h6)
    # if h1, h3, h4, h6 have coverage 1 and h2 and h4 have coverage 2:
    # keep only (h1, h2) and (h2, h3) [delete (h4, h5) and (h5, h6)
    for genotype1 in genotypes_with_two_resolutions:
        (h1, h2, h3, h4, valid) = \
            determine_h1_h2_h3_h4(genotype1, genotype_resolutions,
                                  haplotype_resolution)
        if not valid:
            continue
        for genotype2 in genotypes_with_two_resolutions:
            if genotype1 == genotype2:
                continue
            (h5, h6, h7, h8, valid) = \
                determine_h1_h2_h3_h4(genotype2, genotype_resolutions,
                                      haplotype_resolution)
            if not valid:
                continue
            if (h2 == h6 and h4 == h8) or (h2 == h8 and h4 == h6):
                # delete h3,h4 and h7,h8 (aka the latter 2 pairs of each
                # genotype)
                genotype_resolutions[genotype1].pop()
                genotype_resolutions[genotype2].pop()

    for solutions in genotype_resolutions.itervalues():
        if not solutions:
            print 'empty solutions list!'

    unique_genotypes.sort(key=lambda x: len(genotype_resolutions[x]))
    # perform depth first search of solutions space
    # opt_solution is a map of genotype -> resolution tuple
    opt_solution = {}
    solution_haplotypes = {}
    unique_haplotype_count = 0
    stack = [(unique_genotypes[0],
                  list(genotype_resolutions[unique_genotypes[0]]))]
    while stack:
        (cur_genotype, cur_genotype_resolutions) = stack[-1]
        # pick first resolution off stack.
        cur_resolution = cur_genotype_resolutions[-1]
        # append haplotypes / update haplotype count
        for haplotype in cur_resolution:
            if not haplotype in solution_haplotypes:
                solution_haplotypes[haplotype] = 0
            if solution_haplotypes[haplotype] == 0:
                unique_haplotype_count += 1
            solution_haplotypes[haplotype] += 1
        opt_solution[cur_genotype] = cur_resolution
        # if picking solution creates partial cost <= greedy:
        if unique_haplotype_count <= greedy_solution_count:
            if len(stack) == len(unique_genotypes):
                break
                # save as current best solution
            else:
                # add next level solutions
                next_genotype = unique_genotypes[len(stack)]
                stack.append((next_genotype,
                              list(genotype_resolutions[next_genotype])))
        else:
            # undo action just committed
            for haplotype in cur_resolution:
                solution_haplotypes[haplotype] -= 1
                if solution_haplotypes[haplotype] == 0:
                    unique_haplotype_count -= 1
            # delete action just undone
            cur_genotype_resolutions.pop()
            # if still solutions left for current depth:
            if cur_genotype_resolutions:
                continue
            else:
                # undo previous layer
                stack.pop()
                (cur_genotype, cur_genotype_resolutions) = stack[-1]
                cur_resolution = cur_genotype_resolutions[-1]
                for haplotype in cur_resolution:
                    solution_haplotypes[haplotype] -= 1
                    if solution_haplotypes[haplotype] == 0:
                        unique_haplotype_count -= 1
                cur_genotype_resolutions.pop()
                while not cur_genotype_resolutions:
                    stack.pop()
                    (cur_genotype, cur_genotype_resolutions) = stack[-1]
                    cur_resolution = cur_genotype_resolutions[-1]
                    for haplotype in cur_resolution:
                        solution_haplotypes[haplotype] -= 1
                        if solution_haplotypes[haplotype] == 0:
                            unique_haplotype_count -= 1
                    cur_genotype_resolutions.pop()
                #del opt_solution[cur_genotype]

    final_solution = []
    for genotype in input_genotypes:
        pair = opt_solution[genotype]
        final_solution.append(pair[0])
        final_solution.append(pair[1])
    return unique_haplotype_count, solution_haplotypes, final_solution


def main():

    #mediumHaplotype = [0, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0]
    #longHaplotype = [0, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    #shortHaplotype = [0, 2, 0, 1, 1]
    #(haplotypes, solution) = greedy_with_regret_solver([longHaplotype])
    #print solution
    avg_diff = 0
    runs = 20
    test1 = ['10101',
             '21010',
             '21211',
             '22112',
             '02120']
    (greedy_haplotypes, solution) = greedy_with_regret_solver(test1)
    print greedy_haplotypes
    print solution
    (opt_haplotypes_count, opt_haplotypes, opt_solution) = \
        branch_and_bound(test1, 7)
    for x in range(runs):
        input_genotypes = generate_genotype_input(N=10, M=10, L=50, print_output=False)
        converted_genotypes = convert_input(input_genotypes)
        (greedy_haplotypes, solution) = greedy_with_regret_solver(converted_genotypes)
        print '{} {}'.format(len(greedy_haplotypes), check_solution2(input_genotypes, solution))

        (opt_haplotypes_count, opt_haplotypes, opt_solution) = \
            branch_and_bound(converted_genotypes, len(greedy_haplotypes))
        print '{} {}'.format(opt_haplotypes_count, check_solution2(input_genotypes, opt_solution))
        print '-----'

    print 'avg diff: {}'.format(avg_diff/float(runs))


if __name__ == "__main__":
    main()





