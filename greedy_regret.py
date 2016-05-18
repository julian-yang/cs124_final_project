from genotype_generator import generate_genotype_input
from solutions_checker import check_solution2
from solutions_checker import check_solution
from greedy import haplotype_solves_genotype
from greedy import get_genotype_hash
from greedy import find_haplotype_complement
from greedy import greedy_haplotype_solver

import math

def generate_permutations(length):
    curList = ['']
    tempList = []
    for i in range(length):
        max = int(math.pow(2, i))
        for x in range(max):
            tempList.append(curList[x] + '0')
            tempList.append(curList[x] + '1')
        oldCurList = curList
        curList = tempList
        tempList = oldCurList
        del tempList[:]
    #print curList
    return curList


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


def haplotype_solves_genotype(haplotype, genotype):
    if len(haplotype) != len(genotype):
        return False
    for i in range(len(haplotype)):
        if genotype[i] == '2':
            continue
        elif genotype[i] == '0' and haplotype[i] != '0':
            return False
        elif genotype[i] == '1' and haplotype[i] != '1':
            return False
    return True


# finds the other haplotype h2 for the input haplotype h1 such that h1 and h2 solve genotype
def find_haplotype_complement(haplotype, genotype):
    haplotype_complement = []
    if haplotype == genotype:
        return genotype

    for i in range(len(genotype)):
        if genotype[i] == '2':
            if haplotype[i] == '0':
                haplotype_complement +='1'
            else:
                haplotype_complement += '0'
        else:
            haplotype_complement += genotype[i]
    return haplotype_complement


# Given NxM matrix, where N = # of individuals, M = # of SNPs
# Find haplotype set.
def greedy_with_regret_solver(genotypes, allow_regret=True):
    M = len(genotypes[0])

    # generate all 2^M possible haplotypes
    allHaplotypes = generate_permutations(M)
    haplotype_list = set()

    genotype_dict = {}

    # converted_genotypes is in str format as opposed to a list of ints
    converted_genotypes = []
    for genotype in genotypes:
        converted = get_genotype_hash(genotype)
        converted_genotypes.append(converted)

    unique_genotypes = set(converted_genotypes)
    genotypes_covered = set()

    # maps haplotypes to genotypes it covers
    haplotype_possible_coverage = {}
    for haplotype in allHaplotypes:
        temp_list = set()
        for genotype in unique_genotypes:
            if haplotype_solves_genotype(haplotype, genotype):
                temp_list.add(genotype)
        haplotype_possible_coverage[haplotype] = temp_list

    # maps haplotypes to the genotypes it covers in genotype_dict
    haplotype_coverage = {}
    # history of max_unique haplotypes
    max_unique_haplotypes = []
    removed_count = 0
    while len(genotypes_covered) < len(unique_genotypes):
        haplotype_add_max_unique_coverage = set()
        max_haplotype = ""
        # find best haplotype to add
        for haplotype in allHaplotypes:
            temp_list = haplotype_possible_coverage[haplotype].difference(genotypes_covered)

            if len(temp_list) > len(haplotype_add_max_unique_coverage):
                max_haplotype = haplotype
                haplotype_add_max_unique_coverage = temp_list

        # find the best haplotype to remove
            # we need to remember which genotypes are covered by a given haplotype
            # we need to remember which genotypes are covered collectively by the set.
                # for each genotype covered in set, list all pairs that cover it.
        min_haplotype = ""
        haplotype_delete_min_unique_coverage = set()
        recovered_dict = {}
        covered_without_min_haplotype = set()
        num_min_unique = int(2**M)
        for haplotype in max_unique_haplotypes:
            # calculate coverage of set after removing haplotype.
            covered_without_cur_haplotype = genotypes_covered.difference(haplotype_coverage[haplotype])
            # attempt to recover lost genotypes
            for lost_genotype in haplotype_coverage[haplotype]:
                # get pairs that can solve lost_genotype
                for (h1, h2) in get_solution_pairs(lost_genotype):
                    if h1 == haplotype or h2 == haplotype:
                        continue
                    if (h1 in haplotype_list) and (h2 in haplotype_list):
                        covered_without_cur_haplotype.add(lost_genotype)
                        recovered_dict[lost_genotype] = (h1, h2)
                        break

            temp_list = haplotype_coverage[haplotype].difference(covered_without_cur_haplotype)
            if len(temp_list) < num_min_unique:
                min_haplotype = haplotype
                haplotype_delete_min_unique_coverage = temp_list
                covered_without_min_haplotype = covered_without_cur_haplotype
                num_min_unique = len(temp_list)

            # perform set subtraction: genotypes covered by haplotype - genotypes covered by set
            # find set that is most expensive -> max(1/|set|) -> min(|set|)
        # add haplotype to the set
        # modify all genotypes in list based on haplotype

        # if len(add) <= len(remove)
            # add
        # else
            # remove
        # we want this to be <= because delete_min_unique should be the most cost-efficient haplotypes,
        # so they should usually cover more unique ones than the next haplotype we want to add.
        if (not max_unique_haplotypes) or \
            (len(haplotype_add_max_unique_coverage) <= len(haplotype_delete_min_unique_coverage)):
            for genotype in haplotype_add_max_unique_coverage:
                # unique_genotypes.remove(genotype)
                haplotype_complement = get_genotype_hash(find_haplotype_complement(max_haplotype, genotype))
                genotype_dict[genotype] = (max_haplotype, haplotype_complement)
                haplotype_list.add(haplotype_complement)
                genotypes_covered.add(genotype)
                if haplotype_complement in haplotype_coverage:
                    haplotype_coverage[haplotype_complement].add(genotype)
                else:
                    haplotype_coverage[haplotype_complement] = set([genotype])
            # add haplotype to haplotype list.
            haplotype_coverage[max_haplotype] = haplotype_add_max_unique_coverage
            haplotype_list.add(max_haplotype)
            max_unique_haplotypes.append(max_haplotype)
            allHaplotypes.remove(max_haplotype)
        else:
            for genotype in haplotype_delete_min_unique_coverage:
                # update new coverage
                genotypes_covered = covered_without_min_haplotype

                (h1, h2) = genotype_dict[genotype]
                haplotype_complement = h1
                if h1 == min_haplotype:
                   haplotype_complement = h2
                # figure out if haplotype_complement has any other connections
                # if haplotype_c has another connection, it should be in haplotype_coverage
                haplotype_coverage[haplotype_complement].remove(genotype)
                if len(haplotype_coverage[haplotype_complement]) == 0 and haplotype_complement != min_haplotype:
                    del haplotype_coverage[haplotype_complement]
                    haplotype_list.remove(haplotype_complement)
                    #allHaplotypes.append(haplotype_complement)

                # remove the unique genotypes from the genotype dict mapping
                del genotype_dict[genotype]
            haplotype_list.remove(min_haplotype)
            max_unique_haplotypes.remove(min_haplotype)
            allHaplotypes.append(min_haplotype)
            # add recovered genotypes
            genotype_dict.update(recovered_dict)
            removed_count += 1


    solution = []
    for genotype in genotypes:
        pair = genotype_dict[get_genotype_hash(genotype)]
        solution.append(pair[0])
        solution.append(pair[1])
    print 'removed_count is {}'.format(removed_count)
    return (haplotype_list, solution)


def main():

    #mediumHaplotype = [0, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0]
    #longHaplotype = [0, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    #shortHaplotype = [0, 2, 0, 1, 1]
    #(haplotypes, solution) = greedy_with_regret_solver([longHaplotype])
    #print solution
    avg_diff = 0
    runs = 20
    for x in range(runs):
        input = generate_genotype_input(N=200, M=10, L=1500)
        (haplotypes, solution) = greedy_with_regret_solver(input)
        (greedy_haplotypes, greedy_solution) = greedy_haplotype_solver(input)
        print '{} {}'.format(len(haplotypes), check_solution2(input, solution))
        print '{} {}'.format(len(greedy_haplotypes), check_solution(input, greedy_solution))
        avg_diff += len(haplotypes) - len(greedy_haplotypes)

    print 'avg diff: {}'.format(avg_diff/float(runs))


if __name__ == "__main__":
    main()
