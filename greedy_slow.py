import itertools
import math
from genotype_generator import generate_genotype_input
from solutions_checker import check_solution

def generate_permutations(length):
    curList = [[]]
    tempList = []
    for i in range(length):
        max = int(math.pow(2, i))
        for x in range(max):
            tempList.append(curList[x] + [0])
            tempList.append(curList[x] + [1])

        oldCurList = curList
        curList = tempList
        tempList = oldCurList
        del tempList[:]

    #print curList
    return curList


def haplotype_solves_genotype(haplotype, genotype):
    if len(haplotype) != len(genotype):
        return False
    for i in range(len(haplotype)):
        if genotype[i] == 2:
            continue
        elif genotype[i] == 0 and haplotype[i] != 0:
            return False
        elif genotype[i] == 1 and haplotype[i] != 1:
            return False
    return True

def get_genotype_hash(genotype):
    string_hash = ""
    for i in range(len(genotype)):
        string_hash += str(genotype[i])
    return string_hash


# finds the other haplotype h2 for the input haplotype h1 such that h1 and h2 solve genotype
def find_haplotype_complement(haplotype, genotype):
    haplotype_complement = []
    if haplotype == genotype:
        return genotype

    for i in range(len(genotype)):
        if genotype[i] == 2:
            if haplotype[i] == 0:
                haplotype_complement.append(1)
            else:
                haplotype_complement.append(0)
        else:
            haplotype_complement.append(genotype[i])
    return haplotype_complement


# Given NxM matrix, where N = # of individuals, M = # of SNPs
# Find haplotype set.
def greedy_haplotype_solver(genotypes):
    M = len(genotypes[0])

    # generate all 2^M possible haplotypes
    allHaplotypes = generate_permutations(M)
    haplotype_list = set()

    genotype_dict = {}
    remaining_genotypes = list(genotypes)
    while len(remaining_genotypes) > 0:
        haplotype_match_list = []
        max_haplotype = []
        for haplotype in allHaplotypes:
            temp_list = []
            for genotype in remaining_genotypes:
                if haplotype_solves_genotype(haplotype, genotype):
                    temp_list.append(genotype)
            if len(temp_list) > len(haplotype_match_list):
                max_haplotype = haplotype
                haplotype_match_list = temp_list

        # modify all genotypes in list based on haplotype
        for genotype in haplotype_match_list:
            remaining_genotypes.remove(genotype)
            haplotype_complement = find_haplotype_complement(max_haplotype, genotype)
            genotype_dict[get_genotype_hash(genotype)] = (max_haplotype, haplotype_complement)
            haplotype_list.add(get_genotype_hash(haplotype_complement))
        # add haplotype to haplotype list.
        haplotype_list.add(get_genotype_hash(max_haplotype))
        allHaplotypes.remove(max_haplotype)

    solution = []
    for genotype in genotypes:
        pair = genotype_dict[get_genotype_hash(genotype)]
        solution.append(pair[0])
        solution.append(pair[1])
    return (haplotype_list, solution)


def main():
    mediumHaplotype = [0, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0]
    longHaplotype = [0, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    shortHaplotype = [0, 2, 0, 1, 1]
    solution = greedy_haplotype_solver([mediumHaplotype])
    # print solution

    input = generate_genotype_input(N=10, M=5, L=5)
    (haplotypes, solution) = greedy_haplotype_solver(input)
    print solution
    print len(haplotypes)
    print check_solution(input, solution)


if __name__ == "__main__":
    main()
