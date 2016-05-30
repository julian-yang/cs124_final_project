from read_input import read_input, convert_input
from branch_and_bound import branch_and_bound
from greedy import greedy_solver

input_files = ['data/very_easy_test_reads.txt', 'data/very_easy_training_reads.txt']

def run_solution():
    for file_name in input_files:
        genotypes = convert_input(read_input(file_name))
        (greedy_haplotypes, greedy_solution) = greedy_solver(genotypes)
        (opt_haplotypes_count, opt_solution, total_branches, branches_explored,
         case_1_prunes, case_2_prunes) = \
            branch_and_bound(genotypes, len(greedy_haplotypes),
                             greedy_solution)
        print 'solution for \'{}\':'.format(file_name)
        for line in opt_solution:
            print line
        print ''


def main():
    run_solution()


if __name__ == "__main__":
    main()