import random


# generates an NxM matrix with L 2's and randomly dispersed 0's and 1's
def generate_genotype_input(N = 100, M=15, L=200, print_output=False):
    # we have an N x M matrix.
    Matrix = [[0 for x in range(M)] for y in range(N)]

    # randomly pick L cells to be a 2
    for x in range(L):
        loc = random.randint(0, N*M - 1)
        # print loc
        # print loc / N
        # print loc % M
        added_2 = False
        while not added_2:
            if Matrix[loc / M][loc % M] != 2:
                Matrix[loc / M][loc % M] = 2
                added_2 = True
            else:
                loc = (loc + 1) % (N*M)

    for r in range(N):
        for c in range(M):
            if Matrix[r][c] != 2:
                Matrix[r][c] = random.randint(0, 1)
    if print_output:
        for row in Matrix:
            print row

    return Matrix


def main():
    generate_genotype_input()

if __name__ == "__main__":
    main()

