
def valid_entry(x):
    if x == 0 or x == 1 or x == 2:
        return True
    else:
        return False


def valid_entry2(x):
    x = int(x)
    if x == 0 or x == 1 or x == 2:
        return True
    else:
        return False


def check_solution(input, solution):
    N = len(input)
    for i in range(N):
        h1 = solution[i*2]
        h2 = solution[i*2+1]
        genotype = input[i]
        for j in range(len(genotype)):
            if (not valid_entry(genotype[j])) or (not valid_entry(h1[j])) or (not valid_entry(h2[j])):
                return False

            if genotype[j] == 0 and (h1[j] != 0 or h2[j] != 0):
                return False
            if genotype[j] == 1 and (h1[j] != 1 or h2[j] != 1):
                return False
            if genotype[j] == 2 and (h1[j] + h2[j] != 1):
                return False
    return True


def check_solution2(input, solution):
    N = len(input)
    for i in range(N):
        h1 = solution[i*2]
        h2 = solution[i*2+1]
        genotype = input[i]
        for j in range(len(genotype)):
            if (not valid_entry2(genotype[j])) or (not valid_entry2(h1[j])) or (not valid_entry2(h2[j])):
                return False

            if genotype[j] == '0' and (h1[j] != '0' or h2[j] != '0'):
                return False
            if genotype[j] == '1' and (h1[j] != '1' or h2[j] != '1'):
                return False
            if genotype[j] == '2' and (h1[j] + h2[j] != '1'):
                return False
    return True


def main():
    assert(check_solution([[0, 2, 0, 1]], [[0, 1, 0, 1], [0, 0, 0, 1]]) == True)
    assert(check_solution([[0, 3, 0, 1]], [[0, 1, 0, 1], [0, 0, 0, 1]]) == False)
    assert(check_solution([[0, 2, 0, 1]], [[0, 3, 0, 1], [0, 0, 0, 1]]) == False)
    assert(check_solution([[0, 2, 0, 1]], [[0, 1, 0, 1], [0, 1, 0, 1]]) == False)
    assert(check_solution([[0, 2, 0, 1]], [[1, 1, 0, 1], [0, 1, 0, 1]]) == False)


if __name__ == "__main__":
    main()


