from branch_and_bound import determine_h1_h2_h3_h4

def test_determine_h1_h2_h3_h4():
    h1 = 'h1'
    h2 = 'h2'
    h3 = 'h3'
    h4 = 'h4'
    h5 = 'h5'
    h6 = 'h6'
    haplotype_resolution = {h1: 1, h2: 2, h3: 1, h4: 2, h5: 0, h6: 4}
    genotype_solutions = {'1': [(h1, h2), (h3, h4)],
                          '2': [(h2, h1), (h4, h3)],
                          '3': [(h1, h5), (h3, h4)],
                          '4': [(h1, h2), (h1, h5)]}
    (a1, a2, a3, a4, valid) = \
        determine_h1_h2_h3_h4('1', genotype_solutions, haplotype_resolution)
    assert(valid and h1 == a1 and h2 == a2 and h3 == a3 and h4 == a4)
    (a1, a2, a3, a4, valid) = \
        determine_h1_h2_h3_h4('2', genotype_solutions, haplotype_resolution)
    assert(valid and h1 == a1 and h2 == a2 and h3 == a3 and h4 == a4)
    (a1, a2, a3, a4, valid) = \
        determine_h1_h2_h3_h4('3', genotype_solutions, haplotype_resolution)
    assert(not valid)
    (a1, a2, a3, a4, valid) = \
        determine_h1_h2_h3_h4('4', genotype_solutions, haplotype_resolution)
    assert(not valid)
    print 'test_determine_h1_h2_h3_h4 passed!'


def main():
    test_determine_h1_h2_h3_h4()


if __name__ == "__main__":
    main()
