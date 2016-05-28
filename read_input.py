

def read_input(input_file_name = 'test.txt'):
    with open(input_file_name) as f:
        content = f.readlines()
    content = map(lambda x: x.replace('\n', ''), content)
    return content


def convert_input(genotypes):
    content = map(lambda x: x.replace('1', '3'), genotypes)
    content = map(lambda x: x.replace('2', '1'), content)
    return map(lambda x: x.replace('3', '2'), content)


def main():
    content = read_input('data/very_easy_test_reads.txt')
    for line in content:
        print line

    print '-----'
    content = convert_input(content)
    for line in content:
        print line


if __name__ == "__main__":
    main()