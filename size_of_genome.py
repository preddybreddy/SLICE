def genome_size():
    sizes = []
    with open('gene_sizes.txt') as f:
        line = f.readline()#.split()[1]
        while line != "":
            sizes.append(line.split()[1]) 
        #   print(line.split())
            line = f.readline()
    #print(sizes)

    sizes = map(lambda x: int(x), sizes)
    return sum(sizes)