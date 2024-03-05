import numpy as np
import argparse
from itertools import combinations
from timeit import default_timer
import matplotlib.pyplot as plt


def main():

    parser = argparse.ArgumentParser(description='Runs the brute force method and the greedy algorithm for minimum set cover of whole genomes and RNA reads.')

    parser.add_argument('--k', type = int, default = 100, help='Number of times to test a matrix of size m by n. Default is 100.')
    parser.add_argument('--test-path', type = str, default = None, help='Path to the test file.')
    parser.add_argument('--save', type = int, default = 0, help = '1 if you want to save the plots, 0 otherwise. Default is 0.')
    parser.add_argument('--test-size', type= int, default = 0, help = 'How large to raise m and n in parallel to compare brute force and greedy algorithm. Default is 0.')
    # parser.add_argument('--m', type=int, default=1, help='number of whole genomes')
    # parser.add_argument('--n', type=int, default=1, help='number of reads')
    parser.add_argument('--seed', type=int, default = 2023, help='Seed for randomizing matrix g. Default is 2023.')

    seed = parser.parse_args().seed
    np.random.seed(seed)

    # m = parser.parse_args().m
    # n = parser.parse_args().n
    k = parser.parse_args().k
    test_path = parser.parse_args().test_path
    
    if test_path is not None:
        with open(test_path) as f:
            n = int(f.readline().strip().split()[0])
            m = int(f.readline().strip().split()[0])
            g = []
            for line in f:
                g.append([int(i) for i in line.strip().split()])
            g = np.array(g)
            if g.shape != (m, n):
                raise Exception(f'Matrix dimensions and provided dimensions do not match.\n m = {m}, n = {n}, ' \
                                f'matrix dimensions = {g.shape}. Matrix dimensions should be ({m}, {n})')
        print('Greedy Solution: ', Greedy(m, n, g))
        print('Brute Force Solution: ', BruteForce(m, n, g))

    if parser.parse_args().test_size > 0:
        CompareGreedyAndBruteForce(k, parser.parse_args().test_size, bool(parser.parse_args().save))
    
    if parser.parse_args().test_size == 0 and parser.parse_args().test_path == None:
        print('Run "py main.py -h" to see options.')

def CompareGreedyBruteForceWrapper(m: int, n: int, k: int):
    """
    Simulates both the greedy and brute force methods of finding a minimum set of genomes that covers all reads and 
    returns average runtime and solution size of both the greedy and brute force methods.

    Input:
        m: number of genomes
        n: number of reads
        k: number of times to simulate
    
    Output:
        Average greedy algorithm runtime, Average greedy algorithm solution size, Average brute-force algorithm runtime, Average brute-force algorithm solution size
    """
    GreedyTime = 0.0
    GreedySize = 0.0
    BruteForceTime = 0.0
    BruteForceSize = 0.0
    for _ in range(k):
        g = GenerateMatrix(m, n)
        start = default_timer()
        subsets = Greedy(m, n, g)
        end = default_timer()
        GreedyTime += start - end
        GreedySize += len(subsets)

        start = default_timer()
        subsets = BruteForce(m, n, g)
        end = default_timer()
        BruteForceTime += end - start
        BruteForceSize += len(subsets)
    GreedyTime /= k
    GreedySize /= k
    
    BruteForceTime /= k
    BruteForceSize /= k
    return ((GreedyTime, GreedySize), (BruteForceTime, BruteForceSize))

def CompareGreedyAndBruteForce(k:int , total_num: int, save = False):
    """
    Compares brute force and greedy approaches for m = n = 1 to total_num.
    Then plots the average runtime and solution size of the greedy and brute force methods.
    Input:
        k (int): How many times to test each size of m
        total_num (int): How large to increase m and n.
        save (bool): Whether or not to save the plots

    Output:
        None
    """
    GreedyTime = [0.0 for _ in range(total_num)]
    BruteForceTime = [0.0 for _ in range(total_num)]

    GreedySize = [0.0 for _ in range(total_num)]
    BruteForceSize = [0.0 for _ in range(total_num)]
    absolute_start = default_timer()
    for i in range(total_num):
        print(f'i = %d \t\t %8.5s' % (i, default_timer()-absolute_start), end = '\r')
        m = i+1
        n = i+1
        stats = CompareGreedyBruteForceWrapper(m, n, k)

        GreedyTime[i] = stats[0][0]
        GreedySize[i] = stats[0][1]

        BruteForceTime[i] = stats[1][0]
        BruteForceSize[i] = stats[1][1]

    if save == True:
        np.savetxt(f'record_{total_num}.txt', np.array([GreedyTime, BruteForceTime, GreedySize, BruteForceSize]))
    # runtime_data = np.loadtxt(f'record_{total_num}.txt')
    # GreedyTime = runtime_data[0]
    # BruteForceTime = runtime_data[1]
    # GreedySize = runtime_data[2]
    # BruteForceSize = runtime_data[3]
    print(f'Finished in %8.5s' % (default_timer()-absolute_start))

    fig, ax = plt.subplots(1, 1)
    ax.plot(GreedyTime, color = 'red', label = 'Greedy')
    ax.plot(BruteForceTime, color = 'blue', label = 'Brute Force')
    ax.set_xlabel('n (Number of Reads and Genomes)')
    ax.set_ylabel('Average Runtime (s)')
    ax.legend()
    ax.set_title('Average Runtime vs Number of Reads and Genomes')
    plt.show()
    if save == True:
        fig.savefig(f'RunTime_{total_num}.png')

    fig2, ax2 = plt.subplots(1, 1)
    ax2.plot(GreedySize, color = 'red', label = 'Greedy')
    ax2.plot(BruteForceSize, color = 'blue', label = 'Brute Force')
    ax2.set_xlabel('n (Number of Reads and Genomes)')
    ax2.set_ylabel('Average Solution Size')
    ax2.set_title('Average Solution Size vs Number of Reads and Genomes')
    ax2.legend()
    plt.show()
    if save == True:
        fig2.savefig(f'Size_{total_num}.png')

def GenerateMatrix(m: int, n: int):
    """
        n: The number of reads
        m: The number of whole genomes
    """
    g = np.random.randint(0, 2, (m,n))
    while np.sum(np.where(np.sum(g, axis = 0)> 0, 1, 0)) < n:
        g = np.random.randint(0, 2, (m,n))
    return g

def Greedy(m: int, n: int, g: np.ndarray):
    """
    Greedy approach for solving the minimum set cover problem for read alignments to genomes
    Input:
    m: The number of whole genomes
    n: The number of reads
    g: m x n matrix of 1's and 0's. If genome i contains read j then g[i, j] = 1, otherwise g[i, j] = 0

    Output:
    The rows in g that are part of the minimum set cover as determined by the greedy approach. (Not Exact)
    """
    remaining = np.ones((n))
    collection = np.zeros((m))
    if np.sum(np.where(np.sum(g, axis = 0)> 0, 1, 0)) < n:
        return tuple(range(m))
    
    while np.sum(remaining) > 0:
        bestSet = 0
        bestScore = 0
        for i in np.nonzero(collection == 0)[0]:
            row = g[i,:]
            score = np.sum(np.multiply(remaining, row))
            if score > bestScore:
                bestSet = i
                bestScore = score
        collection[bestSet] = 1
        remaining = np.multiply(remaining, 1-g[bestSet]) # Update the array holding the remaining reads
    return tuple(np.nonzero(collection)[0])

def BruteForce(m: int, n: int, g: np.ndarray):
    """
    Brute Force approach for solving the minimum set cover problem for read alignments to genomes
    Input:
    m: The number of whole genomes
    n: The number of reads
    g: m x n matrix of 1's and 0's. If genome i contains read j then g[i, j] = 1, otherwise g[i, j] = 0

    Output:
    The rows in g that are part of the minimum set cover.
    """
    foundSolution = False
    collection = (range(m))
    for i in range(m):
        for comb in combinations(range(m), i):
            cover = np.zeros((n))
            for r in comb:
                cover = np.add(cover, g[r])
            cover = np.where(cover > 0, 1, 0)
            if np.sum(cover) == n:
                collection = comb
                foundSolution = True
                break
        if foundSolution == True:
            break
    return collection

if __name__ == '__main__':
    main()
