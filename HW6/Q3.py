# import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

gen = np.random.default_rng(2022)

class Graph:
    def __init__(self):
        self.nodes = []
    
    def remove_node(self, node):
        for neighbor in node.neighbors:
            neighbor.remove_neighbor(node)
            node.remove_neighbor(neighbor)
        self.nodes.remove(node)

    def append(self, element):
        self.nodes.append(element)

    def extend(self, iterable):
        self.nodes.extend(iterable)

    def __getitem__(self, index):
        return self.nodes[index]

    def __array__(self):
        return np.array(self.nodes)
    
    def __len__(self):
        return len(self.nodes)

class Node:
    def __init__(self, val: int = 0) -> None:
        self._adjacencyList = []
        self.neighbors = []
        self.value = val
        self.degree = 0
    
    def add_neighbor(self, node_obj, idx = None):
        self.neighbors.append(node_obj)
        if idx is not None:
            self._adjacencyList.append(idx)
        self.degree += 1

    def remove_neighbor(self, node_obj, idx = None, update_adj = False):
        self.neighbors.remove(node_obj)
        if update_adj:
            if idx is not None:
                self._adjacencyList.remove(idx)
            else:
                raise Exception('Cannot remove index from adjacency list of Node object if the index is not given.')
        self.degree -= 1

    def __str__(self):
        return str(self.value)

    def __add__(self, other):
        if isinstance(other, Node):
            return self.value + other.value
        else:
            return self.value + other
        
    def __sub__(self, other):
        if isinstance(other, Node):
            return self.value - other.value
        else:
            return self.value - other
        
    def __mul__(self, other):
        if isinstance(other, Node):
            return self.value * other.value
        else:
            return self.value * other

    def __truediv__(self, other):
        if isinstance(other, Node):
            return self.value/other.value
        else:
            return self.value/other
    
    def __radd__(self, other):
        return self.value+other
    
    def __rsub__(self, other):
        return other - self.value
    
    def __rmul__(self, other):
        return other*self.value
    
    def __rtruediv__(self, other):
        return other/self.value    

def InitializeGraph(num_nodes, min_node_degree):
    """
    Initialize a graph such that there is there is a minimum required degree for all nodes.
    Input:
    num_nods: The number of nodes in the graph.
    min_node_ndegree: The minimum degree a node can have.
    Output:
    g: Graph object
    """
    if min_node_degree >= num_nodes:
        raise Exception('Degree of a node cannot exceed number of nodes - 1 in a graph')
    
    g = Graph()
    for _ in range(num_nodes):
        g.append(Node())

    all_idxs = np.arange(num_nodes)
    for i in range(num_nodes):
        possible_idxs = np.delete(all_idxs,[i] + g[i]._adjacencyList)
        while g[i].degree < min_node_degree:
            idx_1 = gen.integers(0, possible_idxs.shape[0])
            neighbor_idx = possible_idxs[idx_1]
            possible_idxs = np.delete(possible_idxs, idx_1)

            g[i].add_neighbor(g[neighbor_idx], neighbor_idx)
            g[neighbor_idx].add_neighbor(g[i], i)

    return g


def MoranProcess(g, s, init_infections = 1):
    """
    Implement the Moran birth death process on a graph
    Input:
        g: The graph
        s: The selective fitness
        init_infections: initial number of infections
    Output:
        infected_record: the number of infected overtime.
    """
    # Reset Graph:
    for node in g:
        node.value = 0

    num_nodes = len(g)
    start_nodes = gen.integers(0, num_nodes, size = init_infections)
    for start_node in start_nodes:
        g[start_node].value = 1
    
    num_mutants = np.sum(g)
    mutant_count_record = [num_mutants]

    iteration = 0
    # mem_addresses = np.array([id(node) for node in g])
    while num_mutants > 0 and num_mutants < len(g):
        death_rn = gen.integers(num_nodes)
        death_node = g[death_rn]
        
        num_neighbors = len(death_node.neighbors)
        num_mutant_neighbors = np.sum(death_node.neighbors)

        if num_mutant_neighbors > 0:
            total_neighbor_weights = num_neighbors + num_mutant_neighbors*s
            mutant_neighbor_weight = num_mutant_neighbors*(1+s)
            
            birth_rn = gen.random()*total_neighbor_weights
            num_mutants -= death_node
            if birth_rn < mutant_neighbor_weight:
                num_mutants += 1
                death_node.value = 1
            else:
                death_node.value = 0
        else:
            num_mutants -= death_node
            death_node.value = 0
        # curr_sum = 0
        # i = 0
        # while curr_sum < birth_rn:
        #     curr_sum += 1 + death_node.neighbors[i]*s
        #     if curr_sum >= birth_rn:
        #         parent_node = death_node.neighbors[i]
        #         break
        #     i += 1

        # if parent_node.value and num_infected == 1:
        #     print(death_node._adjacencyList[i] , death_rn, death_node._adjacencyList[i])

        # num_mutants += parent_node
        # num_mutants -= death_node
        # death_node.value = parent_node.value
        mutant_count_record.append(num_mutants)

        # if infected_record[-1] - infected_record[-2] > 1:
        #     raise Exception()

        # neighbor_lengths =  np.array([len(node.neighbors) for node in g])
        iteration += 1
    
    return mutant_count_record

def MultiS(g):
    num_nodes = len(g)
    num_simulations = 250
        
    s_list = np.arange(20)/20
    success = np.zeros(s_list.shape[0])
    
    for i, s in enumerate(s_list):
        for j in range(num_simulations):
            print(f'Simulation {j:< 4d}', end = '\r')
            record = MoranProcess(g, s)
            if record[-1] == num_nodes:
                success[i] += 1
        print('')
    
    success /= num_simulations

    fig, ax = plt.subplots(1)
    ax.scatter(s_list, success)
    ax.set_xlabel('s')
    ax.set_ylabel('Fixation Rate')
    ax.set_title(f's vs Fixation Rate. Fully Connected Graph with {num_nodes} nodes.\n{num_simulations} Simulations per s')
    fig.savefig('q3_multi_s.png')
    # plt.show()

def SingleS(g, s = None):
    num_nodes = len(g)
    num_simulations = 250

    if s is None:
        s = 0.0

    record_list = []
    success_rate = 0
    for i in range(num_simulations):
        print(f'Simulation {i:< 4d}', end = '\r')
        record = MoranProcess(g, s)
        record_list.append(np.array(record)/num_nodes)
        if record[-1] == num_nodes:
            success_rate += 1
    success_rate /= num_simulations

    fig, ax = plt.subplots(1)
    for record in record_list:
        ax.plot(record, alpha = 0.5)
    ax.set_xlabel('Generations')
    ax.set_ylabel('Proportion of Population With Mutation')
    ax.set_title(f'Fully Connected Graph with {num_nodes} nodes.\n{num_simulations} Simulations    s = {s}   Fixation Rate = {success_rate*100:4.1f}%')
    fig.savefig(f'q3_s{s:.1f}.png')
    # plt.show()

def MultiN():
    num_simulations = 200
        
    s = 0.0
    n_list = np.arange(1, 21)*10
    success = np.zeros(n_list.shape[0])
    
    for i, num_nodes in enumerate(n_list):
        g = InitializeGraph(num_nodes, num_nodes-1)
        for j in range(num_simulations):
            print(f'Simulation {j:< 4d}', end = '\r')
            record = MoranProcess(g, s)
            if record[-1] == num_nodes:
                success[i] += 1
        print('')
    
    success /= num_simulations

    line_x = (1+np.arange(50))*4
    line_eq = 1/line_x

    fig, ax = plt.subplots(1)
    ax.scatter(n_list, success)
    line1, = ax.plot(line_x, line_eq, color = 'orange', label = '1/n')
    ax.legend(handles= [line1])
    ax.set_xlabel('number of nodes (n)')
    ax.set_ylabel('Fixation Rate')
    ax.set_title(f'Graph Size (fully connected) vs Fixation Rate. \n{num_simulations} Simulations per n    s = {s}')
    fig.savefig('q3_multi_n.png')

def MultiNodeDegree():
    num_simulations = 200
    num_nodes = 100
        
    s = 40
    node_degrees = np.arange(1, 20)/19*95
    success = np.zeros(node_degrees.shape[0])
    
    mean_node_degree = np.zeros(node_degrees.shape)
    for i, min_node_degree in enumerate(node_degrees):
        g = InitializeGraph(num_nodes, min_node_degree)
        mean_node_degree[i] = np.mean([node.degree for node in g])
        for j in range(num_simulations):
            print(f'Simulation {j:< 4d}', end = '\r')
            record = MoranProcess(g, s)
            if record[-1] == num_nodes:
                success[i] += 1
        print('')
    
    success /= num_simulations

    fig, ax = plt.subplots(1)
    ax.scatter(mean_node_degree, success)
    ax.set_xlabel('Degree of Nodes(n)')
    ax.set_ylabel('Fixation Rate')
    ax.set_title(f'Fully Connected Graph with {num_nodes} \n{num_simulations} Simulations per n    s = {s}')
    fig.savefig(f'q3_multi_node_degree_{s}.png')

if __name__ == '__main__':

    # number_of_nodes = 100
    # min_degree = int(1*number_of_nodes)-1
    # g = InitializeGraph(number_of_nodes, min_degree)

    # for val in [0.0, 0.5, 1.0]:
    #     SingleS(g, val)
    # MultiN()
    # MultiS(g)
    MultiNodeDegree()

    