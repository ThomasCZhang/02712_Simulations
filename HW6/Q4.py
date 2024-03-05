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
        self.nextValue = 0

    
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

def InitializeGraphWithMinDegree(num_nodes, min_node_degree):
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

def InitializeGraphByDensity(num_nodes, edge_density):
    """
    Initialize a graph such that there is the density of edges is edge_density.
    Input:
    num_nods: The number of nodes in the graph.
    min_node_ndegree: The minimum degree a node can have.
    Output:
    g: Graph object
    """
    if edge_density <= 0:
        raise Exception('edge density too low')
    if edge_density > 1:
        raise Exception('Edge density cannot be greater than 1')
    
    g = Graph()
    for _ in range(num_nodes):
        g.append(Node())

    for i in np.arange(num_nodes-1):
        for j in np.arange(i, num_nodes):
            if gen.random() < edge_density:
                g[i].add_neighbor(g[j], j)
                g[j].add_neighbor(g[i], i)

    return g


def sirModel(g, beta, gamma, xi, alpha, init_infections = 1):
    """
    Implement the Moran birth death process on a graph
    Input:
        g: The graph
        s: The selective fitness
        beta: rate at which zombies infect susceptible
        gamma: susceptibe death rate
        zeta: dead resurrected as zombies
        alpha: susceptible defeat zombies   
        init_infections: initial number of infections
    Output:
        infected_record: the number of infected overtime.
    """

    assert beta < 1, 'beta must be less than 1'
    assert gamma < 1, 'gamma must be less than 1'
    assert xi < 1, 'zeta must be less than 1'
    assert alpha < 1, 'alpha must be less than 1'

    # Reset Graph:
    for node in g:
        node.value = 0

    num_nodes = len(g)
    start_nodes = gen.integers(0, num_nodes, size = init_infections)
    for start_node in start_nodes:
        g[start_node].value = 1
    
    num_mutants = np.sum(g)
    num_susceptible = num_nodes - num_mutants
    mutant_count_record = [num_mutants]
    susceptible_count_record = [num_susceptible]
    removed_count_record = [0]

    iteration = 0
    # mem_addresses = np.array([id(node) for node in g])
    while num_mutants > 0 and num_susceptible > 0:
        # For each node update the node value
        for i in range(num_nodes):
            # num_neighbors = len(death_node.neighbors)
            neighbor_arr = np.array([node.value for node in g[i].neighbors])
            # Check whether node is zombie, susceptible or removed.       
            if g[i].value == 1: # Mutant
                num_normal_neighbors = np.sum(neighbor_arr == 0)
                p_zombie_lives = np.power(1-alpha, num_normal_neighbors)
                if gen.random() < p_zombie_lives:
                    g[i].nextValue = 1  
                else:
                    g[i].nextValue = -1 # KILL THE ZOMBIE!!!!!
            elif g[i].value == 0: # Susceptible
                num_mutant_neighbors = np.sum(neighbor_arr == 1)
                p_stay = np.power(1-beta, num_mutant_neighbors)*(1-gamma)
                p_zombify = 1 - p_stay - gamma
                rn = gen.random()
                if rn < p_zombify:
                    g[i].nextValue = 1  # ZOOOOOOMBIFYYYYYY BABY
                elif p_zombify <= rn and rn < p_zombify + gamma:
                    # Dies of old age. RIP
                    g[i].nextValue = -1 
                else:
                    g[i].nextValue = 0               
            elif g[i].value == -1: # Removed
                if gen.random() < xi:
                    # The king is back (from the dead) baby!!!!!!
                    g[i].nextValue = 1 # But only as a zombie :(
                else:
                    g[i].nextValue = -1

        num_susceptible = 0
        num_mutants = 0
        num_removed = 0
        for i in range(num_nodes):
            g[i].value = g[i].nextValue
            if g[i].value == 0:
                num_susceptible += 1
            elif g[i].value == 1:
                num_mutants += 1
            else:
                num_removed += 1
        
        mutant_count_record.append(num_mutants)
        susceptible_count_record.append(num_susceptible)
        removed_count_record.append(num_removed)
        iteration += 1
    
    return susceptible_count_record, mutant_count_record, removed_count_record

def Q5():
    graph_density = 0.1
    num_nodes = 100

    beta = 0.1          # rate at which zombies infect susceptible
    gamma = 0.1        # susceptibe death rate (non-zombie involved)
    xi = 0.1         # dead resurrected as zombies
    alpha = 0.1        # susceptible defeat zombies  
    
    num_simulations = 100

    success_rate = 0

    fig, ax = plt.subplots(1)
    fig1, ax1 = plt.subplots(1)
    fig2, ax2 = plt.subplots(1)

    for i in range(num_simulations):
        g = InitializeGraphByDensity(num_nodes, graph_density)
        S, I, R= sirModel(g, beta, gamma, xi, alpha)
        ax.plot(S)
        ax1.plot(I)
        ax2.plot(R)

        if S[-1] == 0:
            success_rate += 1
    
    success_rate /= num_simulations
    
    ax.set_xlabel(f'Time Steps')
    ax.set_ylabel('Susceptible Count')
    ax.set_title(f'Graph Edge Density={graph_density}  Fixation Rate:{success_rate}\n$\\beta=${beta}   $\\gamma=${gamma}    $\\xi=${xi}    $\\alpha=${alpha}')
    ax1.set_xlabel(f'Time Steps')
    ax1.set_ylabel('Zombie Count')
    ax1.set_title(f'Graph Edge Density={graph_density}  Fixation Rate:{success_rate}\n$\\beta=${beta}   $\\gamma=${gamma}    $\\xi=${xi}    $\\alpha=${alpha}') 
    ax2.set_xlabel(f'Time Steps')
    ax2.set_ylabel('Removed Count')
    ax2.set_title(f'Graph Edge Density={graph_density}  Fixation Rate:{success_rate}\n$\\beta=${beta}   $\\gamma=${gamma}    $\\xi=${xi}    $\\alpha=${alpha}')
    fig.savefig(f'q4_S.png')
    fig1.savefig(f'q4_I.png')
    fig2.savefig(f'q4_R.png')

def Q5_density():
    density_list = np.arange(1,20)/19*.95
    num_nodes = 100

    beta = 0.1          # rate at which zombies infect susceptible
    gamma = 0.1        # susceptibe death rate (non-zombie involved)
    xi = 0.1         # dead resurrected as zombies
    alpha = 0.1        # susceptible defeat zombies     
    

    # beta_list = np.arange(20)/20*0.05         
    num_simulations = 100

    success_rate_list = []
    for graph_density in density_list:
        success_rate = 0
        for _ in range(num_simulations):
            g = InitializeGraphByDensity(num_nodes, graph_density)
            S, I, R= sirModel(g, beta, gamma, xi, alpha)
            if S[-1] == 0:
                success_rate += 1
        success_rate /= num_simulations
        success_rate_list.append(success_rate)

    fig, ax = plt.subplots(1)
    ax.scatter(density_list, success_rate_list)
    ax.set_xlabel('Graph Edge Density')
    ax.set_ylabel('Fixation Rate')
    ax.set_title(f'$\\beta=${beta}    $\\gamma=${gamma}    $\\xi=${xi}    $\\alpha=${alpha}')
    fig.savefig('q4_density.png')

def Q5_beta():
    graph_density = 0.1
    num_nodes = 100

    beta_list = np.arange(1, 21)/20*0.2          # rate at which zombies infect susceptible
    gamma = 0.1        # susceptibe death rate (non-zombie involved)
    xi = 0.1         # dead resurrected as zombies
    alpha = 0.1        # susceptible defeat zombies     
    

    # beta_list = np.arange(20)/20*0.05         
    num_simulations = 100

    success_rate_list = []
    for beta in beta_list:
        success_rate = 0
        for _ in range(num_simulations):
            g = InitializeGraphByDensity(num_nodes, graph_density)
            S, I, R= sirModel(g, beta, gamma, xi, alpha)
            if S[-1] == 0:
                success_rate += 1
        success_rate /= num_simulations
        success_rate_list.append(success_rate)

    fig, ax = plt.subplots(1)
    ax.scatter(beta_list, success_rate_list)
    ax.set_xlabel('$\\beta$')
    ax.set_ylabel('Fixation Rate')
    ax.set_title(f'Graph Edge Density={graph_density}\n$\\gamma=${gamma}    $\\xi=${xi}    $\\alpha=${alpha}')
    fig.savefig('q4_beta.png')

def Q5_gamma():
    graph_density = 0.1
    num_nodes = 100

    beta = 0.1          # rate at which zombies infect susceptible
    gamma_list = np.arange(1, 21)/20*0.2        # susceptibe death rate (non-zombie involved)
    xi = 0.1         # dead resurrected as zombies
    alpha = 0.1        # susceptible defeat zombies     
    

    # beta_list = np.arange(20)/20*0.05         
    num_simulations = 100

    success_rate_list = []
    for gamma in gamma_list:
        success_rate = 0
        for _ in range(num_simulations):
            g = InitializeGraphByDensity(num_nodes, graph_density)
            S, I, R= sirModel(g, beta, gamma, xi, alpha)
            if S[-1] == 0:
                success_rate += 1
        success_rate /= num_simulations
        success_rate_list.append(success_rate)

    fig, ax = plt.subplots(1)
    ax.scatter(gamma_list, success_rate_list)
    ax.set_xlabel('$\\gamma$')
    ax.set_ylabel('Fixation Rate')
    ax.set_title(f'Graph Edge Density={graph_density}\n$\\beta=${beta}    $\\xi=${xi}    $\\alpha=${alpha}')
    fig.savefig('q4_gamma.png')

def Q5_xi():
    graph_density = 0.1
    num_nodes = 100

    beta = 0.1          # rate at which zombies infect susceptible
    gamma = 0.1          # susceptibe death rate (non-zombie involved)
    xi_list = np.arange(1, 21)/20*0.2       # dead resurrected as zombies
    alpha = 0.1        # susceptible defeat zombies     
    

    # beta_list = np.arange(20)/20*0.05         
    num_simulations = 100

    success_rate_list = []
    for xi in xi_list:
        success_rate = 0
        for _ in range(num_simulations):
            g = InitializeGraphByDensity(num_nodes, graph_density)
            S, I, R= sirModel(g, beta, gamma, xi, alpha)
            if S[-1] == 0:
                success_rate += 1
        success_rate /= num_simulations
        success_rate_list.append(success_rate)

    fig, ax = plt.subplots(1)
    ax.scatter(xi_list, success_rate_list)
    ax.set_xlabel('$\\xi$')
    ax.set_ylabel('Fixation Rate')
    ax.set_title(f'Graph Edge Density={graph_density}\n$\\beta=${beta}    $\\gamma=${gamma}    $\\alpha=${alpha}')
    fig.savefig('q4_xi.png')

def Q5_alpha():
    graph_density = 0.1
    num_nodes = 100

    beta = 0.1                              # rate at which zombies infect susceptible
    gamma = 0.1                             # susceptibe death rate (non-zombie involved)
    xi = 0.1                                # dead resurrected as zombies
    alpha_list = np.arange(1, 21)/20*0.2    # susceptible defeat zombies     
    

    # beta_list = np.arange(20)/20*0.05         
    num_simulations = 100

    success_rate_list = []
    for alpha in alpha_list:
        success_rate = 0
        for _ in range(num_simulations):
            g = InitializeGraphByDensity(num_nodes, graph_density)
            S, I, R= sirModel(g, beta, gamma, xi, alpha)
            if S[-1] == 0:
                success_rate += 1
        success_rate /= num_simulations
        success_rate_list.append(success_rate)

    fig, ax = plt.subplots(1)
    ax.scatter(alpha_list, success_rate_list)
    ax.set_xlabel('$\\alpha$')
    ax.set_ylabel('Fixation Rate')
    ax.set_title(f'Graph Edge Density={graph_density}\n$\\beta=${beta}    $\\gamma=${gamma}    $\\xi=${xi}')
    fig.savefig('q4_alpha.png')

if __name__ == '__main__':
    Q5()
    pass