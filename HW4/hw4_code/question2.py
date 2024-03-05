import numpy as np
import matplotlib.pyplot as plt
"""
Simulating a SIR continuous time markov model.
"""


def dn1dt(theta, a, c, n1, n2):
    return theta - a*c*n1*n2

def n1null(theta, a, c, n1):
    if n1 != 0:
        return theta/a/c/n1
    else:
        return np.nan

def dn2dt(epsilon, gamma, a, c, n1, n2):
    print(epsilon*a*c*n1*n2-gamma*n2)
    return epsilon*a*c*n1*n2-gamma*n2

def n2null(gamma, epsilon, a, c):
    return gamma/epsilon/a/c


if __name__ == "__main__":
    arrow_scale = 50 # Scale for sizing the vector arrows
    num_ticks = 10
    
    theta = 1000
    epsilon = 0.0005
    gamma = 0.001
    a = 1
    c = 0.01

    # Vector origin location 
    x = np.linspace(0, 1000, num_ticks) 
    y = np.linspace(0, 1000, num_ticks)
    X, Y = np.meshgrid(x, y)
    X = X.flatten()
    Y = Y.flatten()
    assert len(X) == len(Y), 'X and Y lengths are not the same.'
    
    # Directional vectors 
    # U = np.array([dn1dt(theta, a, c, n1, n2) for n1, n2 in zip(X,Y)])   
    # V = np.array([dn2dt(epsilon, gamma, a, c, n1, n2) for n1, n2 in zip(X,Y)])
    U = theta - a*c*X*Y
    V = epsilon*a*c*X*Y-gamma*Y  
    
    fig, ax = plt.subplots(1)
    # Creating plot 
    ax.quiver(X, Y, U, V, color='gray', units='xy', scale=arrow_scale)


    num_ticks = 100
    x = np.linspace(0, 1000, num_ticks) 
    y = np.linspace(0, 1000, num_ticks)

    y_n1null = np.array([n1null(theta, a, c, n1) for n1 in x])
    x_n2null = np.array([n2null(gamma, epsilon, a, c) for _ in x])
    y_n2null_trivial = np.zeros(num_ticks)

    h1, = ax.plot(x, y_n1null, linewidth = 5, linestyle = '--', label = 'n1 null cline')
    h2, = ax.plot(x_n2null, y, color = 'C1', linewidth = 5, linestyle = '--', label = 'n2 null cline')
    ax.plot(x, y_n2null_trivial, color = 'C1', linewidth = 5, linestyle = '--', label = 'n2 null cline')
    fig.legend(handles = [h1, h2], bbox_to_anchor = (1, 0.6))
    ax.set_position((0.125, 0.10, 0.65, 0.8))
    # fig.set_size_inches(10, 7)    

    fig.suptitle('Vector Field') 
    
    # x-lim and y-lim 
    ax.set_xlim(0, 1000) 
    ax.set_ylim(0, 1000) 
    ax.set_xlabel('$n_1$')
    ax.set_ylabel('$n_2$')
    
    # Show plot with grid 
    ax.grid() 
    plt.savefig('question2_with_magnitude.png')


    arrow_scale = 50 # Scale for sizing the vector arrows
    num_ticks = 10
    
    theta = 1000
    epsilon = 0.0005
    gamma = 0.001
    a = 1
    c = 0.01

    # Vector origin location 
    x = np.linspace(0, 1000, num_ticks) 
    y = np.linspace(0, 1000, num_ticks)
    X, Y = np.meshgrid(x, y)
    X = X.flatten()
    Y = Y.flatten()
    assert len(X) == len(Y), 'X and Y lengths are not the same.'
    
    # Directional vectors 
    # U = np.array([dn1dt(theta, a, c, n1, n2) for n1, n2 in zip(X,Y)])   
    # V = np.array([dn2dt(epsilon, gamma, a, c, n1, n2) for n1, n2 in zip(X,Y)])
    U = theta - a*c*X*Y
    V = epsilon*a*c*X*Y-gamma*Y  
    magnitude = np.sqrt(U**2 + V**2)
    U /= magnitude
    V /= magnitude
    
    fig, ax = plt.subplots(1)
    # Creating plot 
    ax.quiver(X, Y, U, V, color='gray', units='xy', scale=0.015)

    num_ticks = 100
    x = np.linspace(0, 1000, num_ticks) 
    y = np.linspace(0, 1000, num_ticks)

    y_n1null = np.array([n1null(theta, a, c, n1) for n1 in x])
    x_n2null = np.array([n2null(gamma, epsilon, a, c) for _ in x])
    y_n2null_trivial = np.zeros(num_ticks)

    h1, = ax.plot(x, y_n1null, linewidth = 5, linestyle = '--', label = 'n1 null cline')
    h2, = ax.plot(x_n2null, y, color = 'C1', linewidth = 5, linestyle = '--', label = 'n2 null cline')
    ax.plot(x, y_n2null_trivial, color = 'C1', linewidth = 5, linestyle = '--', label = 'n2 null cline')
    fig.legend(handles = [h1, h2], bbox_to_anchor = (1, 0.6))
    ax.set_position((0.125, 0.10, 0.65, 0.8))
    # fig.set_size_inches(10, 7)    

    fig.suptitle('Normalized Vector Field') 
    
    # x-lim and y-lim 
    ax.set_xlim(0, 1000) 
    ax.set_ylim(0, 1000) 
    ax.set_xlabel('$n_1$')
    ax.set_ylabel('$n_2$')
    
    # Show plot with grid 
    ax.grid() 

    plt.savefig('question2_normalized.png')