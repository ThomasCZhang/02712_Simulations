import numpy as np
import matplotlib.pyplot as plt

gen = np.random.default_rng(2023)

def ClosedForm(t, n_0, r, K):
    return (n_0*K*np.exp(r*t))/((K-n_0)+n_0*np.exp(r*t))

def Lambda(x, r, K):
    return (1 + r*(1-x/K))*x

def Deterministic(n, r, K):
    x = [n]
    for t in range(100):
        x.append(ClosedForm(t+1, n, r, K))
    return x

def DeterministicForwardEuler(n, r, K):
    x = [n]
    for t in range(100):
        curr_x = x[-1]
        x.append(Lambda(curr_x, r, K))
    return x

def Stochastic(n, r, K):
    x = [n]
    for _ in range(100):
        curr_x = x[-1]
        lambda_1 = Lambda(curr_x, r, K)
        if lambda_1 <= 0:
            break
        x.append(gen.poisson(lambda_1))
    return x

def MakePlot(y1, y2, ax, title):
    ax.plot(y1, label = 'Stochastic')
    ax.plot(y2, label = 'Deterministic')
    ax.legend(loc = 'lower right')
    ax.set_xlabel('Generation')
    ax.set_ylabel('Population')
    ax.set_title(title)

if __name__ == '__main__':
    K = 100
    n = 10

    fig, axes = plt.subplots(2, 1)
    fig.set_size_inches(12, 14)

    r = 0.2
    y1_s1 = Stochastic(n, r, K)
    y2 = DeterministicForwardEuler(n, r, K)
    MakePlot(y1_s1, y2, axes[0], 'Forward Euler Determinisitic vs Stochastic')
    
    y2 = Deterministic(n, r, K)
    MakePlot(y1_s1, y2, axes[1], 'Closed Form Determinisitic vs Stochastic' )
    fig.suptitle('r = 0.2', y = 0.92)
    fig.savefig('question5.png')

    fig, axes = plt.subplots(2, 1)
    fig.set_size_inches(12, 14)

    r = 2.4
    y1_s1 = Stochastic(n, r, K)

    r = 2.7
    y2 = DeterministicForwardEuler(n, r, K)
    MakePlot(y1_s1, y2, axes[0], 'Forward Euler Determinisitic vs Stochastic')
    
    y2 = Deterministic(n, r, K)
    MakePlot(y1_s1, y2, axes[1], 'Closed Form Determinisitic vs Stochastic' )
    fig.suptitle('Stochastic with r = 2.4, Deterministic with r = 2.7', y = 0.92)
    fig.savefig('question5_2.png')

