'''
Moran Model Simulation
'''
import numpy as np
import matplotlib.pyplot as plt

gen = np.random.default_rng(2023)

def MoranModel(A0: int, B0: int , weightA: float, weightB: float):
    """
    A0 (int): Starting number of allele A.
    B0 (int): Starting number of allele B.
    weightA (float): The birth/death rate of allele A.
    weightB (float): The birth/death rate of allele B.

    Returns A, B (list[int]): The record of the A and B allele populations
    """

    A = [A0]
    B = [B0]

    while A[-1] > 0 and B[-1] > 0:
        A.append(A[-1])
        B.append(B[-1])

        # Chance of choosing allele A or B for birth: 
        birth_A = weightA*A[-1] / (weightA*A[-1] + weightB*B[-1])
        # birth_B = weightB*B[-1] / (weightA*A[-1] + weightB*B[-1])

        # For birth
        rn = gen.random()
        if rn < birth_A:
            A[-1] += 1
        else:
            B[-1] += 1

        deathA = A[-1]/(A[-1]+B[-1])
        # deathB = B[-1]/(A[-1]+B[-1])

        # For death
        rn = gen.random()
        if rn < deathA:
            A[-1] -= 1
        else:
            B[-1] -= 1

    success = False
    if B[-1] == 0:
        success = True
    
    return A, B, success

def MoranModelQ5(A0: int, BS: float, c: float, d:float, m: int):
    """
    A0 (int): Starting number of allele A.
    B0 (int): Starting number of allele B.
    weightA (float): The birth/death rate of allele A.
    weightB (float): The birth/death rate of allele B.

    Returns A, B (list[int]): The record of the A and B allele populations
    """

    A = [A0]

    while A[-1] > 0 and A[-1] < m:
        A.append(A[-1])

        # Chance of choosing allele A or B for birth: 
        birth = BS*A[-1]
        death = (c+d)*A[-1]
        # stay = 1 - birth - death

        # For birth
        rn = gen.random()
        if rn < birth:
            A[-1] += 1
        elif birth <= rn and rn < birth + death:
            A[-1] -= 1


    success = False
    if A[-1] == m:
        success = True
    
    return A, success


def Question4():
    N = 1000
    A0 = 1
    B0 = N - 1
    weightA = 0.5
    weightB = 1 - weightA

    NumTrials = 1000

    # fig, ax = plt.subplots(1)
    # success = 0
    # for i in range(NumTrials):
    #     A_list, _, suc = MoranModel(A0, B0, weightA, weightB)
    #     if suc:
    #         success += 1
    #     ax.plot(A_list, alpha = 0.5, linewidth = 1)

    # ax.set_xlabel('Time')
    # ax.set_ylabel('Population of A')
    # ax.set_title(f'Population of A vs Time. {NumTrials} Simulations\nFixation Percent: {success/N:.2%}')
    # fig.savefig('Q4a.png')

    fig, ax = plt.subplots(1)
    N = np.arange(100, 1100, 50)
    for i in N:
        print(i)
        A0 = 1
        B0 = i - 1
        success = 0
        for j in range(NumTrials):
            print(j)
            _, _, suc = MoranModel(A0, B0, weightA, weightB)
            if suc:
                success += 1
        ax.scatter(i, success)
    
    ax.set_xlabel('Population Cap')
    ax.set_ylabel('Fixation Rate')
    ax.set_title(f'Population of A vs Population Cap. {NumTrials} simulations')
    fig.savefig('Q4a2.png')
    # s = [0.001, 0.01, 0.1, 0.5, 1]
    # weightB = 1
    # for i in range(len(s)):
    #     fig, ax = plt.subplots(1)
    #     weightA = weightB + s[i]
    #     success = 0
    #     for j in range(N):
    #         A_list, _, suc = MoranModel(A0, B0, weightA, weightB)
    #         if suc:
    #             success += 1
    #         ax.plot(A_list, alpha = 0.5, linewidth = 1)
    #     ax.set_title(f'Population vs Time. {NumTrials} Simulations\nFixation Percent: {success/N:.2%}  s = {s[i]}')
    #     ax.set_xlabel('Time')
    #     ax.set_ylabel('Population of A')
    #     fig.savefig(f'Q4b_{i}.png')
    #     plt.close('all')

def Question5():
    A0 = 1
    BS = np.arange(1, 11, 1)/1000
    c = 0.000
    d = np.arange(1, 11, 1)/10000
    m = np.arange(100, 1001, 100)

    N = 1000

    # BS vs Fixation Rate
    fr_list = []
    for bs in BS:
        print(bs)
        success = 0
        for i in range(N):
            A_list, suc =MoranModelQ5(A0, bs, c, d[-1], m[0])
            if suc:
                success += 1
        fr_list.append(success/N)
    
    fig, ax = plt.subplots(1)
    ax.scatter(BS, fr_list)
    ax.set_title(f'Fixation Rate vs BS. {N} simulations per BS value.\nc+d={c+d[-1]}  m={m[0]}')
    ax.set_xlabel('BS')
    ax.set_ylabel('fixation rate')
    fig.savefig('Q5_bs.png')

    # c+d vs Fixation Rate
    fr_list = []
    for d_i in d:
        success = 0
        for i in range(N):
            A_list, suc =MoranModelQ5(A0, BS[0], c, d_i, m[0])
            if suc:
                success += 1
        fr_list.append(success/N)
    
    fig, ax = plt.subplots(1)
    ax.scatter(d, fr_list)
    ax.set_title(f'Fixation Rate vs c+d. {N} simulations per c+d value.\nBS={BS[0]}  m={m[0]}')
    ax.set_xlabel('c+d')
    ax.set_ylabel('fixation rate')
    fig.savefig('Q5_cd.png')

    # BS vs Fixation Rate
    fr_list = []
    for m_i in m:
        success = 0
        for i in range(N):
            A_list, suc =MoranModelQ5(A0, BS[0], c, d[-1], m_i)
            if suc:
                success += 1
        fr_list.append(success/N)
    
    fig, ax = plt.subplots(1)
    ax.scatter(m, fr_list)
    ax.set_title(f'Fixation Rate vs m. {N} simulations per m value.\nBS={BS[0]}  c+d={c+d[-1]}')
    ax.set_xlabel('m')
    ax.set_ylabel('fixation rate')
    fig.savefig('Q5_m.png')

    # BS/(c+d) 
    fr_list = []
    ratio =  np.arange(10, 30, 1)/10
    for i in ratio:
        print(i)
        success = 0
        for j in range(N):
            A_list, suc =MoranModelQ5(A0, d[-1]*i, c, d[-1], m[0])
            if suc:
                success += 1
        fr_list.append(success/N)
    
    fig, ax = plt.subplots(1)
    ax.scatter(ratio, fr_list)
    ax.set_title(f'Fixation Rate vs BS/(c+d). {N} simulations per ratio.\nm={m[0]}')
    ax.set_xlabel('BS/(c+d)')
    ax.set_ylabel('fixation rate')
    fig.savefig('Q5_ratio.png')
    

if __name__ == '__main__':
    Question4()
    # Question5()