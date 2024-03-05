import numpy as np
import argparse
import matplotlib.pyplot as plt

gen = np.random.default_rng(2023)

def SI(lambda_1, S, I, N):
    return lambda_1 * S * I/N

def IR(lambda_2, I):
    return lambda_2 * I

def IR_2(lambda_3, S, I, N):
    return lambda_3 * S * I /N

def CTMM_Fight( N, m, lambda_1, lambda_2, lambda_3):
    S = N-m
    I = m
    R = 0

    record = [(S, I, R, 0)]
    while I > 0 and S > 0:
        t1 = np.inf
        t2 = np.inf
        t3 = np.inf

        SI_rate = SI(lambda_1, S, I, N)
        IR_rate1 = IR(lambda_2, I)
        IR_rate2 = IR_2(lambda_3, S, I, N)
        if SI_rate != 0:
            t1 = gen.exponential(1/SI_rate)
        if IR_rate1 != 0:
            t2 = gen.exponential(1/IR_rate1)
        if IR_rate2 != 0:    
            t3 = gen.exponential(1/IR_rate2)        

        if t1 < t2 and t1 < t3:
            S -= 1
            I += 1
        else:
            I -= 1
            R += 1
           
    return S, I, R, record

def SimulateCTMM_Fight(N, m, lambda_1, lambda_2, lambda_3, t):
    success = 0
    for _ in range(t):
        S, _, _ , _ = CTMM_Fight( N, m, lambda_1, lambda_2, lambda_3)
        # S, _, _ , rec = CTMM_Fight( N, m, lambda_1, lambda_2, lambda_3)
        if S == 0:
            success += 1

    # length = len(rec)
    # S = [x[0] for x in rec][:length]
    # I = [x[1] for x in rec][:length]
    # R = [x[2] for x in rec][:length]
    # time = [x[3] for x in rec][:length]
    # time = np.log10(time)

    # fig, ax = plt.subplots(1)
    # ax.plot(time, S, label = 'Susceptible')
    # ax.plot(time, I, label = 'Infected')
    # ax.plot(time, R, label = 'Dead')
    # ax.legend()
    # plt.show()
    return success/t

def ReadInputs(filepath):
    params = []
    with open(filepath) as f:
        for line in f:
            line = line.strip().split()
            params.append({'N': int(line[0]),
                           'm': int(line[1]),
                           'lambda_1': float(line[2]),
                           'lambda_2': float(line[3]),
                           'lambda_3': float(line[4]),
                           't': int(line[5])})
    return params


if __name__ == '__main__':
    N = 10000
    m = 100
    lambda_1 = 1.5
    lambda_2s = [0.001, 0.01, 0.1, 1, 10, 100]
    lambda_3 = 0
    t = 100

    success_rates = []
    for lambda_2 in lambda_2s:
        success_rates.append(SimulateCTMM_Fight(N, m, lambda_1, lambda_2, lambda_3, t))

    fig, ax = plt.subplots(1)
    ax.plot(lambda_2s, success_rates, ls = 'None', marker = 'o')
    ax.set_xlabel('$\lambda_2$')
    ax.set_xscale('log')
    ax.set_ylabel('Success Rate')
    ax.set_title('Success Rate Without Fighting vs $\lambda_2$')
    fig.savefig('BaseSix_noFight.png')

    lambda_3 = 1
    success_rates = []
    for lambda_2 in lambda_2s:
        success_rates.append(SimulateCTMM_Fight(N, m, lambda_1, lambda_2, lambda_3, t))

    fig, ax = plt.subplots(1)
    ax.plot(lambda_2s, success_rates, ls = 'None', marker = 'o')
    ax.set_xlabel('$\lambda_2$')
    ax.set_xscale('log')
    ax.set_ylabel('Success Rate')
    ax.set_title('Success Rate With Fighting vs $\lambda_2$')
    fig.savefig('BaseSix_Fight.png')

    # success_rates = []
    
    # lambda_2s = [0.001, 0.01]
    # lambda_2s.extend(np.around(np.arange(0.1, 0.3, 0.01), 2))
    # lambda_2s.extend([1, 10, 100])
    # for lambda_2 in lambda_2s:
    #     print(lambda_2)
    #     success_rates.append(SimulateCTMM_Fight(N, m, lambda_1, lambda_2, lambda_3, t))

    # fig, ax = plt.subplots(1)
    # ax.plot(lambda_2s, success_rates, ls = 'None', marker = 'o')
    # ax.set_xlabel('$\lambda_2$')
    # ax.set_xscale('log')
    # ax.set_ylabel('Success Rate')
    # ax.set_title('Success Rate Without Fighting vs $\lambda_2$')
    # ax.set_xlim(0.09, 0.3)
    # fig.savefig('Extra_noFight.png')
    

    # success_rates = []
    # lambda_2s = [0.001, 0.01]
    # lambda_2s.extend(np.around(np.arange(0.04, 0.08, 0.002), 3))
    # lambda_2s.extend([0.1, 1, 10, 100])
    # lambda_3 = 1

    # for lambda_2 in lambda_2s:
    #     print(lambda_2)
    #     success_rates.append(SimulateCTMM_Fight(N, m, lambda_1, lambda_2, lambda_3, t))

    # fig, ax = plt.subplots(1)
    # ax.plot(lambda_2s, success_rates, ls = 'None', marker = 'o')
    # ax.set_xlabel('$\lambda_2$')
    # ax.set_xscale('log')
    # ax.set_ylabel('Success Rate')
    # ax.set_title('Success Rate With Fighting vs $\lambda_2$')
    # ax.set_xlim(0.04, 0.1)
    # fig.savefig('Extra_Fight.png')



#### Don't Use this version
# def CTMM_Fight( N, m, lambda_1, lambda_2, lambda_3):
#     S = N-m
#     I = m
#     R = 0

#     SI_rate = SI(lambda_1, S, I, N)
#     IR_rate1 = IR(lambda_2, I)
#     IR_rate2 = IR_2(lambda_3, S, I, N)
    
#     t1 = np.inf
#     t2 = np.inf
#     t3 = np.inf
#     if SI_rate != 0:
#         t1 = gen.exponential(1/SI_rate)
#     if IR_rate1 != 0:
#         t2 = gen.exponential(1/IR_rate1)
#     if IR_rate2 != 0:    
#         t3 = gen.exponential(1/IR_rate2)

#     record = [(S, I, R, 0)]
#     while I > 0:
#         if t1 < t2 and t1 < t3:
#             S -= 1
#             I += 1
#             SI_rate = SI(lambda_1, S, I, N)
#             record.append((S,I,R,t1))
#             if SI_rate != 0:
#                 IR_rate1 = IR(lambda_2, I)
#                 IR_rate2 = IR_2(lambda_3, S, I, N)
#                 if IR_rate1 != 0:
#                     t2 = min(t1+gen.exponential(1/IR_rate1), t2)
#                 if IR_rate2 != 0:
#                     t3 = min(t1+gen.exponential(1/IR_rate2), t3)

#                 t1 += gen.exponential(1/SI_rate)
#             else:
#                 t1 = np.inf
#                 break
#         else:
#             I -= 1
#             R += 1
#             if t2 < t3:
#                 IR_rate1 = IR(lambda_2, I)
#                 record.append((S,I,R,t2))
#                 if IR_rate1 != 0:
#                     SI_rate = SI(lambda_1, S, I, N)
#                     IR_rate2 = IR_2(lambda_3, S, I, N)
#                     if SI_rate != 0:
#                         t1 = min(t2 + gen.exponential(1/SI_rate), t1)
#                     if IR_rate2 != 0:
#                         t3 = min(t2 + gen.exponential(1/IR_rate2), t3)

#                     t2 += gen.exponential(1/IR_rate1)
#                 else:
#                     t2 = np.inf
#             else:
#                 IR_rate2 = IR_2(lambda_3, S, I, N)
#                 record.append((S,I,R,t3))
#                 if IR_rate2 != 0:
#                     SI_rate = SI(lambda_1, S, I, N)
#                     IR_rate1 = IR(lambda_2, I)
#                     if SI_rate != 0:
#                         t1 = min(t3 + gen.exponential(1/SI_rate), t1)
#                     if IR_rate1 != 0:
#                         t2 = min(t3 + gen.exponential(1/IR_rate1), t2)

#                     t3 += gen.exponential(1/IR_rate2)
#                 else:
#                     t3 = np.inf
#     return S, I, R, record
####