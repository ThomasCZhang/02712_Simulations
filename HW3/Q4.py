import numpy as np
import argparse
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file-path', type = str, help='Path to input file')
    parser.add_argument('-x0', type = float, default = 50, help='x0 value')
    parser.add_argument('-y0', type = float, default = 1, help='y0 value')
    parser.add_argument('-a', type = float, default = 1, help='a')
    parser.add_argument('-b', type = float, default = 0.25, help='b')
    parser.add_argument('-c', type = float, default = 0.1, help='c')
    parser.add_argument('-d', type = float, default = 1, help='d')
    parser.add_argument('-dt', type = float, default = 0.01, help='delta t (time step)')
    parser.add_argument('-T', type = float, default = 100, help='Total time to run simulation.')

    args = parser.parse_args()

    x0 = args.x0
    y0 = args.y0
    a = args.a
    b = args.b
    c = args.c
    d = args.d
    delta_t = args.dt
    T = args.T

    print(f'x0: {x0}\ny0: {y0}\na: {a}\nb: {b}\nc: {c}\nd: {d}\ndelta t: {delta_t}\nT: {T}')
    t = np.arange(0, T, delta_t)

    x, y = ForwardEuler(x0, y0, delta_t, T, a, b, c, d)

    PlotGraph(x, y, t,'Forward Euler', 'ForwardEuler')

    x, y = LeapFrog(x0, y0, delta_t, T, a, b, c, d)

    PlotGraph(x, y, t, 'Leap Frog', 'leapfrog')

    x, y = BackwardEuler(x0, y0, delta_t, T, a, b, c, d)

    PlotGraph(x, y, t,'Backward Euler', 'BackwardEuler')
    
    
def PlotGraph(x, y, t, methodName, fileName):
    fig, ax = plt.subplots(1)
    ax.plot(t, x, label = 'Lantern Fly')
    ax.plot(t, y, label = 'Praying Mantis', alpha = 0.5)
    ax.set_xlabel('Time (T)')
    ax.set_ylabel('Population Size')
    ax.legend()
    fig.suptitle(f'{methodName} Simulation \n Lantern Fly vs Praying Mantis Popuplation')
    fig.savefig(f'{fileName}.png')
    # plt.show()

def BackwardEuler(x0, y0, delta_t, T, a, b, c, d):
    num_iters = int(T/delta_t)
    if T/delta_t%1 == 0:
        num_iters = num_iters - 1
    x_arr = np.zeros(num_iters+1)
    y_arr = np.zeros(num_iters+1)
    x_arr[0] = x0
    y_arr[0] = y0
    for i in range(num_iters):
        print(f'{i}', end = "\r")
        x_arr[i+1], y_arr[i+1] = BackwardEulerStep(x_arr[i], y_arr[i], delta_t, a, b, c, d)
    return x_arr, y_arr

def BackwardEulerStep(x, y, delta_t, a, b, c, d):
    "Newton Rathson estimate of the Backward Euler Step."
    p = np.array([[x],[y]])
    for i in range(10):
        p = p - np.linalg.pinv(jacobian_g(p[0, 0], p[1, 0], delta_t, a, b, c, d)) @ g(p[0,0], p[1,0], x, y, delta_t, a, b, c, d)
        p = np.array([[max(p[0, 0], 0)],[max(p[1,0], 0)]])
    return p[0,0], p[1,0]

def g(x, y, x0, y0, delta_t, a, b, c, d):
    arr = np.array([[x - x0 - delta_t*(a*x-b*x*y)],
                    [y - y0 - delta_t*(c*x*y-d*y)]])
    return arr

def jacobian_g(x, y, delta_t, a, b, c, d):
    arr = np.array([[1-delta_t*(a - b*y), delta_t*b*x],
                    [-delta_t*c*y, 1 - delta_t*(c*x-d)]])
    return arr

def LeapFrog(x0, y0, delta_t, T, a, b, c, d):
    num_iters = int(T/delta_t)
    if T/delta_t%1 == 0:
        num_iters = num_iters - 1
    x_arr = np.zeros(num_iters+1)
    y_arr = np.zeros(num_iters+1)
    x_arr[0] = x0
    y_arr[0] = y0
    for i in range(num_iters):
        if i == 0:
            x_arr[i+1], y_arr[i+1] = ForwardEulerStep(x_arr[i], y_arr[i], delta_t, a, b, c, d)
        else:
            x_arr[i+1], y_arr[i+1] = LeapFrogStep(x_arr[i-1], y_arr[i-1], x_arr[i], y_arr[i], delta_t, a, b, c, d)
    return x_arr, y_arr

def LeapFrogStep(u, v, x, y, delta_t, a, b, c, d):
    """
    u, v: values at time i-1
    x, y: value at time i
    """
    x1 = max(0, u + 2*delta_t*dxdt(x, y, a, b))
    y1 = max(0, v + 2*delta_t*dydt(x, y, c, d))
    # x1 = u + 2*delta_t*dxdt(x, y, a, b)
    # y1 = v + 2*delta_t*dydt(x, y, c, d)
    return x1, y1
    
def ForwardEuler(x0, y0, delta_t, T, a, b, c, d):
    num_iters = int(T/delta_t)
    if T/delta_t%1 == 0:
        num_iters = num_iters - 1
    x_arr = np.zeros(num_iters+1)
    y_arr = np.zeros(num_iters+1)
    x_arr[0] = x0
    y_arr[0] = y0
    for i in range(num_iters):
        x_arr[i+1], y_arr[i+1] = ForwardEulerStep(x_arr[i], y_arr[i], delta_t, a, b, c, d)
    return x_arr, y_arr

def ForwardEulerStep(x, y, delta_t, a, b, c, d):
    x1 = max(0, x + delta_t*dxdt(x, y, a, b))
    y1 = max(0, y + delta_t*dydt(x, y, c, d))
    return x1, y1
    
def dxdt(x, y, a, b):
    return a*x - b*x*y

def dydt(x, y, c, d):
    return c*x*y- d*y

if __name__ == '__main__':
    main()