import numpy as np
import argparse

def DistSquared(a, b):
    """
    Input:
    a, b: Two pairs of 2d coordinates.
    Output:
    The distance between a and b squared.
    """
    return (a[0]-b[0])**2 + (a[1]-b[1])**2

def Interpolate(A, x, y, sigma = 5):
    """
    Input:
    A (np.ndarray): numpy array of dimension M x N
    x (float): x coordinate of interpolation
    y (float): y coordinate of interpolation
    sigma (float): scaling factor
    """
    M, N = A.shape
    a_hat = 0
    scale_sum = 0
    for i in range(M):
        for j in range(N):
            scale_factor = np.exp(-sigma*DistSquared([x, y],[i, j])) 
            scale_sum += scale_factor
            a_hat += A[i, j]*scale_factor
    return a_hat/scale_sum


def Objective(A, T, x, y, sigma):
    """
    Input:
    A (np.ndarray): numpy array of dimension M x N
    T (np.ndarray): numpy array of dimension m x n
    x (float): x coordinate of interpolation
    y (float): y coordinate of interpolation
    sigma (float): scaling factor
    """
    m, n = T.shape
    score = 0
    for i in range(m):
        for j in range(n):
            score += (T[i, j]-Interpolate(A, x+i, y+j, sigma))**2
    return score

def Gradient(A, T, x, y, dx, dy, sigma):
    """
    Input:
    A (np.ndarray): numpy array of dimension M x N
    T (np.ndarray): numpy array of dimension m x n
    x (float): x coordinate of interpolation
    y (float): y coordinate of interpolation
    dx (float): step size in x direction
    dy (float): step size in y direction
    sigma (float): scaling factor
    
    Output:
    2 x 1 Jacobian Matrix
    """
    return np.array([[Objective(A, T, x+dx, y, sigma)-Objective(A, T, x-dx, y, sigma)],[Objective(A, T, x,y+dy, sigma)-Objective(A, T, x-dx, y-dy, sigma)]])

def Hessian(A, T, x, y, dx, dy, sigma):
    """
    Input:
    A (np.ndarray): numpy array of dimension M x N
    T (np.ndarray): numpy array of dimension m x n
    x (float): x coordinate of interpolation
    y (float): y coordinate of interpolation
    dx (float): step size in x direction
    dy (float): step size in y direction
    Output:
    2 x 2 Hessian matrix
    """
    f_xx = (Objective(A, T, x+2*dx, y, sigma) - 2*Objective(A, T, x, y, sigma) + Objective(A, T, x-2*dx, y, sigma))/(4*dx**2)
    f_xy = (Objective(A, T, x+dx, y+dy, sigma) - Objective(A, T, x-dx, y+dy, sigma) - Objective(A, T, x+dx, y-dy, sigma) + Objective(A, T, x-dx, y-dy, sigma))/(4*dx*dy)
    f_yy = (Objective(A, T, x, y+2*dy, sigma) - 2*Objective(A, T, x, y, sigma) + Objective(A, T, x, y-2*dy, sigma))/(4*dy**2)
    return np.array([[f_xx, f_xy],[f_xy, f_yy]])

def GetCoordinateUpdate(J, H):
    """
    Input:
    J (np.ndarray): a 2 x 1 jacobian matrix
    H (np.ndarray): a 2 x 2 hessian matrix
    Output:
    A 2 x 1 matrix. First element is change in x, second element is change in y.
    """
    if (H[0, 0]*H[1, 1]-H[0, 1]*H[1, 0]) != 0:
        change_in_x = (H[1, 1]*J[0, 0]-H[0, 1]*J[1, 0])/(H[0, 0]*H[1, 1]-H[0, 1]*H[1, 0])
        change_in_y = (H[0, 0]*J[1, 0]-H[1, 0]*J[0, 0])/(H[0, 0]*H[1, 1]-H[0, 1]*H[1, 0])
    else:
        change_in_x = np.inf
        change_in_y = np.inf
    return np.array([[change_in_x], [change_in_y]])

def Inbounds(x, y, M, N):
    """"
    Checks if X and Y are in the boundries.
    
    Inputs:
    x, y: input coordinates
    M, N: matrix shapes. So N is the max "x" coordinate and M is the max "y" coordinate
    """
    if x < 0 or x > N:
        return False
    
    if y < 0 or y > M:
        return False
    
    return True

def NewtonRathson(A, T, dx, dy, delta_x, delta_y, sigma):
    M, N = A.shape
    m, n = T.shape
    best_x = -1
    best_y = -1
    best_score = np.inf
    num_row_start = int((M-m)/delta_y)+1
    num_col_start = int((N-n)/delta_x)+1
    for i in range(num_row_start):
        for j in range(num_col_start):
            print(f'i: {i}/{num_row_start-1} \t j: {j}/{num_col_start-1: <4}', end = "\r")
            x_new = i*delta_x
            y_new = j*delta_y 
            diverged = False
            for k in range(10):
                # print(f"k: {k}\tx: {x_new}\ty: {y_new}")
                grad = Gradient(A, T, x_new, y_new, dx, dy, sigma)
                # print(f'JACOBIAN: {jacob}')
                hess = Hessian(A, T, x_new, y_new, dx, dy, sigma)
                # print(f'HESSIAN: {hess}')
                change = GetCoordinateUpdate(grad, hess)
                x_new -= change[0, 0]
                y_new -= change[1, 0]
                if not Inbounds(x_new, y_new, M-m, N-n):
                    diverged = True
                    break
            if not diverged:
                score = Objective(A, T, x_new, y_new, sigma)
            else:
                score = np.inf
            if score < best_score and not diverged:
                best_score = score
                best_x = x_new
                best_y = y_new
    print('\n')
    return (best_x, best_y)

def ReadTestFile(path):
    with open(path) as f:
        line = [float(x) for x in f.readline().strip().split()]
        delta_x = line[0]
        delta_y = line[1]
        dx = line[2]
        dy = line[3]
        sigma = line[4]
        line = [int(x) for x in f.readline().strip().split()]
        M = line[0]
        N = line[1]
        A = []
        for i in range(N):
            A.append([float(x) for x in f.readline().strip().split()])
        A = np.array(A)
        line = [int(x) for x in f.readline().strip().split()]
        m = line[0]
        n = line[1]
        T = []
        for i in range(n):
            T.append([float(x) for x in f.readline().strip().split()])
        T = np.array(T)
    
    return (A, T, dx, dy, delta_x, delta_y, sigma)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file-path', type = str, help='Path to input file')
    parser.add_argument('--delta-x', type = float, help='Override delta x value')
    parser.add_argument('--delta-y', type = float, help='Override delta y value')
    parser.add_argument('--dx', type = float, help='Override dx value')
    parser.add_argument('--dy', type = float, help='Override dy value')
    parser.add_argument('--sigma', type = float, help='Override delta sigma value')

    args = parser.parse_args()

    filepath = args.file_path
    if filepath:
        print(filepath)
        A, T, dx, dy, delta_x, delta_y, sigma = ReadTestFile(filepath) 
        if args.delta_x:
            delta_x = args.delta_x
        if args.delta_y:
            delta_y = args.delta_y
        if args.dx:
            dx = args.dx
        if args.dy:
            dy = args.dy
        if args.sigma:
            sigma = sigma
        print(f"A: \n{A}\nT: \n{T}\ndx: {dx}\ndy: {dy}\ndelta_x: {delta_x}\ndelta_y: {delta_y}\nsigma: {sigma}")
    else:
        print('No File Chosen. Use -h for options.')
    
    print(NewtonRathson(A, T, dx, dy, delta_x, delta_y, sigma))


if __name__ == '__main__':
    main()

