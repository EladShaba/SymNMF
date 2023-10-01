import sys
import numpy as np
import pandas as pd
import symnmfmod as symnmf
import math #check if allowed

# print matrix that is given as list
def print_mat(mat_list): 
    for row in mat_list:
         print(','.join([format(row[i], ".4f") for i in range(len(row))]))


#calculate average element of matrice, used to initiliaze H
def avg_mat_entries(mat, n):
    m = 0.0
    for i in range(n):
        for j in range(n):
            m += mat[i][j]
    m = m / (n * n)
    return m


#check if we assume valid input
# get output from user from cmd
def CMD_input():
    if (len(sys.argv) != 4):
        print("An Error Has Occurred")
        sys.exit(1)
    k = sys.argv[1]
    goal = sys.argv[2]
    file_name = sys.argv[3]

    return int(k), goal, file_name


# return the A matrix - Simlarity matrice
def get_sym(data_points):
    sym_mat = symnmf.sym(data_points)
    if None == sym_mat:
        print("An Error Has Occurred")
        sys.exit(1)
    return sym_mat

# return the D matrix - Diagonal Degree Matrix
def get_ddg(data_points):
    ddg_mat = symnmf.ddg(data_points)
    if None == ddg_mat:
        print("An Error Has Occurred")
        sys.exit(1)
    return ddg_mat

# return the W matrix - the normalized similarity matrix - W = D^(-1/2) * A * D^(-1/2)
def get_norm(data_points):
    norm_mat = symnmf.norm(data_points)
    if None == norm_mat:
        print("An Error Has Occurred")
        sys.exit(1)
    return norm_mat


# return the H matrix
def get_symnmf(data_points, n, k, eps, max_iter):
        W_mat = symnmf.norm(data_points)
        if None == W_mat:
            print("An Error Has Occurred")
            sys.exit(1)
        
        #randomly intiliaze H^0
        m = avg_mat_entries(W_mat, n)
        endpoint = 2 * math.sqrt(m / k)
        H_mat = np.random.uniform(0, endpoint, (n, k))
        #calculate H^i until converge in C
        H_mat = symnmf.updateH(H_mat.tolist(), W_mat, k, eps, max_iter)
        
        if None == H_mat:
            print("An Error Has Occurred")
            sys.exit(1)
        
        return H_mat    


if __name__ == "__main__":
    np.random.seed(0)
    k, goal, file_name = CMD_input()
    max_iter = 300
    eps= 0.0001 # 1e-4
    df = pd.read_csv(file_name, header=None)
    n = df.shape[0]
    data_points = df.to_numpy().tolist()
    
    if goal == "sym":
        sym_mat = get_sym(data_points)
        print_mat(sym_mat)

    elif goal == "ddg":
        ddg_mat = get_ddg(data_points)
        print_mat(ddg_mat)

    elif goal == "norm":
        norm_mat = get_norm(data_points)
        print_mat(norm_mat)

    elif goal == "symnmf":
        H_mat = get_symnmf(data_points, n, k, eps, max_iter)
        print_mat(H_mat)

