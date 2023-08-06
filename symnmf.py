import sys
import numpy as np
import pandas as pd
import symnmf as symnmf
import math #check if allowed

def print_mat(mat_list):
    # print matrix that is given as list in format ?
    for row in mat_list:
         print(','.join([format(row[i], ".4f") for i in range(len(row))]))

def avg_mat_entries(mat,n):#calculate average element of matrice, used to initiliaze H
    m=0.0
    for i in range(n):
        for j in range(n):
            m+=mat[i][j]
    m = m/(n*n)
    return m

def CMD_input():#check if we assume valid input
    # get output from user from cmd
    if (len(sys.argv) != 4):
        print("An Error Has Occurred")
        sys.exit(1)
    k = sys.argv[1]
    goal = sys.argv[2]
    file_name = sys.argv[3]

    return k, goal, file_name

if __name__ == "__main__":
    np.random.seed(0)
    k, goal, file_name = CMD_input()
    max_iter = 300
    eps= 0.0001 # 1e-4
    df = pd.read_csv(file_name, header=None)
    n = df.shape[0]
    d = df.shape[1]
    data_points = df.to_numpy().tolist()
    
    if goal == "sym":
        sym_mat = symnmf.sym(data_points, n, d)  # use module
        if None == sym_mat:
            print("An Error Has Occurred")
            sys.exit(1)
        print_mat(sym_mat)

    elif goal == "ddg":
        ddg_mat = symnmf.ddg(data_points, n, d)  # use module
        if None == ddg_mat:
            print("An Error Has Occurred")
            sys.exit(1)
        print_mat(ddg_mat)

    elif goal == "norm":
        norm_mat = symnmf.norm(data_points, n, d)  #use module
        if None == norm_mat:
            print("An Error Has Occurred")
            sys.exit(1)
        print_mat(norm_mat)

    elif goal == "symnmf":
        W_mat = symnmf.norm(data_points,n,d) #use module to get the normalized similarity W
        if None == W_mat:
            print("An Error Has Occurred")
            sys.exit(1)
        m = avg_mat_entries(W_mat,n)
        endpoint = 2 * math.sqrt(m/k)
        H_mat = np.random.uniform(0,endpoint,(n,k)) #randomly intiliaze H^0
        #calculate H^i until converge in C
        symnmf.updateH(H_mat, W_mat, n, k, eps, iter)

