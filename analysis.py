import numpy as np
from sklearn.metrics import silhouette_score
import kmeans as kmeansHW1
import sys
import pandas as pd
import symnmfmod as symnmf
import math


# calculate average element of matrice, used to initiliaze H
def avg_mat_entries(mat,n):
    m = 0.0
    for i in range(n):
        for j in range(n):
            m += mat[i][j]
    m = m / (n*n)
    return m


# returns the H matrice
def getSymnmf(data_points, n, k, eps, max_iter):
        # use module to get the normalized similarity W
        W_mat = symnmf.norm(data_points)
        if None == W_mat:
            print("An Error Has Occurred")
            sys.exit(1)
        
        # randomly intiliaze H^0
        m = avg_mat_entries(W_mat, n)
        endpoint = 2 * math.sqrt(m / k)
        H_mat = np.random.uniform(0, endpoint, (n, k))
        
        #calculate H^i until converge in C
        H_mat = symnmf.updateH(H_mat.tolist(), W_mat, k, eps, max_iter)
        if None == H_mat:
            print("An Error Has Occurred")
            sys.exit(1)
        
        return H_mat


# predicts the labels for the datapoints given an H matrix
def predict(H_mat):
    labels = []
    for i in range(len(H_mat)):
        labels.append(H_mat[i].index(max(H_mat[i])))
    return labels


# print matrix that is given as list
def print_mat(mat_list):
    for row in mat_list:
         print(','.join([format(row[i], ".4f") for i in range(len(row))]))
         
         
if __name__ == "__main__":
    np.random.seed(0)
    
    if (len(sys.argv) != 3):
        print("An Error Has Occurred")
        sys.exit(1)

    k = int(sys.argv[1])
    file_name = sys.argv[2]
    max_iter = 300
    eps= 0.0001 # 1e-4

    # read the data from the file
    df = pd.read_csv(file_name, header=None)
    n = df.shape[0]
    data_points = df.to_numpy().tolist()

    # fit and predict both the kmeans algorithm from HW1 and the symnmf algorithm
    clusterLst = kmeansHW1.fit(file_name, k, max_iter)
    K_labels = kmeansHW1.predict(data_points, clusterLst)
    H_mat = getSymnmf(data_points, n, k, eps, max_iter)
    H_labels = predict(H_mat)

    # Calculate the silhouette scores
    silhouette_score_symnmf = silhouette_score(data_points, H_labels)
    silhouette_score_kmeans = silhouette_score(data_points, K_labels)

    # Print the silhouette scores
    print("nmf: {:.4f}".format(silhouette_score_symnmf))
    print("kmeans: {:.4f}".format(silhouette_score_kmeans))
    