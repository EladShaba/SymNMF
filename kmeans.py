import math
import sys

#many of the functions can be deleted since they are not used in our current program!!!!!!!!!!!!!!!!!!!

class Cluster(object):
    def __init__(self, dataPoint):
        self.centroid = dataPoint.copy()
        self.sumDataPoints = dataPoint.copy()
        self.size = 1

    def distance(self, dataPoint):
        dist = 0
        for i in range(len(self.centroid)):
            dist += (self.centroid[i] - dataPoint[i])**2
        return math.sqrt(dist)

    def addDataPoint(self, dataPoint):
        for i in range(len(self.centroid)):
            self.sumDataPoints[i] += dataPoint[i]
        self.size += 1


# returns the closest cluster and its' index in clusterLst
def minDistCluster(dataPoint, clusterLst):
    min_distance = float("inf")
    min_cluster = 0
    min_cluster_ind = 0
    for i, cluster in enumerate(clusterLst):
        dist = cluster.distance(dataPoint)
        if dist < min_distance:
            min_distance = dist
            min_cluster = cluster
            min_cluster_ind = i

    return min_cluster, min_cluster_ind


# returns true if the distance between the previous centroid and the new centroid is less than epsilon for all centroids 
def deltaCentroids(clusterLst, prevCentroidsLst):
    for i in range(len(clusterLst)):
        if clusterLst[i].distance(prevCentroidsLst[i]) >= 0.0001:
            return False

    return True


# counts the number of lines in the file
def countLines(file):
    N = 0
    for line in file:
        if line != "":
            N+=1
    file.close()
    return N


# check if the user input is legal
def CheckLegal(k,N,iter):
    flag=0
    if(N==0):
        #print error
        print("An Error Has Occurred")
        sys.exit(1)
    try:
        if(k<=1 or k>N or int(k)!=k):
            print("Invalid number of clusters!")
            flag=1
    except:
        print("Invalid number of clusters!")
        flag=1
        
    if(iter<=1 or iter>=1000 or int(iter)!=iter):
        print("Invalid maximum iteration!")
        flag=1
        
    if 1==flag:
        sys.exit(1)


# run the kmeans algorithm on the data in fileName
def fit(fileName, k, iter):

    file = open(fileName, 'r')
    closest_cluster_ind = 0
    
    # create a list of k initial clusters
    clusterLst = []
    for i in range(k):
        line = file.readline()
        dataPointStrLst = line.split(",")
        dataPointLst = [float(num) for num in dataPointStrLst]
        clusterLst.append(Cluster(dataPointLst))

    # run the kmeans algorithm until convergence or maximum number of iterations
    while iter > 0:
        while True:
            line = file.readline()
            if line == "":
                break
            dataPointStrLst = line.split(",")
            dataPointLst = [float(num) for num in dataPointStrLst]
            closest_cluster, closest_cluster_ind = minDistCluster(dataPointLst, clusterLst)
            closest_cluster.addDataPoint(dataPointLst)

        prevCentroidsLst = [cluster.centroid.copy() for cluster in clusterLst]
        for cluster in clusterLst:
            for i in range(len(cluster.sumDataPoints)):
                cluster.centroid[i] = cluster.sumDataPoints[i]/cluster.size
                cluster.sumDataPoints[i] = 0
            cluster.size = 0

        if deltaCentroids(clusterLst, prevCentroidsLst):
            break

        file.close()
        file = open(fileName, 'r')

    file.close()
    
    return clusterLst


# predict the labels for the data-points given a list of clusters
def predict(dataPoints, clusterLst):
    labels = []
    for point in dataPoints:
        closest_cluster, closest_cluster_ind = minDistCluster(point, clusterLst)
        labels.append(closest_cluster_ind)
    
    return labels
    
    
if __name__ == "__main__":
    if len(sys.argv) >3:
        k = int(sys.argv[1])
        iter = int(sys.argv[2])
        fileName = sys.argv[3]
    elif(len(sys.argv)<2):
        print("An Error Has Occurred")
        sys.exit(1)
    else:
        k = int(sys.argv[1])
        iter = 300
        fileName = sys.argv[2]

    # Open the file in read mode
    try:
        file = open(fileName, 'r')
    except:
        print("An Error Has Occurred")
        sys.exit(1)

    
    N = countLines(file)
    CheckLegal(k, N, iter)
    
    clusterLst = fit(fileName, k, iter)
    
    for cluster in clusterLst:
        print(",".join("%.4f" % num for num in cluster.centroid))
    