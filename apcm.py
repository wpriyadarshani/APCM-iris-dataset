
from PIL import Image
import random
import numpy
import pdb

from PIL import Image

import array
import logging

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot

import math

from matplotlib.patches import Rectangle


class Cluster(object):
    def __init__(self):
        self.centroid = None # define centroid for particular cluster


class PCA(object):
    def __init__(self, k = 3, m = 2.0, epsilon = 0.0001, max_fcm_iteration =100, max_apcm_iteration = 30):
        self.pixels = []
        self.k = k
        self.epsilon = epsilon
        self.max_fcm_iteration = max_fcm_iteration
        self.max_apcm_iteration = max_apcm_iteration
        self.clusters = []
        self.max_diff = 10.0
        self.s = 150
        self.m = m
        self.degree_of_membership = []
       


    def openIrishDataset(self):
        #read the file

        file = open("dataset.txt","r")

        for i in range(self.s):
            self.pixels.append(numpy.random.dirichlet(numpy.ones(5), size=1))

        i = 0
        for line in file:
            fields = line.split(",")
            field1 = fields[0]
            field2 = fields[1]
            field3 = fields[2]
            field4 = fields[3]
            field5 = fields[4]

            print field1, field2, field3, field4, field5

            self.pixels[i] = [float(field1), float(field2), float(field3), float(field4), field5]
            i+=1

       


    def run(self):

        for i in range(self.s):
            num_1 = random.randint(1, 2) * 0.1
            num_2 = random.randint(1, 2) * 0.1
            num_3 = 1.0 - (num_1+num_2)
            degreelist = [num_1, num_2, num_3]
            self.degree_of_membership.append(degreelist)

        randomPixels = random.sample(self.pixels, self.k)

        #initialize the clusters
        self.clusters = [None for i in range(self.k)]

        for idx in range(self.k):
            self.clusters[idx] = Cluster()
            self.clusters[idx].centroid = randomPixels[idx]

        #run FCM algorithm
        iterations = 0

        while self.shouldExitFCM(iterations) is False:
            print "FCM Iteration -------------->", iterations  
            self.calculate_centre_vector()
            self.update_membership()
            iterations += 1
        
        iterations = 0

        # APCM
        while self.shouldExitAPCM(iterations) is False:
            print "HELLO I A AM ITERATIONS:", iterations
            self.calculate_centre_vector()
            self.update_degree_of_membershipPCA(iterations)   
            iterations += 1


        return [cluster.centroid for cluster in self.clusters]



    def shouldExitFCM(self, iteration):
        if self.max_diff < self.epsilon:
            return True

        if self.max_fcm_iteration < iteration:
            return True

        return False

    def shouldExitAPCM(self, iterations):
        if iterations <= self.max_apcm_iteration:
            return False
        return True

    def initializeEta(self, idx):
        sum_membership = 0.0
        eta_numerator = 0.0
        eta_k = 1.0

        for i in range(self.s):
            dis = self.calcDistance(self.clusters[idx].centroid, self.pixels[i])

            membership = self.degree_of_membership[i][idx]

            eta_numerator += (membership * dis)

            sum_membership += membership

        eta =eta_numerator / sum_membership
        eta = eta * eta_k
        return eta

    def calcEta(self, idx):
        #mean of the cluster

        sum_numerator = [0,0,0,0]
        count = 0

        cluster_pixel = []

        for i in range(self.s):
            max = 0.0
            highest_index = 0
            for j in range(self.k):
                if (self.degree_of_membership[i][j] > max):
                    max = self.degree_of_membership[i][j]
                    highest_index = j

            if highest_index == idx:
                sum_numerator[0] = self.pixels[i][0] + sum_numerator[0]
                sum_numerator[1] = self.pixels[i][1] + sum_numerator[1] 
                sum_numerator[2] = self.pixels[i][2] + sum_numerator[2] 
                sum_numerator[3] = self.pixels[i][3] + sum_numerator[3] 
                
                count = count+1

                cluster_pixel.append(self.pixels[i])

        mean = [sum_numerator[0]/count, sum_numerator[1]/count, sum_numerator[2]/count, sum_numerator[3]/count]

        # print "\n ______________________> ", mean, count, sum_numerator
        #eta
        sum = 0
        for i in cluster_pixel:
            sum += self.calcDistance(mean,i)


        eta = sum/count

        return eta


    def eta_bar(self, iterations):
        eta_list = []

        if iterations == 0:
            for idx in range(self.k):
                eta_list.append(self.initializeEta(idx))

        else:
            for idx in range(self.k):
                eta_list.append(self.calcEta(idx))

        eta_list.sort()

        print "eta list ", eta_list

        print eta_list[0]

        return eta_list[0]

    def calcGamma(self, idx, iterations):
        eta_bar = self.eta_bar(iterations)
        alpha = 1.0
        eta_j = 0

        if iterations == 0:
            eta_j = self.initializeEta(idx)
        else:
            eta_j = self.calcEta(idx)

        gamma = (eta_bar/alpha) *eta_j

        return gamma,eta_bar, eta_j


    # update the degree of membership for PCA
    def update_degree_of_membershipPCA(self, iterations):
        for idx in range(self.k):
            
            #get gamma for particular cluster
            gamma, eta_bar, eta_j = self.calcGamma(idx, iterations)

            if gamma > 0.0:
                # print "******************* eta", eta
                for i in range(self.s):
                    if (i == 0):
                        print "This is the Update degree centroid number:", idx, self.clusters[idx].centroid

                    dis = pow(self.calcDistance(self.clusters[idx].centroid, self.pixels[i]), 2.0)

                    factor = dis /gamma
                    factor = factor * -1.0

                    updated_membership_degree = math.exp(factor)

                    self.degree_of_membership[i][idx] = updated_membership_degree

    def get_new_value(self, i, j):
        sum = 0.0
        val = 0.0
        p = (2 * (1.0) / (self.m - 1))  # cast to float value or else will round to nearst int
        for k in self.clusters:
            num = self.calcDistance(i, j)
            denom = self.calcDistance(i, k.centroid)
            val = num / denom
            val = pow(val, p)
            sum += val
        return (1.0 / sum)

    def update_membership(self):
        self.max_diff = 0.0

        for idx in range(self.k):
            for i in range(self.s):
                new_uij = self.get_new_value(self.pixels[i], self.clusters[idx].centroid)
                if (i == 0):
                    print "This is the Updatedegree centroid :", idx, self.clusters[idx].centroid
                diff = new_uij - self.degree_of_membership[i][idx]

                if (diff > self.max_diff):
                    self.max_diff = diff
                self.degree_of_membership[i][idx] = new_uij
        return self.max_diff

    def calcDistance(self, a, b):
        dis1 = pow((a[0]-b[0]),2.0)
        dis2 = pow((a[1]-b[1]),2.0)
        dis3 = pow((a[2]-b[2]),2.0)
        dis4 = pow((a[3]-b[3]),2.0)

        sumation = dis1 + dis2 + dis3 + dis4
        result = numpy.sqrt(sumation)
        return result

    # Calculates the centroids using degree of membership and fuzziness.
    def calculate_centre_vector(self):
        for cluster in range(self.k):
            sum_numerator = [0,0,0,0]
            sum_denominator = 0.0
            for i in range(self.s):
                pow_uij= pow(self.degree_of_membership[i][cluster], self.m)
                sum_denominator +=pow_uij
                
                num= (pow_uij * self.pixels[i][0], pow_uij * self.pixels[i][1], pow_uij * self.pixels[i][2], pow_uij * self.pixels[i][3])

                sum_numerator[0] = num[0] + sum_numerator[0]
                sum_numerator[1] = num[1] + sum_numerator[1] 
                sum_numerator[2] = num[2] + sum_numerator[2] 
                sum_numerator[3] = num[3] + sum_numerator[3] 
              

            updatedcluster_center = (sum_numerator[0]/sum_denominator, sum_numerator[1]/sum_denominator, sum_numerator[2]/sum_denominator,sum_numerator[3]/sum_denominator)

            self.clusters[cluster].centroid = updatedcluster_center

    def showclustering(self):

        for i in range(self.s):
            max = 0.0
            highest_index = 0
            # Find the index with highest probability
            for j in range(self.k):
                if (self.degree_of_membership[i][j] > max):
                    max = self.degree_of_membership[i][j]
                    highest_index = j
            # Normalize, set highest prob to 1 rest to zero

            # print self.degree_of_membership[i], max
            for j in range(self.k):

                if (j != highest_index):
                    self.degree_of_membership[i][j] = 0
                else:
                    self.degree_of_membership[i][j] = 1
                    # print j

        f = open('result.txt','w')
        for i in range(self.s):
                # Find the index with highest probability
            for j in range(self.k):
                if self.degree_of_membership[i][j] == 1:
                    # print self.pixels[i], j
                    f.write('%s \t' %self.pixels[i] )
                    f.write('%d \n' %j )

        f.close()


    def analysis(self):

        x=0
        cl1 = 0
        cl2 = 0
        cl3 = 0

        sum = 0

        for i in range(self.s):
            
            x = x +1 

            if x % 50 != 0:
                
                for j in range(self.k):

                    if self.degree_of_membership[i][j] == 1:
                        if j == 0:
                            cl1 = cl1 + 1

                        elif j == 1 :
                            cl2 = cl2 + 1

                        elif j == 2:
                            cl3 = cl3 + 1
            else:

                for j in range(self.k):

                    if self.degree_of_membership[i][j] == 1:
                        if j == 0:
                            cl1 = cl1 + 1

                        elif j == 1 :
                            cl2 = cl2 + 1

                        elif j == 2:
                            cl3 = cl3 + 1

                print "first class ", cl1, i, x
                print "sec class ", cl2
                print "third class ", cl3

                if cl1 < 35:
                    sum += cl1

                if cl2 < 35:
                    sum += cl2

                if cl3 < 35:
                    sum += cl3


                cl1 = 0
                cl2 = 0
                cl3 = 0

                x = 0        

        accuracy =(150.0 - sum) / 150.0 
        accuracy = accuracy * 100.0
        print accuracy
        return accuracy




if __name__ == "__main__":
    p = PCA()
    p.openIrishDataset()
    p.run()
    p.showclustering()
    p.analysis()
    
  
    


