# -*- coding: utf-8 -*-
"""
This version was created on 03 18 2024
author: Amanda Erin Wilson
@author: @amandaerinwilson


Takes csv of individual residual values (observed-expected) across best models
Tests for normalcy
Tests if residual distribution matches a random normal and random laplace distribution with the same mean and standard devation

https://www.statology.org/normality-test-python/
https://docs.scipy.org/doc/scipy/reference/stats.html#continuous-distributions
"""
import csv
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy.stats import shapiro 
from scipy.stats import kstest
import statistics
import scipy.stats as stats



print("start time: " + datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

file_name_input = input("Enter your csv file of your residuals and models: ") or 'residuals_model_selection_03_16_2024_coarse.csv' 
number_of_models = int(input("Enter the number of models you are testing: ") or '7')
data_points = int(input("Enter the number of data points you have: ")or '11')


class Row_values:
    def __init__(self, row):
        self.t1 = row["t1"]
        self.t2 = row["t2"]
        self.expected_probability_ratio = row["expected_pratio"]
        self.b_alt_func = row["b_alt"]
        self.c_alt_func = row["c_alt"]
        self.d_alt_func = row["d_alt"]
        self.f_alt_func = row["f_alt"]
        self.b_dos = row["b_dos"]
        self.c_dos = row["c_dos"]
        self.d_dos = row["d_dos"]
        self.f_dos = row["f_dos"]
        self.b_non = row["b_non"]
        self.c_non = row["c_non"]
        self.d_non = row["d_non"]
        self.f_non = row["f_non"]
        self.alt = row["percent_alt"]
        self.dos = row["percent_dos"]
        self.non = row["percent_non"]
        self.switch_percent = row["percent_switch"]
        self.observed_pratio = row["observed_pratio"]
        self.residual = row["residual"]
        self.model_identifier = row["model_identifier"]
        self.data_point_identifier = row["data_point_identifier"]

def read_csv(filename):
    with open(filename, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            yield Row_values(row)
        
        
def plot_histogram(residual_list):
    plt.hist(residual_list, alpha=0.5, bins = 10)
    #plt.set_xticks(bins)
    #plt.xticks(residual_list, rotation=90)
    plt.xticks(rotation=90, ha = 'left')
    plt.title('Distribution of Residuals')
    plt.xlabel('Binned Residuals')
    plt.ylabel('Count of Residuals')
    plt.show()   

instances = []
for item in read_csv(file_name_input):
    instances.append(item)

residual_list = []
log_residual_list = []
sq_residual_list = []
for each_instance in range(len(instances)):
    resid = float(instances[each_instance].residual)
    residual_list.append(resid)

plot_histogram(residual_list) 

residual_array = np.array(residual_list)
mean = statistics.mean(residual_list)
sd = statistics.stdev(residual_list)   
##############################################################################

dataset1 = np.random.normal(loc=mean, scale=sd, size=100000)

np.random.seed(1)
dataset2 = stats.norm.rvs(mean, sd, size=100000)

dataset3 = stats.laplace.rvs(mean, sd, size=100000)

#create histograms to visualize values in dataset
plt.hist(dataset1, edgecolor='black')
plt.title('Random Normal Distrubtion #1')
plt.show()   

plt.hist(dataset2, edgecolor='black')
plt.title('Random Normal Distrubtion #2')
plt.show()   

plt.hist(dataset3, edgecolor='black')
plt.title('Random Laplace Distrubtion #3')
plt.show()   

###############################################################################

#create Q-Q plot with standard line added to 
fig1, ax1 = plt.subplots()
sm.qqplot(residual_array, line='s', ax = ax1)
ax1.set_title('QQ Plot of Residuals')
plt.show()


print("\n")
#perform Shapiro-Wilk test for normality
shapiro_result = shapiro(residual_array)
print("Shapiro-Wilk Test Result: ")
print(shapiro_result)
print("\n")


#perform Kolmogorov-Smirnov test against given datasets
kstest_result = kstest(residual_array, dataset1)
print("Kolmogorov-Smirnov Test For Given Random Normal Distribution #1: ")
print(kstest_result) 
print("\n")
 

kstest_result = kstest(residual_array, dataset2)
print("Kolmogorov-Smirnov Test For Given Random Normal Distribution #2: ")
print(kstest_result) 
print("\n")

kstest_result = kstest(residual_array, dataset3)
print("Kolmogorov-Smirnov Test For Given Random Laplace Distribution #3: ")
print(kstest_result) 
print("\n")
###############################################################################
print("end time: " + datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
