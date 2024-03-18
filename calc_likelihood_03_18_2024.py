# -*- coding: utf-8 -*-
"""
This version was created on 03 18 2024
author: Amanda Erin Wilson
@author: @amandaerinwilson


Calculate likelihood for each model from a normal distribution
"""

import csv
from datetime import datetime
import statistics
import numpy as np
import scipy.stats as stats
from scipy.stats import norm



print("start time: " + datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

file_name_input = input("Enter your csv file of your residuals and models: ") or 'residuals_model_selection_03_16_2024_coarse.csv' 
#number_of_models = int(input("Enter the number of models you are testing: ") or '7')
data_points = int(input("Enter the number of data points you have: ")or '11')

number_of_models = 7
model_categories_list = ["3_mix_dup", "alt_dos_dup", "alt_non_dup", "non_dos_dup", "alt_non_mut", "ind", "3_mix_mut"]

###############################################################################
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
        self.model_category = row["model_category"]

def read_csv(filename):
    with open(filename, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            yield Row_values(row)

instances = []
for item in read_csv(file_name_input):
    instances.append(item)

all_residual_list = []
for each_instance in range(len(instances)):
    resid = float(instances[each_instance].residual)
    all_residual_list.append(resid)

###############################################################################

three_mix_dup_residual_list = []
alt_dos_dup_residual_list = []
alt_non_dup_residual_list = []
non_dos_dup_residual_list = []
alt_dos_mut_residual_list = []
alt_non_mut_residual_list = []
ind_residual_list = []
three_mix_mut_residual_list = []
for each_instance in range(len(instances)):
    if instances[each_instance].model_identifier == str(0):
        three_mix_dup_residual_list.append(float(instances[each_instance].residual))
    elif instances[each_instance].model_identifier == str(1):
        alt_dos_dup_residual_list.append(float(instances[each_instance].residual))
    elif instances[each_instance].model_identifier == str(2):
        alt_non_dup_residual_list.append(float(instances[each_instance].residual))
    elif instances[each_instance].model_identifier == str(3):
        non_dos_dup_residual_list.append(float(instances[each_instance].residual))
    elif instances[each_instance].model_identifier == str(4):
        alt_non_mut_residual_list.append(float(instances[each_instance].residual))
    elif instances[each_instance].model_identifier == str(5):
        ind_residual_list.append(float(instances[each_instance].residual))
    elif instances[each_instance].model_identifier == str(6):
        three_mix_mut_residual_list.append(float(instances[each_instance].residual))    
    else:
        pass
###############################################################################

mean = statistics.mean(all_residual_list)
sd = statistics.stdev(all_residual_list)
print("Mean: " + str(mean) + "\n")
print("Standard Deviation: " + str(sd) + "\n")

for i in range(number_of_models):
    three_mix_dup_likelihood =  np.prod(norm.cdf(three_mix_dup_residual_list, loc=mean, scale=sd))
    alt_dos_dup_likelihood =    np.prod(norm.cdf(alt_dos_dup_residual_list, loc=mean, scale=sd))
    alt_non_dup_likelihood =    np.prod(norm.cdf(alt_non_dup_residual_list, loc=mean, scale=sd))
    non_dos_dup_likelihood =    np.prod(norm.cdf(non_dos_dup_residual_list, loc=mean, scale=sd))
    alt_non_mut_likelihood =    np.prod(norm.cdf(alt_non_mut_residual_list, loc=mean, scale=sd))
    ind_likelihood =            np.prod(norm.cdf(ind_residual_list, loc=mean, scale=sd))
    three_mix_mut_likelihood =  np.prod(norm.cdf(three_mix_mut_residual_list, loc=mean, scale=sd))
    file_output = 'likelihoods_norm_03_18_2024.csv'

likelihoods = [three_mix_dup_likelihood, alt_dos_dup_likelihood, alt_non_dup_likelihood, non_dos_dup_likelihood, alt_non_mut_likelihood, ind_likelihood, three_mix_mut_likelihood]
print(likelihoods)

# Normalize likelihoods
total_likelihood = np.sum(likelihoods)
normalized_likelihoods = likelihoods / total_likelihood

file = open(file_output , 'w+', newline='')
writer = csv.writer(file, delimiter=',')
write_header_row = ["Model_Category", "Mean", "SD", "Likelihood", "Normalized_Likelihood"]
writer.writerow(write_header_row)  
# Print results
for i in range(number_of_models):
    print("Model: " + str(model_categories_list[i]) + "\n")
    print("Likelihood: " + str(likelihoods[i]) + "\n")
    print("Normalized Likelihood: " + str(normalized_likelihoods[i]) + "\n")
    row_to_write = [model_categories_list[i], mean, sd, likelihoods[i], normalized_likelihoods[i]]
    writer.writerow(row_to_write)
file.close()

###############################################################################
