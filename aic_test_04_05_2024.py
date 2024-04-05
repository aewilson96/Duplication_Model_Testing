# -*- coding: utf-8 -*-
"""
This version was created on 04 05 2024
author: Amanda Erin Wilson
@author: @amandaerinwilson


Calculate AIC statistic for 
"""
import csv
import numpy as np

input_likelihood = "Actual_"
file_name_input_begin = input_likelihood

file_name_input = input("Enter your csv file of your residuals and models: ") or 'likelihoods_norm_03_25_2024.csv'
file_name_output = file_name_input_begin + 'AIC_output_norm_04_05_2024.txt'

# number of elements
number_of_models = int(input("Enter the number of models you are testing: ") or "7")
 
# Enter the list of numbers for each parameter value
num_params = []
default_num_params = [10, 8, 6, 5, 7, 0, 11]
print("Enter the list of numbers for each parameter value one at a time: ")
for i in range(0, number_of_models):
    each_num_param = int(input()or default_num_params[i])
    num_params.append(each_num_param)  
print(num_params)

# Enter the list of names for each model category
model_categories_list = []
default_model_categories_list = ["3_mix_dup", "alt_dos_dup", "alt_non_dup", "non_dos_dup", "alt_non_mut", "ind", "3_mix_mut"]
print("Enter the list of names for each model category one at a time: ")
for i in range(0, number_of_models):
    each_model_cat = input() or default_model_categories_list[i]
    model_categories_list.append(each_model_cat)  
print(model_categories_list)

###############################################################################
#FUNCTIONS
def calculate_aic(likelihood, num_params):
    return 2 * num_params - 2 * np.log(likelihood)

class Row_values:
    def __init__(self, row):
        self.model_category = row["Model_Category"]
        self.mean = row["Mean"]
        self.SD = row["SD"]
        self.likelihood = row["Likelihood"]
        self.norm_likelihood = row["Normalized_Likelihood"]

def read_csv(filename):
    with open(filename, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            yield Row_values(row)

instances = []
for item in read_csv(file_name_input):
    instances.append(item)

likelihood_list = []
norm_likelihood_list = []
for each_instance in range(len(instances)):
    likelihood = float(instances[each_instance].likelihood)
    likelihood_list.append(likelihood)    
    norm_likelihood = float(instances[each_instance].norm_likelihood)
    norm_likelihood_list.append(norm_likelihood)  

# Calculate AIC for each model
if input_likelihood == "Normalized_":
    aic_values = [calculate_aic(norm_likelihood, params) for norm_likelihood, params in zip(norm_likelihood_list, num_params)]
elif input_likelihood == "Actual_":
    aic_values = [calculate_aic(likelihood, params) for likelihood, params in zip(likelihood_list, num_params)]

# Find the index of the minimum AIC value
best_model_index = np.argmin(aic_values)

f = open(file_name_output, 'w+')
# Print AIC values for each model
for i, aic in enumerate(aic_values):
    print(f"Model {model_categories_list[i]}: AIC = {aic}", file = f)

# Print the Models ordered by best
print(f"\nBest Model: Model {model_categories_list[best_model_index]} with AIC = {aic_values[best_model_index]}", file = f)

sorted_array = np.sort(aic_values)
print("\nOrdered List of AIC Values", file = f)
print(sorted_array, file = f)

f.close()