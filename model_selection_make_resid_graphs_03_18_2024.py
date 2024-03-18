# -*- coding: utf-8 -*-
"""
This version was created on 03 18 2024
author: Amanda Erin Wilson
@author: @amandaerinwilson


Takes csv of individual residual values (observed-expected) plotted for each of the models separately as a histogram

"""
import csv
from datetime import datetime
import matplotlib.pyplot as plt


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
        
        
def plot_histogram(residual_list, model_number):
    plt.hist(residual_list)#, bins = 9)
    #plt.set_xticks(bins)
    #plt.xticks(residual_list, rotation=90)
    plt.xticks(rotation=90, ha = 'left')
    plt.title('Distribution of Residuals for Model Number '+ str(model_number))
    plt.xlabel('Binned Residuals')
    plt.ylabel('Count of Residuals')
    plt.show()   

def plot_time_vs_residual(residual_list, time_list, model_number):
    plt.scatter(time_list, residual_list)
    plt.xticks(rotation=90, ha = 'left')
    plt.title('Residuals vs time'+': ' +str(model_number))
    plt.xlabel('Time')
    plt.ylabel('Residuals')
    plt.show()

instances = []
for item in read_csv(file_name_input):
    instances.append(item)
for each_model_number in range(number_of_models):
    residual_list = []
    t1_list = []
    t2_list = []
    for each_instance in range(len(instances)):
        if instances[each_instance].model_identifier == str(each_model_number):
            residual_list.append(float(instances[each_instance].residual))
            t1_list.append(float(instances[each_instance].t1))
            t2_list.append(float(instances[each_instance].t2))
        else:
            continue
    plot_histogram(residual_list, each_model_number)

    
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
