# -*- coding: utf-8 -*-
"""
This version is created on 03 18 2024

author: Amanda Erin Wilson
@author: @amandaerinwilson

taken and modified from gene_dup_oct_2023_submission.py at https://github.com/aewilson96/Gene_Duplicability_Models
Purpose:
    1) Generate the 3D Surface Plots t1 vs t2 vs pratio
    2) Create csv files that contain a 3D and 2D version of the arrays that contain the pratios for various t1 and t2 values.
    3) Generate the survival curves over time for each category Alt_func, Dos, and Non.
    
B-F parameters:
F+d (rate at which fully redundant genes get lost from the genome) to d (d is the rate at which non duplicated genes are lost)
B and C describe the dynamics from how you move from instantaneous rate to the asymptotic rate (shape of the curve - approx of how [process behaves)
For genes in the Non category, where both copies can only be retained by chance, the parameter values for the survival curve are b = 0, c = 1, d > 10, f = anything. 
For genes in the Dos category, that are sensitive to dosage balance effects, the parameter values are b < 0, 0 < c < 1, d = -f for λ(t)0.02 < 0.1. 
For genes in the Alt_func category, that have potential to subfunctionalize or neofunctionalize, the parameter values are b > 0, c > 0, d > 0, f > 0.   

Additional Parameters
Alt_func : percent of starting genome in Alt_func category (Alpha_Alt_func)
Dos : percent of starting genome in Dos category (Alpha_Dos)
Non : percent of starting genome in Non category (Alpha_Non)
switch : percent of gene duplicates retained through ALt_func that switch to Non category in t2. This value is 0 under the gene duplicability model and >0 for mutational opportunity model (beta_switch_mo)


Results:

"""
import math
import csv
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import sys


input_switch = input("Type '0', '0.1', '0.2' for default results or 'other' to manually input parameter values: ")

###########################################################################
#initialize parameters
#n_max = 170
n_max = 100

#number of time points in t1 and t2 
time_points = 99
#time_points = 51
#time_points = 51
#time_points = 81 #may be more clear

file_name_start = '03_18_2024_results' 

#make graph using the pratio OR the log 10 of pratio
p_ratio = 'pratio'
#p_ratio = 'log of pratio'
###############################################################################
#TOP MODELS RESULTS FROM 03/16/2024

if input_switch == '0':
    number_of_combos = 4
    model_categories = ["3mix", "alt_dos", "alt_non", "non_dos"]
    alts = [0.80, 0.90, 0.60, 0.00]
    doses= [0.10, 0.10, 0.00, 0.10]
    nons =  [0.10, 0.00, 0.40, 0.90]
    b_alt_funcs = [10, 5, 5, 35]
    c_alt_funcs = [1, 5, 1, 0.5]
    d_alt_funcs = [5, 0.0005, 5, 50]
    f_alt_funcs = [2, 2, 10, 10]
    b_doses = [-12, -20, -12, -20]
    c_doses = [0.6, 0.6, 0.6, 0.8]
    d_doses = [-0.03, -0.0003, -0.03, -0.03]
    f_nons = [0.01, 0.01, 0.01, 0.01]
    switch = 0
    file_name_end = '_dup'
    hypothesis  = "duplicability"
elif input_switch == '0.1':
    number_of_combos = 1
    model_categories = ["3mix_mut"]
    alts = [0.80]
    doses= [0.10]
    nons =  [0.10]
    b_alt_funcs = [30]
    c_alt_funcs = [3]
    d_alt_funcs = [0.5]
    f_alt_funcs = [5]
    b_doses = [-12]
    c_doses = [0.6]
    d_doses = [-0.03]
    f_nons = [5]
    switch = 0.1
    file_name_end = '_mut1'
    hypothesis  = "mutational opportunity"
elif input_switch == '0.2':
    number_of_combos = 1
    model_categories = ["alt_non"]
    alts = [0.90]
    doses= [0.00]
    nons =  [0.10]
    b_alt_funcs = [5]
    c_alt_funcs = [5]
    d_alt_funcs = [0.5]
    f_alt_funcs = [8]
    b_doses = [-12]
    c_doses = [0.6]
    d_doses = [-0.03]
    f_nons = [0.01]
    switch = 0.2
    file_name_end = '_mut2'
    hypothesis  = "mutational opportunity"
elif input_switch == 'other':
    number_of_combos = 1
    model_categories = input("Enter your model name: ")
    hypothesis = model_categories
    file_name_end =  model_categories
    alt = float(input("Enter the proportion of your starting genome that starts in the Alt_func category: "))
    dos = float(input("Enter the proportion of your starting genome that starts in the Dos category: "))
    non = float(input("Enter the proportion of your starting genome that starts in the Non category: "))
    total = alt + dos + non
    if total != 1:
        print("your proportions do not add up to 1, please try again")
        sys.exit()
    else:
        pass
    switch = float(input("Enter the proportion of the genes retained through Alt_func switch to the Non category in t2: "))
    if switch < 0 or switch > 1:
        print("you've entered an invalid proportion, please try again")
        sys.exit()
    else:
        pass
    b_alt_func = float(input("Enter the b parameter for the Alt_func category: "))
    if b_alt_func <= 0:
        print("you've entered an invalid parameter value, please try again")
        sys.exit()
    else:
        pass
    c_alt_func = float(input("Enter the c parameter for the Alt_func category: "))
    if c_alt_func <= 0:
        print("you've entered an invalid parameter value, please try again")
        sys.exit()
    else:
        pass
    d_alt_func = float(input("Enter the d parameter for the Alt_func category: "))
    if d_alt_func <= 0:
        print("you've entered an invalid parameter value, please try again")
        sys.exit()
    else:
        pass
    f_alt_func = float(input("Enter the f parameter for the Alt_func category: "))
    if f_alt_func <= 0:
        print("you've entered an invalid parameter value, please try again")
        sys.exit()
    else:
        pass
    b_dos = float(input("Enter the b parameter for the Dos category: "))
    if b_dos >= 0:
        print("you've entered an invalid parameter value, please try again")
        sys.exit()
    else:
        pass
    c_dos = float(input("Enter the c parameter for the Dos category: "))
    if c_dos <= 0 or c_dos > 1:
        print("you've entered an invalid parameter value, please try again")
        sys.exit()
    else:
        pass
    d_dos = float(input("Enter the d parameter for the Dos category: "))
    f_non = float(input("Enter the f parameter for the Non category: "))
    if f_non <= 0:
        print("you've entered an invalid parameter value, please try again")
        sys.exit()
    else:
        pass
    number_of_percent_combos = 1
    alts = []
    doses = []
    nons = []
    b_alt_funcs = []
    c_alt_funcs = []
    d_alt_funcs = []
    f_alt_funcs = []
    b_doses= []
    c_doses = []
    d_doses = []
    b_nons = [0]
    c_nons = [1]
    d_nons = [10.01]
    f_nons = []
    alts.append(alt)
    doses.append(dos)
    nons.append(non)
    b_alt_funcs.append(b_alt_func)
    c_alt_funcs.append(c_alt_func)
    d_alt_funcs.append(d_alt_func)
    f_alt_funcs.append(f_alt_func)
    b_doses.append(b_dos)
    c_doses.append(c_dos)
    d_doses.append(d_dos)
    f_nons.append(f_non)
else:
    print("error")
###########################################################################

#Functions

def calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b, c, d, f, time):  
    summation = 0
    for n in range(0,n_max):
        nfac = math.factorial(n)
        beta = (((-b)**n)*(time**((c*n)+1)))/((c*n*nfac) + nfac)
        summation = summation + beta
        # print("summation: " + str(summation))
    survival_probability = math.exp(-d*time - f*summation)
    return survival_probability
    
def calculate_pratio_2d(st1_alt_func, st1_dos, st1_non, st2_alt_func, st2_dos, st2_non, alt_func_percent, dos_percent, non_percent, alt_switch_percent):    
    #probability of (survival in t2 given survived in t1)/(survival in t2 given lost in t1)
    alt_ret =  (2*alt_func_percent*st1_alt_func)
    alt_ret_ret_switch = alt_ret*st2_non *alt_switch_percent
    alt_ret_ret_noswitch = alt_ret*st2_alt_func * (1-alt_switch_percent)
    dos_ret = (2*dos_percent*st1_dos)
    dos_ret_ret = dos_ret * st2_dos   
    non_ret = (2*non_percent*st1_non)
    non_ret_ret = non_ret * st2_non
    alt_noret = ((1-st1_alt_func)*alt_func_percent)
    alt_noret_ret = alt_noret * st2_alt_func
    dos_noret = ((1-st1_dos)*dos_percent)
    dos_noret_ret = dos_noret * st2_dos    
    non_noret = ((1-st1_non)*non_percent)
    non_noret_ret = non_noret * st2_non
    pratio = ((alt_ret_ret_noswitch + alt_ret_ret_switch + dos_ret_ret + non_ret_ret)/(alt_noret_ret + dos_noret_ret + non_noret_ret)) * ((alt_noret + dos_noret + non_noret)/(alt_ret + dos_ret + non_ret))                                 
    return pratio  

def read_csv_file(file_name):
    file_name_full = (file_name_start + file_name +'.csv')
    df = pd.read_csv(file_name_full,
            header=0,
            usecols=["t1", "t2", "pratio","alt_surv_t1", "dos_surv_t1", "non_surv_t1", "alt_surv_t2", "dos_surv_t2", "non_surv_t2", "log of pratio"])    
    #print(df.head())    
    return df
    
def print_3d_graph(percents, file_name, p_ratio, model_category):
    df = read_csv_file(file_name)
    minimum_pratio = df['pratio'].min()
    print("Minimum Pratio: " + str(minimum_pratio))  
    #plot 3D scatter
    fig2 = plt.figure()
    ax2 = plt.axes(projection='3d')
    surf2 = ax2.scatter3D(df['t1'], df['t2'], df[p_ratio], c = df[p_ratio], cmap=cm.cividis)
    fig2.colorbar(surf2)
    ax2.set_title('Pratio for t1 and t2: ' + percents + model_category, fontsize=14)
    ax2.set_xlabel('$t1$', fontsize=12)
    ax2.set_ylabel('$t2$', fontsize=12)
    ax2.set_zlabel(r'probability ratio', fontsize=11)
    ax2.set_xlim([0, 1.01])
    ax2.set_ylim([0, 1.01])
    #ax2.set_zlim([0.8, 1.01])
    ax2.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax2.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    #ax2.view_init(15, 45)
    ax2.view_init(15, 45)
    plt.draw()
    plt.pause(.001)
    
def plot_survival_curves(file_name):
    df = read_csv_file(file_name)
    #Plot survival curves
    fig, t1_plot = plt.subplots()
    t1_plot.plot(df['t1'], df['alt_surv_t1'], color = 'red', label = 'Alt_func')
    t1_plot.plot(df['t1'], df['dos_surv_t1'], color = 'blue', label = 'Dos')
    t1_plot.plot(df['t1'], df['non_surv_t1'], color = 'yellow', label = 'Non')
    t1_plot.legend(loc = 'upper right', shadow = True, fontsize = '12')
    t1_plot.set_title('Survival over Time (t1)', fontsize=14)
    t1_plot.set_xlabel('Time Since Duplication Event', fontsize=12)
    t1_plot.set_ylabel('Proportion Gene Duplicate Copies Surviving', fontsize=12)
    plt.show()   
    
def calculate_and_make_csv(file_name, alt_func_percent, dos_percent, non_percent, alt_switch_percent, b_alt_func, c_alt_func, d_alt_func, f_alt_func, b_dos, c_dos, d_dos, f_dos, b_non, c_non, d_non, f_non):
    file = open(file_name_start+ file_name +'.csv', 'w', newline='')
    writer = csv.writer(file, delimiter=',')
    write_header_row = ["t1", "t2", "pratio", "alt_surv_t1", "dos_surv_t1", "non_surv_t1", "alt_surv_t2", "dos_surv_t2", "non_surv_t2", "log of pratio"]
    writer.writerow(write_header_row)
    each_t1 = 0.01
    for i in range(0, time_points):  
        alt_func_survival_t1 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_alt_func, c_alt_func, d_alt_func, f_alt_func, each_t1)
        dos_survival_t1 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_dos, c_dos, d_dos, f_dos, each_t1)
        non_survival_t1 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_non, c_non, d_non, f_non, each_t1)
        each_t2 = 0.01
        for j in range(0, time_points):
            #print("t1: " + str(each_t1))
            #print("t2: " + str(each_t2))        
            alt_func_survival_t2 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_alt_func, c_alt_func, d_alt_func, f_alt_func, each_t2)
            dos_survival_t2 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_dos, c_dos, d_dos, f_dos, each_t2)
            non_survival_t2 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_non, c_non, d_non, f_non, each_t2)       
            probability_ratio_2d = calculate_pratio_2d(alt_func_survival_t1, dos_survival_t1, non_survival_t1, alt_func_survival_t2, dos_survival_t2, non_survival_t2, alt_func_percent, dos_percent, non_percent, alt_switch_percent)
            #print("Probability ratio: " + str(probability_ratio_2d))
            #log_pratio = math.log10(probability_ratio_2d)
            row_to_write = [each_t1, each_t2, probability_ratio_2d, alt_func_survival_t1, dos_survival_t1, non_survival_t1, alt_func_survival_t2, dos_survival_t2, non_survival_t2]
            writer.writerow(row_to_write)
            each_t2 = each_t2+0.01
        each_t1 = each_t1+0.01                     
    file.close()
#############################################################################

def main(percent_alt_func, percent_dos, percent_non, percent_alt_switch, b_alt_func, c_alt_func, d_alt_func, f_alt_func, b_dos, c_dos, d_dos, f_dos, b_non, c_non, d_non, f_non, model_category):    
    percentages = str(100*switch) +"% switch,  \n" +str(100*percent_alt_func) +"% Alt_func, "+str(100*percent_dos)+"% Dos, "+str(100*percent_non)+"% Non, \n (" + hypothesis + " Hypothesis)"
    percentages_file_name = str(100*percent_alt_func)+'_'+str(100*percent_dos)+'_'+str(100*percent_non)+file_name_end
    print(percentages)
    calculate_and_make_csv(percentages_file_name, percent_alt_func, percent_dos, percent_non, percent_alt_switch,  b_alt_func, c_alt_func, d_alt_func, f_alt_func, b_dos, c_dos, d_dos, f_dos, b_non, c_non, d_non, f_non)
    print_3d_graph(percentages, percentages_file_name, p_ratio, model_category)
    plot_survival_curves(percentages_file_name)


#############################################################################

for i in range (0, number_of_combos):
    alt = alts[i]
    dos = doses[i]
    non = nons[i]
    b_alt_func = b_alt_funcs[i]
    c_alt_func = c_alt_funcs[i]
    d_alt_func = d_alt_funcs[i]
    f_alt_func = f_alt_funcs[i]
    b_dos = b_doses[i]
    c_dos = c_doses[i]
    d_dos = d_doses[i]
    f_dos = -d_dos
    b_non = 0
    c_non = 1
    d_non = 10.01
    f_non = f_nons[i]
    model_category = model_categories[i]
    main(alt, dos, non, switch,  b_alt_func, c_alt_func, d_alt_func, f_alt_func, b_dos, c_dos, d_dos, f_dos, b_non, c_non, d_non, f_non, model_category)

########################################
