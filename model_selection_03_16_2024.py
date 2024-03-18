# -*- coding: utf-8 -*-
"""
This version was created on 03 16 2024
author: Amanda Erin Wilson
@author: @amandaerinwilson


SET MODEL CATEGORIES
1.	3 Mixture Gene Duplicability Model (3mix_dup)
2.	2 Mix (Alt+Dos) Gene Duplicability Model (alt_dos_dup)
3.	2 Mix (Alt+Non) Gene Duplicability Model (alt_non_dup)
4.	2 Mix (Non+Dos) Gene Duplicability Model (non_dos_dup)
5.	2 Mix (Alt+Non) Mutational Opportunity Model (alt_non_mut)
6.	Independence Model (ind)
7.	3 Mixture Mutational Opportunity Model (3mix_mut)

Pieces of code was taken from scripts associated with data and figures for 
Wilson AE, Liberles DA. Expectations of Duplicate Gene Retention Under the Gene Duplicability Hypothesis 
at https://github.com/aewilson96/Gene_Duplicability_Models called "gene_dup_oct_2023_submission"

Purpose: 
    1) calculate expected pratio for each set of time points
    2) calculate residual for each set of time points
    3) compare sum of squared residuals across range of parameter values to get the model that gives the smallest value
    
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

"""

import math
import csv
from datetime import datetime
import sys

print("start time: " + datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

input_model = input("Enter your model category. \n     Pre-set options include '3_mix_dup', 'alt_dos_dup', 'alt_non_dup', 'non_dos_dup', 'alt_non_mut', 'ind', '3_mix_mut'.\n     Or you can input 'other' and each parameter value yourself. \nThe Default is 'ind' : \n") or 'ind'


###########################################################################
# #initialize Data Set
data_set = 11

class Fish:
    def __init__(self, name, t1_mya, t2_mya, observed_pratio):
        self.name = name
        self.t1_mya = t1_mya
        self.t2_mya = t2_mya
        self.observed_pratio = observed_pratio
        self.t1 = self.calculate_t(self.t1_mya)
        self.t2 = self.calculate_t(self.t2_mya) 
    def calculate_t(self, t_mya):
        fish_substitution_per_year = 0.00000000413 #BMC Genomics, Fu et al 2010
        t = (t_mya * 1000000 * fish_substitution_per_year)
        return(t)        

    
class Plants:
    def __init__(self, name, t1_mya, t2_mya, observed_pratio):
        self.name = name
        self.t1_mya = t1_mya
        self.t2_mya = t2_mya
        self.observed_pratio = observed_pratio
        self.t1 = self.calculate_t(t1_mya)
        self.t2 = self.calculate_t(t2_mya)  
    def calculate_t(self, t_mya):
        plant_substitution_per_year = 0.000000006 #PNAS Wolfe et al 1987, PNAS Schultz et al 1999
        t = (t_mya * 1000000 * plant_substitution_per_year)
        return(t)        


atlantic_salmon = Fish("Atlantic Salmon", 240, 80, 0.97)
nick_p_equestrius = Plants("Nick's Phalaenopsis equestrius", 54, 76, 0.94)
nick_p1_p_halli = Plants("Nick's Pair 1 Panicum halli", 8, 99, 0.92) # 110 mya, 107.5mya(t1 = 2.5, t2 = 107.5)OR t1 = 0-35, t2 = 100-120
nick_p1_o_brachyantha = Plants("Nick's Pair 1 Oryza brachyantha", 8, 99, 0.93) #110 mya, 107.5mya (t1 = 2.5, t2 = 107.5)OR t1 = 0-35, t2 = 100-120
nick_p2_p_halli = Plants("Nick's Pair 2 Panicum halli", 25, 99, 0.88)#122.5mya, 107.5mya(t1 = 15, t2 = 107.5) OR t1 = 0-40, t2 = 95-120
nick_p2_o_brachyantha = Plants("Nick's Pair 2 Oryza brachyantha", 25, 99, 0.87)#122.5mya,107.5mya (t1 = 15, t2 = 107.5) OR  t1 = 0-40, t2 = 95-120
nick_p3_p_halli = Plants("Nick's Pair 3 Panicum halli", 17, 107, 0.87)#122.5mya, 110mya (t1 =12.5, t2 =110)OR  t1 = 5-25, t2 = 95-120
nick_p3_o_brachyantha = Plants("Nick's Pair 3 Oryza brachyantha", 17, 107, 0.86)#122.5mya, 110mya(t1 =12.5, t2 =110)OR t1 = 5-25, t2 = 95-120
nick_p3_a_comosus = Plants("Nick's Pair 3 Ananas comosus", 17, 107, 0.87)
nick_p4_e_guineensis = Plants("Nick's Pair 4 Elaeis guineensis", 49, 75, 0.91)
nick_p4_p_dactylifera = Plants("Nick's Pair 4 Phoenix dactylifera", 49, 75, 0.90)

t1_set = [atlantic_salmon.t1, nick_p_equestrius.t1, nick_p1_p_halli.t1, nick_p1_o_brachyantha.t1, nick_p2_p_halli.t1, nick_p2_o_brachyantha.t1, nick_p3_p_halli.t1, nick_p3_o_brachyantha.t1, nick_p3_a_comosus.t1, nick_p4_e_guineensis.t1, nick_p4_p_dactylifera.t1]
t2_set = [atlantic_salmon.t2, nick_p_equestrius.t2, nick_p1_p_halli.t2, nick_p1_o_brachyantha.t2, nick_p2_p_halli.t2, nick_p2_o_brachyantha.t2, nick_p3_p_halli.t2, nick_p3_o_brachyantha.t2, nick_p3_a_comosus.t2, nick_p4_e_guineensis.t2, nick_p4_p_dactylifera.t2]
observed_pratio_set = [atlantic_salmon.observed_pratio, nick_p_equestrius.observed_pratio, nick_p1_p_halli.observed_pratio, nick_p1_o_brachyantha.observed_pratio, nick_p2_p_halli.observed_pratio, nick_p2_o_brachyantha.observed_pratio, nick_p3_p_halli.observed_pratio, nick_p3_o_brachyantha.observed_pratio, nick_p3_a_comosus.observed_pratio, nick_p4_e_guineensis.observed_pratio, nick_p4_p_dactylifera.observed_pratio]


###########################################################################
#INITIALIZE PARAMTER VALUES
if input_model == '3mix_dup':
    #1.	3 Mixture Gene Duplicability Model
    file3_name = 'model_selection_03_16_2024_coarse_3mix_dup_sum_of_squares_minimum'
    number_of_percent_combos = 11
    alts = [0.80, 0.10, 0.10,
            0.60, 0.60, 0.60,
            0.40, 0.40, 0.20,
            0.20, 0.20]
    doses = [0.10, 0.10, 0.80,
             0.10, 0.20, 0.30,
             0.40, 0.20, 0.40,
             0.60, 0.20]
    nons = [0.10, 0.80, 0.10,
            0.30, 0.20, 0.10, 
            0.20, 0.40, 0.40,
            0.20, 0.6]
    switches = [0]
    b_alt_funcs = [5, 10, 15, 30, 35]
    c_alt_funcs = [0.5, 1, 3, 5]
    d_alt_funcs = [50, 5, 0.5, 0.05, 0.0005]
    f_alt_funcs = [0.5, 2, 5, 8, 10]
    b_doses= [-1, -10, -12, -14, -16, -18, -20, -30, -50]
    c_doses = [0.025, 0.05, 0.2, 0.4, 0.6, 0.8]
    d_doses = [-0.03, -0.0003, -0.000003, -0.00000003, -0.0000000003]
    b_nons = [0]
    c_nons = [1]
    d_nons = [10.01]
    f_nons = [0.01, 5, 20, 50]
elif input_model == 'alt_dos_dup':
    #2.	2 Mix (Alt+Dos) Gene Duplicability Model
    file3_name = 'model_selection_03_16_2024_coarse_alt_dos_dup_sum_of_squares_minimum'
    number_of_percent_combos = 5
    alts = [0.9, 0.80, 0.60, 0.20, 0.10]
    doses = [0.1, 0.20, 0.40, 0.80, 0.90]
    nons = [0.0, 0.0, 0.0, 0.0, 0.0]
    switches = [0]
    b_alt_funcs = [5, 10, 15, 30, 35]
    c_alt_funcs = [0.5, 1, 3, 5]
    d_alt_funcs = [50, 5, 0.5, 0.05, 0.0005]
    f_alt_funcs = [0.5, 2, 5, 8, 10]
    b_doses= [-1, -10, -12, -14, -16, -18, -20, -30, -50]
    c_doses = [0.025, 0.05, 0.2, 0.4, 0.6, 0.8]
    d_doses = [-0.03, -0.0003, -0.000003, -0.00000003, -0.0000000003]
    b_nons = [0]
    c_nons = [1]
    d_nons = [10.01]
    f_nons = [0.01]
elif input_model == 'alt_non_dup':
    #3.	2 Mix (Alt+Non) Gene Duplicability Model
    file3_name = 'model_selection_03_16_2024_coarse_alt_non_dup_sum_of_squares_minimum'
    number_of_percent_combos = 5
    alts = [0.9, 0.80, 0.60, 0.20, 0.10]
    doses = [0.0, 0.0, 0.0, 0.0, 0.0]
    nons = [0.1, 0.20, 0.40, 0.80, 0.90]
    switches = [0]
    b_alt_funcs = [5, 10, 15, 30, 35]
    c_alt_funcs = [0.5, 1, 3, 5]
    d_alt_funcs = [50, 5, 0.5, 0.05, 0.0005]
    f_alt_funcs = [0.5, 2, 5, 8, 10]
    b_doses= [-12]
    c_doses = [0.6]
    d_doses = [-0.03]
    b_nons = [0]
    c_nons = [1]
    d_nons = [10.01, 20]
    f_nons = [0.01, 5, 20, 50]
elif input_model == 'non_dos_dup':
    #4.	2 Mix (Non+Dos) Gene Duplicability Model
    file3_name = 'model_selection_03_16_2024_coarse_non_dos_dup_sum_of_squares_minimum'
    number_of_percent_combos = 5
    alts = [0.0, 0.0, 0.0, 0.0, 0.0]
    doses = [0.9, 0.80, 0.60, 0.20, 0.10]
    nons = [0.1, 0.20, 0.40, 0.80, 0.90]
    switches = [0]
    b_alt_funcs = [35]
    c_alt_funcs = [0.5]
    d_alt_funcs = [50]
    f_alt_funcs = [10]
    b_doses= [-1, -10, -12, -14, -16, -18, -20, -30, -50]
    c_doses = [0.025, 0.05, 0.2, 0.4, 0.6, 0.8]
    d_doses = [-0.03, -0.0003, -0.000003, -0.00000003, -0.0000000003]
    b_nons = [0]
    c_nons = [1]
    d_nons = [10.01, 20]
    f_nons = [0.01, 5, 10, 20, 50]
elif input_model == 'alt_non_mut':
    #5.	2 Mix (Alt+Non) Mutational Opportunity Model
    file3_name = 'model_selection_03_16_2024_coarse_alt_non_mut_sum_of_squares_minimum'
    number_of_percent_combos = 5
    alts = [0.9, 0.80, 0.60, 0.20, 0.10]
    doses = [0.0, 0.0, 0.0, 0.0, 0.0]
    nons = [0.1, 0.20, 0.40, 0.80, 0.90]
    switches = [0.1, 0.20, 0.5]
    b_alt_funcs = [5, 10, 15, 30, 35]
    c_alt_funcs = [0.5, 1, 3, 5]
    d_alt_funcs = [50, 5, 0.5, 0.05, 0.0005]
    f_alt_funcs = [0.5, 2, 5, 8, 10]
    b_doses= [-12]
    c_doses = [0.6]
    d_doses = [-0.03]
    b_nons = [0]
    c_nons = [1]
    d_nons = [10.01, 20]
    f_nons = [0.01, 5, 20, 50]
elif input_model == 'ind':
    #6.	Independence Model
    file3_name = 'model_selection_03_16_2024_2024_coarse_ind_sum_of_squares_minimum'
    number_of_percent_combos = 1
    alts = [0.0]
    doses = [0.0]
    nons = [1.0]
    switches = [0]
    b_alt_funcs = [35]
    c_alt_funcs = [0.5]
    d_alt_funcs = [50]
    f_alt_funcs = [10]
    b_doses= [-12]
    c_doses = [0.6]
    d_doses = [-0.03]
    b_nons = [0]
    c_nons = [1]
    d_nons = [10.01]
    f_nons = [0.01]
elif input_model == "3mix_mut":
    #7.	3 Mixture Mutational Opportunity Model
    file3_name = 'model_selection_03_16_2024_coarse_3mix_mut_sum_of_squares_minimum'
    number_of_percent_combos = 11
    alts = [0.80, 0.10, 0.10,
            0.60, 0.60, 0.60,
            0.40, 0.40, 0.20,
            0.20, 0.20]
    doses = [0.10, 0.10, 0.80,
             0.10, 0.20, 0.30,
             0.40, 0.20, 0.40,
             0.60, 0.20]
    nons = [0.10, 0.80, 0.10,
            0.30, 0.20, 0.10, 
            0.20, 0.40, 0.40,
            0.20, 0.6]
    switches = [0.1, 0.20, 0.5]
    b_alt_funcs = [5, 10, 15, 30, 35]
    c_alt_funcs = [0.5, 1, 3, 5]
    d_alt_funcs = [50, 5, 0.5, 0.05, 0.0005]
    f_alt_funcs = [0.5, 2, 5, 8, 10]
    b_doses= [-1, -10, -12, -14, -16, -18, -20, -30, -50]
    c_doses = [0.025, 0.05, 0.2, 0.4, 0.6, 0.8]
    d_doses = [-0.03, -0.0003, -0.000003, -0.00000003, -0.0000000003]
    b_nons = [0]
    c_nons = [1]
    d_nons = [10.01]
    f_nons = [0.01, 5, 20, 50]
elif input_model == 'other':
    file3_name = input("Enter your model name: ")
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
    switches = []
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
    switches.append(switch)
    b_alt_funcs.append(b_alt_func)
    c_alt_funcs.append(c_alt_func)
    d_alt_funcs.append(d_alt_func)
    f_alt_funcs.append(f_alt_func)
    b_doses.append(b_dos)
    c_doses.append(c_dos)
    d_doses.append(d_dos)
    f_nons.append(f_non)


###########################################################################

class Data_Point:
    def __init__(self, alt_value, dos_value, non_value, switch_value, b_alt_func_value, c_alt_func_value, d_alt_func_value, f_alt_func_value, b_dos_value, c_dos_value, d_dos_value, f_dos_value, b_non_value, c_non_value, d_non_value, f_non_value, expected_pratio_value, residual_value, absolute_residual_value):
        self.alt_value = alt_value
        self.dos_value = dos_value
        self.non_value = non_value
        self.switch_value = switch_value
        self.b_alt_func_value = b_alt_func_value
        self.c_alt_func_value = c_alt_func_value
        self.d_alt_func_value = d_alt_func_value
        self.f_alt_func_value = f_alt_func_value
        self.b_dos_value = b_dos_value
        self.c_dos_value = c_dos_value
        self.d_dos_value = d_dos_value
        self.f_dos_value = f_dos_value
        self.b_non_value = b_non_value
        self.c_non_value = c_non_value
        self.d_non_value = d_non_value
        self.f_non_value = f_non_value 
        self.expected_pratio_value = expected_pratio_value
        self.residual_value = residual_value
        self.absolute_residual_value = absolute_residual_value
###########################################################################
        
#Functions

def calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b, c, d, f, time):  
    summation = 0
    n_max = 100
    for n in range(0,n_max):
        nfac = math.factorial(n)
        beta = (((-b)**n)*(time**((c*n)+1)))/((c*n*nfac) + nfac)
        summation = summation + beta
        # print("summation: " + str(summation))
    survival_probability = math.exp(-d*time - f*summation)
    return survival_probability
    
def calculate_pratio_2d(st1_alt_func, st1_dos, st1_non, st2_alt_func, st2_dos, st2_non, alt_func_percent, dos_percent, non_percent, alt_switch_percent):    
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
    
def calculate_expected_pratio(t1, t2, alt_func_percent, dos_percent, non_percent, alt_switch_percent, b_alt_func, c_alt_func, d_alt_func, f_alt_func, b_dos, c_dos, d_dos, f_dos, b_non, c_non, d_non, f_non):
    alt_func_survival_t1 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_alt_func, c_alt_func, d_alt_func, f_alt_func, t1)
    dos_survival_t1 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_dos, c_dos, d_dos, f_dos, t1)
    non_survival_t1 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_non, c_non, d_non, f_non, t1)   
    alt_func_survival_t2 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_alt_func, c_alt_func, d_alt_func, f_alt_func, t2)
    dos_survival_t2 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_dos, c_dos, d_dos, f_dos, t2)
    non_survival_t2 = calculate_probability_of_survival_of_duplicate_gene_copy_by_time(b_non, c_non, d_non, f_non, t2)       
    probability_ratio_2d = calculate_pratio_2d(alt_func_survival_t1, dos_survival_t1, non_survival_t1, alt_func_survival_t2, dos_survival_t2, non_survival_t2, alt_func_percent, dos_percent, non_percent, alt_switch_percent)
    return probability_ratio_2d

def calculate_residual(observed_pratio, expected_pratio):
    residual = observed_pratio - expected_pratio
    return residual    

def main(t1, t2, observed_pratio, alt, dos, non, switch, b_alt_func, c_alt_func, d_alt_func, f_alt_func, b_dos, c_dos, d_dos, f_dos, b_non, c_non, d_non, f_non):
    expected_probability_ratio = calculate_expected_pratio(t1, t2, alt, dos, non, switch, b_alt_func, c_alt_func, d_alt_func, f_alt_func, b_dos, c_dos, d_dos, f_dos, b_non, c_non, d_non, f_non)
    residual = calculate_residual(observed_pratio, expected_probability_ratio)
    return expected_probability_ratio, residual

########################################
#LOOP ACROSS ALL VALID PARAMETER VALUES ACROSS GRID
model_identifier = 0
minimum_sum_of_squares = 1
for each_b_alt_func in range(len(b_alt_funcs)):
    for each_c_alt_func in range(len(c_alt_funcs)):
        for each_d_alt_func in range(len(d_alt_funcs)):
            for each_f_alt_func in range(len(f_alt_funcs)):
                for each_b_dos in range(len(b_doses)):
                    for each_c_dos in range(len(c_doses)):
                        for each_d_dos in range(len(d_doses)):
                            b_dos = b_doses[each_b_dos]
                            c_dos = c_doses[each_c_dos]
                            d_dos = d_doses[each_d_dos]
                            f_dos = -1*d_doses[each_d_dos]
                            llambda = f_dos * math.exp((-b_dos*(0.02**c_dos)) + d_dos)
                            if llambda < 0.1:
                                pass
                            elif llambda >= 0.1:
                                continue
                                print('Dos category does not meet requirements of a hazard less than 0.1 for time = 0.02')
                            else:
                                print("error")  
                                pass
                            for each_d_non in range(len(d_nons)):
                                for each_f_non in range(len(f_nons)):
                                    for each_one in range(number_of_percent_combos):
                                        for each_switch in range(len(switches)):   
                                            sum_of_squares_counter = 0
                                            for each_data_point in range(data_set):
                                                t1 = t1_set[each_data_point]
                                                t2 = t2_set[each_data_point]
                                                observed_pratio = observed_pratio_set[each_data_point]       
                                                alt = alts[each_one]
                                                dos = doses[each_one]
                                                non = nons[each_one]
                                                switch = switches[each_switch]
                                                b_alt_func = b_alt_funcs[each_b_alt_func]
                                                c_alt_func = c_alt_funcs[each_c_alt_func]
                                                d_alt_func = d_alt_funcs[each_d_alt_func]
                                                f_alt_func = f_alt_funcs[each_f_alt_func]
                                                b_dos = b_doses[each_b_dos]
                                                c_dos = c_doses[each_c_dos]
                                                d_dos = d_doses[each_d_dos]
                                                f_dos = -1*d_doses[each_d_dos]
                                                b_non = 0
                                                c_non = 1
                                                d_non = d_nons[each_d_non]
                                                f_non = f_nons[each_f_non]  
                                                main_output = main(t1, t2, observed_pratio, alt, dos, non, switch, b_alt_func, c_alt_func, d_alt_func, f_alt_func, b_dos, c_dos, d_dos, f_dos, b_non, c_non, d_non, f_non)
                                                expected_pratio = main_output[0]
                                                residual = main_output[1]
                                                absolute_value_residual = abs(residual)
                                                sum_of_squares_counter = (absolute_value_residual*absolute_value_residual) + sum_of_squares_counter
                                                data1 = Data_Point(alt, dos, non, switch, b_alt_func, c_alt_func, d_alt_func, f_alt_func, b_dos, c_dos, d_dos, f_dos, b_non, c_non, d_non, f_non, expected_pratio, residual, absolute_value_residual) 
                                                if each_data_point == (data_set-1):
                                                    if sum_of_squares_counter < minimum_sum_of_squares:
                                                        minimum_sum_of_squares = sum_of_squares_counter
                                                        row_to_write3 = [model_identifier, sum_of_squares_counter, data1.b_alt_func_value, data1.c_alt_func_value, data1.d_alt_func_value, data1.f_alt_func_value, data1.b_dos_value, data1.c_dos_value, data1.d_dos_value, data1.f_dos_value, data1.b_non_value, data1.c_non_value, data1.d_non_value, data1.f_non_value, data1.alt_value, data1.dos_value, data1.non_value, data1.switch_value]
                                                    else:
                                                        pass
                                                else:
                                                    pass
                                            model_identifier = model_identifier+1
    print(str(each_b_alt_func))

#PRINT BEST MODEL
file3 = open(file3_name +'.csv', 'w+', newline='')
writer3 = csv.writer(file3, delimiter=',')
write_header_row3 = ["model_number", "sum_of_squared_residuals", "b_alt_func", "c_alt_func", "d_alt_func", "f_alt_func", "b_dos", "c_dos", "d_dos", "f_dos", "b_non", "c_non", "d_non", "f_non", "alt_percent", "dos_percent", "non_percent", "percent_switch"]
writer3.writerow(write_header_row3)    
writer3.writerow(row_to_write3)
file3.close()

print("end time: " + datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
