import pandas as pd
import numpy as np
import free_params
import math
import coping_with_eps_and_resolution
import size_of_genome
import utils

#calculate probability of finding a locus in a NP
#forget about it being untied or tied
def calculate_v1(h, R, b, G):
    r_b = coping_with_eps_and_resolution.radius_of_loci(b, G, R)
    h_eff = h + 2*r_b
    return h_eff/(h_eff + 2*R)

def expected_NP_i_j(c_0,v_0, i, j):
    if (i == 0 and j == 0):
        return math.pow(c_0)
    elif (i == 1 and j == 1):
        return 2*((1-2*v_0+c_0)*c_0 + math.pow(v_0 - c_0, 2))
    elif (i == 2 and j == 2):
        return math.pow(1-2*v_0+c_0, 2)
    elif (i + j == 1):
        return 2 * (v_0 - c_0) * c_0
    elif (i + j == 2):
        return math.pow(v_0 - c_0, 2)
    elif (i + j == 3):
        return 2 * (1-2*v_0+c_0)*(v_0-c_0)
    



def expected_NP_co_segregation_components(c_0, v_0):
    #m_0 = expected fraction of NP with neither A nor B
    m_0 = expected_NP_i_j(c_0, v_0, 0, 0)

    #m_1 = expected fraction of NP with either A or B
    m_1 = 2 * (expected_NP_i_j(c_0, v_0, 1, 0) + expected_NP_i_j(c_0, v_0, 2, 0))

    #m_2 = expected fraction of NP with both A and B
    m_2 = 1 - m_0 - m_1

    return (m_0, m_1, m_2)

# m - co-seg ratio
def solve_for_Pi(m, u_0, t_0, v_0):
    # quadratic expression ax^2 + bx + c = 0
    a = (m * math.pow(t_0, 2)) + (m * math.pow(u_0, 2)) - (2 * m * t_0 * u_0) - math.pow(t_0, 2) - math.pow(u_0, 2) + (2 * t_0 * u_0)
    b = (-2*m*math.pow(u_0, 2)) + (2*m*t_0*u_0) + (2*math.pow(u_0, 2)) - (2*t_0*u_0)
    c = (2*m) + (m*math.pow(u_0, 2)) -1 - math.pow(u_0, 2) + (2*math.pow(v_0, 2))
    temp = (math.pow(b, 2) - (4* a * c))/(2 * a)
    root1 = 0
    root2 = 0
    if temp >= 0:
        root1 = -b + math.pow((math.pow(b, 2) - (4* a * c))/(2 * a), 0.5)
        root2 = -b + math.pow((math.pow(b, 2) + (4* a * c))/(2 * a), 0.5)
    return (root1, root2)


def solve_for_c0(eps, m, v_0):
    c_0_0 = coping_with_eps_and_resolution.detection_eff_coeffs(eps, 0, 0)
    c_1_0 = coping_with_eps_and_resolution.detection_eff_coeffs(eps, 1, 0)
    c_2_0 = coping_with_eps_and_resolution.detection_eff_coeffs(eps, 2, 0)
    
    a = -1*c_0_0 + m*c_0_0 + 4*c_1_0 - m*4*c_1_0 + c_2_0 - m*c_2_0 + 4*m*c_1_0 - 2*m*c_2_0
    b = -4*v_0*c_1_0 + 4*m*v_0*c_1_0 - 2*v_0*c_2_0 + 2*m*v_0*c_2_0 - 4*m*v_0*c_1_0 + 4*m*v_0*c_2_0
    c = 1 + c_2_0*math.pow(v_0, 2) - m - 3*m*c_2_0*math.pow(v_0, 2)

    root_1 = math.sqrt(math.pow(b, 2) - 4*a*c) / (2 * a)
    root_2 = math.sqrt(math.pow(b, 2) + 4*a*c) / (2 * a)
    return (root_1, root_2)

# Compare the experimental value of avg. co-seg ratio for 

# c0 here actually refers to u0 with pi=0
# m = experimental average co-segregation ratio for all locus pairs at distance g
def solve_for_c0_2(v_0, m):
    c_0_0 = coping_with_eps_and_resolution.detection_eff_coeffs(0.83, 0, 0)
    c_1_0 = coping_with_eps_and_resolution.detection_eff_coeffs(0.83, 1, 0)
    c_2_0 = coping_with_eps_and_resolution.detection_eff_coeffs(0.83, 2, 0)
    m2_square_coeff = c_0_0 + 4*c_1_0 + c_2_0
    m1_square_coeff = -4*c_1_0 + 2*c_2_0

    m2_linear_coeff = -4*c_1_0*v_0 - 2*c_2_0*v_0
    m1_linear_coeff = 4*c_1_0*v_0 - 4*v_0*c_2_0

    m2_constant = 1 + c_2_0*math.pow(v_0, 2)
    m1_constant = 2*c_2_0*math.pow(v_0, 2)

    a = m2_square_coeff - m*(m2_square_coeff + m1_square_coeff)
    b = m2_linear_coeff - m*(m2_linear_coeff + m1_linear_coeff)
    c = m2_constant - m*(m2_constant + m1_constant)
    root1 = -10000
    root2 = -10000
    
    if math.pow(b, 2) - 4*a*c >= 0:
        root1 = math.sqrt(math.pow(b, 2) - 4*a*c) / (2 * a)
    if math.pow(b, 2) + 4*a*c >= 0:
        root2 = math.sqrt(math.pow(b, 2) + 4*a*c) / (2 * a)
    return (root1, root2)


def solve_for_pi_2(locus_a_start, locus_b_start, t_0):
    locus_a = df.iloc[start_.index(locus_a_start)]
    locus_b = df.iloc[start_.index(locus_b_start)]
    m_exp_avg = free_params.avg_experimental_co_seg_ratio_with_same_genomic_dist(df_chrom13, locus_a_start, locus_b_start)
    u_0_list = solve_for_c0_2(v_0, m_exp_avg)
    print('u_0 list ', u_0_list)
    try:
        u_0 = next(filter(lambda x: x > 0, u_0_list))
    except(StopIteration):
        u_0 = 0
    c_0 = solve_for_c0_2(v_0, utils.experimental_co_segreation_ratio(locus_a, locus_b))[1]
    pi = (c_0 - u_0)/(t_0 - u_0)
    return pi

if __name__ == "__main__":
    df = pd.read_csv('GAM.csv', index_col=0)

    #select chrom 13
    df_chrom13 = df.loc['chr13']

    start_ = list(df_chrom13['start'].values)
    loci_with_same_genomic_distance = []

   # start_a_val = 3420000
   # start_b_val = 33420000

    locus_a_start = 3750000
   
    locus_b_start = 33750000

    loci_a_start = [3420000, 3690000, 3750000, 3990000, 4020000]
    loci_b_start = [33420000, 33690000, 33750000, 33990000, 34020000]
   

    #sample h,R values in micrometers
    h = 0.22
    R = 4.5
    b = 30000
    G = size_of_genome.genome_size() * 2

    v_1 = calculate_v1(h, R, b, G)
    v_0 = 1 - v_1

    print('v1 ' + str(v_1))
    print('v0 ' + str(v_0))

   
    t_list = free_params.calculate_tied_loci_probability(v_0, v_1)
    for locus_a, locus_b in zip(loci_a_start, loci_b_start):
        print(locus_a, locus_b)
        print(solve_for_pi_2(locus_a, locus_b, t_list[0]))
    #print(solve_for_c0(0.83, m_exp_avg, v_0))
  
    
    #u_list_same_chromosome = free_params.calculate_untied_loci_probability_same_chromosome(v_0, free_params.avg_experimental_co_seg_ratio(df_chrom13, locus_a_start, locus_b_start))
    #print(u_list_same_chromosome)
    #c_0 = 0

    ## TODO need Pi to calculate the actual expected NP that contain 2,1,0 loci

# ---------------------------------------------------------------------------------------------------------------------------------------------# 
    #u_0 = 0.9533351567059093
    #t_0 = t_list[0]
#
    #experimental_co_seg_ration_locus_a_locus_b = utils.experimental_co_segreation_ratio(locus_a, locus_b) 
    #Pi = solve_for_Pi(experimental_co_seg_ration_locus_a_locus_b, u_0, t_0, v_0)
    #print("Pi", Pi)
    #print(Pi)


# TODO: Compute u_0 correctly with average_experimental_co_segregation_ratio - need to find out genomic distance 'g'
# Where do the detection efficiency, genomic resolution play in?
    # detection efficiency - plays in wherever co-segregation ratio is involved
    # genomic resolution - plays in wherever v_1 is involved