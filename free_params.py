import math
import utils
import numpy as np
import pandas as pd


#probability of finding 2, 1, 0 loci in a non-interacting state in an NP
def calculate_untied_loci_probability_different_chromosome(v_0, v_1):
    # return u_0, u_1, u_2
    return (math.pow(v_0, 2), v_0*v_1, math.pow(v_1, 2))

# Compare the experimental value of avg. co-segregation ratio
# for all couples at a distance 'g'
def avg_experimental_co_seg_ratio_with_same_genomic_dist(df, locus_a_start, locus_b_start):
    diff = locus_b_start - locus_a_start
    rows = df.shape[0]
    start_values = list(df['start'])
    loci_pair_with_same_genomic_distance = []
    for i in range(rows):
        temp_locus_a_start = df.iloc[i].iloc[0]
        temp_locus_b_start = temp_locus_a_start + diff
        if temp_locus_b_start in start_values:
            locus_b_start_index = start_values.index(temp_locus_b_start)
            locus_b = df.iloc[locus_b_start_index]
            loci_pair_with_same_genomic_distance.append(utils.experimental_co_segreation_ratio(df.iloc[i], locus_b))
    return np.mean(loci_pair_with_same_genomic_distance)


def calculate_untied_loci_probability_same_chromosome(v_0, experimental_average_co_seg_ratio):
    #compare the exp co-segregation ratio with the theorectical co-segregation ratio
    #u_0 = math.sqrt(expected_co_segregation_ratio + 2 * math.pow(v_0, 2) - 1) / (1 + expected_co_segregation_ratio) 
    a = 1 - (2 * math.pow(v_0, 2)) - experimental_average_co_seg_ratio
    b = experimental_average_co_seg_ratio - 1
    u_0 = math.sqrt(a/b)
    #u_1 = v_0 - u_0
    #u_2 = 1 - u_0 - u_1
    return u_0

# TODO: Check if d < h
def calculate_tied_loci_probability(v_0, v_1):
    # return t_0, t_1, t_2
    return (v_1, 0, v_0)

if __name__ == '__main__':
    df = pd.read_csv('GAM.csv', index_col=0)
    avg_experimental_co_seg_ratio(df, 3420000, 33420000)