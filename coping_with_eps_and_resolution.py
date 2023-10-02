import math
import utils

# i = alpha = given containing A loci
# j = beta - given containing B loci
def expected_NP_i_j_detection_efficiency(N_i_j, eps, i, j):
    final_sum = 0
    for i_ in range(i, 3):
        for j_ in range(j, 3):
            final_sum += math.pow((1-eps), (i_+j_)-(i+j)) * (utils.kronecker_delta_function(i, 1) * utils.kronecker_delta_function(i_, 2) + 1) * (utils.kronecker_delta_function(j, 1) * utils.kronecker_delta_function(j_, 2) + 1) * N_i_j  
    return final_sum

# i = alpha = given containing A loci
# j = beta - given containing B loci
def detection_eff_coeffs(eps, i, j):
    final_sum = 0
    for i_ in range(i, 3):
        for j_ in range(j, 3):
            final_sum += math.pow((1-eps), (i_+j_)-(i+j)) * (utils.kronecker_delta_function(i, 1) * utils.kronecker_delta_function(i_, 2) + 1) * (utils.kronecker_delta_function(j, 1) * utils.kronecker_delta_function(j_, 2) + 1)  
    return math.pow(eps, i+j)*final_sum

# b = resolutions
# G = size of genome
# R = nuclear radius = 4.5 micrometers
def radius_of_loci(b, G, R):
    return math.pow((b/G) * math.pow(R, 3), 1/3)

