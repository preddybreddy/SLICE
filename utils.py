def kronecker_delta_function(i, j):
    if i == j:
        return 1
    return 0

# New method - try to abandon find_experimental_co_seg
def experimental_co_segreation_ratio(locus_a, locus_b):
    if locus_a.size != locus_b.size:
        print('Number of NP\'s not equal - need for data preprocessing')
        return
    m_2 = 0
    m_1 = 0
    for i in range(locus_a.size):
        if locus_a.iloc[i] + locus_b.iloc[i] == 2:
            m_2 += 1
        elif locus_a.iloc[i] + locus_b.iloc[i] == 1:
            m_1 += 1
    denominatior = 1
    if m_2 + m_1 != 0:
        denominatior = m_2 + m_1
    return m_2 / denominatior

def find_experimental_co_seg(df, locus_a_start, locus_b_start, start_):
    row_a = df.iloc[start_.index(locus_a_start)].values
    row_b = df.iloc[start_.index(locus_b_start)].values
    
    exp_m_2 = 0
    exp_m_1 = 0
    for (i, j) in zip(row_a, row_b):
        if i + j == 2:
            exp_m_2 = exp_m_2 + 1
        elif i + j == 1:
            exp_m_1 = exp_m_1 + 1
    
    # TODO: What happens when m_2 + m_1 = 0
    if (exp_m_2 + exp_m_1 == 0):
        return 0
    
    return exp_m_2 / (exp_m_1 + exp_m_2)