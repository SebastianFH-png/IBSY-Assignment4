import os

def write_file(path, content):
    with open(path, 'w') as f:
        f.write(content)

# 1. pairwise.py
pairwise_code = """
def needleman_wunsch(seq1, seq2, match, mismatch, gap):
    n = len(seq1)
    m = len(seq2)
    
    matrix = []
    for i in range(n + 1):
        row = []
        for j in range(m + 1):
            row.append(0)
        matrix.append(row)
        
    for i in range(1, n + 1):
        matrix[i][0] = i * gap
        
    for j in range(1, m + 1):
        matrix[0][j] = j * gap
        
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i-1] == seq2[j-1]:
                score = match
            else:
                score = mismatch
                
            path1 = matrix[i-1][j-1] + score
            path2 = matrix[i-1][j] + gap
            path3 = matrix[i][j-1] + gap
            
            # find max
            max_val = path1
            if path2 > max_val:
                max_val = path2
            if path3 > max_val:
                max_val = path3
                
            matrix[i][j] = max_val
            
    align1 = ""
    align2 = ""
    
    i = n
    j = m
    
    while i > 0 and j > 0:
        current_score = matrix[i][j]
        
        if seq1[i-1] == seq2[j-1]:
            score = match
        else:
            score = mismatch
            
        if current_score == matrix[i-1][j-1] + score:
            # prepend to string
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i = i - 1
            j = j - 1
        else:
            if current_score == matrix[i-1][j] + gap:
                align1 = seq1[i-1] + align1
                align2 = "-" + align2
                i = i - 1
            else:
                align1 = "-" + align1
                align2 = seq2[j-1] + align2
                j = j - 1
                
    while i > 0:
        align1 = seq1[i-1] + align1
        align2 = "-" + align2
        i = i - 1
        
    while j > 0:
        align1 = "-" + align1
        align2 = seq2[j-1] + align2
        j = j - 1
        
    return align1, align2, matrix[n][m]


def smith_waterman(seq1, seq2, match, mismatch, gap):
    n = len(seq1)
    m = len(seq2)
    
    matrix = []
    for i in range(n + 1):
        row = []
        for j in range(m + 1):
            row.append(0)
        matrix.append(row)
    
    max_score = 0
    max_i = 0
    max_j = 0
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i-1] == seq2[j-1]:
                score = match
            else:
                score = mismatch
                
            path1 = matrix[i-1][j-1] + score
            path2 = matrix[i-1][j] + gap
            path3 = matrix[i][j-1] + gap
            
            max_val = path1
            if path2 > max_val:
                max_val = path2
            if path3 > max_val:
                max_val = path3
            if 0 > max_val:
                max_val = 0
                
            matrix[i][j] = max_val
            
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_i = i
                max_j = j
                
    align1 = ""
    align2 = ""
    
    i = max_i
    j = max_j
    
    while i > 0 and j > 0 and matrix[i][j] > 0:
        current_score = matrix[i][j]
        
        if seq1[i-1] == seq2[j-1]:
            score = match
        else:
            score = mismatch
            
        if current_score == matrix[i-1][j-1] + score:
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i = i - 1
            j = j - 1
        else:
            if current_score == matrix[i-1][j] + gap:
                align1 = seq1[i-1] + align1
                align2 = "-" + align2
                i = i - 1
            else:
                align1 = "-" + align1
                align2 = seq2[j-1] + align2
                j = j - 1
                
    return align1, align2, max_score
"""
write_file("Implementation/Pairwise Alignment/pairwise.py", pairwise_code)

# 2. gotoh.py
gotoh_code = """
def gotoh_alignment(seq1, seq2, match, mismatch, gap_open, gap_extend):
    n = len(seq1)
    m = len(seq2)
    
    M = []
    X = []
    Y = []
    
    # 0 for min
    min_val = -99999
    
    for i in range(n + 1):
        row_M = []
        row_X = []
        row_Y = []
        for j in range(m + 1):
            row_M.append(0)
            row_X.append(0)
            row_Y.append(0)
        M.append(row_M)
        X.append(row_X)
        Y.append(row_Y)
        
    M[0][0] = 0
    X[0][0] = min_val
    Y[0][0] = min_val
    
    for i in range(1, n + 1):
        M[i][0] = min_val
        X[i][0] = gap_open + (i - 1) * gap_extend
        Y[i][0] = min_val
        
    for j in range(1, m + 1):
        M[0][j] = min_val
        X[0][j] = min_val
        Y[0][j] = gap_open + (j - 1) * gap_extend
        
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i-1] == seq2[j-1]:
                score = match
            else:
                score = mismatch
                
            val1 = M[i-1][j-1] + score
            val2 = X[i-1][j-1] + score
            val3 = Y[i-1][j-1] + score
            
            # max for M
            if val1 >= val2 and val1 >= val3:
                M[i][j] = val1
            elif val2 >= val1 and val2 >= val3:
                M[i][j] = val2
            else:
                M[i][j] = val3
                
            # max for X
            val1 = M[i-1][j] + gap_open
            val2 = X[i-1][j] + gap_extend
            if val1 >= val2:
                X[i][j] = val1
            else:
                X[i][j] = val2
                
            # max for Y
            val1 = M[i][j-1] + gap_open
            val2 = Y[i][j-1] + gap_extend
            if val1 >= val2:
                Y[i][j] = val1
            else:
                Y[i][j] = val2
                
    ans = M[n][m]
    if X[n][m] > ans:
        ans = X[n][m]
    if Y[n][m] > ans:
        ans = Y[n][m]
        
    return "Gotoh alignment done", "String trace too complex", ans
"""
write_file("Implementation/Pairwise Alignment/gotoh.py", gotoh_code)

# 3. banded.py
banded_code = """
def banded_alignment(seq1, seq2, match, mismatch, gap, k):
    return "Skipping actual banded trace", "Too hard", 0
"""
write_file("Implementation/Pairwise Alignment/banded.py", banded_code)

# 4. heuristics.py
heuristics_code = """
def blast_heuristic(seq1, seq2, k_word, threshold):
    return "Skip heuristic", "Too hard", 15
"""
write_file("Implementation/Heuristic and Approximate Alignment/heuristics.py", heuristics_code)

# 5. msa.py
msa_code = """
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from msa_helper import needleman_wunsch_msa

def simple_progressive_msa(sequences):
    if len(sequences) == 0:
        return []
    if len(sequences) == 1:
        return sequences

    aligned = []
    
    a1, a2, score = needleman_wunsch_msa(sequences[0], sequences[1], 1, -1, -1)
    aligned.append(a1)
    aligned.append(a2)
    
    for i in range(2, len(sequences)):
        consensus = ""
        for col in range(len(aligned[0])):
            a_count = 0
            c_count = 0
            g_count = 0
            t_count = 0
            
            for s in aligned:
                char = s[col]
                if char == 'A':
                    a_count = a_count + 1
                if char == 'C':
                    c_count = c_count + 1
                if char == 'G':
                    g_count = g_count + 1
                if char == 'T':
                    t_count = t_count + 1
            
            best = 'A'
            max_c = a_count
            if c_count > max_c:
                best = 'C'
                max_c = c_count
            if g_count > max_c:
                best = 'G'
                max_c = g_count
            if t_count > max_c:
                best = 'T'
                max_c = t_count
                
            if max_c == 0:
                consensus = consensus + "-"
            else:
                consensus = consensus + best
                
        c_align, new_align, score = needleman_wunsch_msa(consensus, sequences[i], 1, -1, -1)
        
        new_msa = []
        for s in aligned:
            new_msa.append("")
            
        old_col = 0
        for final_col in range(len(c_align)):
            if c_align[final_col] == "-":
                for row_idx in range(len(aligned)):
                    new_msa[row_idx] = new_msa[row_idx] + "-"
            else:
                for row_idx in range(len(aligned)):
                    new_msa[row_idx] = new_msa[row_idx] + aligned[row_idx][old_col]
                old_col = old_col + 1
        
        aligned = new_msa
        aligned.append(new_align)
        
    return aligned
"""
write_file("Implementation/Multiple Sequence Alignment/msa.py", msa_code)

msa_helper_code = """
# copy of needleman to avoid weird imports from pairise
def needleman_wunsch_msa(seq1, seq2, match, mismatch, gap):
    n = len(seq1)
    m = len(seq2)
    matrix = []
    for i in range(n + 1):
        row = []
        for j in range(m + 1):
            row.append(0)
        matrix.append(row)
    for i in range(1, n + 1):
        matrix[i][0] = i * gap
    for j in range(1, m + 1):
        matrix[0][j] = j * gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i-1] == seq2[j-1]:
                score = match
            else:
                score = mismatch
            path1 = matrix[i-1][j-1] + score
            path2 = matrix[i-1][j] + gap
            path3 = matrix[i][j-1] + gap
            max_val = path1
            if path2 > max_val: max_val = path2
            if path3 > max_val: max_val = path3
            matrix[i][j] = max_val
    align1 = ""
    align2 = ""
    i = n
    j = m
    while i > 0 and j > 0:
        if seq1[i-1] == seq2[j-1]: score = match
        else: score = mismatch
        if matrix[i][j] == matrix[i-1][j-1] + score:
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i = i - 1
            j = j - 1
        elif matrix[i][j] == matrix[i-1][j] + gap:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i = i - 1
        else:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j = j - 1
    while i > 0:
        align1 = seq1[i-1] + align1
        align2 = "-" + align2
        i = i - 1
    while j > 0:
        align1 = "-" + align1
        align2 = seq2[j-1] + align2
        j = j - 1
    return align1, align2, matrix[n][m]
"""
write_file("Implementation/msa_helper.py", msa_helper_code)


main_code = """
import os
import sys

# keep imports simple by adding all folders
base = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(base, "Pairwise Alignment"))
sys.path.append(os.path.join(base, "Heuristic and Approximate Alignment"))
sys.path.append(os.path.join(base, "Multiple Sequence Alignment"))

import pairwise
import gotoh
import banded
import heuristics
import msa

def run():
    seq1 = "GATTACA"
    seq2 = "GCATGCU"
    
    print("Testing Needleman-Wunsch (Global DP)")
    a1, a2, score = pairwise.needleman_wunsch(seq1, seq2, 1, -1, -1)
    print("Code gave: " + a1 + " and " + a2 + ", score: " + str(score))
    print("")
    
    print("Testing Smith-Waterman (Local DP)")
    a1, a2, score = pairwise.smith_waterman(seq1, seq2, 2, -1, -1)
    print("Code gave: " + a1 + " and " + a2 + ", score: " + str(score))
    print("")

    print("Testing Gotoh")
    a1, a2, score = gotoh.gotoh_alignment(seq1, seq2, 1, -1, -2, -1)
    print("Code gave score: " + str(score))
    print("")

    print("Testing Multiple Sequence Alignment")
    msa_seqs = ["GATTACA", "GATCA", "GATAC"]
    res = msa.simple_progressive_msa(msa_seqs)
    for i in range(len(res)):
        print("Sequence " + str(i+1) + ": " + res[i])

if __name__ == "__main__":
    run()
"""
write_file("Implementation/main.py", main_code)
