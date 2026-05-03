class BasicSeedAndExtend:
    """
    Implements a simplified BLAST-style heuristic.
    1. Finds an exact matching k-mer (Seed).
    2. Tries to extend it linearly without gaps.
    """
    def __init__(self, kmer_size=3):
        self.kmer_size = kmer_size

    def align(self, seq1, seq2):
        # 1. Build a hash map / dictionary of k-mers from seq1
        kmers = {}
        for i in range(len(seq1) - self.kmer_size + 1):
            kmer = seq1[i:i + self.kmer_size]
            if kmer not in kmers:
                kmers[kmer] = []
            kmers[kmer].append(i)  # Store indices of all occurrences

        best_score = 0
        best_align1 = ""
        best_align2 = ""

        # 2. Iterate through seq2, finding seeds
        for j in range(len(seq2) - self.kmer_size + 1):
            kmer = seq2[j:j + self.kmer_size]
            
            # If we find a seed match
            if kmer in kmers:
                for i in kmers[kmer]:
                    # 3. Extend the match left and right without gaps
                    score, a1, a2 = self.extend(seq1, seq2, i, j)
                    if score > best_score:
                        best_score = score
                        best_align1 = a1
                        best_align2 = a2

        return best_score, best_align1, best_align2

    def extend(self, seq1, seq2, start1, start2):
        # Extend to the right
        i = start1 + self.kmer_size
        j = start2 + self.kmer_size
        score = self.kmer_size  # Starting score for the exact seed match
        
        while i < len(seq1) and j < len(seq2) and seq1[i] == seq2[j]:
            score += 1
            i += 1
            j += 1
            
        # Extend to the left
        i_left = start1 - 1
        j_left = start2 - 1
        while i_left >= 0 and j_left >= 0 and seq1[i_left] == seq2[j_left]:
            score += 1
            i_left -= 1
            j_left -= 1
            
        align1 = seq1[i_left+1:i]
        align2 = seq2[j_left+1:j]

        return score, align1, align2