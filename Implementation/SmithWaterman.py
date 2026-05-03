from Config import Config

class SmithWaterman:
    """
    Implements Local Alignment algorithm.
    Focuses on finding the highest-scoring matching region between sequences.
    """
    def __init__(self):
        self.match = Config.MATCH_SCORE
        self.mismatch = Config.MISMATCH_SCORE
        self.gap = Config.GAP_PENALTY

    def align(self, seq1, seq2):
        rows = len(seq1) + 1
        cols = len(seq2) + 1

        # 1. Initialize matrix with zeros
        dp = []
        for i in range(rows):
            row = []
            for j in range(cols):
                row.append(0)
            dp.append(row)

        max_score = 0
        max_i = 0
        max_j = 0

        # 2. Fill the matrix
        # In Smith-Waterman, negative scores become 0. We start storing from (1,1)
        for i in range(1, rows):
            for j in range(1, cols):
                if seq1[i-1] == seq2[j-1]:
                    score = self.match
                else:
                    score = self.mismatch
                
                diag = dp[i-1][j-1] + score
                up = dp[i-1][j] + self.gap
                left = dp[i][j-1] + self.gap
                
                # Any score going below 0 is "reset" to 0
                dp[i][j] = max(0, diag, up, left)
                
                # Keep track of the highest score anywhere in the matrix
                if dp[i][j] > max_score:
                    max_score = dp[i][j]
                    max_i = i
                    max_j = j

        # 3. Traceback (starts from the cell with the highest score)
        align1 = ""
        align2 = ""
        i = max_i
        j = max_j

        # We stop as soon as we hit a 0 score
        while i > 0 and j > 0 and dp[i][j] > 0:
            current_score = dp[i][j]
            
            if seq1[i-1] == seq2[j-1]:
                score = self.match
            else:
                score = self.mismatch

            if current_score == dp[i-1][j-1] + score:
                align1 = seq1[i-1] + align1
                align2 = seq2[j-1] + align2
                i -= 1
                j -= 1
            elif current_score == dp[i-1][j] + self.gap:
                align1 = seq1[i-1] + align1
                align2 = "-" + align2
                i -= 1
            else:
                align1 = "-" + align1
                align2 = seq2[j-1] + align2
                j -= 1

        return max_score, align1, align2
