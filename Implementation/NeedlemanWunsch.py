from Config import Config

class NeedlemanWunsch:
    """
    Implements the Global Alignment algorithm.
    It aligns sequences end-to-end.
    """
    def __init__(self):
        # We load parameters from our Config class
        self.match = Config.MATCH_SCORE
        self.mismatch = Config.MISMATCH_SCORE
        self.gap = Config.GAP_PENALTY

    def align(self, seq1, seq2):
        rows = len(seq1) + 1
        cols = len(seq2) + 1

        # 1. Initialize matrices (2D arrays in Java)
        # We use nested standard lists because Python doesn't have native 2D arrays built-in
        dp = []
        for i in range(rows):
            row = []
            for j in range(cols):
                row.append(0)
            dp.append(row)

        # 2. Fill the first row and column
        for i in range(rows):
            dp[i][0] = i * self.gap
        for j in range(cols):
            dp[0][j] = j * self.gap

        # 3. Fill the rest of the DP matrix
        for i in range(1, rows):
            for j in range(1, cols):
                # Check for match or mismatch
                if seq1[i-1] == seq2[j-1]:
                    score = self.match
                else:
                    score = self.mismatch
                
                # Calculate the 3 possible paths: diagonal, up, and left
                diag_score = dp[i-1][j-1] + score
                up_score = dp[i-1][j] + self.gap
                left_score = dp[i][j-1] + self.gap
                
                # Take the maximum of all three
                dp[i][j] = max(diag_score, up_score, left_score)

        # 4. Traceback to find the alignment string
        align1 = ""
        align2 = ""
        i = len(seq1)
        j = len(seq2)

        while i > 0 and j > 0:
            current_score = dp[i][j]
            
            if seq1[i-1] == seq2[j-1]:
                score = self.match
            else:
                score = self.mismatch

            # Follow the path back optimally
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

        # Fill any remaining characters
        while i > 0:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1

        while j > 0:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1

        final_score = dp[len(seq1)][len(seq2)]
        return final_score, align1, align2
