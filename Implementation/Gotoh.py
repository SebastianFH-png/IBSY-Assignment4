from Config import Config

class Gotoh:
    """
    Implements the Gotoh Algorithm for global alignment with
    Affine Gap Penalties.
    It uses 3 matrices (M, X, Y) instead of 1.
    """
    def __init__(self):
        self.match = Config.MATCH_SCORE
        self.mismatch = Config.MISMATCH_SCORE
        self.gap_open = Config.GAP_OPEN
        self.gap_extend = Config.GAP_EXTEND

    def align(self, seq1, seq2):
        rows = len(seq1) + 1
        cols = len(seq2) + 1
        INF = float('inf')

        # Create three matrices (M: match, X: gap in seq1, Y: gap in seq2)
        M = [[-INF] * cols for _ in range(rows)]
        X = [[-INF] * cols for _ in range(rows)]
        Y = [[-INF] * cols for _ in range(rows)]

        M[0][0] = 0

        # Initialization
        for i in range(1, rows):
            X[i][0] = self.gap_open + (i - 1) * self.gap_extend
            M[i][0] = -INF # Invalid
            Y[i][0] = -INF

        for j in range(1, cols):
            Y[0][j] = self.gap_open + (j - 1) * self.gap_extend
            M[0][j] = -INF
            X[0][j] = -INF

        # Matrix filling
        for i in range(1, rows):
            for j in range(1, cols):
                # Penalty for gap in seq1 (X)
                X[i][j] = max(M[i-1][j] + self.gap_open, X[i-1][j] + self.gap_extend)
                
                # Penalty for gap in seq2 (Y)
                Y[i][j] = max(M[i][j-1] + self.gap_open, Y[i][j-1] + self.gap_extend)
                
                # Match / Mismatch
                score = self.match if seq1[i-1] == seq2[j-1] else self.mismatch
                M[i][j] = max(M[i-1][j-1], X[i-1][j-1], Y[i-1][j-1]) + score

        final_score = max(M[rows-1][cols-1], X[rows-1][cols-1], Y[rows-1][cols-1])
        # Note: A full traceback for Gotoh is lengthy, so we just return the score for brevity
        return final_score, "Traceback not fully implemented", "Traceback not fully implemented"
