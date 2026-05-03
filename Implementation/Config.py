class Config:
    """
    This class holds the configuration and scoring parameters for our alignment algorithms.
    Just like a static fields class in Java.
    """
    MATCH_SCORE = 2
    MISMATCH_SCORE = -1
    GAP_PENALTY = -2

    # For affine gap penalty algorithms (Gotoh)
    GAP_OPEN = -3
    GAP_EXTEND = -1
