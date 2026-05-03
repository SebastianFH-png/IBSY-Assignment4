import os
import time
import tracemalloc

from SequenceReader import SequenceReader
from Config import Config
from NeedlemanWunsch import NeedlemanWunsch
from SmithWaterman import SmithWaterman
from Gotoh import Gotoh
from BasicSeedAndExtend import BasicSeedAndExtend

class ExperimentRunner:
    """
    Handles Task B from the assignment:
    Experiments - Compare accuracy, runtime, memory, and parameter sensitivity.
    """
    def __init__(self):
        base_dir = os.path.dirname(os.path.abspath(__file__))
        self.datasets_dir = os.path.join(base_dir, "..", "Datasets")
        self.reader = SequenceReader()

    def run_experiments(self):
        print("=== Starting Experiments (Task B) ===\n")
        
        # We will test on 3 datasets as requested
        datasets = [
            "Short related proteins.fasta",
            "Distant homologs.fasta",
            "Escherichia coli str. K-12 substr. MG1655, complete genome.fasta"
        ]

        for ds in datasets:
            file_path = os.path.join(self.datasets_dir, ds)
            print(f"--- Loading Dataset: {ds} ---")
            sequences = self.reader.read_fasta(file_path)
            
            if len(sequences) < 2 and ds != datasets[2]:
                print(f"Not enough sequences found in {ds}\n")
                continue
            
            # Setup sequence 1 and sequence 2 for comparison
            # For the genome (single long sequence), we will just cut two slices from it to simulate comparison
            if len(sequences) == 1:
                seq1_str = sequences[0]["sequence"][0:200]
                seq2_str = sequences[0]["sequence"][500:700]
                print("Single sequence dataset. Using slices from the genome.")
            else:
                # We limit to 200 characters to prevent Memory Errors (O(N^2) complexity with standard DP lists)
                seq1_str = sequences[0]["sequence"][:200]
                seq2_str = sequences[1]["sequence"][:200]
            
            # 1. Benchmark: Runtime & Memory
            self.benchmark_algorithms(seq1_str, seq2_str)
            
        # 2. Parameter Sensitivity (Using the first dataset)
        print("\n=== Parameter Sensitivity Test ===")
        print("Testing how different match/mismatch/gap scores affect Needleman-Wunsch.")
        short_prots_path = os.path.join(self.datasets_dir, "Short related proteins.fasta")
        seqs = self.reader.read_fasta(short_prots_path)
        if len(seqs) >= 2:
            s1 = seqs[0]["sequence"][:100]
            s2 = seqs[1]["sequence"][:100]
            self.test_parameter_sensitivity(s1, s2)

    def benchmark_algorithms(self, seq1, seq2):
        """
        Calculates time taken and memory used for the 4 algorithms.
        """
        algorithms = [
            ("Needleman-Wunsch (Global)", NeedlemanWunsch()),
            ("Smith-Waterman (Local)", SmithWaterman()),
            ("Gotoh (Affine Gap)", Gotoh()),
            ("Basic Seed-and-Extend (Heuristic)", BasicSeedAndExtend(kmer_size=4))
        ]

        # Table header
        print(f"{'Algorithm':<35} | {'Score':<6} | {'Time (ms)':<10} | {'Peak Memory (KB)':<15}")
        print("-" * 75)

        for name, algo in algorithms:
            # Start memory tracking
            tracemalloc.start()
            
            # Start timer
            start_time = time.time()
            
            # Run the alignment
            score, align1, align2 = algo.align(seq1, seq2)
            
            # Stop timer
            end_time = time.time()
            
            # Get memory usage
            current, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            
            time_ms = (end_time - start_time) * 1000
            peak_kb = peak / 1024
            
            print(f"{name:<35} | {score:<6} | {time_ms:<10.2f} | {peak_kb:<15.2f}")
        print()

    def test_parameter_sensitivity(self, seq1, seq2):
        """
        Changes configuration statically to see how algorithms respond.
        """
        # Save original values
        orig_match = Config.MATCH_SCORE
        orig_mismatch = Config.MISMATCH_SCORE
        orig_gap = Config.GAP_PENALTY

        # Parameter Set 1: High gap penalty (forces fewer gaps)
        Config.MATCH_SCORE = 2
        Config.MISMATCH_SCORE = -1
        Config.GAP_PENALTY = -5
        nw1 = NeedlemanWunsch()
        score1, a1_1, a2_1 = nw1.align(seq1, seq2)
        print(f"Set 1 (Match: 2, Mis: -1, Gap: -5) -> Score: {score1}")

        # Parameter Set 2: Free gap penalty (allows lots of gaps)
        Config.MATCH_SCORE = 2
        Config.MISMATCH_SCORE = -1
        Config.GAP_PENALTY = 0
        nw2 = NeedlemanWunsch()
        score2, a1_2, a2_2 = nw2.align(seq1, seq2)
        print(f"Set 2 (Match: 2, Mis: -1, Gap:  0) -> Score: {score2}")
        
        # Parameter Set 3: High mismatch penalty
        Config.MATCH_SCORE = 1
        Config.MISMATCH_SCORE = -5
        Config.GAP_PENALTY = -1
        nw3 = NeedlemanWunsch()
        score3, a1_3, a2_3 = nw3.align(seq1, seq2)
        print(f"Set 3 (Match: 1, Mis: -5, Gap: -1) -> Score: {score3}")

        # Restore original values
        Config.MATCH_SCORE = orig_match
        Config.MISMATCH_SCORE = orig_mismatch
        Config.GAP_PENALTY = orig_gap
