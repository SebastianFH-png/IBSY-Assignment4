import os
from SequenceReader import SequenceReader
from NeedlemanWunsch import NeedlemanWunsch
from SmithWaterman import SmithWaterman
from Gotoh import Gotoh
from BasicSeedAndExtend import BasicSeedAndExtend
from ExperimentRunner import ExperimentRunner

class Main:
    def __init__(self):
 
        base_dir = os.path.dirname(os.path.abspath(__file__)) # The Implementation folder
        self.datasets_dir = os.path.join(base_dir, "..", "Datasets") # Navigates up one layer to Datasets
        self.reader = SequenceReader()

    def run(self):
        print("--- Starting Sequence Alignment Assignment ---")

        # Dynamically build path to the fasta file
        short_proteins_file = os.path.join(self.datasets_dir, "Short related proteins.fasta")
        
        # 1. Read Sequences
        print("Reading file:", short_proteins_file)
        sequences = self.reader.read_fasta(short_proteins_file)
        
        if len(sequences) < 2:
            print("Not enough sequences to align!")
            return

        # Pick the first two sequences
        seq1 = sequences[0]
        seq2 = sequences[1]
        
        # To avoid massive console prints, we take a substring for demonstration
        s1_str = seq1["sequence"][:50]
        s2_str = seq2["sequence"][:50]

        print("\nAligning:")
        print("Seq1 (" + seq1["header"][:30] + "...): " + s1_str + "...")
        print("Seq2 (" + seq2["header"][:30] + "...): " + s2_str + "...")

        # 2. Run Needleman-Wunsch Global Alignment
        print("\n--- Algorithm 1: Needleman-Wunsch (Global) ---")
        nw = NeedlemanWunsch()
        score_nw, align1_nw, align2_nw = nw.align(s1_str, s2_str)
        print("Score: " + str(score_nw))
        print(align1_nw)
        print(align2_nw)

        # 3. Run Smith-Waterman Local Alignment
        print("\n--- Algorithm 2: Smith-Waterman (Local) ---")
        sw = SmithWaterman()
        score_sw, align1_sw, align2_sw = sw.align(s1_str, s2_str)
        print("Score: " + str(score_sw))
        print(align1_sw)
        # 4. Run Gotoh (Affine Gap)
        print("\n--- Algorithm 3: Gotoh (Affine Gap) ---")
        gotoh = Gotoh()
        score_go, align1_go, align2_go = gotoh.align(s1_str, s2_str)
        print("Score: " + str(score_go))

        # 5. Run Basic Seed and Extend (BLAST style)
        print("\n--- Algorithm 4: Basic Seed and Extend (BLAST-style Heuristic) ---")
        blast = BasicSeedAndExtend(kmer_size=4)
        score_bl, align1_bl, align2_bl = blast.align(s1_str, s2_str)
        print("Score: " + str(score_bl))
        print("Matching Region 1: " + align1_bl)
        print("Matching Region 2: " + align2_bl)
        
        print("\nAll 4 algorithms successfully demonstrated! (Task A)")

        print("\n" + "="*50)
        print("Beginning Benchmarking and Analysis (Task B)")
        print("="*50)
        runner = ExperimentRunner()
        runner.run_experiments()

# standard python main hook
if __name__ == "__main__":
    app = Main()
    app.run()
