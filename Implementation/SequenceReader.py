import os

class SequenceReader:
    """
    This class handles the parsing of FASTA files.
    """
    def __init__(self):
        pass

    def read_fasta(self, file_path):
        """
        Reads a FASTA file from the given file_path.
        Returns a list of dictionaries containing 'header' and 'sequence'.
        """
        sequences = []
        current_sequence = ""
        current_header = ""

        # Check if file exists, just like File.exists() in Java
        if not os.path.exists(file_path):
            print("Error: File not found at " + file_path)
            return sequences

        # Open the file for reading ('r')
        with open(file_path, 'r') as file:
            for line in file:
                # Remove newline characters at the end of the line
                line = line.strip()
                
                if line.startswith(">"):
                    # If we already have a header, save the previous sequence
                    if current_header != "":
                        sequences.append({"header": current_header, "sequence": current_sequence})
                    
                    # Start a new sequence
                    current_header = line[1:] # Remove the '>' character
                    current_sequence = ""
                else:
                    # Append the string line to the current sequence
                    current_sequence += line

        # Add the very last sequence in the file to our list
        if current_header != "":
            sequences.append({"header": current_header, "sequence": current_sequence})

        return sequences
