from pathlib import Path
import logging

class deepTMHMMRunner:
    def __init__(self,
                 fasta_path
        ):
        self.fasta_path = fasta_path
    
    def test(self):
        command = f"biolib run DTU/DeepTMHMM --fasta {self.fasta_path}"

def main():
    deeptmhmm_runner = deepTMHMMRunner(fasta_path = "data/processed/test.fa")
    deeptmhmm_runner.test()

if __name__ == "__main__":
    main()
                