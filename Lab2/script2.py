import sys
from collections import defaultdict

# -----------------------------
# Parameters
# -----------------------------
K = 8
T = 3
L = 1000

# -----------------------------
# FASTA reader
# -----------------------------
def read_fasta(filename):
    seq = []
    with open(filename) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)

# -----------------------------
# Main
# -----------------------------
if len(sys.argv) != 2:
    print("Usage: python script2.py genomic.fa")
    sys.exit(1)

fasta_file = sys.argv[1]
genome = read_fasta(fasta_file)
n = len(genome)

clump_kmers = set()

# Slide window of length L
for start in range(0, n - L + 1):
    window = genome[start:start + L]

    kmer_count = defaultdict(int)

    # Count k-mers in window
    for i in range(L - K + 1):
        kmer = window[i:i + K]
        kmer_count[kmer] += 1

        # Early stopping for efficiency
        if kmer_count[kmer] >= T:
            clump_kmers.add(kmer)

# -----------------------------
# Output
# -----------------------------
print(f"(L={L}, k={K}, t={T}) clumps found:\n")
for kmer in sorted(clump_kmers):
    print(kmer)

print(f"\nTotal clump-forming k-mers: {len(clump_kmers)}")

