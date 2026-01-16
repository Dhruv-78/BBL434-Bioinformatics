import sys
from collections import Counter
import matplotlib.pyplot as plt

# -----------------------------
# Parameters
# -----------------------------
K = 8
WINDOW_SIZE = 5000
STEP_SIZE = 500

# -----------------------------
# FASTA reader
# -----------------------------
def read_fasta(filename):
    sequence = []
    with open(filename) as f:
        for line in f:
            if not line.startswith(">"):
                sequence.append(line.strip().upper())
    return "".join(sequence)

# -----------------------------
# k-mer counter
# -----------------------------
def count_kmers(seq, k):
    return Counter(seq[i:i+k] for i in range(len(seq) - k + 1))

# -----------------------------
# Main
# -----------------------------
if len(sys.argv) != 2:
    print("Usage: python script1.py genomic.fa")
    sys.exit(1)

fasta_file = sys.argv[1]

# Read genome
genome = read_fasta(fasta_file)
genome_len = len(genome)

# Global k-mer frequencies
global_kmers = count_kmers(genome, K)
total_global = sum(global_kmers.values())

global_freq = {
    kmer: count / total_global
    for kmer, count in global_kmers.items()
}

positions = []
enrichment_scores = []

# Sliding window analysis
for start in range(0, genome_len - WINDOW_SIZE + 1, STEP_SIZE):
    end = start + WINDOW_SIZE
    window_seq = genome[start:end]

    window_kmers = count_kmers(window_seq, K)
    total_window = sum(window_kmers.values())

    enrichment = 0
    for kmer, obs in window_kmers.items():
        expected = global_freq.get(kmer, 0) * total_window
        if expected > 0:
            enrichment += obs / expected

    enrichment /= len(window_kmers)

    positions.append(start + WINDOW_SIZE // 2)
    enrichment_scores.append(enrichment)

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(12, 5))
plt.plot(positions, enrichment_scores)
plt.xlabel("Genomic position (bp)")
plt.ylabel("Average k-mer enrichment")
plt.title("ORI signal detection using k-mer enrichment (k=8)")
plt.grid(True)
plt.tight_layout()
plt.show()

