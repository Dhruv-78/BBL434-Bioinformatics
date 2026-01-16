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
    seq = []
    with open(filename) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)

# -----------------------------
# k-mer counter
# -----------------------------
def count_kmers(seq, k):
    return Counter(seq[i:i+k] for i in range(len(seq) - k + 1))

# -----------------------------
# Main
# -----------------------------
if len(sys.argv) != 2:
    print("Usage: python script4.py genomic.fa")
    sys.exit(1)

genome = read_fasta(sys.argv[1])
n = len(genome)

# ---- Global k-mer frequencies ----
global_kmers = count_kmers(genome, K)
total_global = sum(global_kmers.values())

global_freq = {
    kmer: count / total_global
    for kmer, count in global_kmers.items()
}

positions = []
enrichment_scores = []
gc_skews = []
cumulative_gc_skew = []

cum_gc = 0

# ---- Sliding window analysis ----
for start in range(0, n - WINDOW_SIZE + 1, STEP_SIZE):
    end = start + WINDOW_SIZE
    window = genome[start:end]

    # ---- k-mer enrichment ----
    window_kmers = count_kmers(window, K)
    total_window = sum(window_kmers.values())

    enrichment = 0
    for kmer, obs in window_kmers.items():
        expected = global_freq.get(kmer, 0) * total_window
        if expected > 0:
            enrichment += obs / expected

    enrichment /= len(window_kmers)

    # ---- GC skew ----
    G = window.count("G")
    C = window.count("C")
    skew = (G - C) / (G + C) if (G + C) > 0 else 0
    cum_gc += skew

    # ---- Store ----
    positions.append(start + WINDOW_SIZE // 2)
    enrichment_scores.append(enrichment)
    gc_skews.append(skew)
    cumulative_gc_skew.append(cum_gc)

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(12, 8))

plt.subplot(2, 1, 1)
plt.plot(positions, enrichment_scores)
plt.ylabel("k-mer enrichment")
plt.title("ORI Detection using k-mer Enrichment and GC Skew")
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(positions, cumulative_gc_skew)
plt.xlabel("Genomic position (bp)")
plt.ylabel("Cumulative GC skew")
plt.grid(True)

plt.tight_layout()
plt.show()

# -----------------------------
# ORI prediction (heuristic)
# -----------------------------
ori_index = cumulative_gc_skew.index(min(cumulative_gc_skew))
ori_position = positions[ori_index]

print(f"\nPredicted ORI position ≈ {ori_position} bp")
