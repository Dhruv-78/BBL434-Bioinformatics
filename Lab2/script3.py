import sys
import matplotlib.pyplot as plt

# -----------------------------
# Parameters
# -----------------------------
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
# Main
# -----------------------------
if len(sys.argv) != 2:
    print("Usage: python script3.py genomic.fa")
    sys.exit(1)

fasta_file = sys.argv[1]
genome = read_fasta(fasta_file)
n = len(genome)

positions = []
gc_skews = []
cumulative_gc_skew = []

cum_sum = 0

# Sliding window GC skew
for start in range(0, n - WINDOW_SIZE + 1, STEP_SIZE):
    end = start + WINDOW_SIZE
    window = genome[start:end]

    G = window.count("G")
    C = window.count("C")

    if G + C == 0:
        skew = 0
    else:
        skew = (G - C) / (G + C)

    cum_sum += skew

    positions.append(start + WINDOW_SIZE // 2)
    gc_skews.append(skew)
    cumulative_gc_skew.append(cum_sum)

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(12, 5))
plt.plot(positions, cumulative_gc_skew)
plt.xlabel("Genomic position (bp)")
plt.ylabel("Cumulative GC skew")
plt.title("Cumulative GC Skew Plot")
plt.grid(True)
plt.tight_layout()
plt.show()
