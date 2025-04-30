# compute_assembly_index.py

import sqlite3
import argparse
import math
import lzma
import networkx as nx
from collections import Counter

# --- Helper Functions ---

def fetch_motifs(conn, species):
    query = """
        SELECT m.consensus FROM motifs m
        WHERE m.species = ?
    """
    cur = conn.cursor()
    cur.execute(query, (species,))
    return [row[0] for row in cur.fetchall()]

# --- Method 1: Entropy-Based ---
def compute_entropy_ai(motifs):
    counts = Counter(motifs)
    total = sum(counts.values())
    probs = [v / total for v in counts.values()]
    entropy = -sum(p * math.log2(p) for p in probs if p > 0)
    return entropy / math.log2(len(motifs)) if motifs else 0

# --- Method 2: Compression-Based ---
def compute_compression_ai(motifs):
    concat = "".join(motifs).encode("utf-8")
    compressed = lzma.compress(concat)
    return len(compressed) / len(concat) if concat else 0

# --- Method 3: Diversity-Based ---
def compute_diversity_ai(motifs):
    return len(set(motifs)) / len(motifs) if motifs else 0

# --- Method 4: Graph Path AI (using substring graph) ---
def build_substring_graph(motifs):
    G = nx.DiGraph()
    for i, m1 in enumerate(motifs):
        G.add_node(i, seq=m1)
    for i, m1 in enumerate(motifs):
        for j, m2 in enumerate(motifs):
            if i != j and m1 in m2:
                G.add_edge(i, j)
    return G

def compute_path_ai(motifs):
    G = build_substring_graph(motifs)
    lengths = []
    for node in G.nodes:
        try:
            path_lengths = nx.single_source_shortest_path_length(G, node)
            lengths.append(min(path_lengths.values()))
        except:
            lengths.append(0)
    return sum(lengths) / len(lengths) if lengths else 0

# --- CLI Interface ---

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", required=True, help="Path to motifs.db")
    parser.add_argument("--species", required=True, help="Species name (as in species.organism)")
    parser.add_argument("--method", required=True, choices=["entropy", "compression", "diversity", "path"], help="AI computation method")
    args = parser.parse_args()

    conn = sqlite3.connect(args.db)
    motifs = fetch_motifs(conn, args.species)

    if args.method == "entropy":
        ai = compute_entropy_ai(motifs)
    elif args.method == "compression":
        ai = compute_compression_ai(motifs)
    elif args.method == "diversity":
        ai = compute_diversity_ai(motifs)
    elif args.method == "path":
        ai = compute_path_ai(motifs)

    print(f"{args.species}\t{args.method}_ai\t{ai:.4f}")

if __name__ == "__main__":
    main()
