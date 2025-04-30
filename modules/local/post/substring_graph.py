# motif_substring_graph.py

import sqlite3
import argparse
import networkx as nx
import os

# --- Load motifs from database ---
def fetch_motifs(conn, species):
    query = """
        SELECT m.motif_id, m.consensus, s.organism
        FROM motifs m
        JOIN species s ON m.species = s.organism
        WHERE s.organism = ?
    """
    cur = conn.cursor()
    cur.execute(query, (species,))
    return cur.fetchall()  # list of tuples (motif_id, consensus, organism)

# --- Build substring graph ---
def build_substring_graph(motifs):
    G = nx.DiGraph()
    for motif_id, seq, species in motifs:
        G.add_node(motif_id, consensus=seq, species=species)

    for i, (id1, s1, _) in enumerate(motifs):
        for j, (id2, s2, _) in enumerate(motifs):
            if i != j and s1 in s2:
                G.add_edge(id1, id2, relationship="substring")

    return G

# --- Export graph ---
def export_graph(G, out_path, fmt):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    if fmt == "gml":
        nx.write_gml(G, out_path)
    elif fmt == "graphml":
        nx.write_graphml(G, out_path)
    elif fmt == "json":
        from networkx.readwrite import json_graph
        with open(out_path, "w") as f:
            json.dump(json_graph.node_link_data(G), f, indent=2)
    elif fmt == "pickle":
        nx.write_gpickle(G, out_path)
    else:
        raise ValueError(f"Unsupported format: {fmt}")

# --- CLI ---
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", required=True, help="Path to motifs.db")
    parser.add_argument("--species", required=True, help="Species name")
    parser.add_argument("--output", required=True, help="Output file (e.g., graphs/axolotl.graphml)")
    parser.add_argument("--format", default="graphml", choices=["gml", "graphml", "json", "pickle"], help="Output graph format")
    args = parser.parse_args()

    conn = sqlite3.connect(args.db)
    motifs = fetch_motifs(conn, args.species)

    print(f"Building substring graph for {args.species} with {len(motifs)} motifs...")
    G = build_substring_graph(motifs)

    print(f"Saving graph to {args.output} as {args.format}...")
    export_graph(G, args.output, args.format)
    print("Done.")

if __name__ == "__main__":
    main()
