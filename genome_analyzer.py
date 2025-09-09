from Bio import SeqIO
import pandas as pd
import os
import subprocess
import json
from typing import Dict, List, Any, Optional
import numpy as np
import itertools


def compute_stats(fasta_file: str) -> Dict[str, Any]:
    """Compute genome statistics from a FASTA file.
    
    Args:
        fasta_file: Path to the FASTA file
        
    Returns:
        Dictionary containing genome statistics
    """
    lengths = []
    gc_count = 0
    total_bases = 0
    num_contigs = 0

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = record.seq.upper()
        seq_len = len(seq)
        lengths.append(seq_len)
        gc_count += seq.count("G") + seq.count("C")
        total_bases += seq_len
        num_contigs += 1

    if not lengths:
        return {"Total length": 0, "Num contigs": 0, "N50": 0, "L90": 0, "GC%": 0}

    lengths.sort(reverse=True)
    gc_percent = (gc_count / total_bases) * 100 if total_bases > 0 else 0

    # N50
    half_length = total_bases / 2
    cumulative = 0
    n50 = 0
    for length in lengths:
        cumulative += length
        if cumulative >= half_length:
            n50 = length
            break

    # L90
    l90_target = total_bases * 0.9
    cumulative = 0
    l90 = 0
    for length in lengths:
        cumulative += length
        l90 += 1
        if cumulative >= l90_target:
            break

    return {
        "Total length": total_bases,
        "Num contigs": num_contigs,
        "N50": n50,
        "L90": l90,
        "GC%": round(gc_percent, 2)
    }


def run_busco(fasta_file: str, output_dir: str = "./temp/busco_results/", lineage: str = "bacteria_odb10") -> Optional[Dict[str, Any]]:
    """Run BUSCO analysis on a FASTA file (Docker-optimized).
    
    Args:
        fasta_file: Path to the FASTA file
        output_dir: Directory to store BUSCO results
        lineage: BUSCO lineage database to use (e.g., "bacteria_odb10", "eukaryota_odb10")
        
    Returns:
        Dictionary containing BUSCO results or None if BUSCO is not available
    """
    try:
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Get filename without extension for BUSCO output
        base_name = os.path.splitext(os.path.basename(fasta_file))[0]
        busco_name = f"busco_{base_name}"
        busco_output_dir = os.path.join(output_dir, busco_name)

        # Run BUSCO command (Docker-optimized)
        cmd = [
            "busco",
            "-i", fasta_file,
            "-o", busco_name,
            "-m", "genome",
            "-f",  # Force overwrite if output exists
            "-l", lineage,
            "--out_path", output_dir,
            "--cpu", "1"  # Limit CPU usage for container
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600
        )
        
        if result.returncode != 0:
            print(f"BUSCO failed for {fasta_file}: {result.stderr}")
            return None

        # Parse BUSCO results
        return parse_busco_results(busco_output_dir, lineage, busco_name)
        
    except subprocess.TimeoutExpired:
        print(f"BUSCO timed out for {fasta_file}")
        return None
    except FileNotFoundError:
        print("BUSCO not found. Please install BUSCO first.")
        return None
    except ModuleNotFoundError:
        print("BUSCO dependencies missing. Please install BUSCO with: conda install -c bioconda busco")
        return None
    except Exception as e:
        print(f"Error running BUSCO: {e}")
        return None


def parse_busco_results(busco_dir: str, lineage: str = "bacteria_odb10", busco_name: str = 'sample') -> Optional[Dict[str, Any]]:
    
    """Parse BUSCO results from the output directory.
    
    Args:
        busco_dir: Directory containing BUSCO results
        lineage: BUSCO lineage database used (for finding the correct summary file)
        
    Returns:
        Dictionary containing parsed BUSCO results
    """
    try:
        # Look for the short summary file
        summary_file = os.path.join(busco_dir, "short_summary.txt")
        if not os.path.exists(summary_file):
            summary_file = os.path.join(busco_dir, f"short_summary.specific.{lineage}.{busco_name}.txt")
        
        if not os.path.exists(summary_file):
            return None
            
        with open(summary_file, 'r') as f:
            lines = f.readlines()
        
        # Parse the summary
        results = {}
        for line in lines:
            line = line.strip()
            if "Complete" in line and "single-copy" in line:
                results["Complete_Single"] = int(line.split()[0])
            elif "Complete" in line and "duplicated" in line:
                results["Complete_Duplicated"] = int(line.split()[0])
            elif "Fragmented" in line:
                results["Fragmented"] = int(line.split()[0])
            elif "Missing" in line:
                results["Missing"] = int(line.split()[0])
            elif "Total" in line and "BUSCO" in line:
                results["Total"] = int(line.split()[0])
        
        # Calculate percentages
        if results.get("Total", 0) > 0:
            results["Complete_Single_Percent"] = round((results.get("Complete_Single", 0) / results["Total"]) * 100, 2)
            results["Complete_Duplicated_Percent"] = round((results.get("Complete_Duplicated", 0) / results["Total"]) * 100, 2)
            results["Fragmented_Percent"] = round((results.get("Fragmented", 0) / results["Total"]) * 100, 2)
            results["Missing_Percent"] = round((results.get("Missing", 0) / results["Total"]) * 100, 2)
            results["Complete_Percent"] = round(((results.get("Complete_Single", 0) + results.get("Complete_Duplicated", 0)) / results["Total"]) * 100, 2)
        
        return results
        
    except Exception as e:
        print(f"Error parsing BUSCO results: {e}")
        return None


def compute_stats_with_busco(fasta_file: str, lineage: str = "bacteria_odb10") -> Dict[str, Any]:
    """Compute genome statistics including BUSCO analysis.
    
    Args:
        fasta_file: Path to the FASTA file
        lineage: BUSCO lineage database to use
        
    Returns:
        Dictionary containing genome statistics and BUSCO results
    """
    # Get basic stats
    basic_stats = compute_stats(fasta_file)
    
    # Add BUSCO results
    busco_results = run_busco(fasta_file, lineage=lineage)
    if busco_results:
        basic_stats.update({
            "BUSCO_Complete": busco_results.get("Complete_Percent", 0),
            "BUSCO_Single": busco_results.get("Complete_Single_Percent", 0),
            "BUSCO_Duplicated": busco_results.get("Complete_Duplicated_Percent", 0),
            "BUSCO_Fragmented": busco_results.get("Fragmented_Percent", 0),
            "BUSCO_Missing": busco_results.get("Missing_Percent", 0),
            "BUSCO_Total": busco_results.get("Total", 0)
        })
    else:
        basic_stats.update({
            "BUSCO_Complete": None,
            "BUSCO_Single": None,
            "BUSCO_Duplicated": None,
            "BUSCO_Fragmented": None,
            "BUSCO_Missing": None,
            "BUSCO_Total": None
        })
    
    return basic_stats


def process_directory(directory: str, include_busco: bool = True, busco_lineage: str = "bacteria_odb10") -> pd.DataFrame:
    """Process all FASTA files in a directory and return statistics as a DataFrame.
    
    Args:
        directory: Path to directory containing FASTA files
        include_busco: Whether to include BUSCO analysis (default: True)
        busco_lineage: BUSCO lineage database to use (default: "bacteria_odb10")
        
    Returns:
        DataFrame with genome statistics for each file
    """
    results = []
    for filename in os.listdir(directory):
        if filename.endswith((".fasta", ".fa", ".fna")):
            filepath = os.path.join(directory, filename)
            if include_busco:
                stats = compute_stats_with_busco(filepath, lineage=busco_lineage)
            else:
                stats = compute_stats(filepath)
            stats["File"] = str(filename.rsplit('.', 1)[0])
            results.append(stats)
    return pd.DataFrame(results)


def run_fastani(query_file, ref_file, output_file, threads=1):
    """
    Run fastANI on two genome files.
    Args:
        query_file: Path to the query genome (FASTA)
        ref_file: Path to the reference genome (FASTA)
        output_file: Path to write fastANI results
        threads: Number of threads to use
    Returns:
        ANI value (float) if successful, None otherwise
    """
    cmd = [
        "fastANI",
        "-q", query_file,
        "-r", ref_file,
        "-o", output_file,
        "-t", str(threads)
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # fastANI output: query, ref, ANI, fragments, matches
        with open(output_file) as f:
            line = f.readline()
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                return float(parts[2])
        return None
    except Exception as e:
        print(f"fastANI failed: {e}")
        return None


def all_vs_all_fastani(file_paths, threads=1):
    """
    Run all-vs-all fastANI comparisons and return ANI matrix.
    Args:
        file_paths: List of genome file paths
        threads: Number of threads
    Returns:
        ANI matrix (numpy array), genome names (list)
    """
    n = len(file_paths)
    ani_matrix = np.zeros((n, n))
    genome_names = [os.path.basename(fp) for fp in file_paths]
    for i, j in itertools.combinations(range(n), 2):
        out_file = f"./temp/fastani_{i}_{j}.txt"
        ani = run_fastani(file_paths[i], file_paths[j], out_file, threads=threads)
        if ani is not None:
            ani_matrix[i, j] = ani
            ani_matrix[j, i] = ani
        else:
            ani_matrix[i, j] = np.nan
            ani_matrix[j, i] = np.nan
    np.fill_diagonal(ani_matrix, 100.0)
    return ani_matrix, genome_names


def compute_distance_matrix(ani_matrix):
    """
    Convert ANI matrix to a distance matrix (1 - ANI/100).
    """
    return 1 - (ani_matrix / 100.0)


def neighbor_joining_tree(distance_matrix, genome_names):
    """
    Build a neighbor-joining tree from a distance matrix.
    Returns a Newick string.
    """
    try:
        from skbio import DistanceMatrix
        from skbio.tree import nj
        dm = DistanceMatrix(distance_matrix, genome_names)
        tree = nj(dm)
        return tree.ascii_art(), str(tree)
    except ImportError:
        return "scikit-bio not installed", ""
   