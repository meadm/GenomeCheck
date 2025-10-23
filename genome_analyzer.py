from Bio import SeqIO
import pandas as pd
import os
import subprocess
import shutil
import json
from typing import Dict, List, Any, Optional, Tuple
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
        return {"Total length (bp)": 0, "Num contigs": 0, "N50": 0, "L90": 0, "GC%": 0}

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
        "Total length (bp)": total_bases,
        "Num contigs": num_contigs,
        "N50": n50,
        "L90": l90,
        "GC%": round(gc_percent, 2)
    }


def run_busco(fasta_file: str, output_dir: str = None, lineage: str = "bacteria_odb10", cpus: int = 1) -> Optional[Dict[str, Any]]:
    """Run BUSCO analysis on a FASTA file (Docker-optimized).
    
    Args:
        fasta_file: Path to the FASTA file
        output_dir: Directory to store BUSCO results
        lineage: BUSCO lineage database to use (e.g., "bacteria_odb10", "eukaryota_odb10")
        
    Returns:
        Dictionary containing BUSCO results or None if BUSCO is not available
    """
    try:
        # Determine session-specific output directory
        if output_dir is None:
            base_tmp = os.environ.get('SESSION_TEMP_DIR', './temp/')
            output_dir = os.path.join(base_tmp, 'busco_results')
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Get filename without extension for BUSCO output
        base_name = os.path.splitext(os.path.basename(fasta_file))[0]
    busco_name = f"busco_{base_name}"
    busco_output_dir = os.path.join(output_dir, busco_name)

        # Resolve BUSCO binary dynamically (works in Docker image where busco is on PATH)
        busco_bin = shutil.which("busco")
        if not busco_bin:
            raise FileNotFoundError("BUSCO binary not found in PATH")

        # Run BUSCO command
        # cpus controls BUSCO's --cpu flag
        cmd = [
            busco_bin,
            "-i", fasta_file,
            "-o", busco_name,
            "-m", "genome",
            "-f",  # Force overwrite if output exists
            "-l", lineage,
            "--out_path", output_dir,
            "--cpu", str(cpus)  # Limit CPU usage for container
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            # Bubble up full context so the caller/UI can render it
            raise RuntimeError(
                f"BUSCO failed (code {result.returncode}).\ncmd: {' '.join(cmd)}\nSTDERR:\n{result.stderr}\nSTDOUT:\n{result.stdout}"
            )

        # Parse BUSCO results
        return parse_busco_results(busco_output_dir, lineage, busco_name)
        
    except subprocess.TimeoutExpired:
        print(f"BUSCO timed out for {fasta_file}")
        return None
    except FileNotFoundError:
        # BUSCO not available in this environment (e.g., Streamlit Cloud)
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


def compute_stats_with_busco(fasta_file: str, lineage: str = "bacteria_odb10", cpus: int = 1) -> Dict[str, Any]:
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
    busco_results = run_busco(fasta_file, lineage=lineage, cpus=cpus)
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


def process_directory(directory: str, include_busco: bool = True, busco_lineage: str = "bacteria_odb10", busco_cpus: int = 1) -> pd.DataFrame:
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
                stats = compute_stats_with_busco(filepath, lineage=busco_lineage, cpus=busco_cpus)
            else:
                stats = compute_stats(filepath)
            stats["File"] = str(filename.rsplit('.', 1)[0])
            results.append(stats)
    return pd.DataFrame(results)


def get_fastani_binary():
    """Find fastANI binary in common locations."""
    # Try common locations in order of preference
    paths = [
        shutil.which("fastANI"),  # PATH (most flexible)
        "/home/adminuser/.conda/bin/fastANI",  # Streamlit Cloud
        "/usr/local/bin/fastANI",  # Docker/conda
        "/opt/conda/bin/fastANI",  # Alternative Docker
        "/usr/bin/fastANI"  # System package
    ]
    for path in paths:
        if path and os.path.exists(path):
            return path
    raise FileNotFoundError("fastANI not found in any common locations. Install fastANI or ensure it's on PATH.")


def run_fastani(query_file, ref_file, output_file, threads=1, timeout=600):
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
    try:
        fastani_bin = get_fastani_binary()
    except FileNotFoundError as e:
        raise RuntimeError(f"fastANI binary not found: {e}")
    
    cmd = [
        fastani_bin,
        "-q", query_file,
        "-r", ref_file,
        "-o", output_file,
        "-t", str(threads)
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=timeout)
        if result.returncode != 0:
            stderr = (result.stderr or "").strip()
            raise RuntimeError(f"fastANI failed (returncode={result.returncode}): {stderr}\ncmd: {' '.join(cmd)}")

        # fastANI output: query, ref, ANI, fragments, matches
        if not os.path.exists(output_file):
            raise RuntimeError(f"fastANI did not create expected output file: {output_file}")

        with open(output_file) as f:
            line = f.readline()
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                try:
                    return float(parts[2])
                except ValueError:
                    raise RuntimeError(f"fastANI produced non-numeric ANI value: {parts[2]}")
        # If we reach here, output file didn't contain expected columns
        raise RuntimeError(f"fastANI output file malformed: {output_file}")
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"fastANI timed out after {timeout} seconds for {query_file} vs {ref_file}")
    except FileNotFoundError:
        raise RuntimeError("fastANI binary not found. Install fastANI or run in an environment that has fastANI available.")


def all_vs_all_fastani(file_paths, threads=1) -> Tuple[np.ndarray, List[str], Dict[Tuple[int,int], str]]:
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
    errors: Dict[Tuple[int,int], str] = {}
    for i, j in itertools.combinations(range(n), 2):
        base_tmp = os.environ.get('SESSION_TEMP_DIR', './temp/')
        out_file = os.path.join(base_tmp, f"fastani_{i}_{j}.txt")
        try:
            ani = run_fastani(file_paths[i], file_paths[j], out_file, threads=threads)
            ani_matrix[i, j] = ani
            ani_matrix[j, i] = ani
        except Exception as e:
            err_msg = str(e)
            errors[(i, j)] = err_msg
            ani_matrix[i, j] = np.nan
            ani_matrix[j, i] = np.nan
    np.fill_diagonal(ani_matrix, 100.0)
    return ani_matrix, genome_names, errors


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
   