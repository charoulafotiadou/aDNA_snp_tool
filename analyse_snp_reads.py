import pandas as pd
import pysam
import logging
from datetime import datetime
import argparse
import sys
from collections import defaultdict

# --- Logging Setup ---
def setup_logging(output_prefix):
    log_filename = f"{output_prefix}_snp_read_analysis_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_filename),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logging.info("Logging setup complete.")

# --- Load Node Definitions ---
def load_node_definitions(node_file):
    try:
        df = pd.read_csv(node_file, sep="\t")
        required_cols = {"Node", "Position", "Reference", "DerivedAllele"}
        if not required_cols.issubset(df.columns):
            raise ValueError(f"Node file missing required columns. Needed: {required_cols}")
        snp_info = {}
        for _, row in df.iterrows():
            pos = int(row["Position"])
            snp_info[pos] = {
                'node': row["Node"],
                'ref': row["Reference"].upper(),
                'derived': row["DerivedAllele"].upper()
            }
        logging.info(f"Loaded SNP definitions for {len(snp_info)} positions")
        return snp_info
    except Exception as e:
        logging.error(f"Error loading node definitions: {str(e)}")
        raise

# --- Analyze Reads Covering SNPs ---
def analyze_snp_reads(bam_file, snp_info):
    try:
        # Open BAM file
        bam = pysam.AlignmentFile(bam_file, "rb")
        ref_name = bam.references[0] 
        print(ref_name)
        logging.info(f"Processing BAM file: {bam_file}, reference: {ref_name}")

        # Initialize data structures
        read_details = []
        summary_dict = defaultdict(lambda: {
            'coverage': 0,
            'forward_reads': 0,
            'reverse_reads': 0,
            'forward_mismatch_count': 0,
            'forward_derived_count': 0,
            'reverse_mismatch_count': 0,
            'reverse_derived_count': 0,
            'forward_deamination_count': 0,
            'reverse_deamination_count': 0
        })

        # Process each SNP position
        for pos, info in snp_info.items():
            ref_pos = pos  # 
            ref_base = info['ref']
            derived_base = info['derived']
            node = info['node']

            # Skip positions based on deamination filter (C or G in ref/derived)
            skip_forward = 'C' in (ref_base, derived_base)
            skip_reverse = 'G' in (ref_base, derived_base)

            # Fetch reads covering this position
            for read in bam.fetch(ref_name, ref_pos - 1, ref_pos): 
                if read.is_unmapped:
                    continue

                # Determine read direction
                direction = 'forward' if not read.is_reverse else 'reverse'
                if (direction == 'forward' and skip_forward) or (direction == 'reverse' and skip_reverse):
                    continue

                # Get read sequence and alignment positions
                read_seq = read.query_sequence
                ref_positions = read.get_reference_positions(full_length=True)

                # Find the read position corresponding to the SNP
                read_pos = None
                for i, ref_pos_align in enumerate(ref_positions):
                    if ref_pos_align == ref_pos - 1:  
                        read_pos = i
                        break

                if read_pos is None:
                    continue 

                read_base = read_seq[read_pos] if read_pos < len(read_seq) else None
                if read_base is None:
                    continue

                # Update summary counts
                summary_dict[pos]['coverage'] += 1
                if direction == 'forward':
                    summary_dict[pos]['forward_reads'] += 1
                    if read_base != ref_base:
                        summary_dict[pos]['forward_mismatch_count'] += 1
                        if read_base == derived_base:
                            summary_dict[pos]['forward_derived_count'] += 1
                        if ref_base == 'C' and read_base == 'T':
                            summary_dict[pos]['forward_deamination_count'] += 1
                else:
                    summary_dict[pos]['reverse_reads'] += 1
                    if read_base != ref_base:
                        summary_dict[pos]['reverse_mismatch_count'] += 1
                        if read_base == derived_base:
                            summary_dict[pos]['reverse_derived_count'] += 1
                        if ref_base == 'G' and read_base == 'A':
                            summary_dict[pos]['reverse_deamination_count'] += 1

                # Collect detailed read information
                read_details.append({
                    'node': node,
                    'position': pos,
                    'reference_base': ref_base,
                    'derived_allele': derived_base,
                    'read_name': read.query_name,
                    'read_base': read_base,
                    'direction': direction,
                    'read_position': read_pos,
                    'is_derived': read_base == derived_base,
                    'is_deamination': (ref_base == 'C' and read_base == 'T' and direction == 'forward') or
                                     (ref_base == 'G' and read_base == 'A' and direction == 'reverse')
                })

        bam.close()
        logging.info("Finished processing BAM file")
        return read_details, summary_dict
    except Exception as e:
        logging.error(f"Error analyzing BAM file: {str(e)}")
        raise

# --- Compute Summary per Node ---
def compute_node_summary(summary_dict, snp_info):
    # Initialize node summary with all nodes from snp_info
    node_summary = {}
    
    # First pass: Initialize all nodes with 0 values
    for pos, info in snp_info.items():
        node = info['node']
        if node not in node_summary:
            node_summary[node] = {
                'total_coverage': 0,
                'derived_count': 0,
                'deaminated_count': 0,
                'matched_positions': 0,
                'total_positions': 0
            }

    # Second pass: Update with actual data where available
    for pos, info in snp_info.items():
        node = info['node']
        summary_data = summary_dict.get(pos, None)
        
        if summary_data:
            total_depth = summary_data['coverage']
            derived_count = summary_data['forward_derived_count'] + summary_data['reverse_derived_count']
            deaminated_count = summary_data['forward_deamination_count'] + summary_data['reverse_deamination_count']
            ref_count = total_depth - (summary_data['forward_mismatch_count'] + summary_data['reverse_mismatch_count'])
            node_summary[node]['total_coverage'] += total_depth
            node_summary[node]['derived_count'] += derived_count
            node_summary[node]['deaminated_count'] += deaminated_count
            
            if derived_count > 0:
                node_summary[node]['matched_positions'] += 1
                
            if derived_count > 0 or ref_count > 0:
                node_summary[node]['total_positions'] +=1

    # Create summary DataFrame
    summary_rows = []
    for node, stats in node_summary.items():
        total_cov = stats['total_coverage']
        derived = stats['derived_count']
        matched_pos = stats['matched_positions']
        total_pos = stats['total_positions']
        
        percent_with_deam = (derived / total_cov) * 100 if total_cov else 0

        summary_rows.append({
            'node': node,
            'matched_positions': matched_pos,
            'total_coverage': total_cov,
            'derived_count': derived,
            'derived_percent': round(percent_with_deam, 2),
            'total_positions': total_pos
        })

    # Sort nodes phylogenetically
    node_order = [
        "A", "B", "C", "D", "E", "F", "G", "H", "I", "J"
    ]
    summary_df = pd.DataFrame(summary_rows)
    summary_df["node"] = pd.Categorical(summary_df["node"], categories=node_order, ordered=True)
    return summary_df.sort_values("node")

# --- Save Results ---
def save_results(read_details, summary_df, output_prefix):
    try:
        # Save detailed read information
        read_details_df = pd.DataFrame(read_details)
        read_details_file = f"{output_prefix}_snp_read_details.csv"
        read_details_df.to_csv(read_details_file, index=False)
        logging.info(f"Saved detailed read information to {read_details_file}")

        # Save summary
        summary_file = f"{output_prefix}_snp_summary.csv"
        summary_df.to_csv(summary_file, index=False)
        logging.info(f"Saved summary to {summary_file}")
    except Exception as e:
        logging.error(f"Error saving results: {str(e)}")
        raise

# --- Main Execution ---
def main():
    parser = argparse.ArgumentParser(description="Detailed analysis of reads covering SNP positions in a BAM file")
    parser.add_argument("bam_file", help="Path to the BAM file")
    parser.add_argument("node_file", help="Path to the node definitions TSV file")
    parser.add_argument("output_prefix", help="Prefix for output files")
    args = parser.parse_args()

    setup_logging(args.output_prefix)

    try:
        logging.info("Starting SNP read analysis...")
        snp_info = load_node_definitions(args.node_file)
        read_details, summary_dict = analyze_snp_reads(args.bam_file, snp_info)
        summary_df = compute_node_summary(summary_dict, snp_info)
        save_results(read_details, summary_df, args.output_prefix)
        logging.info("Analysis completed successfully.")
    except Exception as e:
        logging.error(f"Analysis failed: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
