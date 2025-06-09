from Bio import AlignIO
from collections import Counter
import json
import sys
import logging
from datetime import datetime

def setup_logging():
    log_filename = f"unique_snps_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_filename),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
def load_node_mapping(node_mapping_file):
    with open(node_mapping_file) as f:
        content = f.read().strip()
        if not content:
            raise ValueError("File is empty")
        if content.startswith("node_mapping = "):
            content = content.split('=', 1)[1].strip()
        node_mapping = json.loads(content)
        logging.info(f"Successfully loaded mapping for {len(node_mapping)} nodes")
        return node_mapping

def identify_unique_snps_per_node(alignment_file, node_mapping, outgroup_id):
    alignment = AlignIO.read(alignment_file, "fasta")
    if len(alignment) == 0:
        raise ValueError("Alignment file is empty")
        
    sample_to_seq = {record.id: str(record.seq) for record in alignment}
    logging.info(f"Found {len(sample_to_seq)} sequences in alignment")
    
    if outgroup_id not in sample_to_seq:
        raise ValueError(f"Outgroup '{outgroup_id}' not found in alignment.")

    reference_seq = sample_to_seq[outgroup_id]

    alignment_len = len(reference_seq)
    logging.info(f"Alignment length: {alignment_len} bp")
    
    all_samples = set(sample_to_seq.keys()) - {outgroup_id}
    sorted_nodes = sorted(node_mapping.items(), 
                         key=lambda x: len(x[1]), 
                         reverse=True)
    
    logging.info("Identifying all SNPs for all nodes...")
    all_node_snps = {}
    total_snps = 0
    
    for node, sample_ids in sorted_nodes:
        node_snps = set()
        logging.info(f"Processing node {node} ({len(sample_ids)} samples)")
        
        node_samples = set(sample_ids) - {outgroup_id}  
        external_samples = all_samples - node_samples  
        
        for i in range(alignment_len):
            ref_base = reference_seq[i]
            if ref_base == '-':
                continue
                
            bases = [sample_to_seq[sample][i] for sample in sample_ids 
                     if sample in sample_to_seq and sample != outgroup_id]
            if not bases:
                continue
            
            external_bases = [sample_to_seq[sample][i] for sample in external_samples if sample in sample_to_seq]
               
            base_counts = Counter(b for b in bases if b != '-')
            if not base_counts:
                continue
                
            consensus_base, count = base_counts.most_common(1)[0]
            if consensus_base != ref_base and count == len(bases):
                if all(sample_to_seq[sample][i] == ref_base for sample in external_samples if sample in sample_to_seq):
                    node_snps.add((i+1, ref_base, consensus_base))
        all_node_snps[node] = node_snps
        total_snps += len(node_snps)
        logging.info(f"Found {len(node_snps)} SNPs for node {node}")
    
    logging.info(f"Identified {total_snps} total SNPs across all nodes")
    
    logging.info("Filtering to unique SNPs per node...")
    unique_node_snps = {}
    unique_counts = 0
    
    for i, (node, sample_ids) in enumerate(sorted_nodes):
        parent_nodes = [n for n, _ in sorted_nodes[:i]]
        parent_snps = set()
        for parent in parent_nodes:
            parent_snps.update(all_node_snps[parent])
        unique_snps = all_node_snps[node] - parent_snps
        unique_node_snps[node] = unique_snps
        unique_counts += len(unique_snps)
        
        logging.info(f"Node {node} has {len(unique_snps)} unique SNPs")
    
    logging.info(f"Found {unique_counts} unique SNPs across all nodes")
    return unique_node_snps

def save_snps(node_snps, output_file):
    try:
        logging.info(f"Saving results to {output_file}")
        with open(output_file, "w") as f:
            f.write("Node\tPosition\tReference\tDerivedAllele\n")
            for node, snps in node_snps.items():
                for pos, ref, alt in sorted(snps):
                    f.write(f"{node}\t{pos}\t{ref.upper()}\t{alt.upper()}\n")
        logging.info(f"Successfully saved results")
    except Exception as e:
        logging.error(f"Error saving SNPs: {str(e)}", exc_info=True)
        raise

if __name__ == "__main__":
    setup_logging()
    
    try:
        if len(sys.argv) != 5:
            logging.error("Incorrect number of arguments")
            logging.info("Usage: python script.py <alignment.fasta> <node_mapping.json> <outgroup_id> <output.tsv>")
            sys.exit(1)
            
        alignment_file = sys.argv[1]
        node_mapping_file = sys.argv[2]
        outgroup_id = sys.argv[3]
        output_file = sys.argv[4]

        logging.info("Starting unique SNP identification")
        logging.info(f"Alignment file: {alignment_file}")
        logging.info(f"Node mapping: {node_mapping_file}")
        logging.info(f"Output file: {output_file}")

        # Load and validate node mapping
        node_mapping = load_node_mapping(node_mapping_file)
        
        # Identify unique SNPs
        unique_snps = identify_unique_snps_per_node(alignment_file, node_mapping, outgroup_id)
        
        # Save results
        save_snps(unique_snps, output_file)
        
        logging.info("Analysis completed successfully")
        sys.exit(0)
        
    except Exception as e:
        logging.error(f"Script failed: {str(e)}", exc_info=True)
        sys.exit(1)