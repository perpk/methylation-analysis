import pandas as pd
import pybedtools
import os

def map_distal_enhancers_to_genes(test_bed_path, output_mapped_path):
    """
    Downloads hg19 RefSeq genes and maps intergenic CpGs to their closest target gene.
    """
    print("Downloading hg19 RefSeq gene bounds from UCSC...")
    # Fetch hg19 refGene table. Columns: 2=chrom, 4=txStart, 5=txEnd, 12=name2 (Gene Symbol)
    os.system(f"curl -s 'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz' | zcat | awk -v OFS='\t' '{{print $3, $5, $6, $13}}' > {output_mapped_path}/hg19_refseq.bed")
    
    # Sort the reference file (required for the closest function)
    os.system(f"sort -k1,1 -k2,2n {output_mapped_path}/hg19_refseq.bed > {output_mapped_path}/hg19_refseq_sorted.bed")
    
    # 1. Load your top 500 strict BED file and the reference into pybedtools
    enhancers = pybedtools.BedTool(test_bed_path)
    genes = pybedtools.BedTool(f"{output_mapped_path}/hg19_refseq_sorted.bed")
    
    print("Intersecting regions and calculating distances...")
    # 2. Find the closest gene. 
    # - d=True adds a column with the distance in base pairs.
    # - t='first' ensures we only get one gene per CpG to avoid duplicates.
    mapped = enhancers.closest(genes, d=True, t='first')
    
    # 3. Convert the results back to a Pandas DataFrame
    # Columns: [chr, start, end, name, ref_chr, ref_start, ref_end, Gene_Symbol, Distance]
    mapped_df = mapped.to_dataframe(
        names=['chr', 'start', 'end', 'cpg_id', 'ref_chr', 'ref_start', 'ref_end', 'target_gene', 'distance_bp']
    )
    
    # Filter out anything mapped to weird contigs and drop the temporary coordinate columns
    final_df = mapped_df[['cpg_id', 'chr', 'start', 'target_gene', 'distance_bp']].copy()
    final_df = final_df[final_df['target_gene'] != '.'] # Remove unmapped
    
    final_df.to_csv(f"{output_mapped_path}/mapped_genes.csv", index=False)
    print(f"Successfully mapped {len(final_df)} distal regions to target genes.")
    
    return final_df

def _format_and_export_strict_bed(df, output_path):
    """Helper function to apply strict genomic formatting and export."""
    
    # 1. Clean Chromosome Names (prevents 'chrchr' bugs if 'chr' is already present)
    chroms = 'chr' + df['CHR'].astype(str).str.replace(r'^chr', '', regex=True, case=False)
    
    # 2. Build coordinates (0-based, half-open, ensure no negative starts)
    starts = (df['MAPINFO'].astype(int) - 1).clip(lower=0)
    ends = df['MAPINFO'].astype(int) + 1
    
    # 3. Create raw BED DataFrame (Strict 4 columns: no scores/decimals allowed)
    bed_df = pd.DataFrame({
        'chrom': chroms,
        'start': starts,
        'end': ends,
        'name': df['IlmnID']
    })
    
    # 4. Filter to standard human chromosomes only
    valid_chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    bed_df = bed_df[bed_df['chrom'].isin(valid_chroms)].copy()
    
    # 5. Sort biologically (chr1 -> chr22 -> X -> Y) and by start position
    bed_df['chrom'] = pd.Categorical(bed_df['chrom'], categories=valid_chroms, ordered=True)
    bed_df = bed_df.sort_values(['chrom', 'start']).dropna(subset=['chrom'])
    
    # 6. Export as plain text (no header, no index)
    bed_df.to_csv(output_path, sep='\t', header=False, index=False)
    print(f"Strict BED saved to {output_path} ({len(bed_df)} regions)")
    
    return bed_df

def export_strict_test_bed(drivers_csv_path, manifest_df, top_n=500, output_path="top_500_pd_enhancers_strict.bed"):
    """Extracts top unannotated GAT drivers and exports directly to strict BED."""
    print(f"Generating strict test BED from {drivers_csv_path}...")
    drivers_df = pd.read_csv(drivers_csv_path)
    
    # Filter for unannotated regions
    unannotated = drivers_df[
        drivers_df['UCSC_RefGene_Name'].isna() | 
        (drivers_df['UCSC_RefGene_Name'].astype(str).str.lower() == 'nan')
    ]
    
    # Get top drivers and merge with physical coordinates
    top_cpgs = unannotated.sort_values(by='Attention_Score', ascending=False).head(top_n)
    merged_df = top_cpgs.merge(manifest_df[['IlmnID', 'CHR', 'MAPINFO']], on='IlmnID', how='inner')
    
    return _format_and_export_strict_bed(merged_df, output_path)

def export_strict_background_bed(manifest_df, common_probes, output_path="gat_background_universe_strict.bed"):
    """Extracts the entire GAT probe universe and exports directly to strict BED."""
    print("Generating strict background BED from manifest...")
    
    # Filter manifest to only probes used in the neural network
    bg_df = manifest_df[manifest_df['IlmnID'].isin(common_probes)].copy()
    bg_df = bg_df.dropna(subset=['CHR', 'MAPINFO'])
    
    return _format_and_export_strict_bed(bg_df, output_path)

