# Original Hotspots Pipeline

## Input requirements

- VCF file needs to use the hg19 grch37 coordinate system
- limited to primary chromosomes 1-X/Y
- sorted
- indexed
- using ensemble nomeclature for chromsomes
- only report SNV
<br><br>

# Mutation Counts Pipeline

## Introduction
This pipeline aggregates mutations across multiple datasets at the per-Gene, per-Codon level. The pipeline creates an output table that counts up the 'Num_Samples_with_Mutation' at each GENE+Codon position. Currently, only considers Missense SNV mutations.

Additionally use the information from target BED files to figure out the 'Num_Samples_Targeted' at each GENE+Codon position.

Start with version 'mc_v0.5', the pipeline now consumes clinical data (<b>note</b>: limited to ER, PR, and HER2 receptor status at this point). The output table now provides Pvalue and OddsRatio columns indicating associations of each hotspot with ER, PR, and HER2 receptor groupings.

Finally, the output table also provide following the annotations:

1. ExAC frequency of each GENE+Codon position (computed as max of the ExAC frequencies of corresponding nucleotide positions).

2. Transcript_ID of each GENE+Codon position.

3. Alternate amino-acids generated by the Missense mutations at each GENE+Codon position.

## Input requirements

- Input datasets should be divided up by panel (e.g. DFCI-ONCOPANEL-1, DFCI-ONCOPANEL-2). For each panel, provide an annotated MAF file and corresponding target BED file.

- Annotated MAF file:

	- Annotated with maf2maf v1.6.9 (VEP 85). See [maf2maf.wdl](https://github.com/Intel-HSS/CCC_HotspotUsecase/blob/master/mutation_counts/workflow/maf2maf.wdl) for details.

	    - Pre-annotated MAF file must pass MAF Validation script ([link](https://github.com/Intel-HSS/CCC_HotspotUsecase/blob/master/mutation_counts/scripts/validate_maf.py)).

	- Use GRCh37 Coordinate System (limited to primary chromosomes: 1-22, X, Y, MT).

	- For an example, click [here](https://github.com/Intel-HSS/CCC_HotspotUsecase/blob/master/mutation_counts/examples/inputs/public/metabric.SNV.subset_genes.maf).
- Target BED file:

	- Use GRCh37 Coordinate System (limited to primary chromosomes: 1-22, X, Y, MT).

    - Must pass BED Validation script ([link](https://github.com/Intel-HSS/CCC_HotspotUsecase/blob/master/mutation_counts/scripts/count_mutations.py)).

	- For an example, click [here](https://github.com/Intel-HSS/CCC_HotspotUsecase/blob/master/mutation_counts/examples/inputs/public/metabric.targetedIntervals.r1.bed).

- Clinical Data file:
    - Currently contains data on ER, PR, and HER2 receptor status of each sample.

    - Must pass Clinical Data Validation script ([link](https://github.com/Intel-HSS/CCC_HotspotUsecase/blob/master/mutation_counts/scripts/validate_clinical.py)).

	- For an example, click [here](https://github.com/Intel-HSS/CCC_HotspotUsecase/blob/master/mutation_counts/examples/inputs/public/clinical/metabric.clinical_data.tsv).


## Output

- The pipeline generates a mutation counts table.

     - For an example, click [here](https://github.com/Intel-HSS/CCC_HotspotUsecase/blob/master/mutation_counts/examples/outputs/public/count_table.minimal.txt).
