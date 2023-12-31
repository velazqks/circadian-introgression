Circadian Archaic Introgression Project
=======================================
` Keila S. Velazquez-Arcelay `


FILES
-----
Note: All files are in build hg19 unless otherwise spcified.

Raw data:
- raw_gencode.v29lift37.bed.xz
  - Wrangles the Gencode raw file and converts it to bed.  
- raw_cCREs.liftOver.to.Hg19.bed.xz
  - Candidate cis-regulatory elements (cCREs). Source: https://doi.org/10.1038/s41586-020-2493-4.  
- raw_kuhlwilm2019_filtered.bed.xz
  - Set of archaic variants analyzed by Kuhlwilm and Boeckx 2019: https://doi.org/10.1038/s41598-019-44877-x. Only bi-allelic sites are included. Generated by filter_kuhlwilm2019.py
- raw_kuhlwilm19_archaic_fixed.bed, raw_kuhlwilm19_human_fixed.bed, raw_kuhlwilm19_human_fixed_tony.bed
  - Human- and archaic- specific variants from Kuhlwilm and Boeckx 2019: https://www.nature.com/articles/s41598-019-44877-x. Generated by ../bin/parse_kuhlwilm2019.py. Source: https://figshare.com/articles/dataset/Variants_and_annotations_of_Neandertals/8184038.
- raw_introgressed_browning2018.bed.xz
  - Introgression map published by Browning et al., 2018.
- raw_introgressed_browninig2018.list
  - List or Browning et al., 2018 rsids.
- raw_introgressed_sankararaman2014.bed.xz
  - Introgression map published by Sankararaman et al., 2014
- raw_introgressed_schaefer2021.bed.xz
  - Introgression map published by Schaefer et al., 2021. SARGE
- raw_introgressed_skov2020.bed.xz
  - Introgression map published by Skov et al., 2020. hg38. Lifted to hg19. 
- raw_introgressed_steinruecken2018.bed.xz
  - Introgression map published by Steinrücken et al., 2018
- raw_introgressed_vernot2016.bed.xz
  - Introgression map published by Vernot et al., 2016
- raw_GTEx_v8_browning2018.bed.xz
  - Introgressed variants (Browning et al. 2018) that are eQTLs. Generated by circadian_and_introgressed_eqtls.py
- raw_predixcan_divergently_regulated.tsv
  - Divergently regulated genes in 3 archaic hominins. DR defined as genes with PrediXcan results containing an empirical p-value = 0 in Colbran et al., 2019.
- raw_nealelab_round2_1180.bed.xz
  - Morningness/eveningness (1180) associated variants in the UKBiobank GWAS run by the Neale lab.
- raw_SpliceAI_circadian_SAVs.tsv
  - Variants in circadian genes identified by SpliceAI as splice altering in archaics. General data generated by Colin Brand.
- raw_GTEx_v8_RNASeQCv1.1.9_gene_median_tpm.gct
  - RNASeq median gene-level TPM by GTEx tissue, 2017-06-05 analysis. 
- raw_predixcan_pvals_circadian_Altai.tsv, raw_predixcan_pvals_circadian_Vindija.tsv, raw_predixcan_pvals_circadian_Denisova.tsv
  - Circadian gene p-values from Altai, Vindija, and Denisovan archaic Predixcan divergently regulated genes analysis.
- raw_KGP_metadata.tsv 
  - 1000 Genomes Project population metadata for the 26 populations, including latitudes. Retrieved from the kgp  R package: https://stephenturner.github.io/kgp/
  - NOTE: Added lat_ancestral column to represent the latitude of the country of origin of diaspora populations:
    - GIH: 23.223, Gandhinagar in Gujarat
    - STU: 6.916667 Sri Jayawardenepura Kotte
    - ITU: 16.5131 Amaravati in Andhra Pradesh
    - CEU: 52.372778 Amsterdam (https://doi.org/10.1016/j.cub.2008.07.049)


Circadian gene evidence sources:
- raw_circadian_genes_biosystems.tsv
  - Retrieved from the old NCBI Biosystems database: 
    - NCBI: https://www.ncbi.nlm.nih.gov/Structure/biosystems/docs/biosystems_about.html
    - FTP: https://ftp.ncbi.nih.gov/pub/biosystems/
    - Citation: Geer LY, Marchler-Bauer A, Geer RC, Han L, He J, He S, Liu C, Shi W, Bryant SH. The NCBI BioSystems database. Nucleic Acids Res. 2010 Jan; 38(Database issue):D492-6, https://doi.org/10.1093/nar/gkp858
    - https://www.wikipathways.org/index.php/Pathway:WP3594
- raw_circadian_genes_cgdb_experimental.tsv
  - Retrieved from: http://cgdb.biocuckoo.org/animals_index.php -> Homo sapiens. Includes only genes identified by Experimental sources.
- raw_circadian_genes_go.tsv
  - Human protein-coding genes in the Gene Ontology database annotated with the GO:0007623 ("circadian rhythm") term or terms annotated with relationship "is_a”, “part_of”, “occurs_in”, or “regulates” circadian rhythm.
- raw_circadian_genes_gwas.tsv
  - Retrieved from the GWAS Catalog, accessed 2021-05-22 (or Dec.8). Extracted SNPs in EFO_0008328 and EFO_00004354.
- raw_circadian_genes_mcmahon.tsv
  - List of circadian genes and their pathways curated by Dr. Douglas McMahon
  

Gene and region sets:
- circadian_genes_candidate.tsv
  - Candidate circadian genes by confidence level.
- circadian_genes.bed
  - Circadian gene loci (medium and high confidence).
- circadian_genes_flanking.bed
  - Flanking circadian gene loci (1M up- and downstream; medium and high confidence).
- circadian_genes.list
  - List of 246 circadian GeneID and GeneName.
- circadian_promoter.bed
  - Promoter region for each circadian gene. Promoter region defined as -5kb and 1kb upstream and downstream from the transcription start site.
- circadian_genes_predixcan_dr.tsv
  - Divergently regulated circadian genes, defined as PrediXcan predictions with an empirical p-value = 0.
- circadian_genes_predixcan_pctl.tsv
  - Percentiles indicating where the archaic PrediXcan predictions fall in the human distribution.
- circadian_variants_sav.tsv
  - Archaic-specific variants in circadian genes with evidence of having splice-altering effects in at least 1 of the 4 high-coverage archaic genomes. General SAV dataset generated by Colin Brand using SpliceAI predictions.


Variant sets:
- circadian_variants.bed
  - 1000 Genomes Project variants in circadian genes or flanking these genes by 1Mb.
- circadian_variants_promoter.bed
  - Variants in circadian promoter regions (-5kb upstream and +1kb downstream of the TSS).
- circadian_variants_ccres.bed
  - Candidate cis-regulatory elements (cCREs) in the circadian flanking regions (1Mb). Published in https://doi.org/10.1038/s41586-020-2493-4.
- circadian_variants_hhmc.bed, circadian_variants_ahmc.bed
  - Circadian variants in a set of human- and archaic-specific variants. Source: Kuhlwilm and Boeckx 2019.
- circadian_variants_hhmc_promoters.bed, circadian_variants_ahmc_promoters.bed
  - Human- and archaic-specific variants in circadian promoter regions.
- circadian_variants_hhmc_ccres.bed, circadian_variants_ahmc_ccres.bed
  - Extracts human- and archaic-specific circadian variants that are candidate cis-regulatory elements (cCREs). Set of variants from Kuhlwilm and Boeckx 2019.
- circadian_variants_introgressed.bed.xz
  - Introgressed variants (Browning et al. 2018) in genic, promoter, and cCRE regions.
- circadian_variants_eqtl.bed.xz
  - eQTLs that influence the expression of circadian genes.


Statistical analysis:
- enrichment_circadian_fixed_snps_genic.txt
  - Fisher's exact test results of genic SNP enrichment in circadian gene regions.
- enrichment_circadian_fixed_snps_promoter.txt
  - Fisher's exact test results of promoter SNP enrichment in circadian gene regions.
- enrichment_circadian_fixed_snps_regulatory.txt
  - Fisher's exact test results of regulatory (cCRE) SNP enrichment in circadian gene flanking regions.
- enrichment_circadian_introgressed_gtex.txt
  - Fisher's exact test results of introgressed variants in circadian genes with evidence of being eQTL in GTEx.
- enrichment_introgressed_gwas_or_no_gwas.txt
  - Fisher's exact test results of difference between the circadian or non-circadian introgressed variants associated with at least 1 GWAS phenotype or no phenotype.


Other
- opentargets_input_browning18_hg38.txt.xz
  - List of Browning et al., 2018 introgressed variants in format: CHR_POS_REF_ALT
- opentargets_results_browning2018_signif_hg19.bed.xz
  - A filtered version of opentargets_output_introgressed_variants.tsv.xz. P-values are filtered by the genome-wide significance level (5x10^-8), and some columns have been removed.
- circadian_gene_count_by_predixcan_tissue_prediction.tsv
  - Counts of circadian genes containing  PrediXcan predictions in each tissue.
- circadian_gene_count_by_gtex_tissue_rnaseq_expression.tsv
  - Fraction of circadian genes expressed in each GTEx tissue with an RNASeq value equal or greater to 1.
- table_circadian_genes.tsv
  - Table 1: Circadian gene descriptions and evidence sources.
- table_circadian_genes_sav.xlsx
  - Table 3: Circadian genes containing archaic specific splice altering variants.
- table_circadian_genes_dr.xlsx
  - Table 4: Divergently regulated genes in each of the 3 archaic hominins. 
- chronotype_chr2_r2.tsv
  - R2 values between introgressed morningness variants generated by LDLink
  

Plotting data:
- enrichment_circadian_introgressed_by_tissue.tsv
  - Probability of success of circadian eQTLs in each GTEx tissue, in a background of all introgressed GTEx eQTLs.
- introgressed_morningness_support_by_2.clumped, introgressed_morningness_c_support_by_2.clumped, introgressed_morningness_support_by_all.clumped, introgressed_morningness_c_support_by_all.clumped
  - Introgressed UKBiobank morningness variants LD clumped (R2=0.5 in CEU).
- plotting_morningness_cumulative_fraction_i2.tsv, plotting_morningness_cumulative_fraction_i2c.tsv, plotting_morningness_cumulative_fraction_iall.tsv, plotting_morningness_cumulative_fraction_iallc.tsv
  - Cumulative fraction of UKBiobank morningness variants effect on the propensity of being a morning/evening person, sorted by P-value. i2 = introgressed variants evidenced by browning and at least one other method; i2c = circadian introgressed variants evidenced by browning and at least one other method; iall = introgressed variants evidenced by all the introgression detection methods; iallc = circadian introgressed variants evidenced by all the introgression detection methods.
- plotting_morningness_latitude_cline_browning2018_chr2.bed, plotting_morningness_latitude_cline_browning2018_c_chr2.bed, plotting_morningness_latitude_cline_browning2018_plus1_chr2.bed, plotting_morningness_latitude_cline_browning2018_plus1_c_chr2.bed, plotting_morningness_latitude_cline_browning2018_plusAll_chr2.bed, plotting_morningness_latitude_cline_browning2018_plusAll_c_chr2.bed
  - Circadian introgressed variants associated with chronotype in a UKBiobank GWAS. Population latitudes used for linear regression analysis to identify SNPs in this set following a latitudinal cline in 1000 Genomes Projects populations from European and Asian ancestry.

