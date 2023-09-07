CIRCADIAN INTROGRESSION PROJECT
===============================
` velazqks `

Modern humans and archaic hominins (i.e, Neanderthals and Denisovans) split from an African common ancestor approximately 700,000 years ago. While humans appeared and evolved in the same continent 300,000 years ago, archaic hominins evolved in Eurasia more than 400,000 years ago. Throughout their evolutionary history, these groups evolved at strikingly different latitudes: from near-Equatorial regions with little seasonal pattern fluctuation and constant UV exposure to a Northern temperate zone with marked fluctuation in photoperiod and lower UV exposure. This project shows that these groups maintained an excess of lineage-specific variation in the genomic regions associated with circadian genes. It also shows latitudinal clines in alleles associated with chronotype in present-day humans of Eurasian ancestry.



Directory

    ├──bin/
        ├──circadian_gene_set.py
        ├──circadian_gene_loci.py
        ├──circadian_variants.py
        ├──circadian_promoter.py
        ├──circadian_variants_ccres.py
        ├──parse_kuhlwilm19.py
        ├──circadian_variants_fixed.py
        ├──circadian_variants_fixed_promoter.py
        ├──circadian_variants_fixed_ccres.py
        ├──circadian_variants_introgressed.py
        ├──introgression_maps.py
        ├──circadian_genes_divergently_regulated.py
        ├──circadian_predixcan_pctl.py
        ├──circadian_variants_splice_altering.py
        ├──circadian_gene_count_by_predixcan_tissue_prediction.py
        ├──circadian_gene_count_by_gtex_tissue_rnaseq_expression.py
        ├──enrichment_circadian_fixed_snps_genic.py
        ├──enrichment_circadian_fixed_snps_promoter.py
        ├──enrichment_circadian_fixed_snps_regulatory.py
        ├──enrichment_circadian_introgressed_gtex.py
        ├──enrichment_introgressed_gwas_or_no_gwas.py
        ├──opentargets_generate_input.py
        ├──opentarget.py
        ├──opentarget_run.py
        ├──input_for_morningness_cumulative_fraction.py
        ├──table_circadian_genes.py
        ├──table_circadian_genes_sav.py
        ├──table_circadian_genes_dr.py
        ├──table_gwas_traits_per_snp.py
        ├──input_for_morningness_cumulative_fraction.py
        ├──input_for_morningness_latitude_cline.py
        ├──plot_upsetr_dr_and_sav.ipynb
        ├──plot_predixcan_distribution.ipynb
        ├──plot_enrichment_introgressed_circadian_in_gtex_tissue.ipynb
        ├──plot_dr_and_sav_piechart.ipynb
        ├──plot_introgressed_pleiotropy.ipynb
        ├──plot_morningness_direction_of_effect.ipynb
        ├──plot_morningness_latitude_cline.ipynb
        ├── README.md
    ├──data/
        ├──raw_gencode.v29lift37.bed.xz
        ├──raw_cCREs.liftOver.to.Hg19.bed.xz
        ├──raw_kuhlwilm2019_filtered.bed.xz
        ├──raw_kuhlwilm19_archaic_fixed.bed
        ├──raw_kuhlwilm19_human_fixed.bed
        ├──raw_kuhlwilm19_human_fixed_tony.bed
        ├──raw_introgressed_browning2018.bed.xz
        ├──raw_introgressed_browninig2018.list
        ├──raw_introgressed_sankararaman2014.bed.xz
        ├──raw_introgressed_schaefer2021.bed.xz
        ├──raw_introgressed_skov2020.bed.xz
        ├──raw_introgressed_steinruecken2018.bed.xz
        ├──raw_introgressed_vernot2016.bed.xz
        ├──raw_GTEx_v8_browning2018.bed.xz
        ├──raw_predixcan_divergently_regulated.tsv
        ├──raw_nealelab_round2_1180.bed.xz
        ├──raw_SpliceAI_circadian_SAVs.tsv
        ├──raw_GTEx_v8_RNASeQCv1.1.9_gene_median_tpm.gct
        ├──raw_predixcan_pvals_circadian_Altai.tsv
        ├──raw_predixcan_pvals_circadian_Vindija.tsv
        ├──raw_predixcan_pvals_circadian_Denisova.tsv
        ├──raw_KGP_metadata.tsv
        ├──raw_circadian_genes_biosystems.tsv
        ├──raw_circadian_genes_cgdb_experimental.tsv
        ├──raw_circadian_genes_go.tsv
        ├──raw_circadian_genes_gwas.tsv
        ├──raw_circadian_genes_mcmahon.tsv
        ├──circadian_genes_candidate.tsv
        ├──circadian_genes.bed
        ├──circadian_genes_flanking.bed
        ├──circadian_genes.list
        ├──circadian_promoter.bed
        ├──circadian_genes_predixcan_dr.tsv
        ├──circadian_genes_predixcan_pctl.tsv
        ├──circadian_variants_sav.tsv
        ├──circadian_variants.bed
        ├──circadian_variants_promoter.bed
        ├──circadian_variants_ccres.bed
        ├──circadian_variants_hhmc.bed
        ├──circadian_variants_ahmc.bed
        ├──circadian_variants_hhmc_promoters.bed
        ├──circadian_variants_ahmc_promoters.bed
        ├──circadian_variants_hhmc_ccres.bed
        ├──circadian_variants_ahmc_ccres.bed
        ├──circadian_variants_introgressed.bed.xz
        ├──circadian_variants_eqtl.bed.xz
        ├──enrichment_circadian_fixed_snps_genic.txt
        ├──enrichment_circadian_fixed_snps_promoter.txt
        ├──enrichment_circadian_fixed_snps_regulatory.txt
        ├──enrichment_circadian_introgressed_gtex.txt
        ├──enrichment_introgressed_gwas_or_no_gwas.txt
        ├──opentargets_input_browning18_hg38.txt.xz
        ├──opentargets_results_browning2018_signif_hg19.bed.xz
        ├──circadian_gene_count_by_predixcan_tissue_prediction.tsv
        ├──circadian_gene_count_by_gtex_tissue_rnaseq_expression.tsv
        ├──table_circadian_genes.tsv
        ├──table_circadian_genes_sav.xlsx
        ├──table_circadian_genes_dr.xlsx
        ├──introgressed_morningness_support_by_2.clumped
        ├──introgressed_morningness_c_support_by_2.clumped
        ├──introgressed_morningness_support_by_all.clumped
        ├──introgressed_morningness_c_support_by_all.clumped
        ├──enrichment_circadian_introgressed_by_tissue.tsv
        ├──plotting_morningness_cumulative_fraction_i2.tsv
        ├──plotting_morningness_cumulative_fraction_i2c.tsv
        ├──plotting_morningness_cumulative_fraction_iall.tsv
        ├──plotting_morningness_cumulative_fraction_iallc.tsv
        ├──plotting_morningness_latitude_cline.bed
        ├──plotting_morningness_latitude_cline_browning2018_chr2.bed
        ├──plotting_morningness_latitude_cline_browning2018_c_chr2.bed
        ├──plotting_morningness_latitude_cline_browning2018_plus1_chr2.bed
        ├──plotting_morningness_latitude_cline_browning2018_plus1_c_chr2.bed
        ├──plotting_morningness_latitude_cline_browning2018_plusAll_chr2.bed
        ├──plotting_morningness_latitude_cline_browning2018_plusAll_c_chr2.bed
        ├──chronotype_chr2_r2.tsv
        ├──README.md
    ├──plots/
        ├──Figure5b.pdf

