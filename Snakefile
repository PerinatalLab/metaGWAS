all_allPTD_coh_nms= ['23andme', 'ALSPAC', 'CHOP', 'DECODE', 'EGCUT', 'NFBC1966',  'MOBAGENETICS', 'DBDS', 'BIB', 'DNBCCASES', 'DNBCCONTROLS', 'DNBCPTD', 'HUNT', 'DILT1DGC', 'WTCCC58BC', 'FINNGEN']

all_postTerm_coh_nms= ['ALSPAC', 'EGCUT', 'DECODE', 'NFBC1966', 'MOBAGENETICS', 'DBDS', 'BIB', 'DNBCCASES', 'DNBCCONTROLS', 'HUNT', 'DILT1DGC', 'WTCCC58BC']

all_GAraw_coh_nms= ['23andme', 'ALSPAC', 'CHOP', 'DECODE', 'EFSOCH', 'HAPO', 'NFBC1966', 'STORK', 'STORKGROR', 'MOBAGENETICS', 'DBDS', 'EFSOCH', 'BIB', 'DNBCCASES', 'DNBCCONTROLS', 'HUNT', 'DNBCPTD', 'DILT1DGC', 'WTCCC58BC', 'CCHMC', 'GPN']

all_GAnrm_coh_nms= ['ALSPAC', 'CHOP', 'DECODE', 'EFSOCH', 'HAPO', 'NFBC1966', 'STORK', 'STORKGROR', 'MOBAGENETICS', 'DBDS', 'EFSOCH', 'BIB', 'DNBCCASES', 'DNBCCONTROLS', 'HUNT', 'DNBCPTD', 'DILT1DGC', 'WTCCC58BC']

allPTD_coh_nms= ['23andme', 'ALSPAC', 'CHOP', 'DECODE', 'EGCUT', 'NFBC1966',  'MOBAGENETICS', 'DBDS', 'DNBCGOYACASES', 'DNBCGOYACONTROLS', 'DNBCPTD', 'HUNT', 'DILT1DGC', 'WTCCC58BC', 'FINNGEN', 'PGPIII', 'PGPII']

postTerm_coh_nms= ['ALSPAC', 'EGCUT', 'DECODE', 'NFBC1966', 'MOBAGENETICS', 'DBDS', 'DNBCGOYACASES', 'DNBCGOYACONTROLS', 'HUNT', 'DILT1DGC', 'WTCCC58BC']

GAraw_coh_nms= ['23andme', 'ALSPAC', 'CHOP', 'DECODE',  'NFBC1966', 'STORK', 'STORKGROR', 'MOBAGENETICS', 'DBDS', 'DNBCGOYACASES', 'DNBCGOYACONTROLS', 'HUNT', 'DNBCPTD', 'DILT1DGC', 'WTCCC58BC', 'CCHMC', 'GPN', 'PGPIII', 'PGPII', 'BIB', 'HAPO', 'Viva', 'Gen3G', 'EFSOCH']

GAnrm_coh_nms= ['ALSPAC', 'CHOP', 'DECODE',  'NFBC1966', 'STORK', 'STORKGROR', 'MOBAGENETICS', 'DBDS', 'DNBCGOYACASES', 'DNBCGOYACONTROLS', 'HUNT', 'DNBCPTD', 'DILT1DGC', 'WTCCC58BC', 'PGPIII', 'PGPII', 'BIB', 'HAPO', 'Viva', 'Gen3G', 'EFSOCH']


repr_pheno_nms= ['miscarriage', 'GA_fetal', 'BW_maternal', 'AFB', 'AMenarche', 'AMenopause', 'NLB', 'Testosterone_fem', 'SHBG_fem', 'Oestradiol_fem', 'POP', 'Testosterone_male', 'PCOS', 'endometriosis', 'BW_fetal', 'BW_maternal_effect', 'BW_fetal_effect', 'leiomyoma_uterus', 'Preeclampsia', 'CBAT_fem', 'CBAT_male', 'SHBG_male', 'Ruth_CBAT_female', 'Ruth_CBAT_male', 'Ruth_SHBG_female', 'Ruth_SHBG_male', 'Ruth_Testosterone_female', 'Ruth_Testosterone_male', 'Ruth_oestradiol']

top_coh_PTD_nms= ['23andme', 'DECODE', 'EGCUT', 'MOBAGENETICS', 'DBDS', 'FINNGEN']
top_coh_GAraw_nms= ['23andme', 'DECODE', 'MOBAGENETICS', 'DBDS']
top_coh_GAnrm_nms= ['DECODE', 'MOBAGENETICS', 'DBDS']
top_coh_postTerm_nms= ['DECODE', 'MOBAGENETICS', 'DBDS']

pheno_nms= ['allPTD', 'postTerm', 'GAraw', 'GAnrm']

big5_nms= ['MOBAGENETICS', 'DECODE', 'DBDS', 'HUNT', '23andme']
big5_nms_allPTD= ['23andme', 'MOBAGENETICS', 'DECODE', 'FINNGEN', 'HUNT', 'DBDS', 'EGCUT']

CCHMC_cohort_nms= ['ALSPAC', 'DNBC', 'FIN', 'GPN', 'HAPO']

cell_types_names= ['LED', 'ILC', 'Decidual', 'Macrophage-2', 'Macrophage-3', 'Macrophage-1', 'Monocyte', 'Stromal-1', 'Stromal-2', 'Endothelial-2', 'CD4_T-cell', 'CD8_T-cell', 'Smooth-muscle-cells-1', 'Endothelial-1', 'overall']

auto_CHR_nms= [1, 2, 3, 4, 5, 6, 7, 8,9, 10 ,11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

include: 'scripts/munge_stats/Snakefile'
include: 'scripts/LDscore/Snakefile'
include: 'scripts/repr_traits_PGS/Snakefile'
include: 'scripts/independent/Snakefile'
include: 'scripts/enrichment/Snakefile'
include: 'scripts/meta/Snakefile'
include: 'scripts/reports/Snakefile'
include: 'scripts/gene_based/Snakefile'
include: 'scripts/VEP/Snakefile'
include: 'scripts/LocusZoom/Snakefile'
include: 'scripts/LDscore_reprpheno/Snakefile'
include: 'scripts/top_regions/Snakefile'
include: 'scripts/colocalization/Snakefile'
include: 'scripts/iPSC/Snakefile'
include: 'scripts/other_metas/Snakefile'
include: 'scripts/MR/Snakefile'
include: 'scripts/nonadditive/Snakefile'
#include: 'scripts/ADCY5/Snakefile'
include: 'scripts/stromal_cells/Snakefile'
include: 'scripts/figures/Snakefile'
include: 'scripts/effect_origin/Snakefile'
include: 'scripts/tables/Snakefile'
include: 'scripts/LCV/Snakefile'
include: 'scripts/fetal_SNP/Snakefile'
include: 'scripts/EGG_sumstats/Snakefile'
include: 'scripts/eQTLs/Snakefile'
include: 'scripts/pQTLs/Snakefile'

rule all:
	'Files to collect'
	input:
		expand('/mnt/hdd/common/pol/metaGWAS/sumstats/META/meta_{pheno}_1.txt',pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/reports/QC/{pheno}/meta_QC_{pheno}.pdf', pheno=pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/repr_phenos/LDSC/results/{pheno}_rg_temp', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/repr_phenos/sumstats/{repr_pheno}.txt', repr_pheno= repr_pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/colocalization/{pheno}/pph_allpheno.txt', pheno= pheno_nms), 
		expand('/mnt/hdd/common/pol/metaGWAS/colocalization/{pheno}/results_allpheno.txt', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/sumstats/META/Maternal_GWAMA_{pheno}.txt.gz', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/FINEMAPPING/results/PP_{pheno}.txt', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/gene_based/fastBAT_{pheno}.txt.gene.fastbat', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/LDscore/{pheno}_temp', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/reports/coloc/{pheno}_reproductive_traits.pdf', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/reports/forest/reports/forest_plot_GAraw.pdf'),
		'/mnt/hdd/common/pol/metaGWAS/sumstats/META/other_meta/NOMOBA/noMOBA_Maternal_GAraw.txt.gz',
                '/mnt/hdd/common/pol/metaGWAS/sumstats/META/other_meta/NOMOBA/noMOBA_Maternal_allPTD.txt.gz',
                '/mnt/hdd/common/pol/metaGWAS/sumstats/META/other_meta/NO23andme/no23andme_Maternal_allPTD.txt.gz',
                '/mnt/hdd/common/pol/metaGWAS/sumstats/META/other_meta/NO23andme/no23andme_Maternal_GAraw.txt.gz',
		'/mnt/hdd/common/pol/metaGWAS/reports/other_meta/no23andme_Garaw.pdf',
		'/mnt/hdd/common/pol/metaGWAS/reports/other_meta/noMOBA_allPTD.pdf',
		'/mnt/hdd/common/pol/metaGWAS/reports/other_meta/noMOBA_GAraw.pdf',
		expand('/mnt/hdd/common/pol/metaGWAS/LDscore/part_h2/{pheno}.results', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/LDscore/h2_cts/{cts}/{pheno}.cell_type_results.txt', pheno= pheno_nms, cts= ['chromatin', 'gene_expr']),
		expand('/mnt/hdd/common/pol/metaGWAS/enrichment/GNOMAD/pLI_{pheno}.txt', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/enrichment/HPA/RNA_{pheno}.txt', pheno=pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/enrichment/MacArthur/diseases_{pheno}.txt', pheno=pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/eqtls/coloc/iPSC/results_{pheno}.txt', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/sumstats/nonadditive/META/Maternal_GWAMA_GAraw_{model}.txt.gz', model= ['dom', 'rec']),
		expand('/mnt/hdd/common/pol/metaGWAS/MR/repr_phenos/independent_signals/GAraw/indep_{repr_pheno}.clumped', repr_pheno= repr_pheno_nms),
		'/mnt/hdd/common/pol/metaGWAS/MR/repr_phenos/IVs/all_traits_GAraw.txt',
#		'/mnt/hdd/common/pol/metaGWAS/HESS/step2/allchrs.txt',
#		expand('/mnt/hdd/common/pol/metaGWAS/reports/QC/GAraw/meta_QC_{model}.pdf', model= ['rec', 'dom']),
		'/mnt/hdd/common/pol/metaGWAS/LDscore/individual_cohorts/h2/allcohorts.txt',
		expand('/mnt/hdd/common/pol/metaGWAS/LDscore/h2/{pheno}_h2.log', pheno= pheno_nms),
		'/mnt/hdd/common/pol/metaGWAS/LDscore/individual_cohorts/h2/allcohorts.txt',
		expand('/mnt/hdd/common/pol/metaGWAS/ADCY5/iPSC2/{ext}_GAraw.txt', ext= ['pph', 'results', 'individual_variant']),
		'/mnt/hdd/common/pol/metaGWAS/ADCY5/haplotypes/variant_POS/ADCY5_LD_buddies.txt',
		'/mnt/hdd/common/pol/metaGWAS/topregions/top_variants/GAraw.txt',
		'/mnt/hdd/common/pol/metaGWAS/repr_phenos/LDSC/results/repr_phenos_rg',
		'/mnt/hdd/common/pol/metaGWAS/stromal_cells/rna-seq/diff_expression.txt',
		'/mnt/hdd/common/pol/metaGWAS/figures/GAraw_manhattan.png',
		'/mnt/hdd/common/pol/metaGWAS/figures/gene_based_coloc_eqtl.pdf',
		'/mnt/hdd/common/pol/metaGWAS/figures/spider_BW_coloc_maternal.png',
		'/mnt/hdd/common/pol/metaGWAS/figures/ADCY5_PheWas.pdf',
		'/mnt/hdd/common/pol/metaGWAS/figures/ADCY5_FST_EUR_AFR.tiff',
		'/mnt/hdd/common/pol/metaGWAS/figures/BW_genetic_correlations.tiff',
#		'/mnt/hdd/common/pol/metaGWAS/figures/repr_pheno_genetic_correlations.tiff',
		'/mnt/hdd/common/pol/metaGWAS/figures/partitioned_h2.pdf',
		'/mnt/hdd/common/pol/metaGWAS/figures/MacArthur_enrichment.pdf',
		'/mnt/hdd/common/pol/metaGWAS/figures/ADCY5_effect_direction.pdf',
		#expand('/mnt/hdd/common/pol/metaGWAS/MR/results/{pheno}/MR_repr_phenos.txt', pheno= pheno_nms),
		'/mnt/hdd/common/pol/metaGWAS/MR/results/MVMR_repr_phenos.txt',
		expand('/mnt/hdd/common/pol/metaGWAS/enrichment/stromal_cells/RNA_decidualization_{pheno}.txt', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/effect_origin/META/meta_{haplotype}_1.txt', haplotype= ['h1', 'h2', 'h3']),
		'/mnt/hdd/common/pol/metaGWAS/effect_origin/META/haplotype_based_analysis.txt',
		'/mnt/hdd/common/pol/metaGWAS/figures/GAraw_effect_origin_lm_h2_maternal.tiff',
		'/mnt/hdd/common/pol/metaGWAS/Tables/Table1.txt',
#		'/mnt/hdd/common/pol/metaGWAS/topregions/non_additive/top_variants/GAraw_dom.txt',
		'/mnt/hdd/common/pol/metaGWAS/colocalization/GA/pph_GAraw_allPTD.txt',
		'/mnt/hdd/common/pol/metaGWAS/figures/BW_conditioning.tiff',
		'/mnt/hdd/common/pol/metaGWAS/effect_origin/META/joint/meta_maternal_effect_1.txt',
		'/mnt/hdd/common/pol/metaGWAS/effect_origin/META/joint/meta_fetal_effect_1.txt',
		'/mnt/hdd/common/pol/metaGWAS/colocalization/GA/GW/pph_GAraw_allPTD.txt',
		'/mnt/hdd/common/pol/metaGWAS/topregions/BW/final/BW_maternal_effect.txt',
		'/mnt/hdd/common/pol/metaGWAS/colocalization/BW/GW/pph_GAraw_BW.txt',
		'/mnt/hdd/common/pol/metaGWAS/figures/repr_pheno_coloc_main.tiff',
		'/mnt/hdd/common/pol/metaGWAS/repr_phenos/LDSC/mtCOJO/results/BW_maternal_effect_rg_temp',
		'/mnt/hdd/common/pol/metaGWAS/LDscore/individual_cohorts/h2/allPTD/allcohorts.txt',
		'/mnt/hdd/common/pol/metaGWAS/figures/GAraw_effect_origin_dend.tiff',
		expand('/mnt/hdd/common/pol/metaGWAS/LCV/results/repr_pheno/{pheno}_LCV.txt', pheno= pheno_nms),
		'/mnt/hdd/common/pol/metaGWAS/figures/repr_pheno_LCV.pdf',
		'/mnt/hdd/common/pol/metaGWAS/figures/repr_pheno_rg.tiff',
		expand('/mnt/hdd/common/pol/metaGWAS/enrichment/stromalRNA_{pheno}.txt', pheno= pheno_nms),
		'/mnt/hdd/common/pol/metaGWAS/figures/RNA_enrichment.pdf',
		expand('/mnt/hdd/common/pol/metaGWAS/figures/QQ_plot_{pheno}.png', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/Tables/Table1_extra_{pheno}.txt', pheno= pheno_nms),
		'/mnt/hdd/common/pol/metaGWAS/Tables/Table1_dominant.txt',
		'/mnt/hdd/common/pol/metaGWAS/Tables/effect_origin_table.txt',
		'/mnt/hdd/common/pol/metaGWAS/figures/h2_allpheno.pdf',
		'/mnt/hdd/common/pol/metaGWAS/figures/h2_allcohort.pdf',
		'/mnt/hdd/common/pol/metaGWAS/sumstats/META/NLB/meta_allPTD_single_pregnancy_1.txt',
		'/mnt/hdd/common/pol/metaGWAS/sumstats/META/NLB/meta_allPTD_all_pregnancies_1.txt',
		'/mnt/hdd/common/pol/metaGWAS/repr_phenos/LDSC/NLB/results/rg_temp.txt',
#		expand('/mnt/hdd/common/pol/metaGWAS/reports/QC/allPTD/file_lvl_{allPTD_coh}_allPTD.html', allPTD_coh= allPTD_coh_nms),
#		expand('/mnt/hdd/common/pol/metaGWAS/reports/QC/postTerm/file_lvl_{postTerm_coh}_postTerm.html', postTerm_coh= postTerm_coh_nms),
#		expand('/mnt/hdd/common/pol/metaGWAS/reports/QC/GAraw/file_lvl_{GAraw_coh}_GAraw.html', GAraw_coh= GAraw_coh_nms),
#		expand('/mnt/hdd/common/pol/metaGWAS/reports/QC/allPTD/filtered/file_lvl_{allPTD_coh}_allPTD.html', allPTD_coh= allPTD_coh_nms),
#                expand('/mnt/hdd/common/pol/metaGWAS/reports/QC/postTerm/filtered/file_lvl_{postTerm_coh}_postTerm.html', postTerm_coh= postTerm_coh_nms),
#                expand('/mnt/hdd/common/pol/metaGWAS/reports/QC/GAraw/filtered/file_lvl_{GAraw_coh}_GAraw.html', GAraw_coh= GAraw_coh_nms),
		'/mnt/hdd/common/pol/metaGWAS/ADCY5/PheWas/PAN_UKBB/flag/cleaned_file.txt',
		expand('/mnt/hdd/common/pol/metaGWAS/ADCY5/PheWas/FINNGEN/results/{coloc_out_FIN}_GAraw.txt', coloc_out_FIN= ['pph', 'results', 'individual_variant']),
		expand('/mnt/hdd/common/pol/metaGWAS/repr_phenos/LDSC/mtCOJO/h2/{BW_only}_only_effect_{GA_effect}_h2.log', BW_only= ['BW_maternal', 'BW_fetal'], GA_effect= ['GA']),
		'/mnt/hdd/common/pol/metaGWAS/repr_phenos/LDSC/h2/all_repr_phenoh2.txt',
		'/mnt/hdd/common/pol/metaGWAS/Tables/Maternal_fetal_effects.txt',
		'/mnt/hdd/common/pol/metaGWAS/figures/postTerm_manhattan.png',
		'/mnt/hdd/common/pol/metaGWAS/Tables/ADCY5_PheWas.txt',
		'/mnt/hdd/common/pol/metaGWAS/Tables/Genetic_correlations_males.txt',
#		expand('/mnt/hdd/common/pol/metaGWAS/ADCY5/eQTL_Catalog/results/{coloc_out_eqtl}_GAraw.txt', coloc_out_eqtl= ['pph', 'results', 'individual_variant']),
#		'/mnt/hdd/common/pol/metaGWAS/Tables/ADCY5_eQTL.txt',
		expand('/mnt/hdd/common/pol/metaGWAS/Tables/sample_size_{pheno}.txt', pheno= pheno_nms),
		'/mnt/hdd/common/pol/metaGWAS/figures/GAnrm_manhattan.png',
		'/mnt/hdd/common/pol/metaGWAS/LDScore/big5/RG/results/rp.txt',
		'/mnt/hdd/common/pol/metaGWAS/LDscore/big5/RG/meta/GAraw_allPTD_rg.log',
		expand('/mnt/hdd/common/pol/metaGWAS/figures/MR_GA_BW_{genome}_MT.pdf', genome= ['maternal', 'fetal']),
		'/mnt/hdd/common/pol/metaGWAS/figures/genet_corr_GAraw.pdf',
		'/mnt/hdd/common/pol/metaGWAS/figures/GAraw_effect_origin_ternary.tiff',
		'/mnt/hdd/common/pol/metaGWAS/Tables/genetic_instruments_MR_repr_phenos.txt',
		expand('/mnt/hdd/common/pol/metaGWAS/Tables/2MR_repr_phenos_{pheno}.txt', pheno= pheno_nms),
		'/mnt/hdd/common/pol/metaGWAS/repr_phenos/PGS/IVS/PGS_repr_traits.txt',
		'/mnt/hdd/common/pol/metaGWAS/BW/PGS_fetal_growth.txt',
		'/mnt/hdd/common/pol/metaGWAS/fetal_SNP/META/all_cohort.txt',
		expand('/mnt/hdd/common/pol/metaGWAS/sumstats/EGG/top-10K/EGG_Maternal_GWAMA_{pheno}_top10K.txt', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/sumstats/EGG/no_23andMe/EGG_Maternal_GWAMA_{phenotypes}.txt.gz', phenotypes= ['GAraw', 'postTerm', 'allPTD']),
		'/mnt/hdd/common/pol/metaGWAS/sumstats/META/other_meta/NO23andme_NOMOBA/NO23andme_noMOBA_Maternal_GAraw.txt.gz',
		expand('/mnt/hdd/common/pol/metaGWAS/sumstats/META/suggestive/Maternal_GWAMA_{pheno}.txt', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/het/{pheno}.txt', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/eqtls/coloc/endometrium/pph_{pheno}.txt', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/eqtls/coloc/GTEx/pph_{pheno}_{tissue}.txt', pheno= pheno_nms, tissue= ['Ovary', 'Uterus', 'Vagina']),
		expand('/mnt/hdd/common/pol/metaGWAS/pqtls/coloc/blood/{pheno}/pph_{prot}.txt', pheno= pheno_nms, prot= ['13481_24', '15535_3', '15635_4', '18882_7', '3708_62', '6605_17', '9204_33']),
		expand('/mnt/hdd/common/pol/metaGWAS/enrichment/labour_associated_DEGs/{pheno}.txt', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/LDscore/own_annot/results/all/{pheno}_labor_DEGs.txt', pheno= pheno_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/MR/results/{pheno}/MR_clusters.txt', pheno= pheno_nms),
		'/mnt/hdd/common/pol/metaGWAS/figures/SHBG_GAraw_2SMR.png',
		'/mnt/hdd/common/pol/metaGWAS/figures/cell_type_enrichment.pdf',
		'/mnt/hdd/common/pol/metaGWAS/figures/labor_deg_enrichment_pvalue.pdf',
		'/mnt/hdd/common/pol/metaGWAS/figures/shbg_coloc_GAraw_ternary.tiff',
		'/mnt/hdd/common/pol/metaGWAS/locuszoom/WNT4/plots/GAraw_rs9823520.pdf',
		expand('/mnt/hdd/common/pol/metaGWAS/colocalization/{repr_pheno}/pph_GAraw.txt', repr_pheno= repr_pheno_nms),
		'/mnt/hdd/common/pol/metaGWAS/figures/top_BW_conditioning2.pdf',
		'/mnt/hdd/common/pol/metaGWAS/figures/evo.pdf',
		expand('/mnt/hdd/common/pol/metaGWAS/figures/{prev_locus}_forest.pdf', prev_locus= ['EEFSEC', 'AGTR2', 'WNT4', 'EBF1', 'ADCY5']),
		'/mnt/hdd/common/pol/metaGWAS/figures/GA_BW_PGS_correlations.pdf',
		'/mnt/hdd/common/pol/metaGWAS/figures/GA_PTD_BETA_correlations.pdf',
		expand('/mnt/hdd/common/pol/metaGWAS/repr_phenos/sumstats/{repr_pheno}.txt', repr_pheno= repr_pheno_nms),
		'/mnt/hdd/common/pol/metaGWAS/figures/top_BW_conditioning2_individual.pdf',
		'/mnt/hdd/common/pol/metaGWAS/figures/top_BW_conditioning2_individual_decode.pdf'
