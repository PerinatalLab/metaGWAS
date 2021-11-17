all_allPTD_coh_nms= ['23andme', 'ALSPAC', 'CHOP', 'DECODE', 'EGCUT', 'NFBC1966',  'MOBAGENETICS', 'DBDS', 'BIB', 'DNBCCASES', 'DNBCCONTROLS', 'DNBCPTD', 'HUNT', 'DILT1DGC', 'WTCCC58BC', 'FINNGEN']

all_postTerm_coh_nms= ['ALSPAC', 'EGCUT', 'DECODE', 'NFBC1966', 'MOBAGENETICS', 'DBDS', 'BIB', 'DNBCCASES', 'DNBCCONTROLS', 'HUNT', 'DILT1DGC', 'WTCCC58BC']

all_GAraw_coh_nms= ['23andme', 'ALSPAC', 'CHOP', 'DECODE', 'EFSOCH', 'HAPO', 'NFBC1966', 'STORK', 'STORKGROR', 'MOBAGENETICS', 'DBDS', 'EFSOCH', 'BIB', 'DNBCCASES', 'DNBCCONTROLS', 'HUNT', 'DNBCPTD', 'DILT1DGC', 'WTCCC58BC', 'CCHMC', 'GPN']

all_GAnrm_coh_nms= ['ALSPAC', 'CHOP', 'DECODE', 'EFSOCH', 'HAPO', 'NFBC1966', 'STORK', 'STORKGROR', 'MOBAGENETICS', 'DBDS', 'EFSOCH', 'BIB', 'DNBCCASES', 'DNBCCONTROLS', 'HUNT', 'DNBCPTD', 'DILT1DGC', 'WTCCC58BC']

allPTD_coh_nms= ['23andme', 'ALSPAC', 'CHOP', 'DECODE', 'EGCUT', 'NFBC1966',  'MOBAGENETICS', 'DBDS', 'DNBCGOYACASES', 'DNBCGOYACONTROLS', 'DNBCPTD', 'HUNT', 'DILT1DGC', 'WTCCC58BC', 'FINNGEN', 'PGPIII', 'PGPII']

postTerm_coh_nms= ['ALSPAC', 'EGCUT', 'DECODE', 'NFBC1966', 'MOBAGENETICS', 'DBDS', 'DNBCGOYACASES', 'DNBCGOYACONTROLS', 'HUNT', 'DILT1DGC', 'WTCCC58BC']

GAraw_coh_nms= ['23andme', 'ALSPAC', 'CHOP', 'DECODE',  'NFBC1966', 'STORK', 'STORKGROR', 'MOBAGENETICS', 'DBDS', 'DNBCGOYACASES', 'DNBCGOYACONTROLS', 'HUNT', 'DNBCPTD', 'DILT1DGC', 'WTCCC58BC', 'CCHMC', 'GPN', 'PGPIII', 'PGPII', 'BIB', 'HAPO', 'Viva', 'Gen3G', 'EFSOCH']

GAnrm_coh_nms= ['ALSPAC', 'CHOP', 'DECODE',  'NFBC1966', 'STORK', 'STORKGROR', 'MOBAGENETICS', 'DBDS', 'DNBCGOYACASES', 'DNBCGOYACONTROLS', 'HUNT', 'DNBCPTD', 'DILT1DGC', 'WTCCC58BC', 'PGPIII', 'PGPII', 'BIB', 'HAPO', 'Viva', 'Gen3G', 'EFSOCH']


repr_pheno_nms= ['miscarriage', 'GA_fetal', 'BW_maternal', 'AFB', 'AMenarche', 'AMenopause', 'NLB', 'Testosterone_fem', 'SHBG_fem', 'Oestradiol_fem', 'POP', 'Testosterone_male', 'PCOS', 'endometriosis', 'BW_fetal', 'BW_maternal_effect', 'BW_fetal_effect', 'leiomyoma_uterus', 'Preeclampsia', 'CBAT_fem', 'CBAT_male', 'SHBG_male']

top_coh_PTD_nms= ['23andme', 'DECODE', 'EGCUT', 'MOBAGENETICS', 'DBDS', 'FINNGEN']
top_coh_GAraw_nms= ['23andme', 'DECODE', 'MOBAGENETICS', 'DBDS']
top_coh_GAnrm_nms= ['DECODE', 'MOBAGENETICS', 'DBDS']
top_coh_postTerm_nms= ['DECODE', 'MOBAGENETICS', 'DBDS']

pheno_nms= ['allPTD', 'postTerm', 'GAraw', 'GAnrm']

big5_nms= ['MOBAGENETICS', 'DECODE', 'DBDS', 'HUNT', '23andme']
big5_nms_allPTD= ['23andme', 'MOBAGENETICS', 'DECODE', 'FINNGEN', 'HUNT', 'DBDS', 'EGCUT']

CCHMC_cohort_nms= ['ALSPAC', 'DNBC', 'FIN', 'GPN', 'HAPO']

include: 'scripts/munge_stats/Snakefile'
include: 'scripts/LDscore/Snakefile'

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
include: 'scripts/ADCY5/Snakefile'
include: 'scripts/stromal_cells/Snakefile'
include: 'scripts/figures/Snakefile'
include: 'scripts/effect_origin/Snakefile'
include: 'scripts/tables/Snakefile'
include: 'scripts/LCV/Snakefile'

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
		expand('/mnt/hdd/common/pol/metaGWAS/reports/forest/reports/forest_plot_{pheno}.pdf', pheno= pheno_nms),
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
		'/mnt/hdd/common/pol/metaGWAS/figures/ADCY5_PheWas.tiff',
		'/mnt/hdd/common/pol/metaGWAS/figures/ADCY5_FST_EUR_AFR.tiff',
		'/mnt/hdd/common/pol/metaGWAS/figures/BW_genetic_correlations.tiff',
		'/mnt/hdd/common/pol/metaGWAS/figures/repr_pheno_genetic_correlations.tiff',
		'/mnt/hdd/common/pol/metaGWAS/figures/partitioned_h2.tiff',
		'/mnt/hdd/common/pol/metaGWAS/figures/MacArthur_enrichment.tiff',
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
		'/mnt/hdd/common/pol/metaGWAS/figures/MR_GA_BW_maternal_effect.tiff',
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
		expand('/mnt/hdd/common/pol/metaGWAS/ADCY5/eQTL_Catalog/results/{coloc_out_eqtl}_GAraw.txt', coloc_out_eqtl= ['pph', 'results', 'individual_variant']),
		'/mnt/hdd/common/pol/metaGWAS/Tables/ADCY5_eQTL.txt',
		expand('/mnt/hdd/common/pol/metaGWAS/Tables/sample_size_{pheno}.txt', pheno= pheno_nms),
		'/mnt/hdd/common/pol/metaGWAS/figures/GAnrm_manhattan.png',
		'/mnt/hdd/common/pol/metaGWAS/LDScore/big5/RG/results/rp.txt',
		'/mnt/hdd/common/pol/metaGWAS/LDscore/big5/RG/meta/GAraw_allPTD_rg.log',
		expand('/mnt/hdd/common/pol/metaGWAS/figures/MR_GA_BW_{genome}_effect_haplotype_1.tiff', genome= ['maternal', 'fetal'])
