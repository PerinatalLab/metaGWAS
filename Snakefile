allPTD_coh_nms= ['23andme', 'ALSPAC', 'BIB', 'CHOP', 'DECODE', 'EGCUT', 'MOBA2008', 'NFBC1966', 'HUNT', 'MOBAGENETICS', 'DBDS']

earlyPTD_coh_nms= ['ALSPAC', 'BIB', 'DECODE', 'MOBA2008', 'NFBC1966sub', 'NFBC1966', 'HUNT', 'MOBAGENETICS', 'DBDS']

postTerm_coh_nms= ['ALSPAC', 'BIB', 'DECODE', 'NFBC1966sub', 'NFBC1966', 'HUNT', 'MOBAGENETICS', 'DBDS']

GAraw_coh_nms= ['23andme', 'ALSPAC', 'BIB', 'CHOP', 'DECODE', 'EFSOCH', 'HAPO', 'MOBA2008', 'NFBC1966sub', 'NFBC1966', 'STORK', 'STORKGROR', 'Viva', 'HUNT', 'MOBAGENETICS', 'DBDS']

GAnrm_coh_nms= ['ALSPAC', 'BIB', 'CHOP', 'DECODE', 'EFSOCH', 'HAPO', 'MOBA2008', 'NFBC1966sub', 'NFBC1966', 'STORK', 'STORKGROR', 'Viva', 'HUNT', 'MOBAGENETICS', 'DBDS']

include: 'scripts/validation/Snakefile'
include: 'scripts/munge_stats/Snakefile'
include: 'scripts/LDscore/Snakefile'
include: 'scripts/annotation/Snakefile'

rule all:
	'Files to collect'
	input:
		'/mnt/hdd/common/pol/metaGWAS/annotation/top_hits_annotated.txt',
		'/mnt/hdd/common/pol/metaGWAS/reports/top_hits_annotated.html',
		'/mnt/hdd/common/pol/metaGWAS/reports/independent_top_hits_annotated.html',
		'/mnt/hdd/common/pol/metaGWAS/validation/figures/mht_GA_BW.pdf',
		expand('/mnt/hdd/common/pol/metaGWAS/sumstats/{pheno}/MOBAGENETICS_{pheno}.txt', pheno= ['allPTD', 'earlyPTD', 'postTerm', 'GAraw', 'GAnrm']),
		expand('/mnt/hdd/common/pol/metaGWAS/LDsc/allPTD/{allPTD_coh}.txt.sumstats.gz', allPTD_coh= allPTD_coh_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/LDsc/earlyPTD/{earlyPTD_coh}.txt.sumstats.gz', earlyPTD_coh= earlyPTD_coh_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/LDsc/postTerm/{postTerm_coh}.txt.sumstats.gz', postTerm_coh= postTerm_coh_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/LDsc/GAnrm/{GAnrm_coh}.txt.sumstats.gz', GAnrm_coh= GAnrm_coh_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/LDsc/GAraw/{GAraw_coh}.txt.sumstats.gz', GAraw_coh= GAraw_coh_nms)


	
