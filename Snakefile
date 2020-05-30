allPTD_coh_nms= ['23andme', 'ALSPAC', 'BIB', 'CHOP', 'DECODE', 'EGCUT', 'HARVEST', 'MOBA2008', 'NFBC1966', 'normentfeb', 'normentjan', 'normentjun', 'nomentmay', 'rotterdam1', 'rotterdam2', 'HUNT', 'MOBAGENETICS', 'DBDS']

earlyPTD_coh_nms= ['ALSPAC', 'BIB', 'DECODE', 'MOBA2008', 'NFBC1966sub', 'NFBC1966', 'normentfeb', 'normentjan', 'normentjun', 'nomentmay', 'rotterdam1', 'rotterdam2', 'HUNT', 'MOBAGENETICS', 'DBDS']

postTerm_coh_nms= ['ALSPAC', 'BIB', 'DECODE', 'HARVEST', 'NFBC1966sub', 'NFBC1966', 'normentfeb', 'normentjan', 'normentjun', 'nomentmay', 'rotterdam1', 'rotterdam2', 'HUNT', 'MOBAGENETICS', 'DBDS']

GAraw_coh_nms= ['23andme', 'ALSPAC', 'BIB', 'CHOP', 'DECODE', 'EFSOCH', 'HAPO', 'HARVEST', 'MOBA2008', 'NFBC1966sub', 'NFBC1966', 'normentfeb', 'normentjan', 'normentjun', 'nomentmay', 'rotterdam1', 'rotterdam2', 'STORK', 'STORKGROR', 'Viva', 'HUNT', 'MOBAGENETICS', 'DBDS']

GAnrm_coh_nms= ['ALSPAC', 'BIB', 'CHOP', 'DECODE', 'EFSOCH', 'HAPO', 'HARVEST', 'MOBA2008', 'NFBC1966sub', 'NFBC1966', 'normentfeb', 'normentjan', 'normentjun', 'nomentmay', 'rotterdam1', 'rotterdam2', 'STORK', 'STORKGROR', 'Viva', 'HUNT', 'MOBAGENETICS', 'DBDS']

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
		expand('/mnt/hdd/common/pol/metaGWAS/sumstats/allPTD/ldsc_input_{allPTD_coh}.txt', allPTD_coh= allPTD_coh_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/sumstats/earlyPTD/ldsc_input_{earlyPTD_coh}.txt', earlyPTD_coh= earlyPTD_coh_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/sumstats/postTerm/ldsc_input_{postTerm_coh}.txt', postTerm_coh= postTerm_coh_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/sumstats/GAnrm/ldsc_input_{GAnrm_coh}.txt', GAnrm_coh= GAnrm_coh_nms),
		expand('/mnt/hdd/common/pol/metaGWAS/sumstats/GAraw/ldsc_input_{GAraw_coh}.txt', GAraw_coh= GAraw_coh_nms)


	
