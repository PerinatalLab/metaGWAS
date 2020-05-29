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
		'/mnt/hdd/common/pol/metaGWAS/validation/figures/mht_GA_BW.pdf',
		expand('/mnt/hdd/common/pol/metaGWAS/sumstats/{pheno}/MOBAGENETICS_{pheno}.txt', pheno= ['allPTD', 'earlyPTD', 'postTerm', 'GAraw', 'GAnrm'])


	
