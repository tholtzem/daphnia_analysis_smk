include: "rules/common.smk"


# -----------------------------------------------

rule all:
	input:
		### -----Plotting bcf stats ----- ###
		#expand('bcf/stats/daphnia_{species}_nbrSites.txt', species=['init_10pops']),
		#expand('bcf/stats/daphnia_{species}_100Ksubset.bcf', species=['init_10pops']),
		#expand('bcf/stats/daphnia_{species}_100Ksubset_MQ_DP.tsv', species=['init_10pops']),
		#expand('bcf/stats/daphnia_{species}_100Ksubset_MQ_DP.{ext}', species=['init_10pops'], ext=['pdf', 'filter.list']),
		### ----- Filtering bcf ----- ###
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], ext=['.bcf', '_nbrSites.txt']),
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], ext=['sites80.bcf', 'sites80_nbrSites.txt']),
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], ext=['invariant.bcf', 'invariant_nbrSites.txt']),
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], ext=['variant.bcf', 'variant_nbrSites.txt']),
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], ext=['final.bcf', 'final_nbrSites.txt']),
		### ----- bcf2vcf ----- ###
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], ext=['final.vcf.gz']),
		### ----- select samples ----- ###
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], ext=['final_68.vcf.gz']),
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], ext=['final_68_updated.vcf.gz']),
		### ----- get SNPs with mafs ----- ###
		expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}.{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], maf=['0.01', '0.05'], ext=['vcf.gz']),
		### ----- get SNPs without mafs ----- ###
		expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], ext=['vcf.gz']),
		### ----- admixture ----- ###
		expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}.{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], maf=['0.01', '0.05'], ext=['LD_pruning.done']),
		expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}.{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], maf=['0.01', '0.05'], ext=['LD_pruning.bed.done']),
		expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}.{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], maf=['0.01', '0.05'], ext=['LD_pruning.chromNames.done']),
		expand('analyses/admixture/log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}.LD_pruning_admixture_K{K}_log.out', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], maf=['0.01', '0.05', '0.1'], K=range(2,16)),
		expand( 'analyses/admixture/data/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}.{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], maf=['0.01', '0.05', '0.1'], ext=['png']),
		expand('analyses/admixture/data/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}.{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], maf=['0.01', '0.05', '0.1'], ext=['LD_pruning_admixture.CV.error']),
		expand('analyses/admixture/data/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}_K{K}.{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], maf=['0.01', '0.05', '0.1'], K=range(2,16), ext=['evaladmixOut.corres']),
		#expand('analyses/admixture/data/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}_virids.done', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], maf=['0.01', '0.05']),
		### ----- pixy ----- ###
		expand('analyses/pixy/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], ext=['pixy.done']),
		### ----- svdq---- ###
		expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], ext=['pruned.vcf.gz']),
		expand('analyses/svdq/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], ext=['pruned.min4.nexus']),
		expand('analyses/svdq/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_{ext}', species=['init_10pops'], MQ=['40'], DPmax=['IQR'], ext=['SPECIEStree_SVDQ_boostrapSTD_rooted.nwk']),
		expand('analyses/svdq/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_{ext}', species=['init_10pops'], MQ=['40'], DPmax=['IQR'], ext=['LINEAGEtree_SVDQ_boostrapSTD_rooted.nwk']),
		#expand('analyses/svdq/daphnia_speciesTree_randomSVDQ_boostrapSTD_rooted.nwk')
		#expand('analyses/ngsRelate/daphnia_{species}.{ext}', species=['curvirostris', 'cucullata', 'dentifera', 'galeata', 'lacustris', 'longispina', 'longispinaFIN', 'mendotae', 'umbra', 'zschokkei'], ext=['vcf.gz','alleles_lmiss20_maf0.01.vcf.gz','alleles_lmiss20_maf0.01.ids.txt']),
		#expand('analyses/ngsRelate/daphnia_{species}.{ext}', species=['curvirostris', 'cucullata', 'dentifera', 'galeata', 'longispina', 'longispinaFIN', 'mendotae', 'umbra', 'zschokkei'], ext=['alleles_lmiss20_maf0.01.filtered.stats.txt']),
		#expand('analyses/ngsRelate/plot_ngsRelate_{species}.done', species=['curvirostris', 'cucullata', 'dentifera', 'galeata', 'longispina', 'longispinaFIN', 'mendotae', 'umbra', 'zschokkei']),
		#expand('vcf/stats/daphnia_{species}_{ext}', species=['{species}'], ext=['invariant_nbrSites.txt', 'variant_alleles_nbrSites.txt']),
		#expand('vcf/stats/daphnia_{species}_variant.frq.done'),
		#expand('analyses/scikit_allel/data/daphnia_{species}_variant_maf5NEW.vcf.gz'),
		#expand('analyses/scikit_allel/data/daphnia_{species}_variant_maf5NEW_renamed.vcf.gz'),
		#expand('analyses/ngsRelate/daphnia_{lakes}.{ext}', lakes=['Bielersee', 'Brienzersee', 'Fagerheim', 'Frauenfeld', 'Lago_di_Maggiore', 'Lake_Whitmore', 'LC', 'Mondsee', 'Oneida_Lake', 'Otrovatnet', 'Sihlsee', 'Timmelsjoch'], ext=['vcf.gz', 'alleles_lmiss20_maf0.01.vcf.gz', 'alleles_lmiss20_maf0.01.ids.txt']),
		#expand('analyses/ngsRelate/daphnia_{lakes}.{ext}', lakes=['Bielersee', 'Brienzersee', 'Fagerheim', 'Frauenfeld', 'Lago_di_Maggiore', 'Lake_Whitmore', 'LC', 'Mondsee', 'Oneida_Lake', 'Otrovatnet', 'Sihlsee', 'Timmelsjoch'], ext=['alleles_lmiss20_maf0.01.filtered.stats.txt']),
		#expand('analyses/ngsRelate/plot_ngsRelate_{lakes}.{ext}', lakes=['Bielersee', 'Brienzersee', 'Fagerheim', 'Frauenfeld', 'Lago_di_Maggiore', 'Lake_Whitmore', 'LC', 'Mondsee', 'Oneida_Lake', 'Otrovatnet', 'Sihlsee', 'Timmelsjoch'], ext=['done']),
		#expand('analyses/admixtools2/geno/')
		#expand('LD_decay/daphnia_{species}_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf{maf}.LD_decay.done', species=['9pops'], maf=['0.01'], ext=['LD_decay.done']),
		#expand('LD_decay/daphnia_{species}_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf{maf}.LD_average.done', species=['9pops'], maf=['0.01'], ext=['LD_decay.done']),
		#expand('vcf/angsd/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9088_hf.DP10_vmiss20.sorted.vcf.gz'),
		#expand('LD_decay/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9088_hf.DP10_vmiss20.sorted.LD_decay_{kbWindow}kbWindow.done', kbWindow=['1000', '10000']),
		#expand('LD_decay/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9088_hf.DP10_vmiss20.sorted.LD_average_{kbWindow}kbWindow.done', kbWindow=['1000', '10000'])

# ---------------------------------------------


#include: "rules/hardfilter_allSites.smk"
#include: "rules/genotypefilter_allSites.smk"
include: "rules/get_SNPs.smk"
include: "rules/admixture.smk"
include: "rules/pixy.smk"
include: "rules/svdq.smk"
##include: "rules/scikit_allel.smk"
##include: "rules/ngsRelate.smk"
##include: "rules/LD_decay.smk"
##include: "rules/admixtools2.smk"
##include: "rules/vcf_stats_filtering.smk"
##include: "rules/reformat_VCF.smk"
