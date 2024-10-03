include: "rules/common.smk"


# -----------------------------------------------

rule all:
	input:
		### ----- Plotting raw all-sites bcf stats ----- ###
		#expand('bcf/stats/daphnia_{species}_nbrSites.txt', species=['init_10pops']),
		#expand('bcf/stats/daphnia_{species}_100Ksubset.bcf', species=['init_10pops']),
		#expand('bcf/stats/daphnia_{species}_100Ksubset_MQ_DP.tsv', species=['init_10pops']),
		#expand('bcf/stats/daphnia_{species}_100Ksubset_MQ_DP_2.{ext}', species=['init_10pops'], ext=['pdf', 'filter.list']),
		### ----- Filtering raw all-sites bcf ----- ###
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}{ext}', species=['init_10pops'], MQ=['30', '35'], DPmax=['IQR'], ext=['.bcf', '_nbrSites.txt']),
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_{ext}', species=['init_10pops'], MQ=['30', '35'], DPmax=['IQR'], ext=['sites.bcf', 'sites_nbrSites.txt']),
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_{ext}', species=['init_10pops'], MQ=['30', '35'], DPmax=['IQR'], ext=['invariant.bcf', 'invariant_nbrSites.txt']),
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_{ext}', species=['init_10pops'], MQ=['30', '35'],  DPmax=['IQR'], ext=['variant.bcf', 'variant_nbrSites.txt']),
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_{ext}', species=['init_10pops'], MQ=['30', '35'],  DPmax=['IQR'], ext=['final.bcf', 'final_nbrSites.txt']),
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_{ext}', species=['init_10pops'], MQ=['30', '35'],  DPmax=['IQR'], ext=['final.imiss']),
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_{ext}', species=['init_10pops'], MQ=['30', '35'],  DPmax=['IQR'], ext=['final.imiss.pdf', 'final.imiss.tsv']),
		### ----- select samples and exclude sites with more than 20 % missing genotypes ---- ###
		#expand('bcf/filtered/{DPmax}/daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites_{ext}', species=['init_10pops'], MQ=['30', '35', '40'],  DPmax=['IQR'], ext=['final_68_sites80.bcf']), ## skip for now
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_{ext}', species=['cucullata', 'curvirostris', 'dentifera', 'galeata', 'lacustris', 'longispinaFIN', 'longispina', 'mendotae', 'umbra', 'zschokkei'], MQ=['30', '35'],  DPmax=['IQR'], ext=['final.bcf']),
		### ----- pixy - pi ----- ###
		#expand('analyses/pixy/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs_{bpWindow}bpwindow_{ext}', species=['cucullata', 'curvirostris', 'dentifera', 'galeata', 'lacustris', 'longispinaFIN', 'longispina', 'mendotae', 'umbra', 'zschokkei'], MQ=['30'], DPmax=['IQR'], bpWindow=['100000'], ext=['pi.done']),
		#expand('bcf/merged/{DPmax}/daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs_merged.{ext}', MQ=['30'], DPmax=['IQR'],ext=['vcf.gz']), 
		#expand('analyses/pixy/daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs_merged_{bpWindow}bpwindow_{ext}', MQ=['30'], DPmax=['IQR'], bpWindow=['100000'], ext=['pixy.done']),
		### ----- Plotting raw variant-sites bcf stats ----- ###
		expand('bcf/stats/daphnia_init_10pops_variants_nbrSites.txt'),
		expand('/media/tania/jugglingjay1/DAPHNIA/daphnia_analysis_smk/bcf/concat/daphnia_init_10pops_variants.bcf'),
		expand('bcf/stats/daphnia_init_10pops_variants_100Ksubset.bcf'),
		expand('bcf/stats/daphnia_init_10pops_variants_100Ksubset_MQ_DP.tsv'),
		expand('bcf/stats/daphnia_init_10pops_variants_100Ksubset_MQ_DP.{ext}', ext=['pdf', 'filter.list']),
		### ----- Filtering raw variant-sites bcf ----- ###
		expand('bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ30_DPminmax{DPmax}{ext}', MQ=['30'], DPmax=['IQR', 'HengLi'], ext=['.bcf', '_nbrSites.txt']),
		expand('bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites{ext}', MQ=['30'], DPmax=['IQR', 'HengLi'], ext=['.bcf', '_nbrSites.txt']),
		expand('bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30{ext}', MQ=['30'], DPmax=['IQR', 'HengLi'], ext=['.bcf', '_nbrSites.txt']),
		expand('bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30{ext}', MQ=['30'], DPmax=['IQR', 'HengLi'], ext=['.imiss']),
		expand('bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.{ext}', MQ=['30'], DPmax=['IQR', 'HengLi'], ext=['imiss.pdf', 'imiss.tsv']),
		expand('bcf/variants/{DPmax}/daphnia_{species}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30{ext}', species=['cucullata', 'curvirostris', 'dentifera', 'galeata', 'lacustris', 'longispinaFIN', 'longispina', 'mendotae', 'umbra', 'zschokkei'], MQ=['30'], DPmax=['IQR', 'HengLi'], ext=['.bcf', '_nbrSites.txt']),
		### ----- isec ----- ### 
		expand('bcf/variants/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_{ext}', sites_set=['set1'], MQ=['30'],  DPmax=['IQR'], ext=['isec80.tsv']),
		#expand('bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT80.{ext}', sites_set=['set1', 'set2', 'set3', 'set4'], MQ=['30'],  DPmax=['IQR'], ext=['bcf']),
		expand('bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT100.{ext}', sites_set=['set1', 'set2', 'set2B', 'set3', 'set4'], MQ=['30'],  DPmax=['IQR'], ext=['bcf']),
		### ----- get SNPs without mafs ----- ###
		expand('bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs.{ext}', sites_set=['set1', 'set2', 'set2B', 'set3', 'set4'], MQ=['30'],  DPmax=['IQR'], GT=['100'], ext=['vcf.gz']),
		### ----- get SNPs with mafs ----- ###
		expand('bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.{ext}', sites_set=['set2', 'set2B','set4'], MQ=['30'],  DPmax=['IQR'], GT=['100'], maf=['0.01'], ext=['vcf.gz']),
		### ----- admixture ----- ###
		expand('bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.{ext}', sites_set=['set2', 'set2B','set4'], MQ=['30'],  DPmax=['IQR'], GT=['100'], maf=['0.01'], ext=['LD_pruning.done']),
		expand('bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.{ext}', sites_set=['set2', 'set2B','set4'], MQ=['30'],  DPmax=['IQR'], GT=['100'], maf=['0.01'], ext=['LD_pruning.bed.done']),
		expand('bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.{ext}', sites_set=['set2', 'set2B','set4'], MQ=['30'],  DPmax=['IQR'], GT=['100'], maf=['0.01'], ext=['LD_pruning.chromNames.done']),
		expand('analyses/admixture/log/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.LD_pruning_admixture_K{K}_log.out', sites_set=['set2', 'set2B','set4'], MQ=['30'],  DPmax=['IQR'], GT=['100'], maf=['0.01'], K=range(2,16)),
		expand('analyses/admixture/log/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.LD_pruning.{ext}', sites_set=['set4'], MQ=['30'],  DPmax=['IQR'], GT=['100'], maf=['0.01'], ext=['png']),
		#expand('analyses/admixture/data/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}.{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], maf=['0.01', '0.05', '0.1'], ext=['LD_pruning_admixture.CV.error']),
		#expand('analyses/admixture/data/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}_K{K}.{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], maf=['0.01', '0.05', '0.1'], K=range(2,16), ext=['evaladmixOut.corres']),
		#expand('analyses/admixture/data/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}_virids.done', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR'], maf=['0.01', '0.05']),
		#expand('bcf/filtered/{DPmax}/daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites_final_68_updated_SNPs_maf{maf}.{ext}', species=['init_10pops'], MQ=['30'], DPmax=['IQR'], maf=['0.01', '0.05'], ext=['vcf.gz']),
		### ----- bcf2vcf ----- ### 
		##expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['HengLi'], ext=['final.vcf.gz']),## skip for now
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_{ext}', species=['init_10pops'], MQ=['30'], DPmax=['HengLi'], ext=['final_68.bcf']),## skip for now
		#expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_{ext}', species=['init_10pops'], MQ=['30'], DPmax=['HengLi'], ext=['final_68_updated.bcf']),## skip for now
		### ----- get common variants/invariants ----- ###	
		### ----- pixy ----- ###
		#expand('analyses/pixy/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_{bpWindow}bpwindow_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR', 'HengLi'], bpWindow=['100000'], ext=['pixy.done']),
		### ---- svdq ---- ###
		##expand('bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR', 'HengLi'], ext=['pruned.vcf.gz']),
		##expand('analyses/svdq/SNPs/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['IQR', 'HengLi'], ext=['pruned.min4.nexus']),
		#expand('analyses/svdq/SNPs/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.{ext}', species=['init_10pops'], MQ=['30'], DPmax=['HengLi'], ext=['min4.nexus']),
		#expand('analyses/svdq/SNPs/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}.{ext}', species=['init_10pops'], MQ=['20'], DPmax=['HengLi'], outgroup=['curvirostris2', 'umbra'], ext=['done']),
		## ---- prepare fasttree input ---- ##
		#expand('analyses/fasttree/SNPs/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['HengLi', 'IQR'], ext=['min4.fasta']),
		## ---- run fasttree on leo5 not locally! ---- ##
		#expand('anexpand('analyses/phyml/{sets}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_bionj.{N}w_{window}_{outgroup}.trees.gz', sets=['ldgmcu'], species=['init_10pops'], MQ=['30', '40'], DPmax=['HengLi'], N=['50'], window=['sites'], outgroup=['curvirostris']),
		## ---- PhyML ML ---- ##
		#expand('analyses/phyml/{sets}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_ML.{N}w_{window}_{outgroup}.trees.gz', sets=['ldgmcu'], species=['init_10pops'], MQ=['30', '40'], DPmax=['HengLi', 'IQR'], N=['50'], window=['sites'], outgroup=['curvirostris']),
		### ---- twisst7 ---- ##
		#expand('analyses/twisst/{sets}/fixedtrees_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_{inf}.{N}w_{window}_{pop1}_{pop2}_{pop3}_{pop4}_{pop5}_{pop6}_{outgroup}_{ext}', sets=['ldgmcu'], species=['init_10pops'], MQ=['30', '40'], DPmax=['HengLi', 'IQR'], inf=['ML'], N=['50'], window=['sites'], pop1=['cucullata'], pop2=['galeata'], pop3=['mendotae'], pop4=['longispina'], pop5=['dentifera'], pop6=['umbra'], outgroup=['curvirostris'], ext=['weights.csv.gz', 'topos.nwk']),
		#### ----- Dsuite using notree ---- ###		
		#expand('analyses/Dsuite/Dtrios_test/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['HengLi', 'IQR'], outgroup=['curvirostris', 'umbra', 'curvirostris2'], ext=['Dtrios_test.done']),
		### ----- Dsuite using mito tree ---- ###
		#expand('analyses/Dsuite/mito/refs_mitoPCG_NT_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}_rf.nwk', outgroup=['curvirostris', 'umbra', 'curvirostris2'], ext=['rf.nwk']),
		#expand('analyses/Dsuite/mito/daphnia_popfile_{outgroup}_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['HengLi', 'IQR'], outgroup=['curvirostris', 'umbra', 'curvirostris2'], ext=['Dtrios_mitosvdq_tree.txt']),
		#expand('analyses/Dsuite/mito/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['HengLi', 'IQR'], outgroup=['curvirostris', 'umbra', 'curvirostris2'], ext=['Dtrios_mitosvdq_Fbranch.txt']),
		#expand('analyses/Dsuite/mito/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['HengLi', 'IQR'], outgroup=['curvirostris', 'umbra', 'curvirostris2'], ext=['Dtrios_mitosvdq_Fbranch.png', 'Dtrios_mitosvdq_Fbranch.svg']),
		### ----- Dsuite using SNPtree ---- ###
		#expand('analyses/Dsuite/SNP/daphnia_topology{topo}_{outgroup}_rf.nwk', outgroup=['curvirostris', 'umbra'], topo=range(1,4)),
		#expand('analyses/Dsuite/SNP/daphnia_popfile_{outgroup}_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_topology{topo}_{outgroup}_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['HengLi', 'IQR'], outgroup=['curvirostris', 'umbra'], topo=range(1,4), ext=['Dtrios_svdq_tree.txt']),
		#expand('analyses/Dsuite/SNP/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_topology{topo}_{outgroup}_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['HengLi', 'IQR'], outgroup=['curvirostris', 'umbra'], topo=range(1,4), ext=['Dtrios_svdq_Fbranch.txt']),
		#expand('analyses/Dsuite/SNP/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_topology{topo}_{outgroup}_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['HengLi', 'IQR'], outgroup=['curvirostris', 'umbra'], topo=range(1,4), ext=['Dtrios_svdq_Fbranch.png', 'Dtrios_svdq_Fbranch.svg']),
		#expand('analyses/Dsuite/SNP/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_topology{topo}_{outgroup}_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['HengLi', 'IQR'], outgroup=['curvirostris', 'curvirostris2', 'umbra'], topo=['1']), ext=['Dtrios_svdq_Fbranch.txt']),
		#expand('analyses/Dsuite/SNP/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_topology{topo}_{outgroup}_{ext}', species=['init_10pops'], MQ=['30', '40'], DPmax=['HengLi', 'IQR'], outgroup=['curvirostris2'], topo=['2']), ext=['Dtrios_svdq_Fbranch.txt']),
		#expand('analyses/Dsuite/SNP/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_topology{topo}_{outgroup}_{ext}', species=['init_10pops'], MQ=['40'], DPmax=['HengLi', 'IQR'], outgroup=['curvirostris'], topo=range(2,4), ext=['Dtrios_svdq_Fbranch.png', 'Dtrios_svdq_Fbranch.svg']),
		############################################################################################################################
		##expand('analyses/svdq/SNPs/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_LINEAGEtree_SVDQ_boostrapSTD_rooted_{outgroup}.{ext}', species=['init_10pops'], MQ=['40'], DPmax=['IQR'], outgroup=['umbra'], ext=['log']),
		###expand('analyses/ngsRelate/daphnia_{species}.{ext}', species=['curvirostris', 'cucullata', 'dentifera', 'galeata', 'lacustris', 'longispina', 'longispinaFIN', 'mendotae', 'umbra', 'zschokkei'], ext=['vcf.gz','alleles_lmiss20_maf0.01.vcf.gz','alleles_lmiss20_maf0.01.ids.txt']),
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
		## testing Dsuite and completeness of GTs
		##expand('bcf/isec/{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{ext}', sites_set=['set1', 'set2', 'set3'], MQ=['40'],  DPmax=['IQR'], ext=['isec80.tsv']),
		##expand('analyses/Dsuite/test/daphnia_popfile_{outgroup}_{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec80_{outgroup}_{ext}', MQ=['40'], DPmax=['IQR'], outgroup=['curvirostris2'], sites_set=['set3'], ext=['Dtrios_mitosvdq_tree.txt']),
		##expand('analyses/Dsuite/test/{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec80_{outgroup}_{ext}', MQ=['40'], DPmax=['IQR'], outgroup=['curvirostris2'], sites_set=['set3'], ext=['Dtrios_mitosvdq_Fbranch.txt']),
		##expand('analyses/Dsuite/test/{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec80_{outgroup}_{ext}', MQ=['40'], DPmax=['IQR'], outgroup=['curvirostris2'], sites_set=['set3'], ext=['Dtrios_mitosvdq_Fbranch.png', 'Dtrios_mitosvdq_Fbranch.svg'])
		#
		### ----- isec ----- ### 
		#expand('bcf/isec/{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites_final_{ext}', sites_set=['set1'], MQ=['30'],  DPmax=['IQR'], ext=['isec80.tsv']),

# ---------------------------------------------


#include: "rules/hardfilter_allSites.smk"
#include: "rules/genotypefilter_allSites.smk"
#include: "rules/pixy.smk"
include: "rules/hardfilter_variantSites.smk"
include: "rules/genotypefilter_variantSites.smk"
include: "rules/get_SNPs.smk"
include: "rules/admixture.smk"
#include: "rules/isec_common_allSites.smk"
#include: "rules/svdq.smk"
#include: "rules/fasttree.smk"
#include: "rules/twisst.smk"
#include: "rules/dsuite.smk"
##include: "rules/scikit_allel.smk"
##include: "rules/ngsRelate.smk"
##include: "rules/LD_decay.smk"
##include: "rules/admixtools2.smk"
##include: "rules/vcf_stats_filtering.smk"
##include: "rules/reformat_VCF.smk"
##include: "rules/test_dsuite.smk"
