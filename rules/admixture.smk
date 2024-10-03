rule LD_pruning:
  input:
    vcf = 'bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.vcf.gz'
  output:
    touch('bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.LD_pruning.done')
  log: 'log/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.LD_pruning.log'
  threads: 4
  message: """--- Perform LD-pruning to identify unlinked sites ---"""
  shell:
    """
    prefix=`echo {input.vcf} | sed 's/.vcf.gz//g'`
    plink2 --vcf {input.vcf} --allow-extra-chr --double-id --set-missing-var-ids @:# --indep-pairwise 100 50 0.2 --out $prefix
    """


rule prune_vcf2plink:
  input:
    vcf = 'bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.vcf.gz',
    touched = 'bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.LD_pruning.done'
  output:
    touch('bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.LD_pruning.bed.done')
  log: 'log/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.LD_pruning.bed.log'
  threads: 4
  message: """---Generate a LD-pruned bed file in plink format from vcf ---"""
  shell:
    """
    prefix=`echo {input.vcf} | sed 's/.vcf.gz//g'`
    plink2 --vcf {input.vcf} --allow-extra-chr --double-id --set-missing-var-ids @:# --extract ${{prefix}}.prune.in --make-bed --out $prefix
    """

rule changeChrom_names:
  input:
    'bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.LD_pruning.bed.done'
  output:
    touch('bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.LD_pruning.chromNames.done')
  log: 'log/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.LD_pruning.chromNames.log'
  threads: 4
  message: """--- Change chromosome names (ADMIXTURE does not accept chr names that are not human chromosomes) ---"""
  shell:
    """
    prefix=`echo {input} | sed 's/.LD_pruning.bed.done//g'`
    awk '{{$1="0";print $0}}' ${{prefix}}.bim > ${{prefix}}.bim.tmp &&
    mv ${{prefix}}.bim.tmp ${{prefix}}.bim
    """


rule admixture:
  input:
    'bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.LD_pruning.chromNames.done'
  output:
    'analyses/admixture/log/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.LD_pruning_admixture_K{K}_log.out'
  log: 'log/admixture/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.LD_pruning_admixture_K{K}_log.log'
  threads: 4
  message: """--- Run admixture ---"""
  shell:
    """
    prefix=`echo {input} | sed 's/.LD_pruning.chromNames.done//g'`
    admixture --cv ${{prefix}}.bed {wildcards.K} > {output} 2> {log}
    mv ./*.Q analyses/admixture/data/
    mv ./*.P analyses/admixture/data/
    """


rule plotAdmixture:
  input:
    pops = 'analyses/pixy/list/daphnia_popfile_{sites_set}.tsv'
  output:
    'analyses/admixture/log/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.LD_pruning.png'
  log: 'log/plotAdmixture_daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.LD_pruning.log'
  threads: 4
  message: """--- Plot admixture between 2-15 K ---"""
  shell:
    """
    prefix=`echo analyses/admixture/data/daphnia_{wildcards.sites_set}_variants_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites_DP10_GQ30_GT{wildcards.GT}_SNPs_maf{wildcards.maf}.LD_pruning`
    Rscript ./scripts/plotADMIXTURE.r -p $prefix -i {input.pops} -k 15 -l cucullata,galeata,mendotae,longispina,zschokkei,longispinaFIN,dentifera 2> {log} 
    """

rule CV_error:
  input:
  'analyses/admixture/data/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}.png'
  output:
    'analyses/admixture/data/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}.LD_pruning_admixture.CV.error'
  log: 'log/admixture/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}.LD_pruning_admixture.CV.error.log'
  threads: 4
  message: """--- Identify the best value of k clusters (with lowest cross-validation error) ---"""
  shell:
    """
    awk '/CV/ {{print $3,$4}}' analyses/admixture/log/daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_maf{wildcards.maf}.LD_pruning_admixture_K*_log.out | sed -e 's/(//;s/)//;s/://;s/K=//' | sort -n > {output} 
    """


rule evalAdmix:
  input:
    bed = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}.LD_pruning.bed.done',
    admix = 'analyses/admixture/log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}.LD_pruning_admixture_K{K}_log.out'
  output:
    'analyses/admixture/data/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}_K{K}.evaladmixOut.corres'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}_K{K}.evaladmixOut.corres.log'
  threads: 4
  message: """--- Evaluate admixture ---"""
  shell:
    """
    bed=`echo {input.bed} | sed 's/.LD_pruning.bed.done//g'`
    admix=`echo analyses/admixture/data/daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_maf{wildcards.maf}`
    ~/bio/evalAdmix/evalAdmix -plink ${{bed}} -fname ${{admix}}.{wildcards.K}.P -qname ${{admix}}.{wildcards.K}.Q -o ${{admix}}_K{wildcards.K}.evaladmixOut.corres -P {threads} 2> {log}
    """


rule plotAdmixture_viridis:
  input:
    pops = 'analyses/pixy/list/daphnia_popfile_68.tsv'
  output:
    touch('analyses/admixture/data/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}_virids.done')
  log: 'log/plotAdmixture_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}_virids.log'
  threads: 4
  message: """--- Plot admixture between 2-15 K ---"""
  shell:
    """
    prefix=`echo analyses/admixture/data/daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_maf{wildcards.maf}`
    Rscript ./scripts/plotADMIXTURE_TH.r -p $prefix -i {input.pops} -k 15 -l curvirostris,lacustris,umbra,cucullata,galeata,mendotae,longispina,zschokkei,longispinaFIN,dentifera -o ${{prefix}}_viridis 2> {log} 
    """

