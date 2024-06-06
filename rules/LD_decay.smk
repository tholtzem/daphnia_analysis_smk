
rule subset_VCF:
  input:
    vcf = 'vcf/angsd/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9088_hf.DP10_vmiss20.vcf.gz',
    samples = 'vcf/angsd/subset_sample.list'
  output:
    vcf = 'vcf/angsd/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9088_hf.DP10_vmiss20.sorted.vcf.gz'
  log: 'log/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9088_hf.DP10_vmiss20.sorted.log'
  threads: 4
  message: """--- Sort the order of the samples in the VCF file with natural sort for number within text (handy for downstream-analysis) ---"""
  shell:
    """
    bcftools view -S {input.samples} {input.vcf} -Oz -o {output.vcf} 2> {log}
    """


rule LD_decay:
  input:
    #vcf = 'vcf/filtered/daphnia_{species}_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf{maf}.vcf.gz'
    vcf = 'vcf/angsd/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9088_hf.DP10_vmiss20.sorted.vcf.gz'
  output:
    #touch('LD_decay/daphnia_{species}_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf{maf}.LD_decay.done')
    touch('LD_decay/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9088_hf.DP10_vmiss20.sorted.LD_decay_{kbWindow}kbWindow.done')
  log:
    #'log/daphnia_{species}_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf{maf}.LD_decay.log'
    'log/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9088_hf.DP10_vmiss20.sorted.LD_decay_{kbWindow}kbWindow.log'
  threads: 4
  message: """--- Determine where linkage disequilibrium decays to the genome-wide background (important for choosing the right window size) ---"""
  shell:
    """
    prefix=`echo {input.vcf} | sed 's!.vcf.gz!!g'`
    #out=`echo $prefix | sed 's!vcf/filtered/!LD_decay/!g'`
    out=`echo $prefix | sed 's!vcf/angsd/!LD_decay/!g'`
    plink --vcf {input.vcf} --allow-extra-chr --double-id --set-missing-var-ids @:# --maf 0.01 --geno 0.1 --mind 0.5 --thin 0.1 -r2 gz --ld-window 100 --ld-window-kb {wildcards.kbWindow} -ld-window-r2 0 --make-bed --out ${{out}}.LD_decay_{wildcards.kbWindow}kbWindow 2> {log}
    """


rule average_LD:
  input:
    #'LD_decay/daphnia_{species}_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf{maf}.LD_decay.done'
    'LD_decay/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9088_hf.DP10_vmiss20.sorted.LD_decay_{kbWindow}kbWindow.done'
  output:
    #touch('LD_decay/daphnia_{species}_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf{maf}.LD_average.done')
    touch('LD_decay/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9088_hf.DP10_vmiss20.sorted.LD_average_{kbWindow}kbWindow.done')
  log:
    #'log/daphnia_{species}_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf{maf}.LD_average.log'
    'log/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9088_hf.DP10_vmiss20.sorted.LD_average_{kbWindow}kbWindow.log'
  threads: 4
  conda: "../envs/py2.yaml"
  message: """ --- Calculate average LD across intervals (set distances)  in the genome. Note: execute python2 environment using "snakemake --cores 4 --use-conda" --- """
  shell:
    """
    prefix=`echo {input} | sed 's!.done!!g'`
    echo $prefix
    python2 scripts/ld_decay_calc.py -i ${{prefix}}.ld.gz -o ${{prefix}} 2> {log}
    """

