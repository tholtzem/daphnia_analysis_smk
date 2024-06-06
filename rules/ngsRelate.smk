#rule vcf2Pop:
#  input:
#    'vcf/filtered/daphnia_10pops_vars_Q30_DP10_GQ30_alleles.imissRM.vcf.gz'
#  output:
#    vcf1 = 'analyses/ngsRelate/daphnia_{species}.vcf.gz',
#    vcf2 = 'analyses/ngsRelate/daphnia_{species}.alleles_lmiss20_maf0.01.vcf.gz',
#    ids = 'analyses/ngsRelate/daphnia_{species}.alleles_lmiss20_maf0.01.ids.txt'
#  log: 'log/vcf2Pop_daphnia_{species}.log'
#  message: """ Make vcfs per pop/species """
#  shell:
#    """
#    bcftools view -S list/{wildcards.species}.list -Oz -o {output.vcf1} {input}  &&
#    bcftools filter -e 'AC==0 || AC==AN || F_MISSING > 0.2 || MAF <= 0.01' -Oz -o {output.vcf2} {output.vcf1} &&
#    bcftools query -l {output.vcf2} > {output.ids} 2> {log}
#    """

rule vcf2Lake:
  input:
    #vcf = "vcf/filtered/daphnia_10pops_vars_Q30_DP10_GQ30_alleles.imissRM.vcf.gz",
    lakes = "analyses/ngsRelate/list/{lakes}.list"
  output:
    vcf1 = "analyses/ngsRelate/daphnia_{lakes}.vcf.gz",
    vcf2 = "analyses/ngsRelate/daphnia_{lakes}.alleles_lmiss20_maf0.01.vcf.gz",
    ids = "analyses/ngsRelate/daphnia_{lakes}.alleles_lmiss20_maf0.01.ids.txt"
  log: "log/vcf2Lake_daphnia_{lakes}_maf0.01.log"
  message: """ Make vcfs per lake """
  shell:
    """
    VCF=(vcf/filtered/daphnia_10pops_vars_Q30_DP10_GQ30_alleles.imissRM.vcf.gz)
    bcftools view -S {input.lakes} -Oz -o {output.vcf1} $VCF &&
    bcftools filter -e 'AC==0 || AC==AN || F_MISSING > 0.2 || MAF <= 0.01' -Oz -o {output.vcf2} {output.vcf1} &&
    bcftools query -l {output.vcf2} > {output.ids} 2> {log}
    """

#rule remove_PCAoutliers:
#  input:
#    'vcf/filtered/daphnia_9pops_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf0.01.vcf.gz'
#  output:
#    vcf1 = 'vcf/filtered/daphnia_9pops_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf0.01_PCAoutliersRM.vcf.gz',
#    vcf2 = 'vcf/filtered/daphnia_9pops_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf0.01_PCAoutliersRM.filtered.vcf.gz',
#    ids = 'vcf/filtered/daphnia_9pops_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf0.01_PCAoutliersRM.filtered.ids.txt'
#  log: 'log/PCAoutliersRM.log'
#  message: """ Remove PCA outliers (hybrids) """
#  shell:
#    """
#    bcftools view -s ^GSB1,MS11,PFAEFF2608,VAR1823 -Oz -o {output.vcf1} {input} &&
#    bcftools view -e 'AC==0 || AC==AN || F_MISSING > 0.2' -Oz -o {output.vcf2} {output.vcf1} &&
#    bcftools query -l {output.vcf2} > {output.ids} 2> {log}
#    """

rule ngsRelate:
  input:
    vcf = 'analyses/ngsRelate/daphnia_{lakes}.alleles_lmiss20_maf0.01.vcf.gz',
    ids = 'analyses/ngsRelate/daphnia_{lakes}.alleles_lmiss20_maf0.01.ids.txt'
  output:
    stats = 'analyses/ngsRelate/daphnia_{lakes}.alleles_lmiss20_maf0.01.filtered.stats.txt'
  log: 'log/ngsRelate_{lakes}_maf0.01.log'
  threads: 4
  message:
    """--- Calculating relatedness (etc.) with ngsRelate2 from vcf using called genotypes ---"""
  shell:
    """
    ~/bio/ngsRelate/ngsRelate -h {input.vcf} -p 4 -T GT -c 1 -z {input.ids} -O {output.stats} 2> {log}
    """


rule plot_ngsRelate:
  input:
    stats = 'analyses/ngsRelate/daphnia_{lakes}.alleles_lmiss20_maf0.01.filtered.stats.txt'
  output:
    touch('analyses/ngsRelate/plot_ngsRelate_{lakes}.done')
  log:
    'log/plot_ngsRelate_{lakes}.log'
  threads: 4
  message:
    """--- Calculating relatedness (etc.) with ngsRelate2 from vcf using called genotypes ---"""
  shell:
    """
    Rscript scripts/plot_ngsrelate.R {input} {wildcards.lakes} 2> {log}
    """
