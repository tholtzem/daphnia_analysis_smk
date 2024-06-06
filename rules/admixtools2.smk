rule removeIND_fromvcf:
  input:
    #vcf = "vcf/filtered/daphnia_10pops_vars_Q30_DP10_GQ30_alleles.imissRM.vcf.gz",
    samplelist = "analyses/admixtools/list/sample.list"
  output:
    vcf1 = "analyses/admixtools2/daphnia_{lakes}.vcf.gz",
    vcf2 = "analyses/admixtools2/daphnia_{lakes}.alleles_lmiss20_maf0.01.vcf.gz",
    ids = "analyses/admixtools2/daphnia_{lakes}.alleles_lmiss20_maf0.01.ids.txt"
  log: "log/vcf2Lake_daphnia_{lakes}_maf0.01.log"
  message: """ Make vcfs per lake """
  shell:
    """
    VCF=(vcf/filtered/daphnia_10pops_vars_Q30_DP10_GQ30_alleles.imissRM.vcf.gz)
    bcftools view -S {input.lakes} -Oz -o {output.vcf1} $VCF &&
    bcftools filter -e 'AC==0 || AC==AN || F_MISSING > 0.2 || MAF <= 0.01' -Oz -o {output.vcf2} {output.vcf1} &&
    bcftools query -l {output.vcf2} > {output.ids} 2> {log}
    """






rule get_plink:
  input:
    fam = 'analyses/admixture/data/daphnia_{species}_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf{maf}.fam',
    bim = 'analyses/admixture/data/daphnia_{species}_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf{maf}.fam',
    touched = 'analyses/admixture/data/daphnia_{species}_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf{maf}.bed.done'
  output:
    directory('analyses/admixtools2/geno/')
  log: 'log/admixtools_geno.log'
  threads: 4
  message: """---Generate a LD-pruned bed file in plink format from vcf ---"""
  shell:
    """
    rsync {input.fam} {output}
    rsync {input.bim} {output}
    """
