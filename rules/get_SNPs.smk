rule SNPs:
  input:
    bcf = 'bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}.bcf'
  output:
    vcf = 'bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs.vcf.gz'
  log: 'log/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs.log'
  message: """ Extract SNPs without maf-filter """
  shell:
    """
    bcftools view -v snps {input.bcf} -Oz -o {output.vcf} --threads {threads} 2> {log}
    """


rule SNP_maf:
  input:
    vcf = 'bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs.vcf.gz'
  output:
    vcf = 'bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.vcf.gz'
  log: 'log/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs_maf{maf}.log'
  message: """ For 68 inidividuals we have 136 alleles (100% of 136 = 136 alleles);
               remove variants if its minor allele is found in less than 1.36 (1%) or 6.8 (5 %) GTs """
  shell:
    """
    bcftools view -v snps {input.vcf} | bcftools filter -e 'MAF <= {wildcards.maf}' -Oz -o {output.vcf} --threads {threads} 2> {log}
    """


