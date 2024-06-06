rule SNP_maf:
  input:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated.vcf.gz'
  output:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}.vcf.gz'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_maf{maf}.log'
  message: """ Subset vcf and remove sites with maf < 0.1; i.e. for 68 inidividuals 136 alleles (100% of 136 = 136 alleles), 136*0.1 ~ 13.6 (136*0.05=6.8) """
  shell:
    """
    bcftools view -v snps {input.vcf} | bcftools filter -e 'MAF <= {wildcards.maf}' -Oz -o {output.vcf} --threads {threads} 2> {log}
    """


rule SNPs:
  input:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated.vcf.gz'
  output:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.log'
  message: """ Extract SNPs without maf-filter """
  shell:
    """
    bcftools view -v snps {input.vcf} -Oz -o {output.vcf} --threads {threads} 2> {log}
    """

