rule get_Invariant:
  input:
    'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80.bcf'
  output:
    bcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_invariant.bcf',
    n_sites = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_invariant_nbrSites.txt'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_invariant.log'
  threads: 4
  message: """ --- Select invariants using max allele count equals (-C0) and filter for minimum site depth --- """
  shell:
    """
    bcftools view -C0 -M2 {input} | bcftools filter -i 'FMT/DP>=10' -Ob -o {output.bcf} &&
    bcftools view -H {output.bcf} | wc -l > {output.n_sites} 2> {log}
    """


rule get_Variant:
  input:
    'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80.bcf'
  output:
    bcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_variant.bcf',
    n_sites = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_variant_nbrSites.txt'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_variant.log'
  threads: 4
  message: """--- Select biallelic variant sites (excluding 0|0 genotypes, i.e. AC=0) and filter for minimum site depth, genotype quality, and QUAL filter ---"""
  shell:
    """
    bcftools view -c1 {input} -m 2 -M 2 | bcftools filter -i 'FMT/DP>=10 & GQ>=30 & QUAL>=30' -Ob -o {output.bcf} &&
    bcftools view -H {output.bcf} | wc -l > {output.n_sites} 2> {log}
    """

rule finalBCF:
  input:
    invariant = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_invariant.bcf',
    variant = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_variant.bcf'
  output:
    bcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final.bcf',
    n_sites = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_nbrSites.txt'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final.log'
  threads: 4
  message: """ --- Combine variantVCF and invariantVCF --- """
  shell:
    """
    tabix -p bcf -f {input.invariant} &&
    tabix -p bcf -f {input.variant} &&
    bcftools concat --allow-overlaps {input.invariant} {input.variant} -Ob -o {output.bcf} --threads {threads} &&
    bcftools view -H {output.bcf} | wc -l > {output.n_sites} 2> {log}
    """

rule bcf2vcf:
  input:
    'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final.bcf'
  output:
    'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final.vcf.gz'
  log:
    'log/LC/BCF2VCF_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final.log'
  message: """ --- Convert bcf 2 uncompressed vcf for downstream analysis --- """
  threads: 4
  shell:
    """
    bcftools index -f --csi {input} --threads {threads} &&
    bcftools convert -Oz -o {output} {input} --threads {threads} 2> {log}
    """

rule select_samples:
  input:
    ids = 'analyses/pixy/list/daphnia_samplelist_68.tsv',
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final.vcf.gz'
  output:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68.vcf.gz'
  log:
    'log/select_samples_vcf_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68.log'
  message: """ --- Select individuals to include in the data set using a sample list --- """
  threads: 4
  shell:
    """
    bcftools view -S {input.ids} {input.vcf} -Oz -o {output.vcf} --threads {threads}    &&
    bcftools index -f --tbi {output.vcf} 2> {log}
    """


rule vcf_update_tags:
  input:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68.vcf.gz'
  output:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated.vcf.gz'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated.log'
  message: """ update tags """
  shell:
    """
    bcftools +fill-tags {input.vcf} -Oz -o {output.vcf} -- -t AN,AC,AF,MAF &&
    bcftools index -f --tbi {output.vcf} 2> {log}
    """
