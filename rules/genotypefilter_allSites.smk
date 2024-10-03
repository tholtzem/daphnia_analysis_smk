rule get_Invariant:
  input:
    'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites.bcf'
  output:
    bcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_invariant.bcf',
    n_sites = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_invariant_nbrSites.txt'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_invariant.log'
  threads: 4
  message: """ --- Select invariants using max allele count equals (-C0) and filter for minimum site depth --- """
  shell:
    """
    bcftools view -C0 -M2 {input} | bcftools filter -S . -e 'FMT/DP<=10 | FMT/GQ<=30' -Ob -o {output.bcf} &&
    bcftools view -H {output.bcf} | wc -l > {output.n_sites} 2> {log}
    """


rule get_Variant:
  input:
    'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites.bcf'
  output:
    bcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_variant.bcf',
    n_sites = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_variant_nbrSites.txt'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_variant.log'
  threads: 4
  message: """--- Select biallelic variant sites (excluding 0|0 genotypes, i.e. AC=0) and filter for minimum site depth, genotype quality, and QUAL filter ---"""
  shell:
    """
    bcftools view -c1 {input} -m 2 -M 2 | bcftools filter -S . -e 'FMT/DP<=10 | FMT/GQ<=30 | QUAL<=30' -Ob -o {output.bcf} &&
    bcftools view -H {output.bcf} | wc -l > {output.n_sites} 2> {log}
    """


rule finalBCF:
  input:
    invariant = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_invariant.bcf',
    variant = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_variant.bcf'
  output:
    bcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    n_sites = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final_nbrSites.txt'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final.log'
  threads: 4
  message: """ --- Combine variantVCF and invariantVCF --- """
  shell:
    """
    tabix -p bcf -f {input.invariant} &&
    tabix -p bcf -f {input.variant} &&
    bcftools concat --allow-overlaps {input.invariant} {input.variant} -Ob -o {output.bcf} --threads {threads} &&
    bcftools view -H {output.bcf} | wc -l > {output.n_sites} 2> {log}
    """


rule imiss:
  input:
    'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf'
  output:
    'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final.imiss'
  log: 'log/aphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final.imiss.log'
  threads: 4
  message: """--- Identify individuals with a high amount of missing data ---"""
  shell:
    """
    bcftools stats -s - {input} | grep -E ^PSC | cut -f3,14 > {output} 2> {log}
    """


rule plot_imiss:
  input:
    bcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    args1 = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final.imiss'
  output:
    args3 = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final.imiss.pdf',
    args4 = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final.imiss.tsv'
  log: 'log/plot_imiss_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final.imiss.log'
  threads: 4
  message: """--- Identify individuals with a high amount of missing data ---"""
  shell:
    """
    args2=$(bcftools view -H {input.bcf} | wc -l)
    Rscript scripts/plot_imiss.R {input.args1} $args2 {output.args3} {output.args4} 2> {log}
    """


rule select_samples_68_lmiss80:
  input:
    ids = 'analyses/pixy/list/daphnia_samplelist_68.tsv',
    bcf = 'bcf/filtered/{DPmax}/daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
  output:
    bcf = 'bcf/filtered/{DPmax}/daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites_final_68_sites80.bcf'
  log:
    'log/select_samples_daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites_final_68_sites80.log'
  message: """ --- Select individuals to include in the data set using a sample list --- """
  threads: 4
  shell:
    """
    bcftools view -S {input.ids} {input.bcf}  | bcftools filter -e 'F_MISSING > 0.2 -Ob -o {output.bcf} --threads {threads} &&
    bcftools index -f --csi {output.bcf} 2> {log}
    """


rule subset_species:
  input:
    ids = 'analyses/pixy/list/daphnia_samplelist_68_{species}.tsv',
    bcf = 'bcf/filtered/{DPmax}/daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf'
  output:
    bcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf'
  log:
    'log/subset_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final.log'
  message: """ --- Select individuals to include in the data set using a sample list --- """
  threads: 4
  shell:
    """
    bcftools view -S {input.ids} {input.bcf} -Ob -o {output.bcf} --threads {threads} &&
    bcftools index -f --csi {output.bcf} 2> {log}
    """


#rule bcf2vcf:
#  input:
#    'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final.bcf'
#  output:
#    'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final.vcf.gz'
#  log:
#    'log/LC/BCF2VCF_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final.log'
#  message: """ --- Convert bcf 2 uncompressed vcf for downstream analysis --- """
#  threads: 4
#  shell:
#    """
#    bcftools index -f --csi {input} --threads {threads} 2> {log} &&
#    bcftools convert -Oz -o {output} {input} --threads {threads} 2> {log}
#    """





