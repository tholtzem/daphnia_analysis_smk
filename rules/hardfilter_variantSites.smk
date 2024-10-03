rule count_variantSites:
  input:
    '/media/tania/jugglingjay1/DAPHNIA/daphnia_analysis_smk/bcf/concat/daphnia_init_10pops.bcf'
  output:
    'bcf/stats/daphnia_init_10pops_variants_nbrSites.txt'
  log: 'log/daphnia_init_10pops_variants_nbrSites.log'
  threads: 4
  message: """--- Only count 'variant-sites' from raw BCF (excluding 0|0 genotypes) ---"""
  shell:
    """
    bcftools view -H -c1 {input} | wc -l > {output} 2> {log}
    """


rule get_variantSites:
  input:
    '/media/tania/jugglingjay1/DAPHNIA/daphnia_analysis_smk/bcf/concat/daphnia_init_10pops.bcf'
  output:
    '/media/tania/jugglingjay1/DAPHNIA/daphnia_analysis_smk/bcf/concat/daphnia_init_10pops_variants.bcf'
  log: 'log/daphnia_init_10pops_variants.log'
  threads: 4
  message: """--- Only get 'variant-sites' from raw BCF (excluding 0|0 genotypes) ---"""
  shell:
    """
    bcftools view -c1 {input} -Oz -o {output} --threads {threads} 2> {log}
    """


rule vcf_randomsample_variantSites:
  input:
    bcf = '/media/tania/jugglingjay1/DAPHNIA/daphnia_analysis_smk/bcf/concat/daphnia_init_10pops_variants.bcf',
    sites = 'bcf/stats/daphnia_init_10pops_variants_nbrSites.txt'
  output:
    'bcf/stats/daphnia_init_10pops_variants_100Ksubset.bcf'
  log: 'log/daphnia_init_10pops_variants_100Ksubset.log'
  threads: 4
  message: """ --- Randomly subsample 'variant-sites' bcf --- """
  shell:
   """
   bcftools view {input.bcf} | vcfrandomsample -r 0.004 | bgzip -c > {output} 2> {log}
   """



rule print_MQ_DP_variantSites:
  input:
    'bcf/stats/daphnia_init_10pops_variants_100Ksubset.bcf'
  output:
    'bcf/stats/daphnia_init_10pops_variants_100Ksubset_MQ_DP.tsv'
  log:
    'log/daphnia_init_10pops_variants_100Ksubset_MQ_DP.log'
  threads: 4
  message: """--- For each variant site print MQ, DP values ---"""
  shell:
    """
    bcftools query {input} -f'%MQ\t%DP\n' > {output} 2> {log}
    """


rule plot_INFOINFO_MQ_DP_variantSites:
  input:
    arg1 = 'bcf/stats/daphnia_init_10pops_variants_100Ksubset_MQ_DP.tsv'
  output:
    arg3 = 'bcf/stats/daphnia_init_10pops_variants_100Ksubset_MQ_DP.pdf',
    arg4 = 'bcf/stats/daphnia_init_10pops_variants_100Ksubset_MQ_DP.filter.list'
  log: 'log/daphnia_init_10pops_variants_100Ksubset_MQ_DP.log'
  threads: 4
  message: """--- Plot MQ and DP values for variant sites and write stats to tsv file. ---"""
  shell:
    """
    Rscript scripts/plotting_MQ_DP_variantSitesVCFs.R {input.arg1} init_10pops {output.arg3} {output.arg4} 2> {log}
    """


rule HardFilter_variantSites_MQ30_DPminmaxIQR:
  input:
    stats = 'bcf/stats/daphnia_init_10pops_variants_100Ksubset_MQ_DP.tsv',
    filterlist = 'bcf/stats/daphnia_init_10pops_variants_100Ksubset_MQ_DP.filter.list',
    bcf = '/media/tania/jugglingjay1/DAPHNIA/daphnia_analysis_smk/bcf/concat/daphnia_init_10pops_variants.bcf'
  output:
    bcf = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmaxIQR.bcf',
    n_sites = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmaxIQR_nbrSites.txt'
  log: 'log/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmaxIQR.log'
  threads: 4
  message: """--- Hard filter average site-level: maxdepth and mapping quality, only keep reads with MQ>30 --"""
  shell:
    """
    bcftools filter -i 'MQ>30 && INFO/DP<5317' {input.bcf} -Oz -o {output.bcf} --threads {threads} &&
    bcftools view -H {output.bcf} | wc -l > {output.n_sites} 2> {log}
    """


rule HardFilter_variantSites_MQ30_DPminmaxHengLi:
  input:
    stats = 'bcf/stats/daphnia_init_10pops_variants_100Ksubset_MQ_DP.tsv',
    filterlist = 'bcf/stats/daphnia_init_10pops_variants_100Ksubset_MQ_DP.filter.list',
    bcf = '/media/tania/jugglingjay1/DAPHNIA/daphnia_analysis_smk/bcf/concat/daphnia_init_10pops_variants.bcf'
  output:
    bcf = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmaxHengLi.bcf',
    n_sites = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmaxHengLi_nbrSites.txt'
  log: 'log/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmaxHengLi.log'
  threads: 4
  message: """--- Hard filter average site-level: maxdepth and mapping quality, only keep reads with MQ>30 --"""
  shell:
    """
    bcftools filter -i 'MQ>30 && INFO/DP<2839' {input.bcf} -Oz -o {output.bcf} --threads {threads} &&
    bcftools view -H {output.bcf} | wc -l > {output.n_sites} 2> {log}
    """


rule exclude_indels_from_variants:
  input: 
    bcf = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}.bcf',
  output:
    bcf = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites.bcf',
    n_sites = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_nbrSites.txt'
  log: 'log/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites.log'
  threads: 4
  message: """--- Exclude sites in close proximity to indels or complex events. Exclude indels ---"""
  shell:
    """
    bcftools filter --SnpGap 2:indel,other {input.bcf} | bcftools view --exclude-types indels -Ob -o {output.bcf} &&
    bcftools view -H {output.bcf} | wc -l > {output.n_sites} 2> {log}
    """



