rule count_allsites:
  input:
    'bcf/concat/daphnia_{species}.bcf'
  output:
    'bcf/stats/daphnia_{species}_nbrSites.txt'
  log: 'log/daphnia_{species}_nbrSites.log'
  threads: 4
  message: """--- Count sites 'all-sites' raw BCF ---"""
  shell:
    """
    bcftools view -H {input} | wc -l > {output} 2> {log}
    """


rule vcf_randomsample:
  input:
    bcf = 'bcf/concat/daphnia_{species}.bcf',
    sites = 'bcf/stats/daphnia_{species}_nbrSites.txt'
  output:
    'bcf/stats/daphnia_{species}_100Ksubset.bcf'
  log: 'log/daphnia_{species}_100Ksubset.log'
  threads: 4
  message: """ --- Randomly subsample vcf --- """
  shell:
   """
   bcftools view {input.bcf} | vcfrandomsample -r 0.0007 | bgzip -c > {output} 2> {log}
   """


rule print_MQ_DP:
  input:
    #'bcf/concat/daphnia_{species}.bcf'
    'bcf/stats/daphnia_{species}_100Ksubset.bcf'
  output:
    #'bcf/stats/daphnia_{species}_MQ_DP.tsv'
    'bcf/stats/daphnia_{species}_100Ksubset_MQ_DP.tsv'
  log:
    #'log/daphnia_{species}_MQ_DP.log'
    'log/daphnia_{species}_100Ksubset_MQ_DP.log'
  threads: 4
  message: """--- For each site print MQ, DP values ---"""
  shell:
    """
    bcftools query {input} -f'%MQ\t%DP\n' > {output} 2> {log}
    """


rule plot_INFOINFO_MQ_DP:
  input:
    arg1 = 'bcf/stats/daphnia_{species}_100Ksubset_MQ_DP.tsv'
    #arg1 = 'bcf/stats/daphnia_{species}_MQ_DP.tsv'
  output:
    arg3 = 'bcf/stats/daphnia_{species}_100Ksubset_MQ_DP_2.pdf',
    arg4 = 'bcf/stats/daphnia_{species}_100Ksubset_MQ_DP_2.filter.list'
  log: 'log/daphnia_{species}_subset_MQ_DP_2.log'
  threads: 4
  message: """--- Plot MQ and DP values and write stats to tsv file. ---"""
  shell:
    """
    Rscript scripts/plotting_MQ_DP_AllSitesVCFs2.R {input.arg1} {wildcards.species} {output.arg3} {output.arg4} 2> {log}
    """


rule HardFilterVCF_MQ40_DPminmax:
  input:
    stats = 'bcf/stats/daphnia_{species}_100Ksubset_MQ_DP.tsv',
    filterlist = 'bcf/stats/daphnia_{species}_100Ksubset_MQ_DP.filter.list',
    bcf = 'bcf/concat/daphnia_{species}.bcf'
  output:
    bcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ40_DPminmax{DPmax}.bcf',
    n_sites = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ40_DPminmax{DPmax}_nbrSites.txt'
  log: 'log/daphnia_{species}_MQ40_DPminmax{DPmax}.log'
  threads: 4
  message: """--- Hard filter average site-level: maxdepth and mapping quality, only keep reads with MQ>40 --"""
  shell:
    """
    bcftools filter -i 'MQ>40 && INFO/DP<5528' -Oz -o {output.bcf} {input.bcf} --threads {threads} &&
    bcftools view -H {output.bcf} | wc -l > {output.n_sites} 2> {log}
    """


rule HardFilterVCF_MQ30_DPminmax:
  input:
    stats = 'bcf/stats/daphnia_{species}_100Ksubset_MQ_DP.tsv',
    filterlist = 'bcf/stats/daphnia_{species}_100Ksubset_MQ_DP.filter.list',
    bcf = 'bcf/concat/daphnia_{species}.bcf'
  output:
    bcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ30_DPminmax{DPmax}.bcf',
    n_sites = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ30_DPminmax{DPmax}_nbrSites.txt'
  log: 'log/daphnia_{species}_MQ30_DPminmax{DPmax}.log'
  threads: 4
  message: """--- Hard filter average site-level: maxdepth and mapping quality, only keep reads with MQ>30 ---"""
  shell:
    """
    bcftools filter -i 'MQ>30 && INFO/DP<5528' -Oz -o {output.bcf} {input.bcf} --threads {threads} &&
    bcftools view -H {output.bcf} | wc -l > {output.n_sites} 2> {log}
    """


rule HardFilterVCF_MQ35_DPminmax:
  input:
    stats = 'bcf/stats/daphnia_{species}_100Ksubset_MQ_DP.tsv',
    filterlist = 'bcf/stats/daphnia_{species}_100Ksubset_MQ_DP.filter.list',
    bcf = 'bcf/concat/daphnia_{species}.bcf'
  output:
    bcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ35_DPminmax{DPmax}.bcf',
    n_sites = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ35_DPminmax{DPmax}_nbrSites.txt'
  log: 'log/daphnia_{species}_MQ35_DPminmax{DPmax}.log'
  threads: 4
  message: """--- Hard filter average site-level: maxdepth and mapping quality, only keep reads with MQ>35 ---"""
  shell:
    """
    bcftools filter -i 'MQ>35 && INFO/DP<5528' -Oz -o {output.bcf} {input.bcf} --threads {threads} &&
    bcftools view -H {output.bcf} | wc -l > {output.n_sites} 2> {log}
    """


rule exclude_indels:
  input: 
    bcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}.bcf'
  output:
    bcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites.bcf',
    n_sites = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_nbrSites.txt'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites.log'
  threads: 4
  message: """--- Exclude sites in close proximity to indels or complex events. Exclude indels ---"""
  shell:
    """
    bcftools filter -g 2:indel,other {input.bcf} | bcftools view --exclude-types indels -Ob -o {output.bcf} &&
    bcftools view -H {output.bcf} | wc -l > {output.n_sites} 2> {log}
    """
