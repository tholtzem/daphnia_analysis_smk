rule vcf_4_pi_GT80:
  input:
    bcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf'
  output:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz',
    tbi = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz.tbi'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.log'
  threads: 4
  message: """ --- Filter out sites with less than 80 % complete GTs in each species  --- """
  shell:
    """
    bcftools filter -e 'F_MISSING > 0.2' {input.bcf} -Oz -o {output.vcf} --threads {threads} 2> {log}
    """


rule indexvcf_4_pi_GT80:
  input:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz'
  output:
    tbi = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz.tbi'
  log: 'log/index_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz.tbilog'
  threads: 4
  message: """ --- Index sites with less than 80 % complete GTs in each species  --- """
  shell:
    """
    bcftools index -f --tbi {input.vcf}
    """



rule pixy_pi:
  input:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz',
    tbi = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz.tbi',
    pops = 'analyses/pixy/list/daphnia_popfile_68_{species}.tsv'
  output:
    touch('analyses/pixy/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs_{bpWindow}bpwindow_pi.done')
  log: 'log/pixy/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs_{bpWindow}bpwindow_pi.log'
  threads: 4
  message: """ --- Calculate pi from VCF over entire genome --- """
  shell:
    """
    pixy --stats pi --vcf {input.vcf} --populations {input.pops} --output_folder analyses/pixy/data/ --output_prefix daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites_final_80GTs_{wildcards.bpWindow}bpwindow --window_size {wildcards.bpWindow} --n_cores {threads} 2> {log}
    """


rule vcf_4_pixy_GT80:
  input:
    cucullata = 'bcf/filtered/{DPmax}/daphnia_cucullata_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz',
    curvirostris = 'bcf/filtered/{DPmax}/daphnia_curvirostris_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz',
    dentifera = 'bcf/filtered/{DPmax}/daphnia_dentifera_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz',
    galeata= 'bcf/filtered/{DPmax}/daphnia_galeata_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz',
    lacustris = 'bcf/filtered/{DPmax}/daphnia_lacustris_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz',
    longispinaFIN = 'bcf/filtered/{DPmax}/daphnia_longispinaFIN_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz',
    longispina = 'bcf/filtered/{DPmax}/daphnia_longispina_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz',
    mendotae = 'bcf/filtered/{DPmax}/daphnia_mendotae_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz',
    umbra = 'bcf/filtered/{DPmax}/daphnia_umbra_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz',
    zschokkei = 'bcf/filtered/{DPmax}/daphnia_zschokkei_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs.vcf.gz'
  output:
    vcf = 'bcf/merged/{DPmax}/daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs_merged.vcf.gz',
    tbi = 'bcf/merged/{DPmax}/daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs_merged.vcf.gz.tbi'
  log: 'log/daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs_merged.vcf.log'

  threads: 4
  message: """ --- Merge vcf with 80 % complete GTs in each species 4 pixy --- """
  shell:
    """
    bcftools merge {input.curvirostris} {input.cucullata} {input.dentifera} {input.galeata} {input.lacustris} {input.longispina} {input.longispinaFIN} {input.mendotae} {input.umbra} {input.zschokkei} -Oz -o {output.vcf} --threads {threads} 2> {log} &&
    bcftools index -f --tbi {output.vcf}
    """


rule pixy_dxy_Fst:
  input:
    vcf = 'bcf/merged/{DPmax}/daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs_merged.vcf.gz',
    tbi = 'bcf/merged/{DPmax}/daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs_merged.vcf.gz.tbi',
    pops = 'analyses/pixy/list/daphnia_popfile_68.tsv'
  output:
    touch('analyses/pixy/daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs_merged_{bpWindow}bpwindow_pixy.done')
  log: 'log/pixy/daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites_final_80GTs_merged_{bpWindow}bpwindow_pixy.log'
  threads: 4
  message: """ --- Calculate summary stats (dxy, fst) from VCF over entire genome --- """
  shell:
    """
    pixy --stats pi fst dxy --fst_type hudson --vcf {input.vcf} --populations {input.pops} --output_folder analyses/pixy/data/ --output_prefix daphnia_init_10pops_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites_final_80GTs_merged_{wildcards.bpWindow}bpwindow --window_size {wildcards.bpWindow} --n_cores {threads} 2> {log}
    """


#rule listChroms_pixyplot:
#  input:
#    touched = 'analyses/pixy/daphnia_allChroms_pixy.done',
#    pixy_file = 'analyses/pixy/data/daphnia_dxy.txt'
#  output:
#    'analyses/pixy/list/Chrom.list'
#  log: 'log/daphnia_allChroms_pixy.log'
#  threads: 4
#  message: """ --- Get a list of chromosoms for pixy (version sorting) --- """
#  shell:
#    """
#    cat {input.pixy_file} | cut -f3 | uniq | sort -V > {output} 2> {log}
#    """
#
#
#rule plot_pixy_chroms:
#  input:
#    'analyses/pixy/list/Chrom.list'
#  output:
#    'analyses/pixy/plots/pixy_{chrompixy}.png'
#  log: 'log/plot_pixy_{chrompixy}.log'
#  threads: 4
#  message: """ --- Calculate summary stats (pi, dxy, fst) from VCF over entire genome --- """
#  shell:
#    """
#    Rscript scripts/plot_pixy.R {wildcards.chrompixy} {output} 
#    """

#rule plot_pixy:
#  input:
#    'analyses/pixy/list/Chrom.list'
#  output:
#    'analyses/pixy/plots/pixy_{chrompixy}.png'
#  log: 'log/plot_pixy_{chrompixy}.log'
#  threads: 4
#  message: """ --- Calculate summary stats (pi, dxy, fst) from VCF over entire genome --- """
#  shell:
#    """
#    Rscript scripts/plot_pixy.R dgal1 {output} 
#    """
