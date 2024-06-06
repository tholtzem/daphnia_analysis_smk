
rule pixy:
  input:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated.vcf.gz',
    tbi = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated.vcf.gz.tbi',
    pops = 'analyses/pixy/list/daphnia_popfile_68.tsv'
  output:
    touch('analyses/pixy/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_pixy.done')
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_pixy.log'
  threads: 3
  message: """ --- Calculate summary stats (pi, dxy, fst) from VCF over entire genome --- """
  shell:
    """
    pixy --stats pi fst dxy --vcf {input.vcf} --populations {input.pops} --output_folder analyses/pixy/data/ --output_prefix daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated --window_size 10000 --n_cores {threads} 2> {log}
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
