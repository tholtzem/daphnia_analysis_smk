rule Dtrios:
  input:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_73.vcf.gz',
    species_sets = 'analyses/pixy/list/daphnia_popfile_73.tsv',
    #nwk = 'analyses/Dsuite/.nwk',
    #trios = 'analyses/Dsuite/test_trios.txt'
  output:
    touch('analyses/Dsuite/Dtrios.done')
  log:
  threads: 4
  message: """ -- Calculate Dstats and f4-ratio stats for all possible trios of species --- """
  shell:
    """
    Dsuite=(~/bio/Dsuite/Build/Dsuite)
    $Dsuite Dtrios -c -n daphnia_geneflow {input.vcf} {input.species_sets} --ABBAclustering
    """


rule Dinvestigate:
  input:
    vcf = 'vcf/filtered/daphnia_9pops_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf0.01_PCAoutliersRM.filtered_pruned.vcf.gz',
    species = 'analyses/Dsuite/SETS.txt',
    nwk = 'analyses/Dsuite/.nwk',
    trios = 'analyses/Dsuite/test_trios.txt'
  output:
    touch('analyses/Dsuite/Dtrios.done')
  log:
  threads: 4
  message: """ -- Calculate Dstats and f4-ratio stats for all possible trios of species --- """
  shell:
    """
    Dsuite Dinvestigate [OPTIONS] {input.vcf} {input.species} {input.trios}
    """

rule Fbranch:
  input:
    nwk = 'analyses/Dsuite/.nwk',
    f_G = 'analyses/Dsuite/FVALS_tree.txt'
  output:
    touch('analyses/Dsuite/Dtrios.done')
  log:
  threads: 4
  message: """ -- Calculate Dstats and f4-ratio stats for all possible trios of species --- """
  shell:
    """
    Dsuite Dinvestigate [OPTIONS] {input.vcf} {input.species} {input.tri -c -n daphnia_geneflow"
