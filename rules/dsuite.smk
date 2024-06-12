rule reformat_speciestrees:
  input:
    nwk = 'analyses/svdq/nquartets100000/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}.nwk'
  output:
    nwk = 'analyses/Dsuite/trees/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}_rf.nwk'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}_rf.log'
  threads: 4
  message: """ --- use 'outgroup' instead of the species name for outgroup  --- """
  shell:
    """
    sed 's/{wildcards.outgroup}/Outgroup/' {input.nwk} > {output.nwk}
    """


rule Dtrios_test:
  input:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz',
    species_sets = 'analyses/Dsuite/daphnia_popfile_{outgroup}.tsv'
  output:
    touch('analyses/Dsuite/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_test.done')
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_test.log'
  threads: 4
  message: """ -- Calculate Dstats and f4-ratio stats for all possible trios of species --- """
  shell:
    """
    Dsuite=(~/bio/Dsuite/Build/Dsuite)
    $Dsuite Dtrios -c -n daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_{wildcards.outgroup}_Dtrios_test {input.vcf} {input.species_sets}
    """
 #--ABBAclustering


rule Dtrios_svdqtree:
  input:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz',
    species_sets = 'analyses/Dsuite/daphnia_popfile_{outgroup}.tsv',
    nwk = 'analyses/Dsuite/trees/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}_rf.nwk'
    #trios = 'analyses/Dsuite/test_trios.txt'
  output:
    touch('analyses/Dsuite/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_svdq.done')
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_svdq.log'
  threads: 4
  message: """ -- Calculate Dstats and f4-ratio stats for all possible trios of species --- """
  shell:
    """
    Dsuite=(~/bio/Dsuite/Build/Dsuite)
    $Dsuite Dtrios -c -t {input.nwk} -n daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_{wildcards.outgroup}_Dtrios_svdq {input.vcf} {input.species_sets} --ABBAclustering
    """


rule Dtrios_mitotree:
  input:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz',
    species_sets = 'analyses/Dsuite/daphnia_popfile_{outgroup}.tsv',
    nwk = 'analyses/Dsuite/trees/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}_rf.nwk'
    #trios = 'analyses/Dsuite/test_trios.txt'
  output:
    touch('analyses/Dsuite/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_mito.done')
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_mito.log'
  threads: 4
  message: """ -- Calculate Dstats and f4-ratio stats for all possible trios of species --- """
  shell:
    """
    Dsuite=(~/bio/Dsuite/Build/Dsuite)
    $Dsuite Dtrios -c -t {input.nwk} -n daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_{wildcards.outgroup}_Dtrios_mito {input.vcf} {input.species_sets}
    """

rule Fbranch:
  input:
    nwk = 'analyses/Dsuite/trees/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}_rf.nwk',
    touched = 'analyses/Dsuite/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_svdq.done',
    f_G = 'analyses/Dsuite/daphnia_popfile_{outgroup}_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_svdq_tree.txt'
  output:
    fbranch = 'analyses/Dsuite/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_svdq_Fbranch.txt'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_svdq_Fbranch.txt'
  threads: 4
  message: """ -- Calculate Dstats and f4-ratio stats for all possible trios of species --- """
  shell:
    """
    Dsuite Fbranch {input.nwk} {input.f_G} > {output.fbranch} 2> {log}
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
