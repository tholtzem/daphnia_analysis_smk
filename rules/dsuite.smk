#rule Dtrios:
#  input:
#    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz',
#    species_sets = 'analyses/Dsuite/daphnia_popfile_{outgroup}.tsv'
#  output:
#    touch('analyses/Dsuite/Dtrios_test/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_test.done')
#  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_test.log'
#  threads: 4
#  message: """ -- Calculate Dstats and f4-ratio stats for all possible trios of species --- """
#  shell:
#    """
#    Dsuite=(~/bio/Dsuite/Build/Dsuite)
#    $Dsuite Dtrios -c -n daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_{wildcards.outgroup}_Dtrios_test {input.vcf} {input.species_sets}
#    """

#rule Dinvestigate:
#  input:
#    vcf = 'vcf/filtered/daphnia_9pops_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf0.01_PCAoutliersRM.filtered_pruned.vcf.gz',
#    species = 'analyses/Dsuite/SETS.txt',
#    nwk = 'analyses/Dsuite/.nwk',
#    trios = 'analyses/Dsuite/test_trios.txt'
#  output:
#    touch('analyses/Dsuite/Dtrios.done')
#  log:
#  threads: 4
#  message: """ -- Calculate Dstats and f4-ratio stats for all possible trios of species --- """
#  shell:
#    """
#    Dsuite Dinvestigate [OPTIONS] {input.vcf} {input.species} {input.trios}
#    """
#

rule reformat_speciestrees_mito:
  input:
    nwk = 'analyses/svdq/mito/refs_mitoPCG_NT_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}.nwk'
  output:
    nwk = 'analyses/Dsuite/mito/refs_mitoPCG_NT_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}_rf.nwk'
  log:
    'log/refs_mitoPCG_NT_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}_rf.log'
  threads: 1
  message: """ --- use 'outgroup' instead of the species name for outgroup  --- """
  shell:
    """
    sed 's/{wildcards.outgroup}/Outgroup/' {input.nwk} > {output.nwk}
    """


rule Dtrios_svdqtree_mito:
  input:
    vcf = 'bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT{GT}_SNPs.vcf.gz'
    species_sets = 'analyses/Dsuite/daphnia_popfile_{outgroup}.tsv',
    nwk = 'analyses/Dsuite/mito/refs_mitoPCG_NT_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}_rf.nwk'
  output:
    f_G = 'analyses/Dsuite/mito/daphnia_popfile_{outgroup}_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_mitosvdq_tree.txt'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_mitosvdq.log'
  threads: 4
  message: """ -- Calculate Dstats and f4-ratio stats for all possible trios of species using mito tree --- """
  shell:
    """
    Dsuite=(~/bio/Dsuite/Build/Dsuite)
    $Dsuite Dtrios -c -t {input.nwk} -n daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_{wildcards.outgroup}_Dtrios_mitosvdq {input.vcf} {input.species_sets} &&
    mv analyses/Dsuite/daphnia_popfile_{wildcards.outgroup}_daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_{wildcards.outgroup}_Dtrios_mitosvdq* analyses/Dsuite/mito/
    """


rule Fbranch_mito:
  input:
    nwk = 'analyses/Dsuite/mito/refs_mitoPCG_NT_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}_rf.nwk',
    f_G = 'analyses/Dsuite/mito/daphnia_popfile_{outgroup}_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_mitosvdq_tree.txt'
  output:
    fbranch = 'analyses/Dsuite/mito/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_mitosvdq_Fbranch.txt'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_mitosvdq_Fbranch.txt'
  threads: 1
  message: """ --- Use fbranch to aid the interpretation of the f4-ratio results --- """
  shell:
    """
    Dsuite Fbranch {input.nwk} {input.f_G} > {output.fbranch} 2> {log}
    """


rule plot_Fbranch_mito:
  input:
    fbranch = 'analyses/Dsuite/mito/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_mitosvdq_Fbranch.txt',
    nwk = 'analyses/Dsuite/mito/refs_mitoPCG_NT_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}_rf.nwk'
  output:
    png = 'analyses/Dsuite/mito/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_mitosvdq_Fbranch.png',
    svg = 'analyses/Dsuite/mito/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_mitosvdq_Fbranch.svg'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_{outgroup}_Dtrios_mitosvdq_Fbranchplot.log'
  threads: 1
  message: """ --- Plot the output of Dsuite fbranch  --- """
  shell:
    """
    dtools.py -n analyses/Dsuite/mito/daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_{wildcards.outgroup}_Dtrios_mitosvdq_Fbranch --outgroup {wildcards.outgroup} {input.fbranch} {input.nwk} --tree-label-size 8 2> {log}
    """



rule reformat_speciestrees_SNPs:
  input:
    nwk = 'analyses/Dsuite/SNP/daphnia_topology{topo}_{outgroup}.nwk'
  output:
    nwk = 'analyses/Dsuite/SNP/daphnia_topology{topo}_{outgroup}_rf.nwk'
  log:
    'log/daphnia_topology{topo}_{outgroup}_rf.log'
  threads: 1
  message: """ --- use 'outgroup' instead of the species name for outgroup  --- """
  shell:
    """
    sed 's/{wildcards.outgroup}/Outgroup/' {input.nwk} > {output.nwk}
    """


rule Dtrios_svdqtree_SNPs:
  input:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz',
    species_sets = 'analyses/Dsuite/daphnia_popfile_{outgroup}.tsv',
    nwk = 'analyses/Dsuite/SNP/daphnia_topology{topo}_{outgroup}_rf.nwk'
  output:
    f_G = 'analyses/Dsuite/SNP/daphnia_popfile_{outgroup}_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_topology{topo}_{outgroup}_Dtrios_svdq_tree.txt'
  log: 'log/daphnia_popfile_{outgroup}_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_topology{topo}_{outgroup}_Dtrios_svdq_tree.log'
  threads: 4
  message: """ -- Calculate Dstats and f4-ratio stats for all possible trios of species using SNP tree --- """
  shell:
    """
    Dsuite=(~/bio/Dsuite/Build/Dsuite)
    $Dsuite Dtrios -c -t {input.nwk} -n daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_topology{wildcards.topo}_{wildcards.outgroup}_Dtrios_svdq {input.vcf} {input.species_sets} 2> {log} &&
    mv analyses/Dsuite/daphnia_popfile_{wildcards.outgroup}_daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_topology{wildcards.topo}_{wildcards.outgroup}_Dtrios_svdq* analyses/Dsuite/SNP/ """


rule Fbranch_SNPs:
  input:
    nwk = 'analyses/Dsuite/SNP/daphnia_topology{topo}_{outgroup}_rf.nwk',
    f_G = 'analyses/Dsuite/SNP/daphnia_popfile_{outgroup}_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_topology{topo}_{outgroup}_Dtrios_svdq_tree.txt'
  output:
    fbranch = 'analyses/Dsuite/SNP/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_topology{topo}_{outgroup}_Dtrios_svdq_Fbranch.txt'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_topology{topo}_{outgroup}_Dtrios_svdq_Fbranch.log'
  threads: 1
  message: """ --- Use fbranch to aid the interpretation of the f4-ratio results --- """
  shell:
    """
    Dsuite Fbranch {input.nwk} {input.f_G} > {output.fbranch} 2> {log}
    """


rule plot_Fbranch:
  input:
    fbranch = 'analyses/Dsuite/SNP/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_topology{topo}_{outgroup}_Dtrios_svdq_Fbranch.txt',
    nwk = 'analyses/Dsuite/SNP/daphnia_topology{topo}_{outgroup}_rf.nwk'
  output:
    png = 'analyses/Dsuite/SNP/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_topology{topo}_{outgroup}_Dtrios_svdq_Fbranch.png',
    svg = 'analyses/Dsuite/SNP/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_topology{topo}_{outgroup}_Dtrios_svdq_Fbranch.svg'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_topology{topo}_{outgroup}_Dtrios_svdq_Fbranch.log'
  threads: 1
  message: """ --- Plot the output of Dsuite fbranch  --- """
  shell:
    """
    dtools.py -n analyses/Dsuite/SNP/daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_topology{wildcards.topo}_{wildcards.outgroup}_Dtrios_svdq_Fbranch --outgroup Outgroup {input.fbranch} {input.nwk} --tree-label-size 8 2> {log}
    """



