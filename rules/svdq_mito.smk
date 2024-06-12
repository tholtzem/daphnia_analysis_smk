rule svdq_SPECIEStree:
  input:
    nexus_vcf = 'analyses/svdq/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned.min4.nexus',
    taxpart = 'analyses/svdq/daphnia_taxpartitions.nex',
    batch = 'analyses/svdq/PaupBlock_SPECIEStree_{outgroup}.batch'
  output:
    touch('analyses/svdq/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}.done')
  log:
    'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}.log'
  threads: 4
  message: """ --- Interference of species tree using SVDquartet in Paup* --- """
  shell:
    """
    paup4a169_ubuntu64 -n {input.batch} 2> {log}
    """


rule svdq_LINEAGEtree:
  input:
    nexus_vcf = 'analyses/svdq/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned.min4.nexus',
    taxpart = 'analyses/svdq/daphnia_taxpartitions.nex',
    batch = 'analyses/svdq/PaupBlock_LINEAGEtree.batch'
  output:
    'analyses/svdq/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_LINEAGEtree_SVDQ_boostrapSTD_rooted.done'
  log:
    'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_LINEAGEtree_SVDQ_boostrapSTD_rooted.log'
  threads: 4
  message: """ --- Interference of lineage tree using SVDquartet in Paup* --- """
  shell:
    """
    paup4a169_ubuntu64 -n {input.batch} 2> {log}
    """

