rule prune_SNPdistance:
  input:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz'
  output:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned.vcf.gz'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned.log'
  threads: 4
  message: """ --- Prune SNPs that are closer to each other than a minimum distance of 100 bp --- """
  shell:
    """
    bcftools +prune -w 100bp -n 1 -N rand -Oz -o {output.vcf} {input.vcf}
    """


rule vcf2phylip_svdq:
  input:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned.vcf.gz'
  output:
    nexus_vcf = 'analyses/svdq/SNPs/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned.min4.nexus'
  log:
    'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_nexus.log'
  threads: 2
  message: """ --- Convert vcf2nexus using vcf2phylip.py --- """
  shell:
    """
    ./scripts/vcf2phylip.py --input {input.vcf} --nexus &&
    mv bcf/filtered/{wildcards.DPmax}/daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_pruned.min4.nexus {output} 2> {log}
    """


rule svdq_SPECIEStree:
  input:
    nexus_vcf = 'analyses/svdq/SNPs/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned.min4.nexus',
    taxpart = 'analyses/svdq/SNPs/daphnia_taxpartitions.nex',
    batch = 'analyses/svdq/SNPs/PaupBlock_SPECIEStree_{outgroup}.batch'
  output:
    touch('analyses/svdq/SNPs/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}.done')
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
    nexus_vcf = 'analyses/svdq/SNPs/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned.min4.nexus',
    taxpart = 'analyses/svdq/SNPs/daphnia_taxpartitions.nex',
    batch = 'analyses/svdq/SNPs/PaupBlock_LINEAGEtree.batch'
  output:
    'analyses/svdq/SNPs/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_LINEAGEtree_SVDQ_boostrapSTD_rooted.done'
  log:
    'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_pruned_LINEAGEtree_SVDQ_boostrapSTD_rooted.log'
  threads: 4
  message: """ --- Interference of lineage tree using SVDquartet in Paup* --- """
  shell:
    """
    paup4a169_ubuntu64 -n {input.batch} 2> {log}
    """

