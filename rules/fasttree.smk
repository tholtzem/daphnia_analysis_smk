rule vcf2phylip_fasttree:
  input:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz'
  output:
    fasta = 'analyses/fasttree/SNPs/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.min4.fasta'
  log:
    'log/fasttree_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.min4.log'
  threads: 2
  message: """ --- Convert vcf2fasta and resolve IUPAC ambiguities using vcf2phylip.py --- """
  shell:
    """
    ./scripts/vcf2phylip.py --input {input.vcf} -p -f -r &&
    mv bcf/filtered/{wildcards.DPmax}/daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs.min4.* analyses/fasttree/SNPs/ 2> {log}
    """


rule fasttree:
  input:
    fasta = 'analyses/fasttree/SNPs/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.min4.fasta'
  output:
    nwk = 'analyses/fasttree/SNPs/fasttree_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.min4.nwk'
  log:
    'log/fasttree_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.min4.log'
  threads: 4
  message: """ --- Convert vcf2fasta and resolve IUPAC ambiguities using vcf2phylip.py --- """
  shell:
    """
    fasttreeMP -nt -gtr -fastest -out {output.nwk} -log {log} {input.fasta}
    """
