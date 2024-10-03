rule beagle_phase:
  input:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz'
  output:
    vcf = 'bcf/phased/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phased.vcf.gz'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_phased.vcf.log'
  message: """ Infer phase using BEAGLE """
  shell:
    """
    beagle gt={input.vcf} out=bcf/phased/daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs.phased impute=true nthreads={threads} window=10000 overlap=1000 gprobs=false 2> {log}
    """


rule vcf2geno:
  input:
    vcf = 'bcf/phased/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phased.vcf.gz'
  output:
    geno = 'bcf/geno/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.geno.gz'
  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.geno.log'
  message: """ Create geno file using the script from https://github.com/simonhmartin/genomics_general/tree/master/VCF_processing """
  shell:
    """
    python ~/bio/genomics_general/VCF_processing/parseVCF.py -i {input.vcf} | gzip > {output.geno} 2> {log}
    """


#rule PhyML_NJ:
#  input:
#    geno = 'bcf/geno/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.geno.gz'
#  output:
#    'analyses/phyml/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_bionj.{N}w_{window}.trees.gz'
#  log: 'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_bionj.{N}w_{window}.log'
#  threads: 4
#  message: """ Get neighbour joining trees for snp windows using the script from https://github.com/simonhmartin/genomics_general/tree/master/phylo (to test use --test flag) """
#  shell:
#    """
#    python /home/tania/bio/genomics_general/phylo/phyml_sliding_windows.py --threads {threads} --genoFile {input.geno} --prefix analyses/phyml/daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs.phyml_bionj.{wildcards.N}w_{wildcards.window} --windSize {wildcards.N} --windType {wildcards.window} --model GTR --optimise n --phyml /usr/bin/phyml-mpi
#    """


rule PhyML_NJ:
  input:
    geno = 'bcf/geno/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.geno.gz',
    ind = 'analyses/phyml/daphnia_{sets}_{outgroup}.tsv'
  output:
    'analyses/phyml/{sets}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_bionj.{N}w_{window}_{outgroup}.trees.gz'
  log: 'log/{sets}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_bionj.{N}w_{window}_{outgroup}.log'
  threads: 4
  message: """ Get neighbour joining trees for snp windows using the script from https://github.com/simonhmartin/genomics_general/tree/master/phylo (to test use --test flag) """
  shell:
    """
    IND=$(cat {input.ind} | cut -f1 | paste -sd ',')

    python /home/tania/bio/genomics_general/phylo/phyml_sliding_windows.py --threads {threads} --genoFile {input.geno} --prefix analyses/phyml/{wildcards.sets}/daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs.phyml_bionj.{wildcards.N}w_{wildcards.window}_{wildcards.outgroup} --windSize {wildcards.N} --windType {wildcards.window} --model GTR --optimise n --phyml /usr/bin/phyml-mpi --individuals $IND
    """


rule PhyML_ML:
  input:
    geno = 'bcf/geno/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.geno.gz',
    ind = 'analyses/phyml/daphnia_{sets}_{outgroup}.tsv'
  output:
    'analyses/phyml/{sets}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_ML.{N}w_{window}_{outgroup}.trees.gz'
  log: 'log/{sets}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_ML.{N}w_{window}_{outgroup}.log'
  threads: 4
  message: """ Get ML trees for snp windows using the script from https://github.com/simonhmartin/genomics_general/tree/master/phylo (to test use --test flag) """
  shell:
    """
    IND=$(cat {input.ind} | cut -f1 | paste -sd ',')
    
    python /home/tania/bio/genomics_general/phylo/phyml_sliding_windows.py --threads {threads} --genoFile {input.geno} --prefix analyses/phyml/{wildcards.sets}/daphnia_{wildcards.species}_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs.phyml_ML.{wildcards.N}w_{wildcards.window}_{wildcards.outgroup} --windSize {wildcards.N} --windType {wildcards.window} --model GTR --optimise tl --phyml /usr/bin/phyml-mpi --individuals $IND
    """




#rule twisst:
#  input:
#    #touched = 'analyses/phyml/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_bionj.{N}w_{window}.done',
#    trees = 'analyses/phyml/{sets}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_bionj.{N}w_{window}_{outgroup}.trees.gz',
#    groups = 'analyses/twisst/daphnia_groups_68.tsv'
#  output:
#    weights = 'analyses/twisst/{sets}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_bionj.{N}w_{window}_{pop1}_{pop2}_{pop3}_{outgroup}_weights.csv.gz',
#    topos = 'analyses/twisst/{sets}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_bionj.{N}w_{window}_{pop1}_{pop2}_{pop3}_{outgroup}_topos.nwk'
#  log: 'log/{sets}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_bionj.{N}w_{window}_{pop1}_{pop2}_{pop3}_{outgroup}_twisst.log'
#  threads: 4
#  message: """ Compute topology weightings """
#  shell:
#    """
#    python /home/tania/bio/twisst/twisst.py -t {input.trees} -w {output.weights} --outputTopos {output.topos} --method complete -g {wildcards.pop1} -g {wildcards.pop2} -g {wildcards.pop3} -g {wildcards.outgroup} --outgroup {wildcards.outgroup} --groupsFile {input.groups} 2> {log}
#    """


rule twisst_5pops:
  input:
    trees = 'analyses/phyml/{sets}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_bionj.{N}w_{window}_{outgroup}.trees.gz',
    groups = 'analyses/twisst/daphnia_groups_68.tsv'
  output:
    weights = 'analyses/twisst/{sets}/p5_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_bionj.{N}w_{window}_{pop1}_{pop2}_{pop3}_{pop4}_{outgroup}_weights.csv.gz',
    topos = 'analyses/twisst/{sets}/p5_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_bionj.{N}w_{window}_{pop1}_{pop2}_{pop3}_{pop4}_{outgroup}_topos.nwk'
  log: 'log/{sets}/p5_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_bionj.{N}w_{window}_{pop1}_{pop2}_{pop3}_{pop4}_{outgroup}_twisst.log'
  threads: 1
  message: """ Compute topology weightings """
  shell:
    """
    python /home/tania/bio/twisst/twisst.py -t {input.trees} -w {output.weights} --outputTopos {output.topos} --method complete -g {wildcards.pop1} -g {wildcards.pop2} -g {wildcards.pop3} -g {wildcards.pop4} -g {wildcards.outgroup} --outgroup {wildcards.outgroup} --groupsFile {input.groups} 2> {log}
    """


rule twisst_7pops:
  input:
    trees = 'analyses/phyml/{sets}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_{inf}.{N}w_{window}_{outgroup}.trees.gz',
    groups = 'analyses/twisst/daphnia_groups_68.tsv',
    topos = 'analyses/twisst/daphnia_ldgmcu_input_topos.nwk'
  output:
    weights = 'analyses/twisst/{sets}/fixedtrees_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_{inf}.{N}w_{window}_{pop1}_{pop2}_{pop3}_{pop4}_{pop5}_{pop6}_{outgroup}_weights.csv.gz',
    topos = 'analyses/twisst/{sets}/fixedtrees_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_{inf}.{N}w_{window}_{pop1}_{pop2}_{pop3}_{pop4}_{pop5}_{pop6}_{outgroup}_topos.nwk'
  log: 'log/{sets}/fixedtrees_daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.phyml_{inf}.{N}w_{window}_{pop1}_{pop2}_{pop3}_{pop4}_{pop5}_{pop6}_{outgroup}_twisst.log'
  threads: 1
  message: """ Compute topology weightings, but only look at the topologies you are specifically interested in (option: --inputTopos) """
  shell:
    """
    python /home/tania/bio/twisst/twisst.py -t {input.trees} --inputTopos {input.topos} -w {output.weights} --outputTopos {output.topos} --method complete -g {wildcards.pop1} -g {wildcards.pop2} -g {wildcards.pop3} -g {wildcards.pop4} -g {wildcards.pop5} -g {wildcards.pop6} -g {wildcards.outgroup} --outgroup {wildcards.outgroup} --groupsFile {input.groups} 2> {log}
    """
