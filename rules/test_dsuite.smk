rule split_samples_lmiss80:
  input:
    ids = 'analyses/pixy/list/daphnia_samplelist_68_{species}.tsv',
    vcf = 'bcf/filtered/{DPmax}/daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz'
  output:
    vcf = 'bcf/filtered/{DPmax}/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz'
  log:
    'log/daphnia_{species}_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.log'
  message: """ --- Select individuals to include in the data set using a sample list --- """
  threads: 4
  shell:
    """
    bcftools view -S {input.ids} {input.vcf} | bcftools filter -e 'F_MISSING > 0.2' -Oz -o {output.vcf} --threads {threads} &&
    bcftools index -f --csi {output.bcf} 2> {log}
    """


rule isec_common:
  input:
    curvirostris = 'bcf/filtered/{DPmax}/daphnia_curvirostris_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz',
    cucullata = 'bcf/filtered/{DPmax}/daphnia_cucullata_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz',
    dentifera = 'bcf/filtered/{DPmax}/daphnia_dentifera_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz',
    galeata = 'bcf/filtered/{DPmax}/daphnia_galeata_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz',
    lacustris = 'bcf/filtered/{DPmax}/daphnia_lacustris_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz',
    longispinaFIN = 'bcf/filtered/{DPmax}/daphnia_longispinaFIN_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz',
    longispina = 'bcf/filtered/{DPmax}/daphnia_longispina_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz',
    mendotae = 'bcf/filtered/{DPmax}/daphnia_mendotae_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz',
    umbra = 'bcf/filtered/{DPmax}/daphnia_umbra_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz',
    zschokkei = 'bcf/filtered/{DPmax}/daphnia_zschokkei_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz'
  output:
    set1 = 'bcf/isec/set1_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec100.tsv',
    set2 = 'bcf/isec/set2_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec100.tsv',
    set3 = 'bcf/isec/set3_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec100.tsv'
  log:
    set1 = 'log/set1_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec100.tsv',
    set2 = 'log/set2_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec100.tsv',
    set3 = 'log/set3_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec100.tsv'
  message: """ --- List common intersections (variants) for different sets --- """
  threads: 4
  shell:
    """
    bcftools isec -e 'F_MISSING > 0' -n=10 {input.curvirostris} {input.cucullata} {input.dentifera} {input.galeata} {input.lacustris} {input.longispina} {input.longispinaFIN} {input.mendotae} {input.umbra} {input.zschokkei} > {output.set1} --threads {threads} 2> {log.set1} &&
    bcftools isec -e 'F_MISSING > 0' -n=9 {input.cucullata} {input.dentifera} {input.galeata} {input.lacustris} {input.longispina} {input.longispinaFIN} {input.mendotae} {input.umbra} {input.zschokkei} > {output.set2} --threads {threads} 2> {log.set2} &&
    bcftools isec -e 'F_MISSING > 0' -n=8 {input.curvirostris} {input.cucullata} {input.dentifera} {input.galeata} {input.longispina} {input.longispinaFIN} {input.mendotae} {input.zschokkei} > {output.set3} --threads {threads} 2> {log.set3}
    """


rule reformat_isec_common:
  input:
    sites = 'bcf/isec/{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec80.tsv'
  output:
    temp('bcf/isec/{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec80_rf.tsv')
  log: 'log/{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec80_rf.log'
  threads: 4
  message: """ --- Reformat isec tsv for bcftools  --- """
  shell:
    """
    cat {input.sites} | cut -f1,2 > {output} 2> {log}
    """
 


rule Dtrios_sites:
  input:
    sites = 'bcf/isec/{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec80_rf.tsv',
    vcf = 'bcf/filtered/{DPmax}/daphnia_init_10pops_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs.vcf.gz',
    species_sets = 'analyses/Dsuite/daphnia_popfile_{outgroup}.tsv',
    nwk = 'analyses/Dsuite/mito/refs_mitoPCG_NT_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}_rf.nwk'
  output:
    f_G = 'analyses/Dsuite/test/daphnia_popfile_{outgroup}_{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec80_{outgroup}_Dtrios_mitosvdq_tree.txt'
  log: 'log/test/daphnia_popfile_{outgroup}_{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec80_{outgroup}_Dtrios_mitosvdq_tree.log'
  threads: 4
  message: """ -- Calculate Dstats and f4-ratio stats for all possible trios of species --- """
  shell:
    """
    Dsuite=(~/bio/Dsuite/Build/Dsuite)
    NUMLINES=$(bcftools view -R {input.sites} {input.vcf} | wc -l)
    bcftools view -R {input.sites} {input.vcf} | $Dsuite Dtrios -c -t {input.nwk} -n {wildcards.sites_set}_daphnia_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_isec80_{wildcards.outgroup}_Dtrios_mitosvdq -l $NUMLINES stdin {input.species_sets} &&
    mv analyses/Dsuite/daphnia_popfile_{wildcards.outgroup}_{wildcards.sites_set}_daphnia_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_isec80_{wildcards.outgroup}_Dtrios_mitosvdq* analyses/Dsuite/test/
    """


rule Fbranch_mito:
  input:
    nwk = 'analyses/Dsuite/mito/refs_mitoPCG_NT_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}_rf.nwk',
    f_G = 'analyses/Dsuite/test/daphnia_popfile_{outgroup}_{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec80_{outgroup}_Dtrios_mitosvdq_tree.txt'
  output:
    fbranch = 'analyses/Dsuite/test/{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec80_{outgroup}_Dtrios_mitosvdq_Fbranch.txt'
  log: 'log/{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec80_{outgroup}_Dtrios_mitosvdq_Fbranch.log'
  threads: 1
  message: """ --- Use fbranch to aid the interpretation of the f4-ratio results --- """
  shell:
    """
    Dsuite Fbranch {input.nwk} {input.f_G} > {output.fbranch} 2> {log}
    """


rule plot_Fbranch_mito:
  input:
    fbranch = 'analyses/Dsuite/test/{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec80_{outgroup}_Dtrios_mitosvdq_Fbranch.txt',
    nwk = 'analyses/Dsuite/mito/refs_mitoPCG_NT_SPECIEStree_SVDQ_boostrapSTD_rooted_{outgroup}_rf.nwk'
  output:
    png = 'analyses/Dsuite/test/{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec80_{outgroup}_Dtrios_mitosvdq_Fbranch.png',
    svg = 'analyses/Dsuite/test/{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec80_{outgroup}_Dtrios_mitosvdq_Fbranch.svg'
  log: 'log/{sites_set}_daphnia_MQ{MQ}_DPminmax{DPmax}_sites80_final_68_updated_SNPs_isec80_{outgroup}_Dtrios_mitosvdq_Fbranchplot.log'
  threads: 1
  message: """ --- Plot the output of Dsuite fbranch  --- """
  shell:
    """
    dtools.py -n analyses/Dsuite/test/{wildcards.sites_set}_daphnia_MQ{wildcards.MQ}_DPminmax{wildcards.DPmax}_sites80_final_68_updated_SNPs_isec80_{wildcards.outgroup}_Dtrios_mitosvdq_Fbranch --outgroup Outgroup {input.fbranch} {input.nwk} --tree-label-size 8 2> {log}
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



