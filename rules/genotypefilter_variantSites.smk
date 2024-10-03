
rule filter_GTs_4_variants:
  input:
    bcf = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites.bcf'
  output:
    bcf = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf',
    n_sites = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_nbrSites.txt'
  log: 'log/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.log'
  threads: 4
  message: """--- Select biallelic variants and filter for minimum site depth, genotype quality, and QUAL filter ---"""
  shell:
    """
    bcftools view -m2 -M2 {input.bcf} | bcftools filter -S . -e 'FMT/DP<=10 | FMT/GQ<=30 | QUAL<=30' -Ob -o {output.bcf} &&
    bcftools view -H {output.bcf} | wc -l > {output.n_sites} 2> {log}
    """


rule imiss_4_variants:
  input:
    'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf'
  output:
    'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.imiss'
  log: 'log/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.imiss.log'
  threads: 4
  message: """--- Identify individuals with a high amount of missing data ---"""
  shell:
    """
    bcftools stats -s - {input} | grep -E ^PSC | cut -f3,14 > {output} 2> {log}
    """


rule plot_imiss_4_variants:
  input:
    bcf = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf',
    n_sites = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_nbrSites.txt',
    args1 = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.imiss'
  output:
    args3 = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.imiss.pdf',
    args4 = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.imiss.tsv'
  log: 'log/plot_daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.imisslog'
  threads: 4
  message: """--- Identify individuals with a high amount of missing data ---"""
  shell:
    """
    args2=$(cat {input.n_sites})
    Rscript scripts/plot_imiss.R {input.args1} $args2 {output.args3} {output.args4} 2> {log}
    """


rule subset_species_variants:
  input:
    imiss = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.imiss.pdf',
    ids = 'analyses/pixy/list/daphnia_samplelist_68_{species}.tsv',
    bcf = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf'
  output:
    bcf = 'bcf/variants/{DPmax}/daphnia_{species}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf',
    n_sites = 'bcf/variants/{DPmax}/daphnia_{species}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_nbrSites.txt'
  log:
    'log/subset_daphnia_{species}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.log'
  message: """ --- Subset bcf using a sample list for each species --- """
  threads: 4
  shell:
    """
    bcftools view -S {input.ids} {input.bcf} -Ob -o {output.bcf} --threads {threads} &&
    bcftools view -H {output.bcf} | wc -l > {output.n_sites} &&
    bcftools index -f --csi {output.bcf} 2> {log}
    """


rule isec_common_variants_GT80:
  input:
    cucullata = 'bcf/variants/{DPmax}/daphnia_cucullata_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf',
    curvirostris = 'bcf/variants/{DPmax}/daphnia_curvirostris_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf',
    dentifera = 'bcf/variants/{DPmax}/daphnia_dentifera_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf',
    galeata= 'bcf/variants/{DPmax}/daphnia_galeata_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf',
    lacustris = 'bcf/variants/{DPmax}/daphnia_lacustris_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf',
    longispinaFIN = 'bcf/variants/{DPmax}/daphnia_longispinaFIN_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf',
    longispina = 'bcf/variants/{DPmax}/daphnia_longispina_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf',
    mendotae = 'bcf/variants/{DPmax}/daphnia_mendotae_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf',
    umbra = 'bcf/variants/{DPmax}/daphnia_umbra_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf',
    zschokkei = 'bcf/variants/{DPmax}/daphnia_zschokkei_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf'
  output:
    set1 = 'bcf/variants/daphnia_set1_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_isec80.tsv',
    set2 = 'bcf/variants/daphnia_set2_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_isec80.tsv',
    set3 = 'bcf/variants/daphnia_set3_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_isec80.tsv',
    set4 = 'bcf/variants/daphnia_set4_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_isec80.tsv'
  log:
    set1 = 'log/daphnia_set1_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_isec80.log',
    set2 = 'log/daphnia_set2_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_isec80.log',
    set3 = 'log/daphnia_set3_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_isec80.log',
    set4 = 'log/daphnia_set4_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_isec80.log'
  message: """ --- List common intersections (variants) for different sets, only keep sites with at least 50 % of GTs in each species --- """
  threads: 4
  shell:
    """
    bcftools isec -e 'F_MISSING > 0.2' --threads {threads} -n=10 {input.curvirostris} {input.cucullata} {input.dentifera} {input.galeata} {input.lacustris} {input.longispina} {input.longispinaFIN} {input.mendotae} {input.umbra} {input.zschokkei} | cut -f1,2 > {output.set1} 2> {log.set1} &&
    bcftools isec -e 'F_MISSING > 0.2' --threads {threads} -n=9 {input.cucullata} {input.dentifera} {input.galeata} {input.lacustris} {input.longispina} {input.longispinaFIN} {input.mendotae} {input.umbra} {input.zschokkei} | cut -f1,2 > {output.set2} 2> {log.set2} &&
    bcftools isec -e 'F_MISSING > 0.2' --threads {threads} -n=8 {input.curvirostris} {input.cucullata} {input.dentifera} {input.galeata} {input.longispina} {input.longispinaFIN} {input.mendotae} {input.zschokkei} | cut -f1,2 > {output.set3} 2> {log.set3} &&
    bcftools isec -e 'F_MISSING > 0.2' --threads {threads} -n=7 {input.cucullata} {input.dentifera} {input.galeata} {input.longispina} {input.longispinaFIN} {input.mendotae} {input.zschokkei} | cut -f1,2 > {output.set4} 2> {log.set4}
    """


rule subset_variants_4_bcf_sets_GT80:
  input:
    bcf = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf',
    ids = 'analyses/pixy/list/daphnia_samplelist_68.tsv',
    sites = 'bcf/variants/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_isec80.tsv',
  output:
    bcf = 'bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT80.bcf'
  log: 'log/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT80.log'
  threads: 4
  message: """ --- Subset bcf using different variant sets having at least 80 % complete GTs in each species of the respective set --- """
  shell:
    """
    bcftools index --csi {input.bcf} &&
    bcftools view -S {input.ids} -R {input.sites} {input.bcf} -Ob -o {output.bcf} --threads {threads} 2> {log}
    """


rule subset_variants_4_bcf_sets_GT100:
  input:
    ids = 'analyses/pixy/list/daphnia_samplelist_{sites_set}.tsv',
    bcf = 'bcf/variants/{DPmax}/daphnia_init_10pops_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30.bcf',
  output:
    bcf = 'bcf/variants/{DPmax}/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT100.bcf'
  log:
    'log/daphnia_{sites_set}_variants_MQ{MQ}_DPminmax{DPmax}_sites_DP10_GQ30_GT100.log'
  message: """ --- Select variants with at least 100% complete GTs in each individual (and all species) for different sets -- """
  threads: 4
  shell:
    """
    bcftools view -S {input.ids} {input.bcf}  | bcftools filter -e 'F_MISSING > 0' -Ob -o {output.bcf} --threads {threads} &&
    bcftools index -f --csi {output.bcf} 2> {log}
    """
