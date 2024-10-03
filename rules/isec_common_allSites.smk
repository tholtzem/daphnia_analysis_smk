rule isec_common_GT100:
  input:
    curvirostris = 'bcf/filtered/{DPmax}/daphnia_curvirostris_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    cucullata = 'bcf/filtered/{DPmax}/daphnia_cucullata_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    dentifera = 'bcf/filtered/{DPmax}/daphnia_dentifera_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    galeata = 'bcf/filtered/{DPmax}/daphnia_galeata_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    lacustris = 'bcf/filtered/{DPmax}/daphnia_lacustris_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    longispinaFIN = 'bcf/filtered/{DPmax}/daphnia_longispinaFIN_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    longispina = 'bcf/filtered/{DPmax}/daphnia_longispina_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    mendotae = 'bcf/filtered/{DPmax}/daphnia_mendotae_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    umbra = 'bcf/filtered/{DPmax}/daphnia_umbra_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    zschokkei = 'bcf/filtered/{DPmax}/daphnia_zschokkei_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf'
  output:
    set1 = 'bcf/isec/set1_daphnia_MQ{MQ}_DPminmax{DPmax}_sites_final_isec100.tsv',
    set2 = 'bcf/isec/set2_daphnia_MQ{MQ}_DPminmax{DPmax}_sites_final_isec100.tsv',
    set3 = 'bcf/isec/set3_daphnia_MQ{MQ}_DPminmax{DPmax}_sites_final_isec100.tsv'
  log:
    set1 = 'log/isec/set1_daphnia_MQ{MQ}_DPminmax{DPmax}_sites_final_isec100.log',
    set2 = 'log/isec/set2_daphnia_MQ{MQ}_DPminmax{DPmax}_sites_final_isec100.log',
    set3 = 'log/isec/set3_daphnia_MQ{MQ}_DPminmax{DPmax}_sites_final_isec100.log'
  message: """ --- List common intersections (variants) for different sets, only keep sites with at least 100 % of GTs in each species --- """
  threads: 4
  shell:
    """
    bcftools isec -e 'F_MISSING > 0' -n=10 {input.curvirostris} {input.cucullata} {input.dentifera} {input.galeata} {input.lacustris} {input.longispina} {input.longispinaFIN} {input.mendotae} {input.umbra} {input.zschokkei} > {output.set1} --threads {threads} 2> {log.set1} &&
    bcftools isec -e 'F_MISSING > 0' -n=9 {input.cucullata} {input.dentifera} {input.galeata} {input.lacustris} {input.longispina} {input.longispinaFIN} {input.mendotae} {input.umbra} {input.zschokkei} > {output.set2} --threads {threads} 2> {log.set2} &&
    bcftools isec -e 'F_MISSING > 0' -n=8 {input.curvirostris} {input.cucullata} {input.dentifera} {input.galeata} {input.longispina} {input.longispinaFIN} {input.mendotae} {input.zschokkei} > {output.set3} --threads {threads} 2> {log.set3}
    """


rule isec_common_GT80:
  input:
    curvirostris = 'bcf/filtered/{DPmax}/daphnia_curvirostris_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    cucullata = 'bcf/filtered/{DPmax}/daphnia_cucullata_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    dentifera = 'bcf/filtered/{DPmax}/daphnia_dentifera_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    galeata = 'bcf/filtered/{DPmax}/daphnia_galeata_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    lacustris = 'bcf/filtered/{DPmax}/daphnia_lacustris_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    longispinaFIN = 'bcf/filtered/{DPmax}/daphnia_longispinaFIN_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    longispina = 'bcf/filtered/{DPmax}/daphnia_longispina_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    mendotae = 'bcf/filtered/{DPmax}/daphnia_mendotae_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    umbra = 'bcf/filtered/{DPmax}/daphnia_umbra_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf',
    zschokkei = 'bcf/filtered/{DPmax}/daphnia_zschokkei_MQ{MQ}_DPminmax{DPmax}_sites_final.bcf'
  output:
    set1 = 'bcf/isec/set1_daphnia_MQ{MQ}_DPminmax{DPmax}_sites_final_isec80.tsv',
    #set2 = 'bcf/isec/set2_daphnia_MQ{MQ}_DPminmax{DPmax}_sites_final_isec80.tsv',
    #set3 = 'bcf/isec/set3_daphnia_MQ{MQ}_DPminmax{DPmax}_sites_final_isec80.tsv'
  log:
    set1 = 'log/isec/set1_daphnia_MQ{MQ}_DPminmax{DPmax}_sites_final_isec80.log',
    #set2 = 'log/isec/set2_daphnia_MQ{MQ}_DPminmax{DPmax}_sites_final_isec80.log',
    #set3 = 'log/isec/set3_daphnia_MQ{MQ}_DPminmax{DPmax}_sites_final_isec80.log'
  message: """ --- List common intersections (variants) for different sets, only keep sites with at least 80 % of GTs in each species --- """
  threads: 4
  shell:
    """
    bcftools isec -e 'F_MISSING > 20' -n=10 {input.curvirostris} {input.cucullata} {input.dentifera} {input.galeata} {input.lacustris} {input.longispina} {input.longispinaFIN} {input.mendotae} {input.umbra} {input.zschokkei} > {output.set1} --threads {threads} 2> {log.set1}
    """

    #&& bcftools isec -e 'F_MISSING > 20' -n=9 {input.cucullata} {input.dentifera} {input.galeata} {input.lacustris} {input.longispina} {input.longispinaFIN} {input.mendotae} {input.umbra} {input.zschokkei} > {output.set2} --threads {threads} 2> {log.set2} &&
    #bcftools isec -e 'F_MISSING > 20' -n=8 {input.curvirostris} {input.cucullata} {input.dentifera} {input.galeata} {input.longispina} {input.longispinaFIN} {input.mendotae} {input.zschokkei} > {output.set3} --threads {threads} 2> {log.set3}
