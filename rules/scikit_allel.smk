rule vcf2zarr:
  input:
    #vcf = 'vcf/filtered/daphnia_{species}_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf{maf}.vcf.gz'
  output:
    zarr = 'analyses/scikit_allel/data/daphnia_{species}_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf{maf}_zarr/.zgroup'
  log: 'log/daphnia_{species}_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf{maf}_zarr.log'
  message: """--- Converting vcf into zarr format ---"""
  shell:
    """
    python scripts/vcf2zarr.py 
    """

rule alStats:
  input:
    zarr = 'analyses/scikit_allel/data/daphnia_9pops_vars_Q30_DP10_GQ30_alleles.imissRM_lmiss20_maf0.01_zarr/.zgroup'
  output:
    touch("analyses/scikit_allel/al.done")
  message:
    """--- Creating scikit-allel statistics ---"""
  shell:
    """
    python script/al.py {input.zarr}
    """
