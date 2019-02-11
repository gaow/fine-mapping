# workhorse(s) for finemapping

# Module input
# ============
# $data: full data; or
# $sumstats: summary statistics; or / and
# $ld: LD information

# Module output
# =============
# $fitted: for diagnostics
# $posterior: for inference

fit_caviar: fit_caviar.R + \
             R(posterior = finemap_mcaviar(sumstats[1,,]/sumstats[2,,], 
                                            ld, args, prefix=cache))
  @CONF: R_libs = (dplyr, magrittr)
  sumstats: $sumstats
  ld: $ld_file
  args: "-g 0.001 -c 1", "-g 0.001 -c 2", "-g 0.001 -c 3"
  cache: file(CAVIAR)
  $posterior: posterior

fit_finemap(fit_caviar): fit_finemap.R + \
             R(posterior = finemap_mvar(sumstats[1,,] / sumstats[2,,],
                                        ld, N, k,
                                        args, prefix=cache))
  N: $N
  k: NULL
  data: $data
  args: "--n-causal-max 1", "--n-causal-max 2", "--n-causal-max 3"
  cache: file(FM)

fit_dap: fit_dap.py + Python(posterior = dap_batch(data['X'], data['Y'], cache, args))
  data: $data
  args: "-ld_control 0.20 --all"
  cache: file(DAP)
  $posterior: posterior

fit_dap_z: fit_dap.py + Python(posterior = dap_batch_z(sumstats[0,:,:]/sumstats[1,:,:],
                                                       ld, cache, args))
  sumstats: $sumstats
  ld: $ld_file
  args: "-t 4"
  cache: file(DAP)
  $posterior: posterior

fit_susie: fit_susie.R
  # Prior variance of nonzero effects.
  @CONF: R_libs = susieR@stephenslab/susieR
  maxI: 200
  maxL: 10
  null_weight: 0, 0.5, 0.9, 0.95
  prior_var: 0, 0.1, 0.4
  data: $data
  $posterior: posterior
  $fitted: fitted

fit_susie_auto: fit_susie.R
  @CONF: R_libs = susieR@stephenslab/susieR
  data: $data
  prior_var: "auto"
  $posterior: posterior
  $fitted: fitted

fit_susie01(fit_susie):
  maxL: 1

fit_susie10(fit_susie):
  maxL: 15

fit_susie_uniroot(fit_susie):
  V_method: "uniroot"
  prior_var: 0

fit_susie_em(fit_susie):
  V_method: "EM"
  prior_var: 0
