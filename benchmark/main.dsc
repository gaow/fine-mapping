#!/usr/bin/env dsc

%include modules/data
%include modules/simulate
%include modules/fit
%include modules/evaluate

DSC:
  define:
    data: liter_data, full_data
  run:
    default: data * summarize_ld * lm_less * get_sumstats * (fit_susie * (score_susie, plot_susie), fit_dap * plot_dap, fit_caviar, fit_finemap)
    susie: liter_data * simple_lm * fit_susie * (score_susie, plot_susie)
  exec_path: modules
  global:
    data_file: data/gtex-manifest.txt
