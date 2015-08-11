mkdir figs_eff
#root -l Tnp_trig_eff_correction_ToyModel.C+ -q -b >& func_log.txt
root -l Tnp_trig_eff_correction_ToyModel.C+ -q -b | grep Erf >& func_log.txt
