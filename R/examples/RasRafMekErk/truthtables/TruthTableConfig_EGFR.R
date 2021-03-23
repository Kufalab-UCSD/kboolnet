config <- c(
  file = "../RasRafMekErk_master.xlsx",
	modules = "EGFR-Src",
	inputStimuli = "EGF_[*]--0",
	inputInhibs = "Src_p+_*, EGFR_p+_*, Erk_p+_*",
	outputs = "EGFR_[*]-{p}, EGFR_[*]--EGFR_[*], PI3K1A_[*]-{p}, Ras_[*]--SOS_[*]",
	out = "./EGFR"
)

