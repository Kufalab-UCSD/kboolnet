config <- c(
  file = "../RasRafMekErk_master.xlsx",
	modules = "PI3K-Akt",
	inputStimuli = c("EGF_[*]--0", "MutantPI3K"),
	inputInhibs = "PI3K*_syn_PIP3",
	outputs = c("Akt_[*]-{p}, GSK3b_[*]-{p}", "Erk_[*]-{p}", "Ras_[*]-{gtp}", "EIF4EBP1_[*]-{p}", "Rheb_[*]-{gtp}"),
	out = "./PI3K"
)

