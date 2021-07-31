config <- list(
  file = "../RasRafMekErk_master.xlsx",
	modules = "S6K",
	inputStimuli = c("Erk_p+_S6K_[CTD(T421S424)]", "mTORC1_p+_S6K_[HM(T389)]", "PDK1_p+_S6K_[Tloop(T229)]"),
	inputInhibs = "S6K_p+_*",
	outputs = c("[pS6KT421S424]", "[pIRS1S636S639]", "[prpS6S235S236]", "[pGSK3bS9]"),
	out = "./S6K"
)

