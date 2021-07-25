config <- list(
    driveFile = "https://docs.google.com/spreadsheets/d/1I0RkJC8NvFx2VlLw7GbQny8U1ldpjFXhlSIqeVgFoCs/edit?usp=sharing",
    runs = list(
      one = list(
        toggle = c('Ligand'),
	      off = c('Ligand2_[*]--0')
      ),
      two = list(
        toggle = c('Ligand2'),
	      on = c('Ligand_[*]--0')
      )
    ),
    out = "./VerifyModel_out"
)

