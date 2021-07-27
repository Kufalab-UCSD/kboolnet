config <- list(
    file = "./inconsistent_example.xlsx",
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

