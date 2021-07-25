config <- list(
    # file = "./CXCR4_toy.xlsx",
    driveFile = "https://docs.google.com/spreadsheets/d/1QbbEsxxu9FIMx1kHuXCPrKEAzPpU4xxKI5QxHTA5gEg/edit#gid=981789074",
    runs = list(
      one = list(
        toggle = c('CXCL12')
      ),
      two = list(
        initial = "./VerifyModel_out/modules_initial_vals.csv",
        on = c('PPase'),
        off = c('PI3K_p+_PIP2_*'),
        toggle = c('CXCL12')
      )
    ),
    out = "./VerifyModel_out"
)

