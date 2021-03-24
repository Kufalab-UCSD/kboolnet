cleanFiles <- function(files) {
  for (file in files) {
    if (file.exists(file)) {
      file.remove(file)
    }
  }
}

callExtractModules <- function(inFile, outFile, modules="", quality=0, args=c()) {
  cleanFiles(outFile)
  path <- paste0(system.file(package="kboolnet"), "/python/extract_modules.py")
  suppressWarnings(stderr <- system2(command = "python3", args = c(path, "--file", inFile,
                                                                   "--modules", paste0('"', paste0(modules, collapse=","), '"'),
                                                                   "--output", outFile, "--quality", quality, args), stderr = TRUE, stdout = ""))

  if (any(grepl("Error", stderr, ignore.case = TRUE)) | !file.exists(outFile)) {
    cat(paste(stderr, "\n"))
    stop("Error during module extraction. Please run extract_modules.py on its own with the -v DEBUG flag.")
  }
}


callRxncon2Reg <- function(inFile, outFile, args=c()) {
  cleanFiles(outFile)
  path <- paste0(system.file(package="kboolnet"), "/python/rxncon2regulatorygraph.py")
  suppressWarnings(stderr <- system2(command = "python3", args = c(path, inFile,
                                                                   "--output", outFile, args), stderr = TRUE, stdout = ""))

  if (any(grepl("Error", stderr, ignore.case = TRUE)) | !file.exists(outFile)) {
    cat(paste(stderr, "\n"))
    stop("Error during module plotting. Please run rxncon2regulatorygraph.py on its own with the -v DEBUG flag.")
  }
}

callRxncon2BNG <- function(inFile, outFile, args=c()) {
  cleanFiles(paste0(outFile, c(".bngl", ".xml", ".rnf")))
  path <- paste0(system.file(package="kboolnet"), "/python/rxncon2bngl.py")
  suppressWarnings(stderr <- system2("python3", args = c(path, inFile, "--output",
                                                         outFile), stderr = TRUE, stdout = ""))
  if (any(grepl("Error", stderr, ignore.case = TRUE))) {
    cat(paste(stderr, "\n"))
    stop("Error during BNG file generation. Please run rxncon2bngl.py on its own with the -v DEBUG flag.")
  }
}

callRxncon2Boolnet <- function(inFile, outFile, args=c()) {
  cleanFiles(outFile)
  path <- paste0(system.file(package="kboolnet"), "/python/rxncon2boolnet.py")
  suppressWarnings(stderr <- system2("python3", args = c(path, inFile, "--output",
                                                       outFile, args), stderr = TRUE, stdout = ""))
  if (any(grepl("Error", stderr, ignore.case = TRUE))) {
    cat(paste(stderr, "\n"))
    stop("Error during BoolNet file generation. Please run rxncon2boolnet.py on its own with the -v DEBUG flag.")
  }
}

callReactionMapping <- function(inFile, outFile, args=c()) {
  cleanFiles(outFile)
  path <- paste0(system.file(package="kboolnet"), "/python/reaction_mapping.py")
  suppressWarnings(stderr <- system2("python3", args = c(path, inFile, "--output",
                                                       outFile, args), stderr = TRUE, stdout = ""))
  if (any(grepl("Error", stderr, ignore.case = TRUE))) {
    cat(paste(stderr, "\n"))
    stop("Error during reaction mapping generation. Please run reaction_mapping.py on its own with the -v DEBUG flag.")
  }

}
