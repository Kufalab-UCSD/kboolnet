cleanFiles <- function(files) {
  for (file in files) {
    if (file.exists(file)) {
      file.remove(file)
    }
  }
}

addQuotes <- function(str) {
  noQuotes <- gsub("^\"", "", str)
  noQuotes <- gsub("\"$", "", noQuotes)

  return(paste0("\"", noQuotes, "\""))
}

callPython <- function(args) {
  config <- loadPackageConfig()
  pythonCommand <- config$value[config$setting == "pythonCommand"]
  return(suppressWarnings(system2(command = pythonCommand, args = args, stderr = TRUE, stdout = "")))
}

loadPackageConfig <- function() {
  configFile <- paste0(system.file(package="kboolnet"), "/config.csv")
  if (!file.exists(configFile)) {
    stop("Package config file could not be found. Please run setupKboolnet() in an interactive terminal first!")
  }

  return(read.csv(file = configFile))
}

callExtractModules <- function(inFile, outFile, modules="", quality=0, args=c()) {
  cleanFiles(outFile)
  path <- paste0(system.file(package="kboolnet"), "/python/extract_modules.py")
  stderr <- callPython(c(addQuotes(path), "--file", addQuotes(inFile), "--modules", paste0('"', paste0(modules, collapse=","), '"'),
                         "--output", addQuotes(outFile), "--quality", quality, args))

  if (any(grepl("Error", stderr, ignore.case = TRUE)) | !file.exists(outFile)) {
    cat(paste(stderr, "\n"))
    stop("Error during module extraction. Please run extract_modules.py on its own with the -v DEBUG flag.")
  }
}

callRxncon2Reg <- function(inFile, outFile, args=c()) {
  cleanFiles(outFile)
  config <- loadPackageConfig()
  path <- paste0(config$value[config$setting == "rxnconDir"], "/rxncon2regulatory.py")
  stderr <-  callPython(c(addQuotes(path), addQuotes(inFile), "--output", addQuotes(outFile), args))


  if (any(grepl("Error", stderr, ignore.case = TRUE))) {
    cat(paste(stderr, "\n"))
    stop("Error during module plotting. Please run rxncon2regulatorygraph.py on its own with the -v DEBUG flag.")
  }
}

callRxncon2BNG <- function(inFile, outFile, args=c()) {
  cleanFiles(addQuotes(paste0(outFile, c(".bngl", ".xml", ".rnf"))))
  config <- loadPackageConfig()
  path <- paste0(config$value[config$setting == "rxnconDir"], "/rxncon2bngl.py")
  stderr <- callPython(c(addQuotes(path), addQuotes(inFile), "--output", addQuotes(outFile)))

  if (any(grepl("Error", stderr, ignore.case = TRUE))) {
    cat(paste(stderr, "\n"))
    stop("Error during BNG file generation. Please run rxncon2bngl.py on its own with the -v DEBUG flag.")
  }
}

callRxncon2Boolnet <- function(inFile, outFile, args=c()) {
  cleanFiles(outFile)
  config <- loadPackageConfig()
  path <- paste0(config$value[config$setting == "rxnconDir"], "/rxncon2boolnet.py")
  stderr <- callPython(c(addQuotes(path), addQuotes(inFile), "--output", addQuotes(outFile), args))

  if (any(grepl("Error", stderr, ignore.case = TRUE))) {
    cat(paste(stderr, "\n"))
    stop("Error during BoolNet file generation. Please run rxncon2boolnet.py on its own with the -v DEBUG flag.")
  }
}

callNFSim <- function(inFile, args=c()) {
  config <- loadPackageConfig()
  BNGDir <- config$value[config$setting == "BNGDir"]
  if (BNGDir == "") {
    stop("BioNetGen install directory has not been set. Run setupKboolnet() to fix this.")
  }

  path <- paste0(BNGDir, "/bin/NFsim")
  suppressWarnings(stderr <- system2(addQuotes(path), args = c("-rnf", inFile)))

  if (any(grepl("Error", stderr, ignore.case = TRUE))) {
    cat(paste(stderr, "\n"))
    stop("Error during NF file generation. Please run NFsim.")
  }
}

callBNG <- function(inFile, outFile, args=c()) {
  config <- loadPackageConfig()
  BNGDir <- config$value[config$setting == "BNGDir"]
  if (BNGDir == "") {
    stop("BioNetGen install directory has not been set. Run setupKboolnet() to fix this.")
  }

  path <- paste0(BNGDir, "/BNG2.pl")
  suppressWarnings(stderr <- system2(path, args = c(addQuotes(inFile), "--outdir", addQuotes(outFile)), stderr = TRUE, stdout = ""))

  if (any(grepl("(Error|ABORT)", stderr, ignore.case = TRUE))) {
    cat(paste(stderr, "\n"))
    stop("Error during BNG simulation. Please run BNG2.pl on its own.")
  }
}

callReactionMapping <- function(inFile, outFile, args=c()) {
  cleanFiles(outFile)
  path <- paste0(system.file(package="kboolnet"), "/python/reaction_mapping.py")
  stderr <- callPython(c(path, addQuotes(inFile), "--output", addQuotes(outFile), args))

  if (any(grepl("Error", stderr, ignore.case = TRUE))) {
    cat(paste(stderr, "\n"))
    stop("Error during reaction mapping generation. Please run reaction_mapping.py on its own with the -v DEBUG flag.")
  }
}

