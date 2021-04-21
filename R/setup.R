getYesNo <- function(prompt, default = TRUE) {
  res <- trimws(tolower(readline(prompt)))
  while (TRUE) {
    if (res %in% c("y", "ye", "yes")) {
      return(TRUE)
    } else if (res %in% c("n", "no")) {
      return(FALSE)
    } else if (res == "") {
      return(default)
    }
  }
}


setupKboolnet <- function() {
  # Load in config file if it exists
  configFile <- paste0(system.file(package="kboolnet"), "/config.csv")
  if (file.exists(configFile)) {
    config <- read.csv(file = configFile)
    oldPythonCommand <- config$value[config$setting == "pythonCommand"]
  } else {
    config <- NA
    oldPythonCommand <- "python3"
  }

  newPythonCommand <- readline(paste0("Set python 3 command name (", oldPythonCommand, "): "))
  if (trimws(newPythonCommand) != "") {
    newPythonCommand <- normalizePath(newPythonCommand, mustWork = FALSE)
  } else {
    newPythonCommand <- oldPythonCommand
  }

  cat("Checking if python3, pip, and rxnxon are installed...\n")
  if (system2(newPythonCommand, args = "--version", stdout = FALSE) != 0) {
    stop("Could not run python3. Please install before continuing!")
  }

  if (system2(newPythonCommand, args = c("-m", "pip"), stdout = FALSE) != 0) {
    stop("Could not run pip. Please install before continuing!")
  }

  checkRxnconInstallScript <- paste0(system.file(package="kboolnet"), "/python/check_rxncon_install.py")
  if (system2(newPythonCommand, args = c(addQuotes(checkRxnconInstallScript))) != 0) {
    stop("Could not load rxncon module. Please install with pip before continuing!")
  }

  if (!any(is.na(config))) {
    defaults <- config
  } else {
    cat("Detecting sensible default values...\n")
    defaults <- sensibleDefaults(newPythonCommand)
  }

  oldRxnconDir <- defaults$value[defaults$setting == "rxnconDir"]
  oldBNGDir <- defaults$value[defaults$setting == "BNGDir"]
  oldInstallDir <- defaults$value[defaults$setting == "installDir"]
  installed <- as.logical(defaults$value[defaults$setting == "installed"])

  cat("Manual configuration entry. Leave response blank to keep default/previous settings.\n")
  # If already installed, ask if we want to reinstall
  reinstall <- FALSE
  if (installed) {
    cat(paste0("rxncon alredy installed in directory ", oldInstallDir, ""), "\n")
    reinstall <- getYesNo("Would you like to reinstall the kboolnet scripts? (y/N): ", default = FALSE)
  }

  # If installation necessary, prompt for new dir and do the installation
  newInstallDir <- oldInstallDir
  if (!installed | reinstall) {
    while (TRUE) {
      newInstallDir <- readline(paste0("Set directory in which to install kboolnet scripts to (", oldInstallDir, "): "))
      if (trimws(newInstallDir) != "") {
        newInstallDir <- normalizePath(newInstallDir)
      } else {
        newInstallDir <- oldInstallDir
      }

      verify <- getYesNo(paste0("Install kboolnet scripts to dir ", newInstallDir, "? (Y/n): "), default = TRUE)
      if (!verify) {
        next
      }

      if (length(list.files(newInstallDir)) > 0) {
        cat(paste0(newInstallDir, " is non-empty. May cause issues.", "\n"))
      }

      dir.create(newInstallDir, recursive = TRUE)
      if (!dir.exists(newInstallDir)) {
        cat(tail(warnings, 1), "\n")
        cat(paste0("Unable to write to directory ", newInstallDir, ", please try a different directory.\n"))
      }
      file.copy(list.files(paste0(system.file(package="kboolnet"), "/scripts"), full.names = TRUE), newInstallDir, recursive = TRUE)
      if (!file.exists(paste0(newInstallDir, "/VerifyModel.R"))) {
        cat(tail(warnings, 1), "\n")
        cat(paste0("Unable to write scripts to directory ", newInstallDir, ", please try a different directory.\n"))
      } else {
        break
      }
    }
  }

  # Prompt for rxncon install directory
  while (TRUE) {
    newRxnconDir <- readline(paste0("Set existing rxncon2____.py script directory (", oldRxnconDir, "): "))
    if (trimws(newRxnconDir) != "") {
      newRxnconDir <- normalizePath(newRxnconDir)
    } else {
      newRxnconDir <- oldRxnconDir
    }

    if (file.exists(paste0(newRxnconDir, "/rxncon2regulatorygraph.py"))) {
      break
    }
    cat(paste0("Could not detect rxncon scripts within directory ", newRxnconDir, ". Please try again."))
  }

  # Prompt for BNG install directory
  newBNGDir <- oldBNGDir
  if (oldBNGDir == "") {
    installBNG <- getYesNo("Do you wish to set the existing BioNetGen install directory? (Y/n): ", default = TRUE)
  } else {
    installBNG <- getYesNo(paste0("Do you wish to change the existing BioNetGen install directory? Currently set to ", oldBNGDir, " (y/N): "), default = FALSE)
  }
  while (installBNG) {
    newBNGDir <- readline("Set existing BioNetGen install directory: ")
    newBNGDir <- normalizePath(newBNGDir)

    if (file.exists(paste0(newBNGDir, "/BNG2.pl"))) {
      break
    } else {
      cat(paste0("Could not detect BioNetGen within directory ", newBNGDir, ". Please try again."))
    }
  }

  if (newBNGDir == "") {
    cat("No BioNetGen install directory provided. Scripts depending on BioNetGen/NFsim will not work, run setup again to properly set this directory.\n")
  }

  newDefaults <- data.frame(setting = c("rxnconDir", "BNGDir", "installDir", "installed", "pythonCommand"),
                            value = c(newRxnconDir, newBNGDir, newInstallDir, TRUE, newPythonCommand))
  write.csv(newDefaults, file = configFile, row.names = FALSE)
  cat("Configuration successfully set.\n")

}

list.dirs.depth.n <- function(p, n) {
  res <- list.dirs(p, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.depth.n(res, n-1)
    c(res, add)
  } else {
    res
  }
}

sensibleDefaults <- function(pythonCommand) {
  res <- data.frame(setting = character(), value = character())

  # Detect rxncon script install dir
  pipOutput <- system2(pythonCommand, args = c("-m", "pip", "show", "rxncon", "-f"), stdout = TRUE)
  rxnconBase <- pipOutput[grepl("Location: ", pipOutput)]
  rxnconBase <- gsub("Location: ", "", rxnconBase)
  rxnconScriptsRel <- trimws(pipOutput[grepl("rxncon2regulatorygraph.py", pipOutput)])
  rxnconScripts <- normalizePath(paste0(rxnconBase, "/", rxnconScriptsRel))
  rxnconScripts <- gsub("rxncon2regulatorygraph.py", "", rxnconScripts)

  if (dir.exists(rxnconScripts)) {
    res[nrow(res) + 1,] <- c("rxnconDir", rxnconScripts)
  }

  # # Try and find BioNetGen somwhere in the home directory
  # homeFiles <- list.dirs.depth.n(Home(), 2)
  # bngDir <- homeFiles[grepl("BioNetGen", homeFiles)]
  # bngDir <- bngDir[which.min(sapply(bngDir, length))]
  #
  # if (dir.exists(bngDir)) {
  #   res[nrow(res) + 1,] <- c("BNGDir", bngDir)
  # }

  # Set BNGdir as empty for now
  res[nrow(res) + 1,] <- c("BNGDir", "")

  # Set arbitrary dir as install directory
  res[nrow(res) + 1,] <- c("installDir", suppressWarnings(normalizePath(paste0(Home(), "/kboolnet_scripts"))))

  # Set installed as false
  res[nrow(res) + 1,] <- c("installed", FALSE)

  return(res)
}

Home <- function() {
  # Returns a string with the user's home directory
  #
  # Serves as a replacement for "~", and works both in RStudio and RScript.
  #
  # Returns:
  #   On Windows returns C:/Users/<username>, where <username> is the current user's username.

  if (.Platform$OS.type == "windows") {
    return(normalizePath(file.path(Sys.getenv("HOMEDRIVE"), Sys.getenv("HOMEPATH")), winslash = .Platform$file.sep))
  } else {
    return(normalizePath("~"))
  }
}

