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


setupKboolnet <- function(forceDefaults = FALSE) {
  if (!forceDefaults) {
    cat("Manual configuration entry. Leave response blank to keep default/previous settings, which are indicated in parentheses.\n")
  }

  # Load in config file if it exists
  configFile <- paste0(rappdirs::user_config_dir(appname="kboolnet"), "/config.csv")
  if (file.exists(configFile)) {
    config <- read.csv(file = configFile)
    pythonCommand <- config$value[config$setting == "pythonCommand"]
  } else {
    config <- NA
    pythonCommand <- "python3"
  }

  if (!forceDefaults) {
    while (TRUE) {
      newPythonCommand <- readline(paste0("Enter python3 command name (", pythonCommand, "): "))
      if (trimws(newPythonCommand) != "") {
        newPythonCommand <- normalizePath(newPythonCommand, mustWork = FALSE)
        if (system2(newPythonCommand, args = "--version", stdout = FALSE) != 0) {
          cat(paste0("Could not run python3 using command: ", newPythonCommand))
        } else {
          pythonCommand <- newPythonCommand
          break
        }
      } else {
        break
      }
    }
  }

  cat("Checking if python3, pip, and rxncon are installed...\n")
  if (system2(pythonCommand, args = "--version", stdout = FALSE) != 0) {
    stop("Could not run python3. Please install before continuing!")
  }

  pythonVer <- system2(pythonCommand, args = "--version", stdout=TRUE)
  pythonLocation <- system2(pythonCommand, args = "-c \"import sys; print(sys.executable)\"", stdout=TRUE)
  cat("Running", pythonVer, "from", pythonLocation, "\n")

  if (system2(pythonCommand, args = c("-m", "pip"), stdout = FALSE) != 0) {
    stop("Could not run pip. Please install before continuing!")
  }

  checkRxnconInstallScript <- paste0(system.file(package="kboolnet"), "/python/check_rxncon_install.py")
  if (system2(pythonCommand, args = c(addQuotes(checkRxnconInstallScript))) != 0) {
    stop("Could not load rxncon module. Please install with pip before continuing!")
  }

  if (!any(is.na(config))) {
    defaults <- config
  } else {
    cat("Detecting sensible default values...\n")
    defaults <- sensibleDefaults(pythonCommand)
  }

  rxnconDir <- defaults$value[defaults$setting == "rxnconDir"]
  BNGDir <- defaults$value[defaults$setting == "BNGDir"]
  installDir <- defaults$value[defaults$setting == "installDir"]
  installed <- as.logical(defaults$value[defaults$setting == "installed"])

  # If already installed, ask if we want to reinstall
  reinstall <- FALSE
  if (installed) {
    cat(paste0("kboolnet scripts already installed in directory ", installDir, ""), "\n")
    if (!forceDefaults) {
      reinstall <- getYesNo("Would you like to reinstall the kboolnet scripts? (y/N): ", default = FALSE)
    }
  }

  # If installation necessary, prompt for new dir and do the installation
  if (!installed | reinstall) {
    while (TRUE) {
      if (!forceDefaults) {
        newInstallDir <- readline(paste0("Enter directory in which to install kboolnet scripts to (", installDir, "): "))
        if (trimws(newInstallDir) != "") {
          installDir <- normalizePath(newInstallDir, mustWork = FALSE)
        }

        verify <- getYesNo(paste0("Install kboolnet scripts to dir ", installDir, "? (Y/n): "), default = TRUE)
        if (!verify) {
          next
        }
      }

      if (length(list.files(installDir)) > 0) {
        cat(paste0(installDir, " is non-empty. May cause issues.", "\n"))
      }

      if (!dir.exists(installDir)) {
        dir.create(installDir, recursive = TRUE)
      }

      if (!dir.exists(installDir)) {
        cat(tail(warnings, 1), "\n")
        cat(paste0("Unable to write to directory ", installDir, ", please try a different directory.\n"))
      }

      file.copy(list.files(paste0(system.file(package="kboolnet"), "/scripts"), full.names = TRUE), installDir, recursive = TRUE)
      if (!file.exists(paste0(installDir, "/VerifyModel.R"))) {
        cat(tail(warnings, 1), "\n")
        if (!forceDefaults) {
          cat(paste0("Unable to write scripts to directory ", installDir, ", please try a different directory.\n"))
        } else {
          stop(paste0("Unable to write scripts to default directory ", installDir))
        }
      } else {
        break
      }
    }
  }

  # Prompt for rxncon install directory
  while (TRUE) {
    if (!forceDefaults) {
      newRxnconDir <- readline(paste0("Enter existing rxncon2____.py script directory (", rxnconDir, "): "))
      if (trimws(newRxnconDir) != "") {
        rxnconDir <- normalizePath(newRxnconDir, mustWork = FALSE)
      }
    }

    if (file.exists(paste0(rxnconDir, "/rxncon2regulatorygraph.py"))) {
      break
    } else if (!forceDefaults) {
      cat(paste0("Could not detect rxncon scripts within directory ", rxnconDir, ". Please try again."))
    } else {
      stop(paste0("Could not detect rxncon scripts within directory ", rxnconDir))
    }
  }

  # Prompt for BNG install directory
  installBNG <- FALSE
  if (!forceDefaults) {
    if (BNGDir == "") {
      installBNG <- getYesNo("Do you wish to enter the existing BioNetGen install directory? (y/N): ", default = FALSE)
    } else {
      installBNG <- getYesNo(paste0("Do you wish to change the existing BioNetGen install directory? Currently set to ", BNGDir, " (y/N): "), default = FALSE)
    }
  }
  while (installBNG) {
    newBNGDir <- readline("Enter existing BioNetGen install directory: ")
    newBNGDir <- normalizePath(newBNGDir, mustWork = FALSE)

    if (file.exists(paste0(newBNGDir, "/BNG2.pl"))) {
      break
    } else {
      cat(paste0("Could not detect BioNetGen within directory ", newBNGDir, ". Please try again."))
    }
  }

  if (BNGDir == "") {
    cat("No BioNetGen install directory provided. Scripts depending on BioNetGen/NFsim will not work, run setup again to properly set this directory.\n")
  }

  newDefaults <- data.frame(setting = c("rxnconDir", "BNGDir", "installDir", "installed", "pythonCommand"),
                            value = c(rxnconDir, BNGDir, installDir, TRUE, pythonCommand))
  configDir <- rappdirs::user_config_dir(appname="kboolnet")
  if (!dir.exists(configDir)) {
    dir.create(configDir, recursive = TRUE)
  }
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

