# Loads in a config file based on an opt list and returns the opt list
loadConfig <- function(opt) {
  source(opt$config)

  if (!exists("config")) {
    stop("No config object found in config file")
  }

  config <- config[!is.na(config)] # Discard NA values

  # Keep only config values that were not passed as command line options
  config <- config[!(names(config) %in% names(opt))]

  # Merge command line args and config values
  opt <- c(opt, config)

  return(opt)
}

# Merges an opt with a list of defaults
setDefaults <- function(opt, defaults) {
  default <- default[!(names(default) %in% names(opt))]
  opt     <- c(opt, default)
  return(opt)
}

