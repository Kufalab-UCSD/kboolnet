# Takes an array of node names and returns an array of the same nodes, but with domains removed
removeDomains <- function(nodes, noDuplicates = TRUE) {
  mapping <- data.frame(orig = nodes, noDomains = gsub("_\\[.*?\\]", "", nodes))

  if (!noDuplicates) {
    return(mapping$noDomains)
  }

  duplicates <- duplicated(mapping$noDomains) | duplicated(mapping$noDomains, fromLast = TRUE)
  mapping$noDomains[duplicates] <- mapping$orig[duplicates]

  return(mapping$noDomains)
}

# Removes domains from an XGMML file and writes the output to another XGMML file
removeDomainsXGMML <- function(inputFile, outputFile, noDuplicates = TRUE) {
  orig <- read.delim(inputFile, header = F, quote = "")[,1]

  nodes <- gsub(".* label=", "", orig)
  nodes <- gsub(">.*", "", nodes)
  nodes <- nodes[!grepl("<", nodes)]
  nodes <- gsub(" ", "", nodes)
  nodes <- nodes[nodes != ""]

  mapping <- data.frame(orig = nodes, noDomains = removeDomains(nodes, noDuplicates))
  mapping <- mapping[mapping$orig != mapping$noDomains,]

  noDomains <- orig
  for (i in 1:nrow(mapping)) {
    old <- paste0("label=", mapping$orig[i])
    new <- paste0("label=", mapping$noDomains[i])

    # We have to properly escape the regexes
    old <- gsub('([.|()\\^{}+$?]|\\[|\\])', '\\\\\\1', old)
    new <- gsub('([.|()\\^{}+$?]|\\[|\\])', '\\\\\\1', new)

    noDomains <- gsub(old, new, noDomains)
  }

  write(noDomains, file = outputFile)
}
