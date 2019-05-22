########################################################
# extractModules.R
# Adrian C
#
# Takes master reaction list excel file and filters by
# module tag and quality. Writes output as Excel file.
#
# Depends on openxlsx, dplyr
#######################################################

extractModules <- function(inPath, outPath, modules = character(0), minQuality = 0) {
  # Load workbook
  wb      <- loadWorkbook(inPath)

  # Delete all comments, erorrs will happen otherwise
  for(sheet in names(wb)){
    removeComment(wb, sheet, 1:1000, 1:1000)
  }

  # Read reaction and contingency lists
  rxnList <- readWorkbook(xlsxFile = wb, sheet = "ReactionList", colNames = FALSE)
  conList <- readWorkbook(xlsxFile = wb, sheet = "ContingencyList", colNames = FALSE)

  # Rename columns in which module and quality info is saved for reaction sheet
  rxnHdrRow <- grep("!UID:Reaction", rxnList[,1], fixed = T)
  colnames(rxnList)[grep("!Quality", rxnList[rxnHdrRow,], fixed = T)] <- "quality"
  colnames(rxnList)[grep("!Module", rxnList[rxnHdrRow,], fixed = T)]  <- "module"
  
  # Do same for contingency sheet
  conHdrRow <- grep("!UID:Contingency", conList[,1], fixed = T)
  colnames(conList)[grep("!Quality", conList[conHdrRow,], fixed = T)] <- "quality"
  colnames(conList)[grep("!Module", conList[conHdrRow,], fixed = T)]  <- "module"
  
  # Select all rule rows
  rxnFiltered <- rxnList[rxnHdrRow+1:nrow(rxnList),]
  conFiltered <- conList[conHdrRow+1:nrow(conList),]
  
  # Filter by modules present in modules argument (only if modules provided)
  if(length(modules) > 0) {
    print("doot")
    # Find which rows are in selected modules
    for (i in 1:nrow(rxnFiltered)) {
      rxnFiltered$hasModule[i] <- any(modules %in% trimws(strsplit(rxnFiltered$module[i], ",")[[1]]))
    }
    for (i in 1:nrow(conFiltered)) {
      conFiltered$hasModule[i] <- any(modules %in% trimws(strsplit(conFiltered$module[i], ",")[[1]]))
    }
    
    # Save only those rows with selected modules
    rxnFiltered <- rxnFiltered %>% filter(hasModule == T)
    conFiltered <- conFiltered %>% filter(hasModule == T)
    
    # Remove hasModule column
    rxnFiltered <- rxnFiltered %>% subset(select = -hasModule)
    conFiltered <- conFiltered %>% subset(select = -hasModule)
  }
  
  # Filter by quality greater or equal to minQuality
  rxnFiltered <- rxnFiltered %>% filter(quality >= minQuality | is.na(quality))
  conFiltered <- conFiltered %>% filter(quality >= minQuality | is.na(quality))
  
  # Delete unfiltered data from workbook
  deleteData(wb, "ReactionList", 1:ncol(rxnList), rxnHdrRow+1:nrow(rxnList), gridExpand = T)
  deleteData(wb, "ContingencyList", 1:ncol(conList), conHdrRow+1:nrow(conList), gridExpand = T)
  
  # Insert new data into workbook
  writeData(wb, "ReactionList", rxnFiltered, startCol = 1, startRow = rxnHdrRow+1, colNames = F)
  writeData(wb, "ContingencyList", conFiltered, startCol = 1, startRow = conHdrRow+1, colNames = F)
  
  # Write workbook to file
  saveWorkbook(wb, outPath, overwrite = T)
}