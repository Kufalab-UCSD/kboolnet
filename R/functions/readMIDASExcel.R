#############################################################
# readMIDASExcel.R
# Adrian C
#
# Reimplementation of CellNOptR's readMIDAS function to fit
# the needs of the Kufareva lab. See specification for modified
# MIDAS format here:
# https://docs.google.com/document/d/1i0h_a7LnEhCerYy-C6f4xoKwpqw_3P_CBR3CzKGyB_4/edit?usp=sharing
#
# Args:
#   -MIDASfile: MIDAS file to be read
#
# Dependencies: openxlsx, dplyr
############################################################

############## This function taken from makeCNOlist in CellNOptR
compareNA <- function(v1,v2) {
    # This function returns TRUE wherever elements are the same, including NA's,
    # and false everywhere else.
    same <- (v1 == v2)  |  (is.na(v1) & is.na(v2))
    same[is.na(same)] <- FALSE
    return(same)
}

readMIDASExcel <- function(MIDASfile) {
  # Load the Excel workbook
  wb <- loadWorkbook(MIDASfile)
  
  ############# Treatment definition parsing ###################
  treatmentDefs <- read.xlsx(wb, sheet = "TreatmentDefs", colNames = F)
  
  # Remove empty and non-MIDAS rows and columns
  treatmentDefs <- treatmentDefs[rowSums(is.na(treatmentDefs)) != ncol(treatmentDefs),]
  treatmentDefs <- treatmentDefs[,colSums(is.na(treatmentDefs)) != nrow(treatmentDefs)]
  treatmentDefs <- treatmentDefs[,grep("^!(Name)|(Type)|(Nodes)|(Components)", treatmentDefs[1,])]
  
  # Set column names
  colnames(treatmentDefs)[grep("^!Name", treatmentDefs[1,])] <- "name"
  colnames(treatmentDefs)[grep("^!Type", treatmentDefs[1,])] <- "type"
  colnames(treatmentDefs)[grep("^!Nodes", treatmentDefs[1,])] <- "nodes"
  colnames(treatmentDefs)[grep("^!Components", treatmentDefs[1,])] <- "components"
  treatmentDefs <- treatmentDefs[2:nrow(treatmentDefs), ,drop=F]
  
  # Replace aliases of treatment types with their correct names
  treatmentDefs$type[grepl("^(Inhibitor|I|Inhib)", ignore.case = T)] <- "Inhibitor"
  treatmentDefs$type[grepl("^(Stimulus|S|Stimulus)", ignore.case = T)] <- "Stimulus"
  treatmentDefs$type[grepl("^(Knockout|Knockdown|KO)", ignore.case = T)] <- "KO"
  
  # Validation of treatment types
  if (any(!grepl("(Inhibitor|Stimulus|KO", treatmentDefs$type))) {
    stop("Invalid treatment type detected in TreatmentDefs sheet")
  }
  
  ########### Experimental data sheet loading ################
  # Load all sheets (except for TreatmentDef)
  sheetNames <- names(wb)[names(wb) != "TreatmentDef"]
  sheets <- lapply(sheetNames, read.xlsx, xlsxFile = wb, colNames = F)
  
  # Extract sub-tables from all sheets
  tables <- list()
  for (i in 1:length(sheets)) {
    headerRows <- grep("TH", sheets[[i]][,1]) # Get row numbers of sub-tables' header rows in each sheet
    
    # If more than one TH row present, iterate through and get all the sub-tables
    if (length(headerRows) > 1) {
      for (j in 1:(length(headerRows)-1)) {
        tables <- c(tables, list(sheets[[i]][(headerRows[j]+1):(headerRows[j+1]-1),])) # Add sub-table to tables list
        names(tables)[length(tables)] <- strsplit(sheets[[i]][headerRows[j],1], ":")[[1]][2] # Set table name based on TH row
      }
    }
    tables <- c(tables, list(sheets[[i]][(headerRows[length(headerRows)]+1):nrow(sheets[[i]]),])) # Do the same for the last sub-table
    names(tables)[length(tables)] <- strsplit(sheets[[i]][headerRows[length(headerRows)],1], ":")[[1]][2] # Set table name based on TH row
  }
  cat("Found ", length(tables), " tables across ", length(sheets), " sheets.", "\n")
  
  
  ################ Table cleanup/validation #################
  # Remove empty (i.e. all NA) and non-MIDAS rows and columns
  tables <- lapply(tables, function(x) {
    x <- x[rowSums(is.na(x)) != ncol(x),]
    x <- x[,colSums(is.na(x)) != nrow(x)]
    x <- x[,grep("(^TH:)|(^TR:)|(^DA$)|(^DV:)", x[1,])]
    return(x)
  })
  
  # Set header row as column names
  tables <- lapply(tables, function(x) {
    names(x) <- x[1,]
    x <- x[-1,]
  })
  
  # Do some validation of data rows
  for (i in 1:length(tables)) {
    # Get all the rows of each type
    TRcol <- grep("^TR:", names(tables[[i]]))
    DAcol <- grep("^DA$", names(tables[[i]]))
    DVcol <- grep("^DV:", names(tables[[i]]))
    
    if (length(TRcol) < 1) stop(paste0("Table ", names(tables)[i], " is missing treatment (TR) columns!"))
    for (j in TRcol) {
      if (any(sapply(tables[[i]][,TRcol[j]], function(x) !(x == 1 | x == 0)))) { # If any value in row doesn't equal 0 or 1
        stop(paste0("Row ", names(tables[[i]])[TRcol[j]], " in table ", names(tables)[i], " has incorrect values (must be 0 or 1)."))
      }
    }
    
    if (length(DAcol) == 0) stop(paste0("Table ", names(tables)[i], " is missing a timepoint (DA) column!"))
    if (length(DAcol) > 1) stop(paste0("Table ", names(tables)[i], " has too many timepoint (DA) columns! There should be only one
                                        DA column per table."))
    if (any(is.na(tables[[i]][,DAcol[1]]))) stop(paste0("DA row in table ", names(tables)[i], " has NA values."))
    
    if (length(DVcol) < 1) stop(paste0("Table ", names(tables)[i], " is missing data value (DV) columns!"))
    
    # Try and coerce all columns to numeric values
    for (j in colnames(tables[[i]])) {
      tables[[i]][,j] <- as.numeric(tables[[i]][,j])
    }
  }
  
  ################ Table merging ################
  # Strip all notes fields from TR and DA columns. This ensures correct table merging
  for (i in 1:length(tables)) {
    col <- grep("^(TR|DV):", names(tables[[i]])) # Get TR and DV column numbers
    stripped <- sapply(strsplit(names(tables[[i]])[col]), function(x) paste0(x[1:2], collapse=":")) # Keep only name field
    names(tables[[i]][col]) <- stripped
  }
  
  # Merge all tables into a single table
  combined <- tables[[1]]
  if (length(tables) > 1) {
    for (i in 2:length(tables)) {
      newCols <- colnames(tables[[i]])[!(colnames(tables[[i]]) %in% colnames(combined))] # Get names of all columns missing from combined
      combined[,newCols] <- 0 # Add said rows to combined
      
      newCols <- colnames(combined)[!(colnames(combined) %in% colnames(tables[[i]]))] # Get names of all columns missing from tables[[i]]
      tables[[i]][,newCols] <- 0 # Add said rows to tables[[i]]
      
      combined <- rbind(combined, tables[[i]]) # Merge combined and tables[[i]]
    }
  }
  
  ########### Get column numbers and names ##########
  # Get different column type names
  TRcol <- grep("^TR", colnames(combined))
  DVcol <- grep("^DV", colnames(combined))
  DAcol <- grep("^DA", colnames(combined))
  
  # Get TR columns' node, names, and type attributes 
  TRheaders <- strsplit(colnames(combined)[TRcol], ":")
  TRnodes <- character()
  TRtypes <- character()
  TRnames <- character()
  for (i in 1:length(TRheaders)) {
    TRnodes[i] <- TRheaders[[i]][2]
    
    # If TR doesn't have enough fields, throw an error
    if (length(TRheaders[[i]]) < 3) {
      stop(paste0("Treatment ", paste(TRheaders[[i]], sep="", collapse=":"), " is missing fields."))
    }
    
    # If name field not provided or blank, set node field as name
    if (length(TRheaders[[i]]) == 3) {
      TRnames[i] <- TRheaders[[i]][2]
    } else if (TRheaders[[i]][4] == "") {
      TRnames[i] <- TRheaders[[i]][2]
    } else {
      TRnames[i] <- TRheaders[[i]][4]
    }
    
    # If node field is blank, throw an error
    if (TRheaders[[i]][2] == "") {
      stop(paste0("Treatment ", paste(TRheaders[[i]], sep="", collapse=":"), " is missing the node field."))
    }
    
    # Match treatment types
    type <- tolower(TRheaders[[i]][3])
    if (type == "inhibitor" | type == "i" | type == "inhib") {
      TRtypes[i] <- "Inhibitor"
    } else if (type == "stimulus" | type == "stim" | type == "s") {
      TRtypes[i] <- "Stimulus"
    } else if (type == "ko" | type == "knockout") {
      TRtypes[i] <- "KO"
    } else {
      stop(paste0("\"", TRheaders[[i]][3], "\" is not a valid treatment type for treatment ", paste(TRheaders[[i]], sep="", collapse=":")))
    }
  }
  
  # Get signal names
  DVheaders <- strsplit(colnames(combined)[DVcol], ":")
  namesSignals <- character()
  # If name field not provided, set node field as name
  for (i in 1:length(DVheaders)) {
    if (length(DVheaders[[i]]) == 2) {
      namesSignals[i] <- DVheaders[[i]][2]
    } else {
      namesSignals[i] <- DVheaders[[i]][3]
    }
  }

  ################# Duplicated condition cleanup ####################
  #### This code taken from makeCNOlist from CellNOptR
  # Get duplicate rows
  duplCond <- as.matrix(combined[,c(TRcol, DAcol)])
  duplRows <- which(duplicated(duplCond) == TRUE)

  # creates a variance matrix
  variances = combined * 0
  while(length(duplRows) != 0){
    # the all(x == ) is buggy in the case of NA hence the compareNA function.
    #dupIndex<-apply(duplCond,MARGIN=1,function(x) all(x == duplCond[duplRows[1],]))
    dupIndex<-apply(duplCond,MARGIN=1,function(x)all(compareNA(x,duplCond[duplRows[1],])))

    dupIndex<-which(dupIndex == TRUE)
    dupMatrix<-combined[dupIndex,]
    #compute the new row as the average across duplicate rows
    newRow<-colMeans(dupMatrix, na.rm=TRUE)
    # variance for these rows
    newVariance = apply(dupMatrix, MARGIN=2, FUN=var, na.rm=T)

    #replace the first duplicated row by the summarised one
    combined[dupIndex[1],]<-newRow

    # same for the variance
    variances[dupIndex[1],]<-newVariance

    #remove the other summarised rows
    combined<-combined[-dupIndex[2:length(dupIndex)],]
    variances<-variances[-dupIndex[2:length(dupIndex)],]

    duplCond<-as.matrix(combined[,c(TRcol,DAcol)])
    duplRows<-which(duplicated(duplCond) == TRUE)
  }
  
  ################ Timepoint extraction ###########
  # create cues matrix, remove duplicate conditions, and order
  cues <- combined[,TRcol]
  cues <- distinct(cues)
  cuesStr <- as.vector(apply(cues, 1, paste, collapse=""))
  cues <- cues[order(cuesStr),]
  
  # get the timepoints
  timeSignals <- sort(unique(combined[,DAcol]))
  
  # build valueSignals and valueVariance lists
  valueSignals <- list()
  valueVariances <- list()
  for (i in 1:length(timeSignals)) {
    # Get data and variance rows for given timepoint
    timeData <- combined[combined[,DAcol] == timeSignals[i],]
    timeVariances <- variances[combined[,DAcol] == timeSignals[i],]
    
    # Check that all cue combinations exist for timepoint and add any that are missing
    if (nrow(timeData) < nrow(cues)) { # If cue rows are missing, add them
      timeCues <- rbind(timeData[,TRcol], cues) # First, combine all cues df with cues for only this timepoint
      duplCues <- duplicated(timeCues) | duplicated(timeCues, fromLast=T) # Then, find all duplicated rows in the timeCues df
      missingTimeCues <- timeCues[!duplCues,] # Remove those rows from timeCues, leaving only those cues which were not present for this timepoint
      missingCols <- colnames(timeData)[!(colnames(timeData) %in% colnames(missingTimeCues))] # Check which columns are missing
      missingTimeCues[,missingCols] <- NA # Add missing columns to missingTimeCues
      timeData <- rbind(timeData, missingTimeCues)
    }
    
    # Order the data and variances
    timeCuesStr <- as.vector(apply(timeData[,TRcol], 1, paste, collapse="")) # Collapse cues as string
    timeData <- timeData[order(timeCuesStr),]
    timeVariances <- timeVariances[order(timeCuesStr),]
    
    # Strip col and row names
    colnames(timeData) <- NULL
    rownames(timeVariances) <- NULL
    colnames(timeData) <- NULL
    rownames(timeVariances) <- NULL
    
    # Get DV data and convert to matrix
    timeData <- timeData[,DVcol, drop=F] # drop=F prevents single columns from becoming vectors
    timeVariances <- timeVariances[,DVcol, drop=F]
    timeData <- as.matrix(timeData, nrow=nrow(timeData), ncol=ncol(timeData))
    timeVariances <- as.matrix(timeVariances, nrow=nrow(timeVariances), ncol=ncol(timeVariances))
    
    # Strip dimnames
    dimnames(timeData) <- NULL
    dimnames(timeVariances) <- NULL
    
    # Save to list
    valueSignals[[i]] <- timeData
    valueVariances[[i]] <- timeVariances
  }
  
  # Make cues a matrix
  cues <- as.matrix(cues, nrow=nrow(cues), ncol=ncol(cues))
  dimnames(cues) <- NULL

  return(list(
    namesCues=TRnames,
    nodesCues=TRnodes,
    typesCues=TRtypes,
    valueCues=cues,
    namesSignals=namesSignals,
    timeSignals=timeSignals,
    valueSignals=valueSignals,
    valueVariances=valueVariances
  ))
}
