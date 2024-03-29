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
# Dependencies: openxlsx
############################################################

############## This function taken from makeCNOlist in CellNOptR
compareNA <- function(v1,v2) {
    # This function returns TRUE wherever elements are the same, including NA's,
    # and false everywhere else.
    same <- (v1 == v2)  |  (is.na(v1) & is.na(v2))
    same[is.na(same)] <- FALSE
    return(same)
}

readMIDASExcel <- function(MIDASfile, cell_types = "all") {
  if (length(cell_types) == 1 & any(cell_types == "all")) {
    cell_types <- c()
  }

  # Load the Excel workbook
  wb <- openxlsx::loadWorkbook(MIDASfile)

  ############# Treatment definition parsing ###################
  treatmentDefs <- openxlsx::read.xlsx(wb, sheet = "TreatmentDefs", colNames = F)

  # Remove empty and non-MIDAS rows and columns
  treatmentDefs <- treatmentDefs[rowSums(is.na(treatmentDefs)) != ncol(treatmentDefs),]
  treatmentDefs <- treatmentDefs[,colSums(is.na(treatmentDefs)) != nrow(treatmentDefs)]
  treatmentDefs <- treatmentDefs[,grep("^!(Name|Type|Nodes)", treatmentDefs[1,])]

  # Set column names
  colnames(treatmentDefs)[grep("^!Name", treatmentDefs[1,])] <- "name"
  colnames(treatmentDefs)[grep("^!Type", treatmentDefs[1,])] <- "type"
  colnames(treatmentDefs)[grep("^!Nodes", treatmentDefs[1,])] <- "nodes"
  treatmentDefs <- treatmentDefs[2:nrow(treatmentDefs), ,drop=F]

  # Replace aliases of treatment types with their correct names
  treatmentDefs$type[grepl("^(Inhibitor|I|Inhib)", treatmentDefs$type, ignore.case = T)] <- "Inhibitor"
  treatmentDefs$type[grepl("^(Stimulus|S|Stimulus)", treatmentDefs$type, ignore.case = T)] <- "Stimulus"
  treatmentDefs$type[grepl("^(Knockout|Knockdown|KO)", treatmentDefs$type, ignore.case = T)] <- "KO"
  treatmentDefs$type[grepl("^(Mutant|Mut|M)", treatmentDefs$type, ignore.case = T)] <- "Mutant"

  # Validation of treatment types
  if (any(!grepl("(Inhibitor|Stimulus|KO|Mutant)", treatmentDefs$type))) {
    stop("Invalid treatment type detected in TreatmentDefs sheet")
  }

  # Turn comma separated lists into actual lists
  for (i in 1:nrow(treatmentDefs)) {
    if (!is.na(treatmentDefs$nodes[i])) treatmentDefs$nodes[i] <- list(trimws(strsplit(treatmentDefs$nodes[i][[1]], ",")[[1]]))
  }

  # Escape all regex characters and expand * to .*? for regexes
  treatmentDefs$regex <- list(nrow(treatmentDefs))
  for (i in 1:nrow(treatmentDefs)) {
    for (j in 1:length(treatmentDefs$nodes[[i]])) {
      treatmentDefs$regex[[i]][j] <- gsub('([.|()\\^{}+$?]|\\[|\\])', '\\\\\\1', treatmentDefs$nodes[[i]][j])
      treatmentDefs$regex[[i]][j] <- gsub('\\*', '\\.\\*\\?', treatmentDefs$regex[[i]][j])
    }
  }

  ########### Experimental data sheet loading ################
  # Load all sheets (except for TreatmentDef)
  sheetNames <- names(wb)[names(wb) != "TreatmentDefs"]
  sheets <- lapply(sheetNames, openxlsx::read.xlsx, xlsxFile = wb, colNames = F)

  # Extract sub-tables from all sheets
  tables <- list()
  for (i in 1:length(sheets)) {
    headerRows <- grep("TH", sheets[[i]][,1]) # Get row numbers of sub-tables' header rows in each sheet

    # If no TH present, throw warning
    if (length(headerRows) == 0) {
      stop("Sheet ", sheetNames[i], " has no TH cells, skipping.")
      next
    }

    # If more than one TH row present, iterate through and get all the sub-tables
    if (length(headerRows) > 1) {
      for (j in 1:(length(headerRows)-1)) {
        # Skip tables marked with :Ignore
        if (grepl(":Ignore$", sheets[[i]][headerRows[j],1])) {
          cat("Skipping table", sheets[[i]][headerRows[j],1], "\n")
          next
        }

        tables <- c(tables, list(sheets[[i]][(headerRows[j]+1):(headerRows[j+1]-1),])) # Add sub-table to tables list
        names(tables)[length(tables)] <- sheets[[i]][headerRows[j],1] # Set table name based on TH row
      }
    }

    if (!grepl(":Ignore$", sheets[[i]][headerRows[length(headerRows)],1])) {
      tables <- c(tables, list(sheets[[i]][(headerRows[length(headerRows)]+1):nrow(sheets[[i]]),])) # Do the same for the last sub-table (if it isnt ignored)
      names(tables)[length(tables)] <- sheets[[i]][headerRows[length(headerRows)],1] # Set table name based on TH row
    } else {
      cat("Skipping table", sheets[[i]][headerRows[length(headerRows)],1], "\n")
    }
  }
  cat("Found", length(tables), "tables across", length(sheets), "sheets.", "\n")


  ################ Table cleanup/validation #################
  # Discard tables not of correct cell type
  table_cell_types <- c()
  for (name in names(tables)) {
    s <- strsplit(name, ":")[[1]]
    if (length(s) > 2) {
      table_cell_types <- c(table_cell_types, s[3])
    } else {
      table_cell_types <- c(table_cell_types, NA)
    }
  }
  if (length(cell_types)> 0) {
    tables <- tables[which(table_cell_types %in% cell_types)]
    if (length(tables) == 0) {
      cat("No tables of cell type(s)", paste0(cell_types, collapse = ", "))
    }
  }

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
      tables[[i]][,TRcol[j]] <- as.numeric(tables[[i]][,TRcol[j]])
      if (any(sapply(tables[[i]][,TRcol[j]], function(x) !(x == 1 | x == 0) | is.na(x)))) { # If any value in row doesn't equal 0 or 1)
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
    stripped <- sapply(strsplit(names(tables[[i]])[col], ":"), function(x) paste0(x[1:2], collapse=":")) # Keep only name field

    # Try to resolve duplicate TRs/DVs
    if (any(grepl("DV:", stripped[duplicated(stripped)]))) {
      stop("Table ", names(tables)[i], " contains duplicate DV columns!")
    } else if (any(grepl("TR:", stripped[duplicated(stripped)]))) {
      warning("Table ", names(tables)[i], " contains duplicate TR columns! Attempting to merge (bitwise OR).")
      while(anyDuplicated(stripped)) {
        duplName <- stripped[duplicated(stripped)][1] # Get a duplicated TR name

        duplCols <- as.matrix(tables[[i]][,!grepl(paste0(duplName, "(:|$)"), names(tables[[i]]))]) # Get all the duplicated columns and remove them from original table
        class(duplCols) <- "numeric"
        tables[[i]] <- tables[[i]][,!grepl(paste0(duplName, "(:|$)"), names(tables[[i]]))]

        tables[[i]][,duplName] <- as.numeric(rowSums(duplCols) > 0) # Do a bitwise or across rows

        col <- grep("^(TR|DV):", names(tables[[i]])) # Update col and stripped
        stripped <- sapply(strsplit(names(tables[[i]])[col], ":"), function(x) paste0(x[1:2], collapse=":"))
      }
    }


    names(tables[[i]])[col] <- stripped
  }

  # Merge all tables into a single table
  combined <- tables[[1]]
  if (length(tables) > 1) {
    for (i in 2:length(tables)) {
      newCols <- colnames(tables[[i]])[!(colnames(tables[[i]]) %in% colnames(combined))] # Get names of all columns missing from combined
      combined[,newCols] <- NA # Add said rows to combined

      newCols <- colnames(combined)[!(colnames(combined) %in% colnames(tables[[i]]))] # Get names of all columns missing from tables[[i]]
      tables[[i]][,newCols] <- NA # Add said rows to tables[[i]]

      combined <- rbind(combined, tables[[i]]) # Merge combined and tables[[i]]
    }
  }

  ########### Get column numbers and names ##########
  # Get different column type names
  TRcol <- grep("^TR", colnames(combined))
  DVcol <- grep("^DV", colnames(combined))
  DAcol <- grep("^DA", colnames(combined))

  # Get TR and DV names
  TRnames <- sapply(strsplit(colnames(combined)[TRcol], ":"), "[[", 2)
  DVnames <- sapply(strsplit(colnames(combined)[DVcol], ":"), "[[", 2)

  # Escape DVnames for regex matching
  DVregex <- gsub('([.|()\\^{}+$?]|\\[|\\])', '\\\\\\1', DVnames)
  DVregex <- gsub('\\*', '\\.\\*\\?', DVregex)

  # Make sure TR names exist in TreatmentDefs
  for (i in 1:length(TRnames)) {
    if(!(TRnames[i] %in% treatmentDefs$name)) {
      stop("Treatment ", TRnames[i], " does not exist in TreatmentDefs.")
    }
  }

  # Set all treatments that are NA after merging to 0
  combined[,TRcol][is.na(combined[,TRcol])] <- 0

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
    # if all the values in a column were NA, set new value to NA
    newRow[colSums(is.na(dupMatrix)) == nrow(dupMatrix)] <- NA

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

    # Save to list
    valueSignals[[i]] <- timeData
    valueVariances[[i]] <- timeVariances
  }

  # Turn valueSignals and valueVariances into 3D arrays
  dimensions <- c(nrow(cues), length(DVnames), length(timeSignals))
  namesDimensions <- list(1:nrow(cues), DVnames, timeSignals)
  valueSignals <- array(unlist(valueSignals), dim=dimensions, dimnames=namesDimensions)
  valueVariances <- array(unlist(valueVariances), dim=dimensions, dimnames=namesDimensions)

  # Make cues a matrix
  cues <- as.matrix(cues, nrow=nrow(cues), ncol=ncol(cues))
  rownames(cues) <- 1:nrow(cues)

  return(list(
    treatmentDefs=treatmentDefs,
    namesCues=TRnames,
    valueCues=cues,
    namesSignals=DVnames,
    regexSignals=DVregex,
    timeSignals=timeSignals,
    valueSignals=valueSignals,
    valueVariances=valueVariances
  ))
}
