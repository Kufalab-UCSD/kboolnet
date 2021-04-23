##################################################
# driveDownload.R
# Adrian C
#
# Function that adds a bit of extra functionality
# to the googledrive package's download function.
#
# Dependencies: googledrive
##################################################

driveDownload <- function(driveFile, out, type=NULL) {
  if (grepl("^https?:\\/\\/", driveFile)) { # If URL provided
    # Deactivate authentication and try and access file
    googledrive::drive_deauth()

    # Try accessing file without authentication
    gDriveID <- googledrive::as_id(tryCatch({
      googledrive::drive_get(id = googledrive::as_id(driveFile))$id[1]
    }, error = function(e) {
      # If access fails, try again with authentication enabled
      repeat {
        cat("Public Google Drive file not found. Trying again with authentication...", "\n")
        googledrive::drive_auth() # Activate authentication

        return(tryCatch({ # Try to find file again
          googledrive::drive_get(id = googledrive::as_id(driveFile))$id[1]

        }, error = function(e){ # If it fails again
          stop("Google Drive file not found. Either the file does not exist or you do not have permission to view it.")
        }))
      }
    }))
  } else { # Search for the file in Drive if name provided
    # Activate Google Drive authentication, required for searching in a user's drive
    drive_auth_config(active = TRUE)

    gDriveID <- googledrive::as_id(googledrive::drive_find(pattern = driveFile, type = type)$id[1])

    if(is.na(gDriveID)) { # If file does not exist
      stop("File not found in Google Drive.")
    }
  }

  # Download file
  googledrive::drive_download(gDriveID, path = out, overwrite = T) # Download the file
}
