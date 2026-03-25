extract_selected_files <- function(zip_path, members) {
  # Create a temporary directory for extraction
  temp_dir <- tempfile("unzip_dir_")
  dir.create(temp_dir)

  # List all files in the archive
  all_files <- utils::unzip(zip_path, list = TRUE)$Name

  # Find files matching any of the glob patterns in 'members'
  selected <- unique(unlist(
    lapply(members, function(pattern) {
      regex <- glob2rx(pattern)
      grep(regex, all_files, value = TRUE)
    })
  ))

  # Extract only the selected files
  utils::unzip(zip_path, files = selected, exdir = temp_dir)

  # Return the path to the extracted folder
  file.path(temp_dir)
}
