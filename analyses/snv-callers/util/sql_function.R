# Function for manipulating SQL file
#
# C. Savonen for ALSF - CCDL
# 2019
#
################################################################################
save_to_sql <- function(df, 
                        sql_file = opt$sql_file, 
                        tbl_name, 
                        overwrite_it = opt$overwrite) {
  # Saves a data.frame as an SQL file
  #
  # Args:
  #   df: the data.frame you would like to save to the SQLite file. 
  #   sql_file: path to sql file where the table is saved. 
  #   tbl_name: The name of the table to extract
  #   overwrite_it: should it be overwritten in the SQL file? 
  
  # Save to SQL Lite file
  con <- DBI::dbConnect(RSQLite::SQLite(), sql_file)
  
  # Read in and save data to 
  dplyr::copy_to(con, 
                 df, 
                 name = tbl_name, 
                 temporary = FALSE, 
                 overwrite = overwrite_it)
  
  # Disconnect 
  DBI::dbDisconnect(con)
}

sql_to_df <- function(sql_file = opt$sql_file, 
                      tbl_name) {
  # Converts a SQL table to a data.frame
  #
  # Args:
  #   sql_file: path to sql file where the table is saved. 
  #   tbl_name: The name of the table to extract
  
  # Save to SQL Lite file
  con <- DBI::dbConnect(RSQLite::SQLite(), sql_file)
  
  # Let's make the SQL table its own object
  sql_table <- dplyr::tbl(con, tbl_name)
  
  # Extract MAF file from the 
  df <- sql_table %>% 
    dplyr::collect()
  
  # Disconnect 
  DBI::dbDisconnect(con)
  
  return(df)
}
