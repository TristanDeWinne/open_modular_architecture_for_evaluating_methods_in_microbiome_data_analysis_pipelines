# This script contains the code for handling the log files such as writing 
# to logs

write_to_log <- function(messages, log_file, folder="log", print_to_console=TRUE)
{
  ##############################################################################
  # messages(string):
  #   Message to write away to log file.
  # log_file (string):
  #   Name of file to write to.
  # folder (string):
  #   Name of folder where the log file is located.
  # print_to_console (boolean):
  #   A flag to set if the message should be printed to the console as well.
  ##############################################################################
  # Output (): 
  #   This function does not return an output.
  ##############################################################################
  
  # Check for valid entries
  if(!all(is.character(messages)))
  {
    stop("Messages should be a string or a vector of strings.")
  }
  
  # Check if file already exists, if not create it.
  full_file_name <- paste(folder, log_file, sep = "/")
  if(!file.exists(full_file_name))
  {
    file.create(full_file_name)
  }
  
  # Open file for writing
  file_conn <- file(full_file_name, open = "a")
  
  # Write each line in messages to the log file
  for(i in 1:length(messages))
  {
    # Add timestamp to message
    message <- paste(paste0("[",Sys.time(),"]"),
                     messages[i],
                     sep = "\t")
    
    # Append file with message
    writeLines(text = message,
               con = file_conn)
    
    # Print file to console if flag is set
    if(print_to_console)
    {
      cat(message, fill = TRUE)
    }
  }
  # Close file
  close(file_conn)
  
}
