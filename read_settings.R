# This R script contains the functions and code to read in the simulation 
# settings (contained in the folder simulation_settings)

# Load used packages
if(! "stringr"  %in% (.packages())) {library(stringr)}

# Helper functions
# Helper function to recognize if string encodes a matrix
is.matrix_setting <- function(val.str)
{
  ##############################################################################
  # val.str (string):
  #   A string obtained from the settings file. Contains the value of the 
  #   setting. Matrixes are nested arrays in the file and arrays are encoded
  #   with "[]". Thus matrixes are encoded as "[[],...,[]]"
  ##############################################################################
  # Output (boolean): 
  #   A boolean object denoting if the setting value is a matrix.
  ##############################################################################
  
  # If it is a matrix it should start with "[[" and end with "]]"
  if(startsWith(val.str, "[[") & endsWith(val.str, "]]"))
  {
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

# Helper function to recognize if string encodes a vector
is.vector_setting <- function(val.str)
{
  ##############################################################################
  # val.str (string):
  #   A string obtained from the settings file. Contains the value of the 
  #   setting. Vectors start with "[" and end with "]"
  ##############################################################################
  # Output (boolean): 
  #   A boolean object denoting if the setting value is a vector
  ##############################################################################
  
  # If it is a vector it should start with "[" and end with "]"
  if( grepl("^\\[(?!\\[)", val.str, perl = TRUE) & 
      grepl("(?<!\\])\\]$", val.str, perl = TRUE))
  {
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

# Helper function to recognize if string encodes a boolean
is.boolean_setting <- function(val.str)
{
  ##############################################################################
  # val.str (string):
  #   A string obtained from the settings file. Contains the value of the 
  #   setting. Boolean values are either FALSE or TRUE.
  ##############################################################################
  # Output (boolean): 
  #   A boolean object denoting if the setting value is a boolean
  ##############################################################################
  
  # If it is a boolean it should either be true or false, the encoded will allow
  # for case insensitivity
  if(toupper(val.str) %in% c("TRUE", "FALSE", "T", "F"))
  {
    return(TRUE)
  }
  else
  {
    return(FALSE)
  }
}

# Helper function to recognize if string encodes a single numeric value
is.numeric_setting <- function(val.str)
{
  ##############################################################################
  # val.str (string):
  #   A string obtained from the settings file. Contains the value of the 
  #   setting. Numeric values contain "1-9", "-" and/or "."
  ##############################################################################
  # Output (boolean): 
  #   A boolean object denoting if the setting value is a single numeric
  ##############################################################################
  
  if(!is.na(suppressWarnings(as.numeric(val.str))))
  {
    return(TRUE)
  }
  else
  {
    return(FALSE)
  }
}

# Helper function for converting text matrix setting to matrix
setting_to_matrix <- function(setting.str)
{
  ##############################################################################
  # setting.str (string):
  #   A string to be converted into a matrix with encoding [[a,b],[c,d]]
  ##############################################################################
  # Output (matrix): 
  #   The matrix containing the value of the setting.
  ##############################################################################
  
  # Split setting string into arrays
  setting_arrays <-  unlist(stringr::str_split(setting.str, pattern = "\\],\\["))
  
  # Clean the arrays
  setting_arrays <- str_remove_all(setting_arrays, "\\[|\\]")
  
  # Split according to ","
  setting_arrays <- stringr::str_split(setting_arrays, ",")
  
  # Determine type of data
  # Is numeric?
  if(all(grepl("-*\\d+\\.*\\d*", unlist(setting_arrays))))
  {
    setting.mat <- matrix(as.numeric(unlist(setting_arrays)),
                           nrow = length(setting_arrays),
                           byrow = TRUE)
  }
  # Is boolean?
  else if(all(sapply(unlist(setting_arrays), is.boolean_setting)))
  {
    setting.mat <- matrix(as.logical(toupper(unlist(setting_arrays))),
                          nrow = length(setting_arrays),
                          byrow = TRUE)
  }
  # If not boolean or numeric then it must be character
  else
  {
    setting.mat <- matrix(unlist(setting_arrays),
                          nrow = length(setting_arrays),
                          byrow = TRUE)
  }
  
  # Return setting in matrix form
  return(setting.mat)
  
}

read_settings <- function(settings_file=NULL)
{
  ##############################################################################
  # settings_file (string):
  #   An object containing the benchmark file to read in and set up.
  #   For information on the formats used in the file see:
  #   the example_benchmark_setting.txt file located in the documentation folder
  #   For even more information on this see the section "User guide" in the thesis
  #   contained in the folder documentation.
  ##############################################################################
  # Output (list): 
  #   A list containing both the simulation and method settings.
  ##############################################################################
  
  if(!is.character(settings_file) | !endsWith(settings_file, ".txt")){
    stop("Input settings file should be a string ending with .txt")
  }
  
  # Initialize settings list
  settings.list <- list()
  
  # Open file for reading
  conn <- file(settings_file, open="r")
  
  # Read through the file line per line and turn the contents into usable settings
  # Do not keep lines that contain no character i.e. white lines
  linn <- readLines(conn) # Get all lines in the settings file
  linn <- linn[linn != ""] # Remove empty lines
  
  for(i in 1:length(linn))
  {
    # Read line i
    line <- linn[i]
    
    # Split line into chunks each containing a setting
    line_split <- unlist(strsplit(line, split = "\t"))
    
    
    # Run through each setting chunck and transform it into the appropriate type
    # And add it to the collector of the settings
    for(j in 1:length(line_split))
    {
      # Read setting and split according to "=", the separating symbol between 
      # name and value
      setting <- unlist(strsplit(line_split[j], "="))
      
      # Each setting has a name and a value
      setting_name <- setting[1]
      setting_value <- str_replace_all(setting[2], fixed(" "), "") # Remove all white space
      
      # Determine what type of setting it is and add it in the appropriate type 
      # to the collection of settings
      if(j == 1)
      {
        settings.list[[setting_name]] <- list("name" = setting_value)
      }
      # Matrix setting?
      else if(is.matrix_setting(setting_value))
      {
        # Matrix format example [[1,2],[3,4]]
        settings.list[[i]][[setting_name]] <- setting_to_matrix(setting_value) 
      }
      # Boolean setting?
      else if(is.boolean_setting(setting_value))
      {
        settings.list[[i]][[setting_name]] <- as.logical(toupper(setting_value))
      }
      # Numeric setting?
      else if(is.numeric_setting(setting_value))
      {
        settings.list[[i]][[setting_name]] <- as.numeric(setting_value)
      }
      # Vector setting?
      else if(is.vector_setting(setting_value))
      {
        # Vector format example [1,2,3,4]
        settings.list[[i]][[setting_name]] <- 
          as.vector(setting_to_matrix(setting_value))
      }
      # String setting?
      else
      {
        settings.list[[i]][[setting_name]] <- setting_value
      }
    }
  }
  
  # Close file 
  close(conn)
  
  # Return settings list
  return(settings.list)
}
