#--------------------------------------------------------------------------------
# setup.R         Setup project
#--------------------------------------------------------------------------------

# 1. set computer --------------------------------------------------------------
user <- "david" 



# 2. Set main data paths -------------------------------------------------------
if(user == "david") main_dir <- "C:/Users/david/SML Dropbox/gitdata/eDNA_Biodiversity"


# set main working directory
setwd(main_dir)


# 3. Create data paths --------------------------------------------------------- 
input_data <- paste(main_dir, "input", sep="/")
if (!dir.exists(input_data)) dir.create(input_data, recursive = TRUE)

temp_data <- paste(main_dir, "temp", sep="/")
if (!dir.exists(temp_data)) dir.create(temp_data, recursive = TRUE)

output_data <- paste(main_dir, "output", sep="/")
if (!dir.exists(output_data)) dir.create(output_data, recursive = TRUE)
