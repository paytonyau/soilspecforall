# Library -----
if(!require(dplyr)){install.packages('dplyr'); require(dplyr)}
if(!require(tidyr)){install.packages('tidyr'); require(tidyr)}
if(!require(readr)){install.packages('readr'); require(readr)}
if(!require(stringr)){install.packages('stringr'); require(stringr)}
if(!require(tibble)){install.packages('tibble'); require(tibble)}


# Setting the folder path ----
dir.p = 'C:/BSSS_project_lowcost_nir/database'
#C:\BSSS_project_lowcost_nir\database

# Importing all .csv files and merging them ----
getfiles = list.files(path = paste0(dir.p, "/raw_data/spectra_usa_ge/isc/"),
                      pattern="_r.csv", 
                      full.names = T)


# Function for getting the files and creating dataframe
combineCSVFiles = function(filePaths){
  combinedData = data.frame() # empty dataframe
  
  for (filePath in filePaths){
    
    data = read_csv(filePath, id = "path", skip = 28) # read each CSV file
    
    combinedData = rbind(combinedData, data) # Combine the data
  }
  return(combinedData)
}

# Using the function
combinedDF = combineCSVFiles(getfiles)

# Checking out the final DF
print(combinedDF)
View(combinedDF)
# Creating ID column using named ID during recording spec
##combinedDF$id = sub("\\..*$", "", basename(combinedDF$path))

combinedDF <- combinedDF %>%
  dplyr::mutate(id = str_extract(sub("\\..*$", "", basename(combinedDF$path)),
                                 "(?<=_)(\\d+)(?=\\_)"))

combinedDF$id = paste("isc", combinedDF$id, sep = "_") # Check if it is necessary or not

combinedDF = combinedDF %>%
  dplyr::select(c('id'), everything())

combinedDF.ok = combinedDF[, -2]
colnames(combinedDF.ok)[2:3] = c("wavelength", "reflectance")

View(combinedDF.ok)
combinedDF.ok[combinedDF.ok$id == c("00000", "00180"),]

# Create a new column by extracting the numbers after the last comma
df <- combinedDF.ok %>%
  mutate(NewColumn = paste0(",", sub(".+,(\\d+)$", "\\1", `reflectance`)))

#df$NewColumn = substr(combinedDF.ok$reflectance, 12, 22)
df[df$id == "00180" & df$wavelength == "1056",]
df[df$id == "00180" & df$wavelength == "1372",]
df[df$id == "00180" & df$wavelength == "1662",]
df[df$id == "00180" & df$wavelength == "1700",]
df[df$id == "00180",]
str(df)

#df$NewColumn = gsub(",0,", "0.", df$NewColumn)

df$NewColumn = as.numeric(gsub(",", ".", df$NewColumn))

#df$NewColumn = as.numeric(df$NewColumn)

df.ok = df[, -3]
df.ok[1:10,]

colnames(df.ok)[3] = "reflectance"

#str(df.ok)
#df.ok$NewColumn = as.numeric(df.ok$NewColumn)

#df.ok[df.ok$id == "00180", ]
#df.ok[df.ok$id == "00080", ]

#combined.df.ok = df.ok[, -3]
#combined.df.ok[1:5, ]

# Reshaping DF from long to wide format
combined.df.ok = df.ok

colnames(combined.df.ok)
combinedDF.wide = tidyr::spread(combined.df.ok,
                                key = "id",
                                value = "reflectance")

row_names <- colnames(combinedDF.wide)[-1]

transposed_tibble <- as_tibble(t(as.matrix(combinedDF.wide)))

colnames(transposed_tibble) = transposed_tibble[1,]
transposed_tibble[2:5, c(1:4, ncol(transposed_tibble))]

combined.df.org = transposed_tibble[-1, ]
combined.df.org[1:5, c(1:5, ncol(combined.df.org))]
combined.df.org$id = row_names

combined.df.org = combined.df.org %>%
  select(id, everything())

combined.df.org[1:5, c(1:5, ncol(combined.df.org))]

# Exporting dataframe
write.csv(combined.df.org, paste0(dir.p, "/raw_data/spectra_usa_ge/isc/combinedSpec.csv"))

# Saving R environment
save.image(paste0(dir.p, "/r_code/isc_importing_and_cleaning.RData"))


