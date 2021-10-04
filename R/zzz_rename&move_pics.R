# TOOLS & DATA
  require(here)
  source(here::here('R/tools.R'))

# April_MAY_2021 
  d = data.table(
    f = c(list.files(path = here::here('April_May_2021/Photos_PRE_sample_ID_MAY/'), pattern = '.jpg', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here('April_May_2021/Photos_PRE_sample_ID_MAY/'), pattern = '.jpg', recursive = TRUE, full.names = FALSE))
    )
  d = d[!grepl('metadata', f2, fixed = TRUE)]
  d[ , sample_ID := substring(f2,6,8)]
  d[, file_name :=sub(".*/", "", f2)]
  d[nchar(file_name)<32]
  d[sample_ID %in% c('400','401'), file_name := paste(substring(file_name,1,8), '2021-04-19', substring(file_name,10), sep =" ")]
  d[sample_ID %in% c('402'), file_name := paste(substring(file_name,1,8), '2021-04-20', substring(file_name,10), sep =" ")]
  
  d[, new_name := paste0("Ruff ", sample_ID, "_", substring(file_name, 6))]
  d[substring(new_name, nchar(new_name)-7, nchar(new_name)-7) == "-", file_name := gsub("-", "-0", file_name)]
  
  # rename & copy
  for(i in 1:nrow(d)){
    #i = 1 
    file.copy(from = d$f[i], to = glue('all_photos/',d$new_name[i])) 
  }

# June_2021 
  d = data.table(
    f = c(list.files(path = here::here('June_2021/Photos_PRE_sample_ID_no'), pattern = '.jpg', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here('June_2021/Photos_PRE_sample_ID_no'), pattern = '.jpg', recursive = TRUE, full.names = FALSE))
  )
  d = d[!grepl('metadata', f2, fixed = TRUE)]
  d[ , sample_ID := substring(f2,6,8)]
  d[, file_name :=sub(".*/", "", f2)]
  #d[nchar(file_name)<32]
  
  d[, new_name := paste0("Ruff ", sample_ID, "_", substring(file_name, 6))]
  
  # rename & copy
  for(i in 1:nrow(d)){
    #i = 1 
    file.copy(from = d$f[i], to = glue('all_photos/',d$new_name[i])) 
  }
  