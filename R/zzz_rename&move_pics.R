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
  #d[nchar(file_name)<32]
  d[sample_ID %in% c('400'), file_name := paste(substring(file_name,1,8), '2021-04-19', substring(file_name,10))]
  d[sample_ID %in% c('401'), file_name := paste(substring(file_name,1,9), '2021-04-19', substring(file_name,11))]
  d[sample_ID %in% c('402'), file_name := paste(substring(file_name,1,8), '2021-04-20', substring(file_name,11))]
  
  #paste(substring('Ruff 320 Snap-591.jpg',1,8), '2021-04-20', substring('Ruff 320 Snap-591.jpg',10))

  #dd = d[sample_ID %in% c('400')]
  d[, new_name := paste0("Ruff ", sample_ID, "_", substring(file_name, 6))]
  d[substring(new_name, nchar(new_name)-8, nchar(new_name)-7) == "p-", new_name := gsub("p-", "p-0", new_name)]
  
  #d$new_name
  
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
  #d[substring(new_name, nchar(new_name)-8, nchar(new_name)-7) == "p-", new_name := gsub("p-", "p-0", new_name)]
  
  
  # rename & copy
  for(i in 1:nrow(d)){
    #i = 1 
    file.copy(from = d$f[i], to = glue('all_photos/',d$new_name[i])) 
  }

 # test
  d = data.table(
    f = c(list.files(path = here::here('all_photos'), pattern = '.jpg', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here('all_photos'), pattern = '.jpg', recursive = TRUE, full.names = FALSE))
  )
  d[ , sample_ID := substring(f2,6,8)]

  nrow(d)
  unique(d$sammple_ID)

  