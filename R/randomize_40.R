# TOOLS & DATA
  require(here)
  source(here::here('R/tools.R'))

# load data  
  d = data.table(
      f = c(list.files(path = here::here('/Pictures/test_40/original/'), pattern = '.jpg', recursive = TRUE, full.names = TRUE)),
      f2 = c(list.files(path = here::here('/Pictures/test_40/original/'), pattern = '.jpg', recursive = TRUE, full.names = FALSE)),
      id = as.character(sample.int(n = 40, size = 40))
      )
  d[nchar(id) == 1, id := paste0(0,id)]


# rename, move & save
  for(i in 1:40){
    #i = 1 
    file.copy(from = d$f[i], to = glue('Pictures/test_40/renamed/',d$id[i], '.jpg')) 
  }

  fwrite(d, file = 'Data/test_40_randomized.csv')