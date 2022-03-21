# TOOLS & DATA
  require(here)
  source(here::here('R/tools.R'))

  # CHECK WHICH package has %do%
  packages = c('anytime','data.table', 'DataEntry.validation', 'DT', 'foreach', 'ggplot2', 'ggthemes', 'glue','googledrive', 'googlesheets4', 'grid', 'htmlTable', 'lattice', 'lubridate', 'magrittr', 'maptools', 'openxlsx','plyr','raster','readxl','stringr','zoo')
> sapply(packages, function(x) suppressPackageStartupMessages(using(x)) )
 
# DATA 
  m1 = data.table(
    f = c(list.files(path = here::here('Measurements/'), pattern = '.erT', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here('Measurements/'), pattern = '.erT', recursive = TRUE, full.names = FALSE))
    )
  m1[, pic := substr(m1$f2, 20, 22)]
  m2 = data.table(
    f = c(list.files(path = here::here('Measurements/'), pattern = '.erA', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here('Measurements/'), pattern = '.erA', recursive = TRUE, full.names = FALSE))
    )
  m2[, pic := substr(m2$f2, 20, 22)]

  d = data.table(
    f = c(list.files(path = here::here('Measurements/'), pattern = '.csv', recursive = TRUE, full.names = TRUE)),
    f2 = c(list.files(path = here::here('Measurements/'), pattern = '.csv', recursive = TRUE, full.names = FALSE))
    )
  d[ , measured := substring(d$f2,21,39)]

  b = foreach(j = 1:nrow(d), .combine = rbind) %do% {
       #j =1
       ff = d[j, ]
       
       x = fread(ff$f, stringsAsFactors = FALSE) #skip = skip, col.names = varnames, colClasses = colclass, 
       x[, measured := ff$measured]
       x[, pic := substr(x$'File Name', 1, 3)]
       x[, part := rep(c('Acrosome','Nucleus','Midpiece','Tail'),50)]
       x[pic %in% m1$pic & part == 'Tail', manip := 'erT']
       x[pic %in% m2$pic & part == 'Acrosome', manip := 'erA']
       return(x)
   }