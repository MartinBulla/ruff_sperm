# TOOLS & DATA
  require(here)
  source(here::here('R/tools.R'))

  v = data.table(read_excel(here::here('Data/ruff_sperm_Vancouver_2018.xlsx'), sheet = 1))
  v = v[!is.na(sample_ID)]
  x = v[!duplicated(bird_ID)]
  
  s = data.table(read_excel(here::here('Data/ruff_males_Seewiesen.xlsx'), sheet = 1))

  d = data.table(read_excel(here::here('Data/ruff_vas-deferens.xlsx'), sheet = 1))

# EXPLORE
  nrow(v)
  
  nrow(x)
  table(x$morph)

  nrow(s[Ind_ID%in%unique(v$bird_ID)])
  table(s$Morph)


  sv = x[bird_ID%in%s$Ind_ID]
  table(sv$morph)


  nrow(d[ID%in%unique(v$bird_ID)])

  Seewiesen
   N = 100 males
   F 8
   I 64
   S 28
  
  Vancouver sample
   N = 86 males
   F 15
   I 51
   S 20

  Vancouver same Seewiesen
    N = 40 males
    F 4
    I 26
    S 10


  