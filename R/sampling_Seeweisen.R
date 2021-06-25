  require(here)
  source(here::here('R/tools.R'))
  require(stringi)
  require(arm)
  require(effects)
  require(RColorBrewer)
  require("PerformanceAnalytics")

  library(ggpubr)
  
  require(rptR) 
  require(magrittr)


  m = data.table(read_excel('Data/ruff_males_Seewiesen.xlsx', sheet = 1))#, range = "A1:G161"))
  m[, hatch_year:=as.numeric(substr(Hatch_date,1,4)) ]
  m[, age := 2021-hatch_year]
  d = data.table(read_excel('Data/sampling_2021-05-03_05.xlsx', sheet = 1))#, range = "A1:G161"))
  unique(d$sperm)
  unique(d$recording)
  x = d[sperm %in% c('sperm', "super","yes","ok", 'some', 'few','moderate'), DB_ID]
  r = d[recording %in% c('good', "super","yes","ok", 'great'), DB_ID]

# needed
  b = unique(d$DB_ID)
  mx = data.table(DB_ID = b[!b%in%x])
  mx$morph = m$Morph[match(mx$DB_ID, m$Ind_ID)]
  mx
  nrow(mx)
  write.csv(file = 'Data/to_sample.csv', mx, row.names = FALSE)

  n = m[!Ind_ID%in%x]
  nrow(n) 
  summary(factor(n$Morph)) 

# have
  h =   m[Ind_ID%in%x]
  summary(factor(h$Morph))
  nrow(h)
  h[order(Ind_ID)]

# exploration
 m$for_morpho = ifelse(m$Ind_ID%in%x, 'yes', 'no')
 m$for_velo  = ifelse(m$Ind_ID%in%r, 'yes', 'no')
 table(m$for_morpho)
 table(m$for_morpho, m$Morph)

 g1 = ggplot(m, aes(x = for_morpho, y = age, fill = Morph)) + geom_boxplot() + ylim(0,15)
 g2 = ggplot(m, aes(x = for_morpho, y = age, fill = Morph)) + geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", dotsize = 0.6) + ylim(0,15)
 #ggplot(m[nchar(Ind_ID)<5], aes(x = for_morpho, y = age, fill = Morph)) + geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", dotsize = 0.4) + ylim(0,15)

 ggsave('Output/sperm_sampling.png', arrangeGrob(g1,g2), width = 3.5, height = 4.5)

 table(m$for_velo)
 table(m$for_velo, m$Morph)
 #g1 = ggplot(m, aes(x = for_velo, y = age, fill = Morph)) + geom_boxplot() + ylim(0,15)
 #g2 = ggplot(m, aes(x = for_velo, y = age, fill = Morph)) + geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", dotsize = 0.6) + ylim(0,15)
 #ggplot(m[nchar(Ind_ID)<5], aes(x = for_morpho, y = age, fill = Morph)) + geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", dotsize = 0.4) + ylim(0,15)

 ggsave('Output/sperm_sampling.png', arrangeGrob(g1,g2), width = 3.5, height = 4.5)



mV = m[nchar(Ind_ID)<5]
table(mV$Morph, mV$for_morpho)
table(mV$for_morpho)

mO = m[nchar(Ind_ID)>4]
table(mO$Morph, mO$for_morpho)
table(mO$for_morpho)
