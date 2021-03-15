# TOOLS & DATA
  require(here)
  source(here::here('R/tools.R'))
  library(ggpubr)
  require(arm)
  require(rptR) 
  require(magrittr)
  require(viridis)
  require(gridExtra)
  require(stringi)

  # function for R output based on sim
  R_out = function(name = "define", model = m, nsim = 5000){
   bsim <- sim(model, n.sim=nsim)  
   l=data.frame(summary(model)$varcor)
   l = l[is.na(l$var2),]
   l$var1 = ifelse(is.na(l$var1),"",l$var1)
   l$pred = paste(l$grp,l$var1)

   q50={}
   q025={}
   q975={}
   pred={}
   
   # variance of random effects
   for (ran in names(bsim@ranef)) {
     #ran =names(bsim@ranef)[1]
     ran_type = l$var1[l$grp == ran]
     for(i in ran_type){
        # i = ran_type[2]
      q50=c(q50,quantile(apply(bsim@ranef[[ran]][,,i], 1, var), prob=c(0.5)))
      q025=c(q025,quantile(apply(bsim@ranef[[ran]][,,i], 1, var), prob=c(0.025)))
      q975=c(q975,quantile(apply(bsim@ranef[[ran]][,,i], 1, var), prob=c(0.975)))
      pred= c(pred,paste(ran, i))
      }
     }
   # residual variance
   q50=c(q50,quantile(bsim@sigma^2, prob=c(0.5)))
   q025=c(q025,quantile(bsim@sigma^2, prob=c(0.025)))
   q975=c(q975,quantile(bsim@sigma^2, prob=c(0.975)))
   pred= c(pred,'Residual')

   ci = c(round(100*q025/sum(q025))[1], round(100*q975/sum(q975))[1])
   ci = ci[order(ci)]
   
   ri=data.table(model = name, repeatability=paste0(round(100*q50/sum(q50)),'%')[1], CI = paste0(paste(ci[1], ci[2], sep ="-"), '%'))
   
   
   return(ri)
   }

# load data  
  d = data.table(read_excel('Data/VD_measurements 2021_26_01.xlsx', sheet = 1))#, range = "A1:G161"))
  d[,Length_um :=as.character(Length_um) ]
  stri_sub(d$Length_um,nchar(d$Length_um)-2, nchar(d$Length_um)-3) = "."
  d[,Length_um :=as.numeric(Length_um) ]
  
  #d[Part%in%c('Acrosome','Head','Midpeace','Tail'), Total_sperm_adding_parts := sum(Length_um), by = Sperm_ID]


  a = d[Part%in%c('Acrosome','Head','Midpiece','Tail')]
  a[, Total_sperm_adding_parts := sum(Length_um), by = Sperm_ID]
  a = a[Part == 'Head']

  d1 = d[Part == 'Total sperm']
 names(d1)[8] = 'Total_sperm'
 
  a = merge(d1,a[,.(Sperm_ID, Total_sperm_adding_parts)], by = "Sperm_ID")

  g1 = ggplot(a, aes(x = Total_sperm_adding_parts, y = Total_sperm)) + 
    stat_smooth(method = 'lm')+geom_point()+
    stat_cor(method="pearson",size = 2) +
    geom_abline(b = 1, col = 'red', lty = 3) +
    xlim(c(485, 565)) +  ylim(c(485, 565)) + 
    ggtitle('All')+
    theme_bw()

  
  a = d[Part%in%c('Head','Midpiece','Tail')]
  a[, Total_sperm_adding_HMT := sum(Length_um), by = Sperm_ID]
  a = a[Part == 'Head']

  d1 = d[Part == 'Total sperm_no_ Ac']
 names(d1)[8] = 'Total_sperm_no_Ac'
 
  a = merge(d1,a[,.(Sperm_ID, Total_sperm_adding_HMT)], by = "Sperm_ID")

  g2 = ggplot(a, aes(x = Total_sperm_adding_HMT, y = Total_sperm_no_Ac)) + 
    stat_smooth(method = 'lm')+geom_point()+
    stat_cor(method="pearson",size = 2) +
    geom_abline(b = 1, col = 'red', lty = 3) +
     xlim(c(485, 565)) +  ylim(c(485, 565)) + 
    ggtitle('no_acro')+
    theme_bw()


  a = d[Part%in%c('Acrosome','Total sperm_no_ Ac')]
  a[, Total_sperm_adding_ac := sum(Length_um), by = Sperm_ID]
  a = a[Part == 'Acrosome']

  d1 = d[Part == 'Total sperm']
  names(d1)[8] = 'Total_sperm'
 
  a = merge(d1,a[,.(Sperm_ID, Total_sperm_adding_ac)], by = "Sperm_ID")

  g3 = ggplot(a, aes(x = Total_sperm_adding_ac, y = Total_sperm)) + 
    stat_smooth(method = 'lm')+geom_point()+
    stat_cor(method="pearson",size = 2) +
    geom_abline(b = 1, col = 'red', lty = 3) +
     xlim(c(485, 565)) +  ylim(c(485, 565)) + 
     ggtitle('no_acro + acro')+
    theme_bw()

  grid.arrange(g1,g2,g3)

  ggsave('Output/adding_measurements_error.png',arrangeGrob(g1,g2,g3,nrow = 3))  