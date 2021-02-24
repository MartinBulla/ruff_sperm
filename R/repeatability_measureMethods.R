# TOOLS & DATA
  require(here)
  source(here::here('R/tools.R'))
  library(ggpubr)
  require(arm)
  require(rptR) 
  require(magrittr)
  require(viridis)
  require(gridExtra)

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
  d = data.table(read_excel('Data/Measurement Validation.xlsx', sheet = 1, range = "A1:G161"))
  d[,date_ := as.factor(Date)]

  dt = d[part %in% c('tail'), ]
  dtm = d[part %in% c('tail+midpiece'), .(Date,date_, sperm_ID, photo_ID, length) ]
  dtm[, tm := length]
  dtm$length =NULL
  dx = merge(dt,dtm)
  dx[, part := 'midpiece_derived']
  dx[,length:=tm-length ]
  dx$tm = NULL
  d = rbind(d,dx)

  d=d[order(Date,method,sperm_ID, part)]

  df = d[method == 'free']
  da = df[, .(sperm_ID, part, length)][1:50] 
  setnames(da, old = names(da), new = c('sperm_ID','part', 'attempt_A'))
  db = df[, .(sperm_ID, part, length)][51:nrow(df)]
  setnames(db, old = names(db), new = c('sperm_ID','part', 'attempt_B'))
  dfAB = merge(da,db)

  ds = d[method == 'segm']
  da = ds[, .(sperm_ID, part, length)][1:50] 
  setnames(da, old = names(da), new = c('sperm_ID','part', 'attempt_A'))
  db = ds[, .(sperm_ID, part, length)][51:nrow(df)]
  setnames(db, old = names(db), new = c('sperm_ID','part', 'attempt_B'))
  dsAB = merge(da,db)

  dfAB[, method := 'free']  
  dsAB[, method := 'segm']  
  dAB = rbind(dsAB, dfAB) 

# plot 
  ggplot(d, aes(x = method, y = length, fill = method )) + geom_boxplot()+facet_wrap(~part, scales = 'free')  + scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5, guide = FALSE) 
  ggsave('Output/measurements_boxplots.png')

  # free
   ggplot(dfAB, aes(x = attempt_A, y = attempt_B)) + 
    stat_smooth(method = 'lm')+geom_point()+
    stat_cor(method="pearson",size = 2) +
    geom_abline(b = 1, col = 'red', lty = 3) + 
    facet_wrap(~part, scales = 'free') + theme_bw()

  # seg
   ggplot(dsAB, aes(x = attempt_A, y = attempt_B)) + 
    stat_smooth(method = 'lm')+geom_point()+
    stat_cor(method="pearson",size = 2) +
    geom_abline(b = 1, col = 'red', lty = 3) + 
    facet_wrap(~part, scales = 'free') + theme_bw() +

  # both
   ggplot(dAB, aes(x = attempt_A, y = attempt_B, fill = method, col = method)) + 
    stat_smooth(method = 'lm')+geom_point()+
    stat_cor(aes(col = method), method="pearson",size = 2) +
    geom_abline(b = 1, col = 'red', lty = 3) + 
    facet_wrap(~part, scales = 'free') + theme_bw() +
    scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
    scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) 
    ggsave('Output/measurements_cor.png')
 
# repeatability
  lfsim = list()
  lfrpt = list()
  lssim = list()
  lsrpt = list()
  for(i in unique(d$part)){
    part_ = i
    # part_ = "head"
    dd = df[part == part_]
    m = lmer(length ~ 1+(1|sperm_ID), dd)
    Rf = R_out(part_)
    lfsim[[i]] = Rf[, method_CI:='arm package']

    R = rpt(length ~ (1 | sperm_ID), grname = "sperm_ID", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
    RR = data.table(merge(data.frame(name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
    RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
    lfrpt[[i]] =  RR[, method_CI := 'rpt package']

    dds = ds[part == part_]
    ms = lmer(length ~ 1+(1|sperm_ID), dds)
    Rs = R_out(name = part_, model = ms)
    lssim[[i]] = Rs[, method_CI:='arm package']

    Rss= rpt(length ~ (1 | sperm_ID), grname = "sperm_ID", data = dds, datatype = "Gaussian")#, nboot = 0, npermut = 0)
    RssR = data.table(merge(data.frame(name =part_), paste0(round(Rss$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
    RssR[, CI := paste0(paste(round(Rss$CI_emp*100)[1], round(Rss$CI_emp*100)[2], sep = "-"), '%')] 
    lsrpt[[i]] =  RssR[, method_CI := 'rpt package']
    print(i)
  }
  fsim = do.call(rbind,lfsim)
  frpt = do.call(rbind,lfrpt)
  ssim = do.call(rbind,lssim)
  srpt = do.call(rbind,lsrpt)

  fsim[ , measurement := 'free']
  ssim[ , measurement := 'segm']
  x = rbind(fsim,ssim)
  names(x)[1] = "part"
  x[, pred:= as.numeric(substr(repeatability,1,2))]
  x[, lwr:= as.numeric(substr(CI,1,2))]
  x[, upr:= as.numeric(substr(CI,4,5))]

  frpt[ , measurement := 'free']
  srpt[ , measurement := 'segm']
  y = rbind(frpt,srpt)
  names(y)[1] = "part"
  y[, pred:= as.numeric(substr(Repeatability,1,2))]
  y[nchar(CI) == 5, CI := paste0(0,CI) ]
  y[, lwr:= as.numeric(substr(CI,1,2))]
  y[, upr:= as.numeric(substr(CI,4,5))]
# plot estimated reputabilities
  # from sim
    g1 = ggplot(x, aes(x = part, y = pred, col = measurement)) +
      geom_errorbar(aes(ymin = lwr, ymax = upr, col = measurement), width = 0.1, position = position_dodge(width = 0.25) ) +
      ggtitle ("Sim based")+
      geom_point(position = position_dodge(width = 0.25)) +
      scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
      scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
      labs(x = NULL, y = "Repeatability [%]")+
      ylim(c(0,100))+
      theme_bw() +
      theme(plot.title = element_text(size=9))
  
  # from rpt 
    g2 = ggplot(y, aes(x = part, y = pred, col = measurement)) +
      geom_errorbar(aes(ymin = lwr, ymax = upr, col = measurement), width = 0.1, position = position_dodge(width = 0.25) ) +
      geom_point(position = position_dodge(width = 0.25)) +
      ggtitle ("Rpt based") + 
      scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
      scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
      labs(x = NULL, y = "Repeatability [%]")+
      ylim(c(0,100))+
      theme_bw()+
      theme(plot.title = element_text(size=9))

  grid.arrange(g1,g2)

  ggsave('Output/measurement_R.png',arrangeGrob(g1,g2,nrow = 2))