# TOOLS & DATA
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

  round_ = 3 # number of decimal places to round model coefficients
  nsim = 5000 # number of simulations to extract estimates and 95%CrI
  ax_lines = "grey60" # defines color of the axis lines
  colors <- c("#999999", "#E69F00", "#56B4E9") #viridis(3)

  # function 
  # for Repeatability output based on sim
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
  # mode output
   m_out = function(model = m, type = "mixed", 
      name = "define", dep = "define", fam = 'Gaussian',
      round_ = 3, nsim = 5000, aic = TRUE, save_sim = FALSE, N = NA, back_tran = FALSE, perc_ = 1){
        # perc_ 1 = proportion or 100%
      bsim = sim(model, n.sim=nsim)  
      
      if(save_sim!=FALSE){save(bsim, file = paste0(save_sim, name,'.RData'))}
     
      if(type != "mixed"){
       v = apply(bsim@coef, 2, quantile, prob=c(0.5))
       ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 

       if(back_tran == TRUE & fam == "binomial"){
        v = perc_*plogis(v)
        ci = perc_*plogis(ci)
       }
      if(back_tran == TRUE & fam == "binomial_logExp"){
            v = perc_*(1-plogis(v))
            ci = perc_*(1-plogis(ci))
            ci = rbind(ci[2,],ci[1,])
           }

       if(back_tran == TRUE & fam == "Poisson"){
        v = exp(v)
        ci = exp(ci)
       }

       oi=data.frame(type='fixed',effect=rownames(coef(summary(model))),estimate=v, lwr=ci[1,], upr=ci[2,])
        rownames(oi) = NULL
        oi$estimate_r=round(oi$estimate,round_)
        oi$lwr_r=round(oi$lwr,round_)
        oi$upr_r=round(oi$upr,round_)
        if(perc_ == 100){
         oi$estimate_r = paste0(oi$estimate_r,"%")
         oi$lwr_r = paste0(oi$lwr_r,"%")
         oi$upr_r = paste0(oi$upr_r,"%")
        }
       x=data.table(oi[c('type',"effect", "estimate_r","lwr_r",'upr_r')]) 
     
      }else{
       v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
       ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 

       if(back_tran == TRUE & fam == "binomial"){
        v = perc_*plogis(v)
        ci = perc_*plogis(ci)
       }
      if(back_tran == TRUE & fam == "binomial_logExp"){
            v = perc_*(1-plogis(v))
            ci = perc_*(1-plogis(ci))
            ci = rbind(ci[2,],ci[1,])
           }

       if(back_tran == TRUE & fam == "Poisson"){
        v = exp(v)
        ci = exp(ci)
       }

       oi=data.frame(type='fixed',effect=rownames(coef(summary(model))),estimate=v, lwr=ci[1,], upr=ci[2,])
          rownames(oi) = NULL
          oi$estimate_r=round(oi$estimate,round_)
          oi$lwr_r=round(oi$lwr,round_)
          oi$upr_r=round(oi$upr,round_)
          if(perc_ == 100){
           oi$estimate_r = paste0(oi$estimate_r,"%")
           oi$lwr_r = paste0(oi$lwr_r,"%")
           oi$upr_r = paste0(oi$upr_r,"%")
          }
       oii=oi[c('type',"effect", "estimate_r","lwr_r",'upr_r')] 
      
       l=data.frame(summary(model)$varcor)
       l = l[is.na(l$var2),]
       l$var1 = ifelse(is.na(l$var1),"",l$var1)
       l$pred = paste(l$grp,l$var1)

       q050={}
       q025={}
       q975={}
       pred={}
       
       # variance of random effects
       for (ran in names(bsim@ranef)) {
         ran_type = l$var1[l$grp == ran]
         for(i in ran_type){
          q050=c(q050,quantile(apply(bsim@ranef[[ran]][,,ran_type], 1, var), prob=c(0.5)))
          q025=c(q025,quantile(apply(bsim@ranef[[ran]][,,ran_type], 1, var), prob=c(0.025)))
          q975=c(q975,quantile(apply(bsim@ranef[[ran]][,,ran_type], 1, var), prob=c(0.975)))
          pred= c(pred,paste(ran, i))
          }
         }
       # residual variance
       q050=c(q050,quantile(bsim@sigma^2, prob=c(0.5)))
       q025=c(q025,quantile(bsim@sigma^2, prob=c(0.025)))
       q975=c(q975,quantile(bsim@sigma^2, prob=c(0.975)))
       pred= c(pred,'Residual')

       ri=data.frame(model = name,type='random %',effect=pred, estimate_r=round(100*q050/sum(q050)), lwr_r=round(100*q025/sum(q025)), upr_r=round(100*q975/sum(q975)))
         rx = ri[ri$effect == 'Residual',]
         if(rx$lwr_r>rx$upr_r){ri$lwr_r[ri$effect == 'Residual'] = rx$upr_r; ri$upr_r[ri$effect == 'Residual'] = rx$lwr_r}
         ri$estimate_r = paste0(ri$estimate_r,'%')
         ri$lwr_r = paste0(ri$lwr_r,'%')
         ri$upr_r = paste0(ri$upr_r,'%')
      
      x = data.table(rbind(oii,ri))
      }
      
      x[1, model := name]                                                                
      x[1, response := dep]                                                                
      x[1, error_structure := fam]                                                                
      x[1, N := N]                                                                

      x=x[ , c('model', 'response', 'error_structure', 'N', 'type',"effect", "estimate_r","lwr_r",'upr_r')] 

      if (aic == TRUE){   
          x[1, AIC := AIC(update(model,REML = FALSE))] 
          }
      if (aic == "AICc"){
          aicc = AICc(model)
          x[1, AICc := aicc] 
      }
      if(type == "mixed"){
        x[1, R2_mar := invisible({capture.output({r2_nakagawa(model)$R2_marginal})})]
        x[1, R2_con := invisible({capture.output({r2_nakagawa(model)$R2_conditional})})]
       }
      x[is.na(x)] = ""
      return(x)
    } 
# load data  
  v = data.table(read_excel('Data/ruff_vas-deferens.xlsx', sheet = 1))#, range = "A1:G161"))
  vv = v[,.(ID,Morph)]

  r = data.table(read_excel('Data/ruff_sperm_2021_test Maggie.xlsx', sheet = 1))#, range = "A1:G161"))
  r[,Length_µm :=as.character(Length_µm) ]
  stri_sub(r$Length_µm,nchar(r$Length_µm)-2, nchar(r$Length_µm)-3) = "."
  r[,Length_µm :=as.numeric(Length_µm) ]
  r[Part == 'Mid piece', Part := 'Midpiece']
  r[Part == 'Flagellum', Part := 'Tail']
  r[Part == 'Head', Part := 'Nucleus']

  d = data.table(read_excel('Data/VD measurements feb_2021.xlsx', sheet = 1))#, range = "A1:G161"))
  d[,Length_µm :=as.character(Length_µm) ]
  stri_sub(d$Length_µm,nchar(d$Length_µm)-2, nchar(d$Length_µm)-3) = "."
  d[,Length_µm :=as.numeric(Length_µm) ]
  d[Part == 'Mid piece', Part := 'Midpiece']
  d[Part == 'Flagellum', Part := 'Tail']
  d[Part == 'Head', Part := 'Nucleus']
  
  d1 = d[,.(Bird_ID, Sample_ID, Sperm_ID, Part, Length_µm)]
  
  d2 = d[, sum(Length_µm), by = list(Bird_ID, Sample_ID, Sperm_ID)]
  names(d2)[4] = 'Length_µm'
  d2[,Part := 'Total']

  d3 = d[Part%in%c('Midpiece', 'Tail'), sum(Length_µm), by = list(Bird_ID, Sample_ID, Sperm_ID)]
  names(d3)[4] = 'Length_µm'
  d3[,Part := 'Flagellum']

  d4 = d[Part%in%c('Acrosome', 'Nucleus'), sum(Length_µm), by = list(Bird_ID, Sample_ID, Sperm_ID)]
  names(d4)[4] = 'Length_µm'
  d4[,Part := 'Head']

  b = merge(d1,d2, all = TRUE)
  b = merge(b,d3, all = TRUE)
  b = merge(b,d4, all = TRUE)
  
  a = b[, mean(Length_µm), by = list(Bird_ID, Part)]
  names(a)[3] = 'Length_avg'
  a1 = data.table(Bird_ID = unique(a$Bird_ID), Blind = c(rep(c('Independent','Satellite', 'Faeder'), 4), 'Independent'))
  a = merge(a,a1)

  a = merge(a, vv, by.x = 'Bird_ID', by.y = 'ID')
  a[Morph == 'Res', Morph := 'Independent']
  a[Morph == 'Sat', Morph := 'Satellite']
  a[Morph == 'Faed', Morph := 'Faeder']

  ax = unique(a[ ,.(Bird_ID, Blind,Morph)])
  ax$Part = 'Midpiece_rel'
  ax$Length_avg = a[Part == 'Midpiece',.(Length_avg)]/a[Part == 'Total',.(Length_avg)]
  a = rbind(a, ax[,.(Bird_ID, Part, Length_avg, Blind, Morph)])
  
  ax = unique(a[ ,.(Bird_ID, Blind,Morph)])
  ax$Part = 'Flagellum_rel'
  ax$Length_avg = a[Part == 'Flagellum',.(Length_avg)]/a[Part == 'Total',.(Length_avg)]
  a = rbind(a, ax[,.(Bird_ID, Part, Length_avg, Blind, Morph)])


  aw = reshape(a[,Blind := NULL], idvar = c('Bird_ID','Morph'), timevar = "Part", direction = "wide")
  names(aw) = c("Bird_ID","Morph",unique(a$Part))


  a[, Part := factor(Part, levels=c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel"))] 
  a[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 

  ax = unique(b[ ,.(Bird_ID, Sample_ID, Sperm_ID)])
  ax$Part = 'Midpiece_rel'
  ax$Length_µm = b[Part == 'Midpiece',.(Length_µm)]/b[Part == 'Total',.(Length_µm)]
  b = merge(b, ax, all = TRUE)

  ax = unique(b[ ,.(Bird_ID, Sample_ID, Sperm_ID)])
  ax$Part = 'Flagellum_rel'
  ax$Length_µm = b[Part == 'Flagellum',.(Length_µm)]/b[Part == 'Total',.(Length_µm)]
  b = merge(b, ax, all = TRUE)

  b = merge(b, vv, by.x = 'Bird_ID', by.y = 'ID')
  b[Morph == 'Res', Morph := 'Independent']
  b[Morph == 'Sat', Morph := 'Satellite']
  b[Morph == 'Faed', Morph := 'Faeder']

  b$Sample_ID = NULL
  bw = reshape(b, idvar = c('Bird_ID',  'Sperm_ID', 'Morph'), timevar = "Part", direction = "wide")
  names(bw) = c("Bird_ID","Sperm_ID","Morph",unique(b$Part))

  b[, Part := factor(Part, levels=c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel"))] 
  b[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 

  # dataset for correlations
    h = b[Part == 'Head']
    names(h)[4]="Head_µm"
    h[, Acrosome_µm := b[Part == 'Acrosome',.(Length_µm)]]
    h[, Nucleus_µm := b[Part == 'Nucleus',.(Length_µm)]]
    h[, Midpiece_µm := b[Part == 'Midpiece',.(Length_µm)]]
    h[, Tail_µm := b[Part == 'Tail',.(Length_µm)]]
    h[, Flagellum_µm := b[Part == 'Flagellum',.(Length_µm)]]
    h[, Total_µm := b[Part == 'Total',.(Length_µm)]]
    h[, MidpieceRel_µm := b[Part == 'Midpiece_rel',.(Length_µm)]]
    h[, FlagellumRel_µm := b[Part == 'Flagellum_rel',.(Length_µm)]]

    ha = a[Part == 'Head']
    names(ha)[3]="Head_µm"
    ha[, Acrosome_µm := a[Part == 'Acrosome',.(Length_avg)]]
    ha[, Nucleus_µm := a[Part == 'Nucleus',.(Length_avg)]]
    ha[, Midpiece_µm := a[Part == 'Midpiece',.(Length_avg)]]
    ha[, Tail_µm := a[Part == 'Tail',.(Length_avg)]]
    ha[, Flagellum_µm := a[Part == 'Flagellum',.(Length_avg)]]
    ha[, Total_µm := a[Part == 'Total',.(Length_avg)]]
    ha[, MidpieceRel_µm := a[Part == 'Midpiece_rel',.(Length_avg)]]
    ha[, FlagellumRel_µm := a[Part == 'Flagellum_rel',.(Length_avg)]]

# explore
  b[, order_ := mean(Length_µm), by = Bird_ID]

  b_ = b[Part =='Total']
  b_[,Bird_ID := reorder(Bird_ID, Length_µm, mean)]
  b[, Bird_ID := factor(Bird_ID, levels = levels(b_$Bird_ID))]
 
  g = ggplot(b, aes(x = reorder(as.factor(Bird_ID),Length_µm,mean), y = Length_µm)) + facet_wrap(~Part, nrow = 7, scales = 'free') + #x = reorder(as.factor(Bird_ID),Length_µm,mean)
    geom_boxplot(aes(col = Morph)) + 
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(), col = 'red', fill ='red')+
    scale_colour_manual(values = colors) + 
    #scale_color_viridis(discrete=TRUE)+
    xlab('Male ID') +
    theme_bw() +
    theme(axis.text.x = element_blank())
    #$, legend.position = "none")

  ggsave('Output/VD_within_male_boxplots_ordered.png', g, width = 20, height = 15, units = 'cm')
  
  
  chart.Correlation(bw[, c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total')], histogram=TRUE, pch=19)
  dev.copy(png,'Output/VD_corr_single.png')
  dev.off()
  
  chart.Correlation(aw[, c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total')], histogram=TRUE, pch=19)
  dev.copy(png,'Output/VD_corr_avg.png')
  dev.off()

  chart.Correlation(bw[Morph == 'Independent', c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total')], histogram=TRUE, pch=19)
  dev.copy(png,'Output/VD_corr_single_I.png')
  dev.off()
  chart.Correlation(bw[Morph == 'Satellite', c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total')], histogram=TRUE, pch=19)
  dev.copy(png,'Output/VD_corr_single_S.png')
  dev.off()
  chart.Correlation(bw[Morph == 'Faeder', c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total')], histogram=TRUE, pch=19)
  dev.copy(png,'Output/VD_corr_single_F.png')
  dev.off()

# repeatability within male
  # estimate
    lfsim = list()
    lfrpt = list()
    for(i in c('Acrosome','Flagellum','Head','Midpiece','Nucleus','Tail','Total')){
      part_ = i
      # part_ = "Acrosome"
      dd = b[Part == part_]
      m = lmer(Length_µm ~ 1+(1|Bird_ID), dd)
      Rf = R_out(part_)
      lfsim[[i]] = Rf[, method_CI:='arm package']

      R = rpt(Length_µm ~ (1 | Bird_ID), grname = "Bird_ID", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
      RR = data.table(merge(data.frame(name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
      RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
      lfrpt[[i]] =  RR[, method_CI := 'rpt package']
      print(i)
    }
    
    x = do.call(rbind,lfsim)
    names(x)[1] = "part"
    y = do.call(rbind,lfrpt)

    x[, pred:= as.numeric(substr(repeatability,1,2))]
    x[, lwr:= as.numeric(substr(CI,1,2))]
    x[, upr:= as.numeric(substr(CI,4,5))]

    y[, pred:= as.numeric(substr(Repeatability,1,2))]
    y[nchar(CI) == 5, CI := paste0(0,CI) ]
    y[, lwr:= as.numeric(substr(CI,1,2))]
    y[, upr:= as.numeric(substr(CI,4,5))]
    names(y)[2] = tolower( names(y)[2])
    xy = rbind(x,y)
    xy[, part := factor(part, levels=c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total"))] 
  # plot
    g1 = 
      ggplot(xy, aes(x = part, y = pred, col = method_CI)) +
        geom_errorbar(aes(ymin = lwr, ymax = upr, col = method_CI), width = 0.1, position = position_dodge(width = 0.25) ) +
        #ggtitle ("Sim based")+
        geom_point(position = position_dodge(width = 0.25)) +
        scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
        scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
        labs(x = NULL, y = "Repeatability [%]")+
        ylim(c(0,100))+
        coord_flip()+
        theme_bw() +
        theme(plot.title = element_text(size=9))
  
    ggsave('Output/VD_within-male_Repeatability.png',g1, width = 10, height =7, units = 'cm')
# repeatability within morph
  # estimate
    lfsim = list()
    lfrpt = list()
    for(i in c('Acrosome','Flagellum','Head','Midpiece','Nucleus','Tail','Total')){
      part_ = i
      # part_ = "Acrosome"
      dd = b[Part == part_]
      m = lmer(Length_µm ~ 1+(1|Morph), dd)
      Rf = R_out(part_)
      lfsim[[i]] = Rf[, method_CI:='arm package']

      R = rpt(Length_µm ~ (1 | Morph), grname = "Morph", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
      RR = data.table(merge(data.frame(name =part_), paste0(round(R$R*100),'%'))) %>% setnames(new = c('part', 'Repeatability'))
      RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
      lfrpt[[i]] =  RR[, method_CI := 'rpt package']
      print(i)
    }
    
    x = do.call(rbind,lfsim)
    names(x)[1] = "part"
    y = do.call(rbind,lfrpt)

    x[, pred:= as.numeric(sub("\\%.*", "", repeatability))]
    x[, lwr:= as.numeric(sub("\\-.*", "", CI))]
    x[, upr:= as.numeric(sub("\\%.*", "",sub(".*-", "", CI)))]
   
    y[, pred:= as.numeric(sub("\\%.*", "", Repeatability))]
    y[, lwr:= as.numeric(sub("\\-.*", "", CI))]
    y[, upr:= as.numeric(sub("\\%.*", "",sub(".*-", "", CI)))]
    #y[nchar(CI) == 5, CI := paste0(0,CI) ]
    #y[, lwr:= as.numeric(substr(CI,1,2))]
    #y[, upr:= as.numeric(substr(CI,4,5))]
    names(y)[2] = tolower( names(y)[2])
    xy = rbind(x,y)
    xy[, part := factor(part, levels=c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total"))] 
  # plot
    g = 
      ggplot(xy, aes(x = part, y = pred, col = method_CI)) +
        geom_errorbar(aes(ymin = lwr, ymax = upr, col = method_CI), width = 0.1, position = position_dodge(width = 0.25) ) +
        #ggtitle ("Sim based")+
        geom_point(position = position_dodge(width = 0.25)) +
        scale_color_viridis(discrete=TRUE, begin=0, end = 0.5, guide = guide_legend(reverse = TRUE)) +
        #scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
        labs(x = NULL, y = "Repeatability [%]")+
        ylim(c(0,100))+
        coord_flip()+
        theme_bw() +
        theme(plot.title = element_text(size=9))
  
    ggsave('Output/VD_within-morph_Repeatability.png',g, width = 10, height =7, units = 'cm')

# MORPH differences BLIND
    ggplot(b, aes(x = Blind, y = Length_avg)) +
    geom_boxplot() + 
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(), col = 'darkgrey', aes(fill =Blind))+
    scale_fill_viridis(discrete=TRUE)+
    facet_wrap(~Part, scales = 'free', nrow = 3)+
    ylab('Length [µm]') +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none")    
# MORPH differences TRUE
  # distributions 
      g1 =  # dummy to extract variables for median calculation
      ggplot(b, aes(x = Morph, y = Length_µm)) +
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 2)+
      geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
      stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
      scale_fill_viridis(discrete=TRUE)+
      facet_wrap(~Part, scales = 'free_y', nrow = 3)+
      ggtitle('Single measurements') +
      ylab('Length [µm]') +
      theme_bw() +
      guides(x =  guide_axis(angle = -45)) +
      theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) #panel.spacing.y = unit(0, "mm")) #, 
 
      g2 =  
      ggplot(a, aes(x = Morph, y = Length_avg)) +
      geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 3)+
      scale_fill_viridis(discrete=TRUE)+
      stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
      facet_wrap(~Part, scales = 'free_y', nrow = 3)+
      ggtitle('Means per male') +
      ylab('Length [µm]') +
      theme_bw() +
      guides(x =  guide_axis(angle = -45)) +
      theme(legend.position = "none",
        plot.title = element_text(size=9)
        ) #panel.spacing.y = unit(0, "mm")) #axis.title.x = element_blank(), axis.text.x = element_blank(), 

      grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), size = "last"))
      #grid.arrange(g1,g2)
      gg1 <- ggplotGrob(g1)
      gg2 <- ggplotGrob(g2) 
      ggsave('Output/VD_Morp_boxplots.png',rbind(gg1,gg2, size = "last"), width = 7*1.5, height =15*1.5, units = 'cm')  

  # CV
   cv_ =  b[, cv(Length_µm), by = list(Part, Morph)]
   cv_[ , Morph123 := as.numeric(Morph)]
   names(cv_) [3]='CV' 
   g=
   ggplot(cv_, aes(x =Morph123 , y = CV, col = Part, fill = Part)) + geom_point() + geom_line() + scale_x_continuous(br eaks =c(1,2,3), labels = c('Independent','Satellite','Faeder')) +
        theme_bw() +
        theme(axis.title.x = element_blank())
   ggsave('Output/VD_CV_overall.png',g, width = 10, height =7, units = 'cm')

   cv_ =  b[, cv(Length_µm), by = list(Bird_ID, Part, Morph)]
   cv_[ , Morph123 := as.numeric(Morph)]
   names(cv_) [4]='CV' 
   g=
   ggplot(cv_, aes(x =Morph , y = CV, col = Part, fill = Part)) + 
   #geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(), col = 'grey',  dotsize = 0.5)+
   geom_boxplot(fill = NA) + 
   theme_bw() +
   theme(axis.title.x = element_blank())
   ggsave('Output/VD_CV_male_alternative.png',g, width = 10, height =7, units = 'cm')
   
   g=
   ggplot(cv_, aes(x =Morph , y = CV)) + 
   geom_boxplot(fill = NA) + 
   geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(), col = 'red',  fill = 'grey', dotsize = 1)+
   facet_wrap(~Part, nrow = 3, scales = "free_y") +
   guides(x =  guide_axis(angle = -45)) +
   theme_bw() +
   theme(axis.title.x = element_blank())
   ggsave('Output/VD_CV_male.png',g, width = 10, height =10, units = 'cm')
    
  # correlations
    # median with 95%
    bs = b[, quantile(Length_µm, prob = c(0.5)), by = list(Part,Morph)]
    names(bs)[3] = 'median'
    bs$lwr = b[, quantile(Length_µm, prob = c(0.025)), by = list(Part,Morph)]$V1
    bs$upr = b[, quantile(Length_µm, prob = c(0.975)), by = list(Part,Morph)]$V1

    # mean +/-sd
    bs = b[, mean(Length_µm), by = list(Part,Morph)]
    names(bs)[3] = 'median'
    bs$lwr = b[, mean(Length_µm)-sd(Length_µm), by = list(Part,Morph)]$V1
    bs$upr = b[, mean(Length_µm)+sd(Length_µm), by = list(Part,Morph)]$V1

    bs_h = bs[Part == 'Head']
    bs_h[,Midpiece_µm:= bs[Part == 'Midpiece',.(median)]]

    bs_m = bs[Part == 'Midpiece']
    bs_m[,Head_µm:= bs[Part == 'Head',.(median)]]
    
    g1 =
    ggplot(h, aes(x = Head_µm, y = Midpiece_µm, col = Morph)) +
      geom_point()+
      geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      geom_segment(data = bs_h, aes(x = lwr, xend = upr, y =Midpiece_µm, yend =Midpiece_µm), col = "black" ) +
      geom_segment(data = bs_m, aes(x = Head_µm, xend = Head_µm, y =lwr, yend =upr), col = "black" )+
      scale_colour_manual(values = colors)+
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      annotate(geom="text", x=32, y=27.5, label="Independent", size = 2, col = colors[1], hjust = 0) +
      annotate(geom="text", x=32, y=27, label="Satellite", size = 2, col = colors[2], hjust = 0) +
      annotate(geom="text", x=32, y=26.5, label="Faeder", size = 2, col = colors[3], hjust = 0) +
      annotate(geom="text", x=32, y=26, label="Mean +/-SD", size = 2, hjust = 0) +
      #scale_colour_viridis(discrete=TRUE)+
      ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 

    bs_h_ = bs[Part == 'Head']
    bs_h_[,Flagellum_µm:= bs[Part == 'Flagellum',.(median)]]

    bs_f = bs[Part == 'Flagellum']
    bs_f[,Head_µm:= bs[Part == 'Head',.(median)]]
      
    g2 =
    ggplot(h, aes(x = Head_µm, y = Flagellum_µm, col = Morph)) +
      geom_point()+
      geom_point(data = bs_h_, aes(x = median, y = Flagellum_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      geom_segment(data = bs_h_, aes(x = lwr, xend = upr, y =Flagellum_µm, yend =Flagellum_µm), col = "black" ) +
      geom_segment(data = bs_f, aes(x = Head_µm, xend = Head_µm, y =lwr, yend =upr), col = "black" )+
      scale_colour_manual(values = colors)+
      #scale_colour_viridis(discrete=TRUE)+
      theme_bw() +
      theme(legend.position = "none",
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )   

    bs_n = bs[Part == 'Nucleus']
    bs_n[,Acrosome_µm:= bs[Part == 'Acrosome',.(median)]]

    bs_a = bs[Part == 'Acrosome']
    bs_a[,Nucleus_µm:= bs[Part == 'Nucleus',.(median)]] 
    
    g3 =
    ggplot(h, aes(x = Nucleus_µm, y = Acrosome_µm, col = Morph)) +
      geom_point()+
      geom_point(data = bs_n, aes(x = median, y = Acrosome_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      geom_segment(data = bs_n, aes(x = lwr, xend = upr, y =Acrosome_µm, yend =Acrosome_µm), col = "black" ) +
      geom_segment(data = bs_a, aes(x = Nucleus_µm, xend = Nucleus_µm, y =lwr, yend =upr), col = "black" )+
      scale_colour_manual(values = colors)+
      #scale_colour_viridis(discrete=TRUE)+
      theme_bw() +
      theme(legend.position = "none",
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )  

    as = a[, quantile(Length_avg, prob = c(0.5)), by = list(Part,Morph)]
    names(as)[3] = 'median'
    as$lwr = a[, quantile(Length_avg, prob = c(0.025)), by = list(Part,Morph)]$V1
    as$upr = a[, quantile(Length_avg, prob = c(0.975)), by = list(Part,Morph)]$V1

    # mean +/-sd
    as = a[, mean(Length_avg), by = list(Part,Morph)]
    names(as)[3] = 'median'
    as$lwr = a[, mean(Length_avg)-sd(Length_avg), by = list(Part,Morph)]$V1
    as$upr = a[, mean(Length_avg)+sd(Length_avg), by = list(Part,Morph)]$V1

    as_h = as[Part == 'Head']
    as_h[,Midpiece_µm:= as[Part == 'Midpiece',.(median)]]

    as_m = as[Part == 'Midpiece']
    as_m[,Head_µm:= as[Part == 'Head',.(median)]]

    g1a =
    ggplot(ha, aes(x = Head_µm, y = Midpiece_µm, col = Morph)) +
      geom_point()+
      geom_point(data = as_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      geom_segment(data = as_h, aes(x = lwr, xend = upr, y =Midpiece_µm, yend =Midpiece_µm), col = "black" ) +
      geom_segment(data = as_m, aes(x = Head_µm, xend = Head_µm, y =lwr, yend =upr), col = "black" )+
      scale_colour_manual(values = colors)+
      #scale_colour_viridis(discrete=TRUE)+
      ggtitle('Means per male') +
      theme_bw() +
      theme(legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        plot.title = element_text(size=9)
        ) 

    as_h_ = as[Part == 'Head']
    as_h_[,Flagellum_µm:= as[Part == 'Flagellum',.(median)]]

    as_f = as[Part == 'Flagellum']
    as_f[,Head_µm:= as[Part == 'Head',.(median)]]

    g2a =
    ggplot(ha, aes(x = Head_µm, y = Flagellum_µm, col = Morph)) +
      geom_point()+
      geom_point(data = as_h_, aes(x = median, y = Flagellum_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      geom_segment(data = as_h_, aes(x = lwr, xend = upr, y =Flagellum_µm, yend =Flagellum_µm), col = "black" ) +
      geom_segment(data = as_f, aes(x = Head_µm, xend = Head_µm, y =lwr, yend =upr), col = "black" )+
      scale_colour_manual(values = colors)+
      #scale_colour_viridis(discrete=TRUE)+
      theme_bw() +
      theme(legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        plot.title = element_text(size=9)
        )   

    as_n = as[Part == 'Nucleus']
    as_n[,Acrosome_µm:= as[Part == 'Acrosome',.(median)]]

    as_a = as[Part == 'Acrosome']
    as_a[,Nucleus_µm:= as[Part == 'Nucleus',.(median)]] 
    g3a =
    ggplot(ha, aes(x = Nucleus_µm, y = Acrosome_µm, col = Morph)) +
      geom_point()+
      geom_point(data = as_n, aes(x = median, y = Acrosome_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      geom_segment(data = as_n, aes(x = lwr, xend = upr, y =Acrosome_µm, yend =Acrosome_µm), col = "black" ) +
      geom_segment(data = as_a, aes(x = Nucleus_µm, xend = Nucleus_µm, y =lwr, yend =upr), col = "black" )+
      scale_colour_manual(values = colors)+
      #scale_colour_viridis(discrete=TRUE)+
      theme_bw() +
      theme(legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        plot.title = element_text(size=9)
        )  

    grid.draw(cbind(rbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "first"),rbind(ggplotGrob(g1a), ggplotGrob(g2a), ggplotGrob(g3a), size = "first"), size = "first"))
      #grid.arrange(g1,g2)
      gg1 <- ggplotGrob(g1)
      gg2 <- ggplotGrob(g2) 
      gg3 <- ggplotGrob(g3) 
      gg1a <- ggplotGrob(g1a)
      gg2a <- ggplotGrob(g2a) 
      gg3a <- ggplotGrob(g3a) 
      ggsave('Output/VD_Morp_corWithMeans.png',cbind(rbind(gg1,gg2, gg3, size = "first"), rbind(gg1a,gg2a, gg3a, size = "first"), size = "first"), width = 7*1.5, height =15, units = 'cm') 
  
  # correlations 3D
    bw[,morph123 :=ifelse(Morph == 'Independent', 1, ifelse(Morph == 'Satellite', 2,3))]             

    # lattice working
      cloud(Tail ~ Head * Midpiece, data=bw, group = Morph, col = colors,
      key = list(text = list(c('Independent', 'Satellite','Faeder'),col=colors,cex=0.6)),
      pch = 16, cex = 1.5, alpha = 0.75, 
      xlab=list(cex=0.7),
      ylab=list(cex=0.7),
      zlab = list(cex=0.7)
      )

      p = cloud(Tail ~ Head * Midpiece, data=bw, group = Morph, col = colors,
              key = list(text = list(c('Independent', 'Satellite','Faeder'),col=colors,cex=0.6)),
              pch = 16, cex = 0.8, alpha = 0.75, 
              xlab=list(cex=0.7),
              ylab=list(cex=0.7),
              zlab = list(cex=0.7) 
              )

      npanel <- c(4, 2)
      rotx <- c(-50, -80)
      rotz <- seq(30, 300, length = npanel[1]+1)  
    
      png('Output/VD_correltations_3D.png',width = 20, height = 10, units = "cm", res = 600)
      update(p[rep(1, prod(npanel))], layout = npanel,
            panel = function(..., screen) {
                panel.cloud(..., screen = list(z = rotz[current.column()],
                                                   x = rotx[current.row()]))
            })
      dev.off()

    # plotly
      plot_ly(data = bw, x=~Head, y=~Midpiece, z=~Tail, color=~Morph, colors = colors)

    # poorly working
      colors_ <- colors[bw$morph123]
      par(las = 1, cex.axis = 0.6, cex.lab = 0.8, cex.main = 0.8)
      s3d=scatterplot3d(bw$Head, bw$Midpiece, bw$Tail, pch = 16, type="h", 
                  color=colors_, grid=TRUE, box=FALSE,
                  #xlab = "",
                  #ylab = "",
                  #zlab = "",
                  #x.ticklabs=c("short","","","","","long"),
                  #y.ticklabs=c("short","","","","long",""),
                  #z.ticklabs=c("slow","","","","","fast"),
                  mar = c(3, 3, 0, 1.5)
                  )     
      text(x = 30, y = 20, "Head", srt = 90,xpd = TRUE, cex = 0.8)
      text(x = 30, y = -0.5, "Midpiece", srt = 0,xpd = TRUE, cex = 0.8)
      text(x = 24, y = 80, "Tail", srt = 0, cex = 0.8)
     
      text(x = 7.5, y = 1, "Tail", srt = 0, cex = 0.8)
      text(x = 2.5, y = -0.5, "Mipiece", srt = 0,xpd = TRUE, cex = 0.8)
      text(x = -0.5, y = 2.5, "Velocity", srt = 90,xpd = TRUE, cex = 0.8)

      legend("bottom", legend = levels(factor(bw$Morph,levels = c('Independent','Satellite','Faeder'))),
          col =  colors, pch = 16,xpd = TRUE, horiz = TRUE,inset = -0.03, bty = "n", cex = 0.7)
    

  # models differences in Morphs  
    # averages     
      # prepare estimates and predictions
        l = list()
        lp =list()
        for(i in unique(a$Part)){
          #i ='Nucleus'
          m = lm(scale(Length_avg) ~ Morph, a[Part == i])
          #summary(m)
          #plot(allEffects(m))
          bsim = sim(m, n.sim=nsim) 
          v = apply(bsim@coef, 2, quantile, prob=c(0.5))
          ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
          l[[i]]=data.frame(response=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

          # get predictions
          m = lm(Length_avg ~ Morph, a[Part == i])
          bsim = sim(m, n.sim=nsim) 
          v = apply(bsim@coef, 2, quantile, prob=c(0.5))
          newD=data.frame(Morph = unique(a$Morph)) # values to predict for
          X <- model.matrix(~ Morph,data=newD) # exactly the model which was used has to be specified here 
          newD$Length_avg <-(X%*%v) 
          predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
          for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                          predmatrix[predmatrix < 0] <- 0
                          newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                          newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                          #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
          newD$Part=i
          lp[[i]] = data.table(newD)

          print(i)     
          }          
        
        ll = data.table(do.call(rbind,l) ) 
        ll[, response := factor(response, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
        ll[, effect := factor(effect, levels=rev(c("(Intercept)", "MorphSatellite", "MorphFaeder")))]
        
        llp = data.table(do.call(rbind,lp) ) 
        llp[, Part := factor(Part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
        llp[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
      # plot effect sizes     
        g1 = 
        ggplot(ll, aes(y = effect, x = estimate, col = response)) +
          geom_vline(xintercept = 0, col = "grey", lty =3)+
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = response), width = 0.1, position = position_dodge(width = 0.25) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.25)) +
          scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          labs(y = NULL ,x = "Standardized effect size")+
          #ylim(c(0,100))+
          #coord_flip()+
          theme_bw() +
          theme(plot.title = element_text(size=9))
  
        ggsave('Output/VD_lm_effectSizes.png',g1, width = 10, height =14, units = 'cm')    
      # prepare plot with distributions with predictions
        g1 = 
          ggplot(a, aes(x = Morph, y = Length_avg)) +
          geom_dotplot(binaxis = 'y', stackdir = 'center',
                        position = position_dodge(), col = 'grey', fill = 'lightgrey', dotsize = 1.5)+
          geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2) + 
          geom_errorbar(data = llp, aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(width = 0.25), col = 'red' ) +
          geom_point(data = llp, aes(x = Morph, y =Length_avg), position = position_dodge(width = 0.25), col = 'red') +
          #stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
          scale_color_viridis(discrete=TRUE)+
          facet_wrap(~Part, scales = 'free_y', nrow = 3)+
          ggtitle('Means per male') +
          ylab('Length [µm]') +
          theme_bw() +
          guides(x =  guide_axis(angle = -45)) +
          theme(legend.position = "none",
            axis.title.x = element_blank(), 
            #axis.text.x = element_blank(),
            plot.title = element_text(size=9)
            ) #panel.spacing.y = unit(0, "mm")) #,    
 
    # single values 
      # prepare estimates and predictions
        ls = list()
        lps = list()
        for(i in unique(a$Part)){
          #i ='Nucleus'
          m = lmer(scale(Length_µm) ~ Morph + (1|Bird_ID), b[Part == i])
          #summary(m)
          #plot(allEffects(m))
          bsim = sim(m, n.sim=5000) 
          v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
          ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
          ls[[i]]=data.frame(response=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])
         
          # get predictions
          m = lmer(Length_µm ~ Morph + (1|Bird_ID), b[Part == i])
          bsim = sim(m, n.sim=nsim) 
          v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
          newD=data.frame(Morph = unique(a$Morph)) # values to predict for
          X <- model.matrix(~ Morph,data=newD) # exactly the model which was used has to be specified here 
          newD$Length_µm <-(X%*%v) 
          predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
          for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                          predmatrix[predmatrix < 0] <- 0
                          newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                          newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                          #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
          newD$Part=i
          lps[[i]] = data.table(newD)

          print(i)         
          }  
        
        lls = data.table(do.call(rbind,ls) ) 
        lls[, response := factor(response, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
        lls[, effect := factor(effect, levels=rev(c("(Intercept)", "MorphSatellite", "MorphFaeder")))]
  
        llps = data.table(do.call(rbind,lps) ) 
        llps[, Part := factor(Part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
        llps[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
      # plot effect sizes     
        g = 
        ggplot(lls, aes(y = effect, x = estimate, col = response)) +
          geom_vline(xintercept = 0, col = "grey", lty =3)+
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = response), width = 0.1, position = position_dodge(width = 0.25) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.25)) +
          scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          labs(y = NULL ,x = "Standardized effect size")+
          #ylim(c(0,100))+
          #coord_flip()+
          theme_bw() +
          theme(plot.title = element_text(size=9))
    
        ggsave('Output/VD_lmer_effectSizes.png',g, width = 10, height =14, units = 'cm')
      # prepare plot for distributions with predictions
        g2 = 
          ggplot(b, aes(x = Morph, y = Length_µm)) +
          geom_dotplot(binaxis = 'y', stackdir = 'center',
                        position = position_dodge(), col = 'grey', fill = 'lightgrey', dotsize = 1)+
          geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2) + 
          geom_errorbar(data = llps, aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(width = 0.25), col = 'red' ) +
          geom_point(data = llps, aes(x = Morph, y =Length_µm), position = position_dodge(width = 0.25), col = 'red') +
          #stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
          scale_color_viridis(discrete=TRUE)+
          facet_wrap(~Part, scales = 'free_y', nrow = 3)+
          ggtitle('Single measurements') +
          ylab('Length [µm]') +
          theme_bw() +
          guides(x =  guide_axis(angle = -45)) +
          theme(legend.position = "none",
            axis.title.x = element_blank(), 
            axis.text.x = element_blank(),
            plot.title = element_text(size=9)
            ) #panel.spacing.y = unit(0, "mm")) #,     

    # plot predictions and boxplots in one plot
      grid.draw(rbind(ggplotGrob(g2), ggplotGrob(g1), size = "last"))
      #grid.arrange(g1,g2)
      gg1 <- ggplotGrob(g1)
      gg2 <- ggplotGrob(g2) 
      ggsave('Output/VD_predictions+boxPlots_boxCol.png',rbind(gg2,gg1, size = "last"), width = 7*1.5, height =15*1.5, units = 'cm')    

    # plot effects sizes in one plot
      ll[,model := paste0(response, ', average sperm')]
      lls[,model := paste0(response, ', single sperm')]
      ll[,data_ := 'average sperm']
      lls[,data_ := 'single sperm']
      ba = rbind(ll,lls)
      #ba[, data_ := factor(data_, levels=c("single sperm", "average sperm"))] 
      ba = ba[!effect=='(Intercept)']
      ba[effect=='MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      ba[effect=='MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      ba[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)"))] 
      g = 
        ggplot(ba, aes(y = effect, x = estimate, col = response, shape = data_)) +
          geom_vline(xintercept = 0, col = "grey", lty =3)+
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = response), width = 0, position = position_dodge(width =0.5) ) +
          #ggtitle ("Relative to Independent")+
          geom_point(position = position_dodge(width = 0.5)) +
          scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          #scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          scale_shape(guide = guide_legend(reverse = TRUE)) +
          scale_x_continuous(limits = c(-3, 3), expand = c(0, 0), breaks = seq(-3,3, by = 1), labels = seq(-3,3, by = 1)) +
          labs(y = NULL ,x = "Standardized effect size", fill = 'Response', col = 'Response', shape = 'Data')+
          #coord_flip()+
          theme_bw() +
          theme( #legend.position ="right",
                legend.title=element_text(size=7), 
                legend.text=element_text(size=6),
                ##legend.spacing.y = unit(0.1, 'cm'), 
                legend.key.height= unit(0.5,"line"),
                #plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r =0.5, unit =  "pt"),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = ax_lines, size = 0.25),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x= element_line( colour = ax_lines, size = 0.25),
                #axis.text.x = element_text()
                axis.ticks.length = unit(1, "pt"),
                axis.text.x=element_text(, size = 6),
                axis.text.y=element_text(colour="black", size = 7),
                axis.title=element_text(size=7)
                )

      ggsave('Output/VD_ALL_effectSizes.png',g, width = 10*1.5, height =5*1.5, units = 'cm')

# stained vs unstained
  bs = b[Bird_ID =='329']
  bs[,Part := as.character(Part)]
  bs$Bird_ID = NULL
  #ds$stained = 'yes'

  u = data.table(read_excel('Data/VD measurements feb_2021.xlsx', sheet = '329 unst'))#, range = "A1:G161"))
    u[,Length_µm :=as.character(Length_µm) ]
    stri_sub(u$Length_µm,nchar(u$Length_µm)-2, nchar(u$Length_µm)-3) = "."
    u[,Length_µm :=as.numeric(Length_µm) ]
    u[Part == 'Mid piece', Part := 'Midpiece']
    u[Part == 'Flagellum', Part := 'Tail']
    u[Part == 'Head', Part := 'Nucleus']

    u1 = u[,.(Bird_ID, Sperm_ID, Part, Length_µm)]
    
    u3 = u[Part%in%c('Midpiece', 'Tail'), sum(Length_µm), by = list(Bird_ID, Sperm_ID)]
    names(u3)[3] = 'Length_µm'
    u3[,Part := 'Flagellum']

    uu = merge(u1,u3, all = TRUE)
    uu = merge(uu, vv, by.x = 'Bird_ID', by.y = 'ID')
    uu[Morph == 'Res', Morph := 'Independent']
    uu[Morph == 'Sat', Morph := 'Satellite']
    uu[Morph == 'Faed', Morph := 'Faeder']
    #uu$stained = 'no'
    names(uu)[4] = 'Length_µm_unstained'
    #u[, Part := factor(Part, levels=c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel"))] 
    uu[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 
    uu$Bird_ID = NULL
  su = merge(bs,uu)


  ggplot(su[Part!='Acrosome'], aes(x = Length_µm, y = Length_µm_unstained)) +
      stat_smooth(method = 'lm') +
      geom_point() + 
      stat_cor(method="pearson",size = 2) +
      geom_abline(b = 1, col = 'red', lty = 3) + 
      facet_wrap(~Part, scales = 'free', nrow = 3)+
      #ggtitle('Single measurements') +
      #ylab('Length [µm]') +
      theme_bw() +
      guides(x =  guide_axis(angle = -45)) +
      theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) #panel.spacing.y = unit(0, "mm")) #, 
  ggsave('Output/VD_stainedVSunstained.png')

# End