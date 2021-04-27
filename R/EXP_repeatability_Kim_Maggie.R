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

# Prepare data
  d = data.table(read_excel('Data/VD measurements feb_2021.xlsx', sheet = 1))#, range = "A1:G161"))
  d[,Length_µm :=as.character(Length_µm) ]
  stri_sub(d$Length_µm,nchar(d$Length_µm)-2, nchar(d$Length_µm)-3) = "."
  d[,Length_µm :=as.numeric(Length_µm) ]
  d[Part == 'Mid piece', Part := 'Midpiece']
  d[Part == 'Flagellum', Part := 'Tail']
  d[Part == 'Head', Part := 'Nucleus']

  d1 = d[,.(Bird_ID, Sample_ID,Photo_ID, Sperm_ID, Line_ID, Part, Length_µm)]
  
  d2 = d[, sum(Length_µm), by = list(Bird_ID, Sample_ID, Photo_ID, Sperm_ID)]
  names(d2)[length(names(d2))] = 'Length_µm'
  d2[,Part := 'Total']
  d2[,Line_ID := 100]

  d3 = d[Part%in%c('Midpiece', 'Tail'), sum(Length_µm), by = list(Bird_ID, Sample_ID,Photo_ID, Sperm_ID)]
  names(d3)[length(names(d3))] = 'Length_µm'
  d3[,Part := 'Flagellum']
  d3[,Line_ID := 100]

  d4 = d[Part%in%c('Acrosome', 'Nucleus'), sum(Length_µm), by = list(Bird_ID, Sample_ID, Photo_ID,Sperm_ID)]
  names(d4)[length(names(d4))] = 'Length_µm'
  d4[,Part := 'Head']
  d4[,Line_ID := 100]

  b = merge(d1,d2, all = TRUE)
  b = merge(b,d3, all = TRUE)
  b = merge(b,d4, all = TRUE)

  r = data.table(read_excel('Data/ruff_sperm_2021_test Maggie.xlsx', sheet = 1))#, range = "A1:G161"))
  r[,Length_µm :=as.character(Length_µm) ]
  stri_sub(r$Length_µm,nchar(r$Length_µm)-2, nchar(r$Length_µm)-3) = "."
  r[,Length_µm :=as.numeric(Length_µm) ]
  r[Part == 'Mid piece', Part := 'Midpiece']
  r[Part == 'Flagellum', Part := 'Tail']
  r[Part == 'Head', Part := 'Nucleus']

  d1 = r[,.(Bird_ID, Sample_ID, Photo_ID, Sperm_ID, Line_ID, Part, Length_µm)]
  
  d2 = r[, sum(Length_µm), by = list(Bird_ID, Sample_ID,Photo_ID, Sperm_ID)]
   names(d2)[length(names(d2))] = 'Length_µm'
  d2[,Part := 'Total']
  d2[,Line_ID := 100]

  d3 = r[Part%in%c('Midpiece', 'Tail'), sum(Length_µm), by = list(Bird_ID, Sample_ID,Photo_ID, Sperm_ID)]
  names(d3)[length(names(d3))] = 'Length_µm'
  d3[,Part := 'Flagellum']
  d3[,Line_ID := 100]

  d4 = r[Part%in%c('Acrosome', 'Nucleus'), sum(Length_µm), by = list(Bird_ID, Sample_ID, Photo_ID,Sperm_ID)]
  names(d4)[length(names(d4))] = 'Length_µm'
  d4[,Part := 'Head']
  d4[,Line_ID := 100]

  rr = merge(d1,d2, all = TRUE)
  rr = merge(rr,d3, all = TRUE)
  rr = merge(rr,d4, all = TRUE)

  dr = b[paste(Bird_ID, Photo_ID, Sperm_ID, Line_ID, Part) %in% paste(rr$Bird_ID, rr$Photo_ID, rr$Sperm_ID, rr$Line_ID, rr$Part)]
  dr[, who := 'Kim']
  #r[!paste(Photo_ID, Sperm_ID, Line_ID, Part)%in%paste(d$Photo_ID, d$Sperm_ID, d$Line_ID, d$Part)]
  r1 = rr[paste(Bird_ID,Photo_ID, Sperm_ID, Line_ID, Part)%in%paste(b$Bird_ID,b$Photo_ID, b$Sperm_ID, b$Line_ID, b$Part)]
  r1[, who := 'Maggie']

  rr_ = rbind(dr,r1)
  rr_[,ID := paste(Photo_ID, Sperm_ID, Line_ID, Part)]

  dr[ , Kim := Length_µm]
  dr$Maggie = r1$Length_µm[match(paste(dr$Bird_ID, dr$Photo_ID, dr$Sperm_ID, dr$Line_ID, dr$Part), paste(r1$Bird_ID, r1$Photo_ID, r1$Sperm_ID, r1$Line_ID, r1$Part))]
 
# plot correlations
  g1 = ggplot(dr[!Part%in%c('Head','Flagellum','Total')], aes(x = Kim, y = Maggie)) +
  facet_wrap(~Part, scales = "free") +  
  stat_smooth(method = 'lm')+geom_point()+
  stat_cor(method="pearson",size = 2) +
  geom_abline(b = 1, col = 'red', lty = 3) + 
  #xlim(c(485, 565)) +  ylim(c(485, 565)) + 
  #ggtitle('All')+
  theme_bw()
  ggsave('Output/VD_between_Kim_Maggie_correlation_BASIC.png',g1, width = 10, height =10, units = 'cm')

  g2 = ggplot(dr[Part%in%c('Head','Flagellum','Total')], aes(x = Kim, y = Maggie)) +
  facet_wrap(~Part, scales = "free", ncol =2) +  
  stat_smooth(method = 'lm')+geom_point()+
  stat_cor(method="pearson",size = 2) +
  geom_abline(b = 1, col = 'red', lty = 3) + 
  #xlim(c(485, 565)) +  ylim(c(485, 565)) + 
  #ggtitle('All')+
  theme_bw()
  ggsave('Output/VD_between_Kim_Maggie_correlation_Composite.png',g2, width = 10, height =10, units = 'cm')
# plot correlations for 370 only
  g1 = ggplot(dr[!Part%in%c('Head','Flagellum','Total') & Bird_ID == 370], aes(x = Kim, y = Maggie)) +
  facet_wrap(~Part, scales = "free") +  
  stat_smooth(method = 'lm')+geom_point()+
  stat_cor(method="pearson",size = 2) +
  geom_abline(b = 1, col = 'red', lty = 3) + 
  #xlim(c(485, 565)) +  ylim(c(485, 565)) + 
  #ggtitle('All')+
  theme_bw()
  ggsave('Output/VD_between_Kim_Maggie_correlation_BASIC_370.png',g1, width = 10, height =10, units = 'cm')

  g2 = ggplot(dr[Part%in%c('Head','Flagellum','Total') & Bird_ID == 370], aes(x = Kim, y = Maggie)) +
  facet_wrap(~Part, scales = "free", ncol =2) +  
  stat_smooth(method = 'lm')+geom_point()+
  stat_cor(method="pearson",size = 2) +
  geom_abline(b = 1, col = 'red', lty = 3) + 
  #xlim(c(485, 565)) +  ylim(c(485, 565)) + 
  #ggtitle('All')+
  theme_bw()
  ggsave('Output/VD_between_Kim_Maggie_correlation_Composite_370.png',g2, width = 10, height =10, units = 'cm')

# plot repeatability
  # estimate BASIC 
    lfsim = list()
    lfrpt = list()
    for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
      part_ = i
      # part_ = "Acrosome"
      dd = rr_[Part == part_]
      m = lmer(Length_µm ~ 1+(1|ID), dd)
      Rf = R_out(part_)
      lfsim[[i]] = Rf[, method_CI:='arm package']

      R = rpt(Length_µm ~ (1 | ID), grname = "ID", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
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
  # plot BASIC
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
  
    ggsave('Output/VD_between_Kim_Maggie_Repeatability_BASIC.png',g1, width = 10, height =7, units = 'cm')
  # estimate COMPOSITE
    lfsim = list()
    lfrpt = list()
    for(i in c("Head", "Flagellum","Total")){
      part_ = i
      # part_ = "Acrosome"
      dd = rr_[Part == part_]
      m = lmer(Length_µm ~ 1+(1|ID), dd)
      Rf = R_out(part_)
      lfsim[[i]] = Rf[, method_CI:='arm package']

      R = rpt(Length_µm ~ (1 | ID), grname = "ID", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
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
  # plot COMPOSITE
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
  
    ggsave('Output/VD_between_Kim_Maggie_Repeatability_COMPOSITE.png',g1, width = 10, height =7, units = 'cm')

# plot repeatability for 370 only
  # estimate BASIC 
    lfsim = list()
    lfrpt = list()
    for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
      part_ = i
      # part_ = "Acrosome"
      dd = rr_[Part == part_ & Bird_ID == 370]
      m = lmer(Length_µm ~ 1+(1|ID), dd)
      Rf = R_out(part_)
      lfsim[[i]] = Rf[, method_CI:='arm package']

      R = rpt(Length_µm ~ (1 | ID), grname = "ID", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
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
  # plot BASIC
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
  
    ggsave('Output/VD_between_Kim_Maggie_Repeatability_BASIC_370.png',g1, width = 10, height =7, units = 'cm')
  # estimate COMPOSITE
    lfsim = list()
    lfrpt = list()
    for(i in c("Head", "Flagellum","Total")){
      part_ = i
      # part_ = "Acrosome"
      dd = rr_[Part == part_ & Bird_ID == 370]
      m = lmer(Length_µm ~ 1+(1|ID), dd)
      Rf = R_out(part_)
      lfsim[[i]] = Rf[, method_CI:='arm package']

      R = rpt(Length_µm ~ (1 | ID), grname = "ID", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
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
  # plot COMPOSITE
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
  
    ggsave('Output/VD_between_Kim_Maggie_Repeatability_COMPOSITE_370.png',g1, width = 10, height =7, units = 'cm')    

# END