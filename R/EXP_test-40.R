# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(stringi)
  require(arm)
  require(effects)
  require(RColorBrewer)
  require("PerformanceAnalytics")
  require('foreach')
  library(ggpubr)
  
  require(rptR) 
  require(magrittr)

  round_ = 3 # number of decimal places to round model coefficients
  nsim = 5000 # number of simulations to extract estimates and 95%CrI
  ax_lines = "grey60" # defines color of the axis lines
  colors <- c("#999999", "#E69F00", "#56B4E9") #viridis(3)
  # functions
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
 
# new protocol - Kim Martin Maggie
  # DATA
    f = data.table(
      f = c(list.files(path = 'Data/', pattern = 'test40', recursive = TRUE, full.names = TRUE)),
      f2 = c(list.files(path = 'Data/', pattern = 'test40', recursive = TRUE, full.names = FALSE))
      )   

    ff = f[substr(f2,8,9) == 'KT']
    d1 = foreach(j = 1:nrow(ff), .combine = rbind) %do% {
      #k = ff[1,]
      k = ff[j,]
      x = fread(file = k[,f])
      x[, who := substr(k[,f2], 8, 9)]
      #x[, datetime_ := substr(j, 27, 45)]
      x[, id := substr(x$'File Name', 1, 2)]
      x[, pk := 1:nrow(x)]
      #x[, manip := substr(x$'File Name', 4, nchar(x$'File Name')-7-(nchar(pk)))]
      
      #table(x$id)
      xx = foreach(i = unique(x$id), .combine = rbind) %do% {
        #i = '02'
        xi = x[id == i]
        if(i == '29'){xi = xi[!Pixels<20]}
        if(nrow(xi) == 4){
          xi[Pixels == max(Pixels), part := 'Tail']
          xi[Pixels != max(Pixels), part := c('Acrosome','Nucleus','Midpiece')]
        }

        if(nrow(xi) == 3){
          xi[Pixels == max(Pixels), part := 'Tail']
          xi[Pixels != max(Pixels), part := c('Nucleus','Midpiece')]
        }
        
        if(nrow(xi) == 1){
          xi[, part := 'Acrosome']
        }
        return(xi)
        }    
      return(xx)
      }
    
    ff = f[substr(f2,8,9) == 'MB']
    d2 = foreach(j = 1:nrow(ff), .combine = rbind) %do% {
      #k = ff[1,]
      k = ff[j,]
      x = fread(file = k[,f])
      x[, who := substr(k[,f2], 8, 9)]
      #x[, datetime_ := substr(j, 27, 45)]
      x[, id := substr(x$'File Name', 1, 2)]
      x[, pk := 1:nrow(x)]
      #x[, manip := substr(x$'File Name', 4, nchar(x$'File Name')-7-(nchar(pk)))]
      
      #table(x$id)
      xx = foreach(i = unique(x$id), .combine = rbind) %do% {
        #i = '02'
        xi = x[id == i]
        xi[Pixels == max(Pixels), part := 'Tail']
        xi[Pixels != max(Pixels), part := c('Acrosome','Nucleus','Midpiece')]
        return(xi)
        }    
      return(xx)
      }
      
    ff = f[substr(f2,8,9) == 'MC']
    d3 = foreach(j = 1:nrow(ff), .combine = rbind) %do% {
      #k = ff[1,]
      k = ff[j,]
      x = fread(file = k[,f])
      x[, who := substr(k[,f2], 8, 9)]
      #x[, datetime_ := substr(j, 27, 45)]
      x[, id := substr(x$'File Name', 1, 2)]
      x[, pk := 1:nrow(x)]
      #x[, manip := substr(x$'File Name', 4, nchar(x$'File Name')-7-(nchar(pk)))]
      
      #table(x$id)
      xx = foreach(i = unique(x$id), .combine = rbind) %do% {
        #i = '02'
        xi = x[id == i]
        xi[Pixels == max(Pixels), part := 'Tail']
        xi[Pixels != max(Pixels), part := c('Acrosome','Nucleus','Midpiece')]
        return(xi)
        }    
      return(xx)
      }
        
    d = rbind(d1,d2,d3)  
    d[, id_part := paste(id, part)]
    b = d[id_part %in% d[duplicated(id_part), id_part]] # only sperm ids and parts measured twice

     # composite parts
         b[, id_who := paste(id, who)]
         bt = data.table(ddply(b,.(who, id, id_who), summarise, part = 'Total', Pixels = sum(Pixels)))
         bt = bt[Pixels>1000]
         bh = data.table(ddply(b[part %in% c('Acrosome','Nucleus')],.(who, id, id_who), summarise, part = 'Head', Pixels = sum(Pixels)))
         bf = data.table(ddply(b[part %in% c('Midpiece','Tail')],.(who, id, id_who), summarise, part = 'Flagellum', Pixels = sum(Pixels)))
         bf = bf[Pixels>1000]
         bc = rbind(bh,bf,bt)

     # wide format
         bw = reshape(b[,.(who,id,part,Pixels)], idvar = c('id','part'), timevar = 'who', direction = "wide")

         bcw = reshape(bc[,.(who,id,part,Pixels)], idvar = c('id','part'), timevar = 'who', direction = "wide")    
  # ANAL
    # plot correlations KT-MB
        g1 = ggplot(bw, aes(x = Pixels.KT, y = Pixels.MB)) +
        facet_wrap(~part, scales = "free") +  
        stat_smooth(method = 'lm')+geom_point()+
        stat_cor(method="pearson",size = 2) +
        geom_abline(b = 1, col = 'red', lty = 3) + 
        #xlim(c(485, 565)) +  ylim(c(485, 565)) + 
        #ggtitle('All')+
        theme_bw()
        ggsave('Output//test-40_cor-KT-MB.png',g1, width = 10, height =10, units = 'cm')
        
        g2 = ggplot(bcw, aes(x = Pixels.KT, y = Pixels.MB)) +
        facet_wrap(~part, scales = "free", nrow = 2) +  
        stat_smooth(method = 'lm')+geom_point()+
        stat_cor(method="pearson",size = 2) +
        geom_abline(b = 1, col = 'red', lty = 3) + 
        #xlim(c(485, 565)) +  ylim(c(485, 565)) + 
        #ggtitle('All')+
        theme_bw()
        ggsave('Output//test-40_cor-KT-MB_composite.png',g2, width = 10, height =10, units = 'cm')
    # plot correlations KT-MC
        g5 = ggplot(bw, aes(x = Pixels.KT, y = Pixels.MC)) +
        facet_wrap(~part, scales = "free") +  
        stat_smooth(method = 'lm')+geom_point()+
        stat_cor(method="pearson",size = 2) +
        geom_abline(b = 1, col = 'red', lty = 3) + 
        #xlim(c(485, 565)) +  ylim(c(485, 565)) + 
        #ggtitle('All')+
        theme_bw()
        ggsave('Output//test-40_cor-KT-MC.png',g5, width = 10, height =10, units = 'cm')
        
        g6 = ggplot(bcw, aes(x = Pixels.KT, y = Pixels.MC)) +
        facet_wrap(~part, scales = "free", nrow = 2) +  
        stat_smooth(method = 'lm')+geom_point()+
        stat_cor(method="pearson",size = 2) +
        geom_abline(b = 1, col = 'red', lty = 3) + 
        #xlim(c(485, 565)) +  ylim(c(485, 565)) + 
        #ggtitle('All')+
        theme_bw()
        ggsave('Output//test-40_cor-KT-MC_composite.png',g6, width = 10, height =10, units = 'cm')
    # plot correlations MB-MC
        g3 = ggplot(bw, aes(x = Pixels.MB, y = Pixels.MC)) +
        facet_wrap(~part, scales = "free") +  
        stat_smooth(method = 'lm')+geom_point()+
        stat_cor(method="pearson",size = 2) +
        geom_abline(b = 1, col = 'red', lty = 3) + 
        #xlim(c(485, 565)) +  ylim(c(485, 565)) + 
        #ggtitle('All')+
        theme_bw()
        ggsave('Output//test-40_cor-MB-MC.png',g3, width = 10, height =10, units = 'cm')
        
        g4 = ggplot(bcw, aes(x = Pixels.MB, y = Pixels.MC)) +
        facet_wrap(~part, scales = "free", nrow = 2) +  
        stat_smooth(method = 'lm')+geom_point()+
        stat_cor(method="pearson",size = 2) +
        geom_abline(b = 1, col = 'red', lty = 3) + 
        #xlim(c(485, 565)) +  ylim(c(485, 565)) + 
        #ggtitle('All')+
        theme_bw()
        ggsave('Output//test-40_cor-MB-MC_composite.png',g4, width = 10, height =10, units = 'cm')
    # plot correlations all together (run previous first)
      grid.draw(rbind(cbind(ggplotGrob(g1), ggplotGrob(g2)), cbind(ggplotGrob(g3), ggplotGrob(g4)),cbind(ggplotGrob(g5), ggplotGrob(g6))))

      ggsave('Output/test-40_cor-ALL.png',rbind(cbind(ggplotGrob(g1), ggplotGrob(g2)), cbind(ggplotGrob(g3), ggplotGrob(g4)),cbind(ggplotGrob(g5), ggplotGrob(g6))), width = 6.5*2.5, height =8*2.5, units = 'cm')  
    # check differences
        bwa = bw[part == 'Acrosome']
        bwa[, diff_KT_MB := Pixels.KT - Pixels.MB]
        bwa[, diff_KT_MC := Pixels.KT - Pixels.MC]
        bwa[, diff_MB_MC := Pixels.MB - Pixels.MC]
        bwa[order(diff_KT_MB)]
        bwa[order(diff_KT_MC)]
        bwa[order(diff_MB_MC)]

    # plot repeatability all
          # estimate BASIC 
              lfsim = list()
              lfrpt = list()
              for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
                part_ = i
                # part_ = "Acrosome"
                dd = b[part == part_]
                m = lmer(Pixels ~ 1+(1|id), dd)
                Rf = R_out(part_)
                lfsim[[i]] = Rf[, method_CI:='arm package']

                R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
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
              xy[, part := factor(part, levels=c("Acrosome", "Nucleus", "Midpiece","Tail"))] 
              xy[nchar(CI) == 7 & upr == 10, upr := 100]
          # plot BASIC
              g1 = 
                ggplot(xy, aes(x = part, y = pred, col = method_CI)) +
                  geom_errorbar(aes(ymin = lwr, ymax = upr, col = method_CI), width = 0.1, position = position_dodge(width = 0.25) ) +
                  #ggtitle ("Sim based")+
                  geom_point(position = position_dodge(width = 0.25)) +
                  scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
                  scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
                  geom_hline(yintercept = 99, col = 'red')+
                  labs(x = NULL, y = "Repeatability [%]")+
                  ylim(c(60,100))+
                  scale_y_continuous(breaks = sort(c(seq(60, 100, length.out=6), 99)))+
                  coord_flip()+
                  theme_bw() +
                  theme(plot.title = element_text(size=9))
            
              ggsave('Output/test-40_Repeatability-MB-KM-MC_BASIC.png',g1, width = 10, height =7, units = 'cm')
          # estimate COMPOSITE
              lfsim = list()
              lfrpt = list()
              for(i in c("Head", "Flagellum","Total")){
                part_ = i
                # part_ = "Flagellum"
                dd = bc[part == part_]
                m = lmer(Pixels ~ 1+(1|id), dd)
                Rf = R_out(part_)
                lfsim[[i]] = Rf[, method_CI:='arm package']

                R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
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
              xy2 = rbind(x,y)
              xy2[, part := factor(part, levels=c("Head","Flagellum","Total"))] 
              xy2[nchar(CI) == 7 & upr == 10, upr := 100]
          # plot COMPOSITE
              g1 = 
                ggplot(xy2, aes(x = part, y = pred, col = method_CI)) +
                  geom_errorbar(aes(ymin = lwr, ymax = upr, col = method_CI), width = 0.1, position = position_dodge(width = 0.25) ) +
                  #ggtitle ("Sim based")+
                  geom_point(position = position_dodge(width = 0.25)) +
                  scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
                  scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
                  labs(x = NULL, y = "Repeatability [%]")+
                  ylim(c(min(xy2$lwr),max(xy2$upr)))+
                  coord_flip()+
                  theme_bw() +
                  theme(plot.title = element_text(size=9))
            
              ggsave('Output/test-40_Repeatability-MB-KM-MC_COMPOSITE.png',g1, width = 10, height =7, units = 'cm')
    
    # plot repeatability MB, KT
          # estimate BASIC 
              lfsim = list()
              lfrpt = list()
              for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
                part_ = i
                # part_ = "Acrosome"
                dd = b[part == part_ & who %in% c('KT','MB')]
                m = lmer(Pixels ~ 1+(1|id), dd)
                Rf = R_out(part_)
                lfsim[[i]] = Rf[, method_CI:='arm package']

                R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
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
              xy[, part := factor(part, levels=c("Acrosome", "Nucleus", "Midpiece","Tail"))] 
              xy[nchar(CI) == 7 & upr == 10, upr := 100]
          # plot BASIC
              g1 = 
                ggplot(xy, aes(x = part, y = pred, col = method_CI)) +
                  geom_errorbar(aes(ymin = lwr, ymax = upr, col = method_CI), width = 0.1, position = position_dodge(width = 0.25) ) +
                  #ggtitle ("Sim based")+
                  geom_point(position = position_dodge(width = 0.25)) +
                  scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
                  scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
                  geom_hline(yintercept = 99, col = 'red')+
                  labs(x = NULL, y = "Repeatability [%]")+
                  ylim(c(60,100))+
                  scale_y_continuous(breaks = sort(c(seq(60, 100, length.out=6), 99)))+
                  coord_flip()+
                  theme_bw() +
                  theme(plot.title = element_text(size=9))
            
              ggsave('Output/test-40_Repeatability-MB-KM_BASIC.png',g1, width = 10, height =7, units = 'cm')
          # estimate COMPOSITE
              lfsim = list()
              lfrpt = list()
              for(i in c("Head", "Flagellum","Total")){
                part_ = i
                # part_ = "Flagellum"
                dd = bc[part == part_ * who %in% c('KT','MB')]
                m = lmer(Pixels ~ 1+(1|id), dd)
                Rf = R_out(part_)
                lfsim[[i]] = Rf[, method_CI:='arm package']

                R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
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
              xy2 = rbind(x,y)
              xy2[, part := factor(part, levels=c("Head","Flagellum","Total"))] 
              xy2[nchar(CI) == 7 & upr == 10, upr := 100]
          # plot COMPOSITE
              g1 = 
                ggplot(xy2, aes(x = part, y = pred, col = method_CI)) +
                  geom_errorbar(aes(ymin = lwr, ymax = upr, col = method_CI), width = 0.1, position = position_dodge(width = 0.25) ) +
                  #ggtitle ("Sim based")+
                  geom_point(position = position_dodge(width = 0.25)) +
                  scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
                  scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
                  labs(x = NULL, y = "Repeatability [%]")+
                  ylim(c(min(xy2$lwr),max(xy2$upr)))+
                  coord_flip()+
                  theme_bw() +
                  theme(plot.title = element_text(size=9))
            
              ggsave('Output/test-40_Repeatability-MB-KM_COMPOSITE.png',g1, width = 10, height =7, units = 'cm')
    # plot repeatability MB, MC
          # estimate BASIC 
              lfsim = list()
              lfrpt = list()
              for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
                part_ = i
                # part_ = "Acrosome"
                dd = b[part == part_ & who %in% c('MC','MB')]
                m = lmer(Pixels ~ 1+(1|id), dd)
                Rf = R_out(part_)
                lfsim[[i]] = Rf[, method_CI:='arm package']

                R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
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
              xy[, part := factor(part, levels=c("Acrosome", "Nucleus", "Midpiece","Tail"))] 
              xy[nchar(CI) == 7 & upr == 10, upr := 100]
          # plot BASIC
              g1 = 
                ggplot(xy, aes(x = part, y = pred, col = method_CI)) +
                  geom_errorbar(aes(ymin = lwr, ymax = upr, col = method_CI), width = 0.1, position = position_dodge(width = 0.25) ) +
                  #ggtitle ("Sim based")+
                  geom_point(position = position_dodge(width = 0.25)) +
                  scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
                  scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
                  geom_hline(yintercept = 99, col = 'red')+
                  labs(x = NULL, y = "Repeatability [%]")+
                  ylim(c(60,100))+
                  scale_y_continuous(breaks = sort(c(seq(60, 100, length.out=6), 99)))+
                  coord_flip()+
                  theme_bw() +
                  theme(plot.title = element_text(size=9))
            
              ggsave('Output/test-40_Repeatability-MB-MC_BASIC.png',g1, width = 10, height =7, units = 'cm')
          # estimate COMPOSITE
              lfsim = list()
              lfrpt = list()
              for(i in c("Head", "Flagellum","Total")){
                part_ = i
                # part_ = "Flagellum"
                dd = bc[part == part_ & who %in% c('MC','MB')]
                m = lmer(Pixels ~ 1+(1|id), dd)
                Rf = R_out(part_)
                lfsim[[i]] = Rf[, method_CI:='arm package']

                R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
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
              xy2 = rbind(x,y)
              xy2[, part := factor(part, levels=c("Head","Flagellum","Total"))] 
              xy2[nchar(CI) == 7 & upr == 10, upr := 100]
          # plot COMPOSITE
              g1 = 
                ggplot(xy2, aes(x = part, y = pred, col = method_CI)) +
                  geom_errorbar(aes(ymin = lwr, ymax = upr, col = method_CI), width = 0.1, position = position_dodge(width = 0.25) ) +
                  #ggtitle ("Sim based")+
                  geom_point(position = position_dodge(width = 0.25)) +
                  scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
                  scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
                  labs(x = NULL, y = "Repeatability [%]")+
                  ylim(c(min(xy2$lwr),max(xy2$upr)))+
                  coord_flip()+
                  theme_bw() +
                  theme(plot.title = element_text(size=9))
            
              ggsave('Output/test-40_Repeatability-MB-MC_COMPOSITE.png',g1, width = 10, height =7, units = 'cm')
    # plot repeatability KT, MC
          # estimate BASIC 
              lfsim = list()
              lfrpt = list()
              for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
                part_ = i
                # part_ = "Acrosome"
                dd = b[part == part_ & who %in% c('MC','KT')]
                m = lmer(Pixels ~ 1+(1|id), dd)
                Rf = R_out(part_)
                lfsim[[i]] = Rf[, method_CI:='arm package']

                R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
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
              xy[, part := factor(part, levels=c("Acrosome", "Nucleus", "Midpiece","Tail"))] 
              xy[nchar(CI) == 7 & upr == 10, upr := 100]
          # plot BASIC
              g1 = 
                ggplot(xy, aes(x = part, y = pred, col = method_CI)) +
                  geom_errorbar(aes(ymin = lwr, ymax = upr, col = method_CI), width = 0.1, position = position_dodge(width = 0.25) ) +
                  #ggtitle ("Sim based")+
                  geom_point(position = position_dodge(width = 0.25)) +
                  scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
                  scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
                  geom_hline(yintercept = 99, col = 'red')+
                  labs(x = NULL, y = "Repeatability [%]")+
                  ylim(c(60,100))+
                  scale_y_continuous(breaks = sort(c(seq(60, 100, length.out=6), 99)))+
                  coord_flip()+
                  theme_bw() +
                  theme(plot.title = element_text(size=9))
            
              ggsave('Output/test-40_Repeatability-KT-MC_BASIC.png',g1, width = 10, height =7, units = 'cm')
          # estimate COMPOSITE
              lfsim = list()
              lfrpt = list()
              for(i in c("Head", "Flagellum","Total")){
                part_ = i
                # part_ = "Flagellum"
                dd = bc[part == part_ & who %in% c('MC','KT')]
                m = lmer(Pixels ~ 1+(1|id), dd)
                Rf = R_out(part_)
                lfsim[[i]] = Rf[, method_CI:='arm package']

                R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
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
              xy2 = rbind(x,y)
              xy2[, part := factor(part, levels=c("Head","Flagellum","Total"))] 
              xy2[nchar(CI) == 7 & upr == 10, upr := 100]
          # plot COMPOSITE
              g1 = 
                ggplot(xy2, aes(x = part, y = pred, col = method_CI)) +
                  geom_errorbar(aes(ymin = lwr, ymax = upr, col = method_CI), width = 0.1, position = position_dodge(width = 0.25) ) +
                  #ggtitle ("Sim based")+
                  geom_point(position = position_dodge(width = 0.25)) +
                  scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
                  scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
                  labs(x = NULL, y = "Repeatability [%]")+
                  ylim(c(min(xy2$lwr),max(xy2$upr)))+
                  coord_flip()+
                  theme_bw() +
                  theme(plot.title = element_text(size=9))
            
              ggsave('Output/test-40_Repeatability-KT-MC_COMPOSITE.png',g1, width = 10, height =7, units = 'cm')
    
    # plot combined
      xyMBMC[,comparison := 'Martin - Maggie']
xyKTMC[,comparison := 'Kim - Maggie']
xyKTMB[,comparison := 'Kim - Martin']
xy[,comparison := 'all three']
xy_ = rbind(xy, xyKTMC, xyKTMB, xyMBMC)
xy_[pred==10,pred:=100]  

xy2MBMC[,comparison := 'Martin - Maggie']
xy2KTMC[,comparison := 'Kim - Maggie']
xy2KTMB[,comparison := 'Kim - Martin']
xy2_ = rbind(xy2MBMC,xy2KTMC,xy2KTMB)


ggplot(xy_, aes(x = part, y = pred, col = comparison, shape = method_CI)) +
          geom_errorbar(aes(ymin = lwr, ymax = upr, col = comparison), width = 0.1, position = position_dodge(width = 0.7) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.7)) +
          #scale_fill_brewer(palette = "Set1", guide = guide_legend(reverse = TRUE)) +
          #scale_color_brewer(palette = "Set1", guide = guide_legend(reverse = TRUE))  +
          scale_color_viridis(discrete=TRUE, guide = guide_legend(reverse = TRUE))  +
          scale_fill_viridis(discrete=TRUE, guide = guide_legend(reverse = TRUE)) + 
          scale_shape(guide = guide_legend(reverse = TRUE)) + 
          geom_hline(yintercept = 99, col = 'red')+
          labs(x = NULL, y = "Repeatability [%]")+
          ylim(c(min(xy_$lwr),max(xy_$upr)))+
          scale_y_continuous(breaks = c(65, 70, 85, 100))+
          geom_text(y =99, x =0.5, label = '99', col = 'red', size = 2.7)+
          coord_flip()+
          theme_bw() +
          theme(plot.title = element_text(size=9),
              axis.line = element_blank(),
              #axis.line = element_line(colour="grey70", size=0.25),
              axis.title = element_text(size=7, colour="grey30"),
              axis.title.y = element_text(vjust=3.5),
              axis.title.x = element_text(vjust=1),
              axis.text = element_text(size=6),#, vjust = 0.5, hjust=1),# margin=units(0.5,"mm")),
              axis.ticks.length=unit(0.5,"mm"),
              axis.ticks = element_line(colour = "grey70", size = 0.1),
              #axis.ticks.margin, 
              #panel.border=element_line(colour="grey70", size=0.25),
              #panel.grid = element_blank(),
              legend.text=element_text(size=6),
              legend.title=element_text(size=6),
              legend.key = element_rect(colour = NA, fill = NA),
              legend.key.height= unit(0.5,"line"),
              legend.key.width = unit(0.25, "cm"),
              legend.margin = margin(0,0,0,0, unit="cm"),
              legend.box.margin = margin(l = -6), #legend.justification = c(-1,0),
              legend.background = element_blank()
              )


# old protocol - within Martin 
  # DATA
    dir = list.dirs('Pictures/test_40')
    dir = dir[grepl('test_40/Results', dir)]
    d = foreach(j = dir, .combine = rbind) %do% {
      #j = dir[1]
      x = fread(paste0(j,'/results.csv'))
      x[, who := substr(j, 48, 49)]
      x[, datetime_ := substr(j, 27, 45)]
      x[, id := substr(x$'File Name', 1, 2)]
      x[, pk := 1:nrow(x)]
      x[, manip := substr(x$'File Name', 4, nchar(x$'File Name')-7-(nchar(pk)))]
      xx = foreach(i = unique(x$id), .combine = rbind) %do% {
        #i = '01'
        xi = x[id == i]
        if(nrow(xi) == 4){
          xi[Pixels == max(Pixels), part := 'Tail']
          xi[Pixels != max(Pixels), part := c('Acrosome','Nucleus','Midpiece')]
        }else if(nrow(xi) == 3){
          xi[, part := c('Acrosome','Nucleus','Midpiece')]
        }else if(nrow(xi) == 2){
          xi[, part := c('Acrosome','Nucleus')]
        }else {
          xi[, part := 'Tail']
        }
        return(xi)
        }    
      return(xx)
      }
     d[datetime_ == '2021-08-10 20-50-34', attempt := 1]
     d[datetime_ != '2021-08-10 20-50-34', attempt := 2]
     d[, id_part := paste(id, part)]

     b = d[id_part %in% d[duplicated(id_part), id_part]] # only sperm ids and parts measured twice

     # composite parts
         b[, id_attempt := paste(id, attempt)]
         bt = data.table(ddply(b,.(who, id, attempt,id_attempt), summarise, part = 'Total', Pixels = sum(Pixels)))
         bt = bt[Pixels>1000]
         bh = data.table(ddply(b[part %in% c('Acrosome','Nucleus')],.(who, id, attempt,id_attempt), summarise, part = 'Head', Pixels = sum(Pixels)))
         bf = data.table(ddply(b[part %in% c('Midpiece','Tail')],.(who, id, attempt, id_attempt), summarise, part = 'Flagellum', Pixels = sum(Pixels)))
         bf = bf[Pixels>1000]
         bc = rbind(bh,bf,bt)

     # wide format
         bw = reshape(b[,.(who,id,attempt,part,Pixels)], idvar = c('who','id','part'), timevar = 'attempt', direction = "wide")

         bcw = reshape(bc[,.(who,id,attempt,part,Pixels)], idvar = c('who','id','part'), timevar = 'attempt', direction = "wide")    
  # ANAL
    # plot correlations
        g1 = ggplot(bw, aes(x = Pixels.1, y = Pixels.2)) +
        facet_wrap(~part, scales = "free") +  
        stat_smooth(method = 'lm')+geom_point()+
        stat_cor(method="pearson",size = 2) +
        geom_abline(b = 1, col = 'red', lty = 3) + 
        #xlim(c(485, 565)) +  ylim(c(485, 565)) + 
        #ggtitle('All')+
        theme_bw()
        ggsave('Output//test-40_cor-MB.png',g1, width = 10, height =10, units = 'cm')
        
        g2 = ggplot(bcw, aes(x = Pixels.1, y = Pixels.2)) +
        facet_wrap(~part, scales = "free", nrow = 2) +  
        stat_smooth(method = 'lm')+geom_point()+
        stat_cor(method="pearson",size = 2) +
        geom_abline(b = 1, col = 'red', lty = 3) + 
        #xlim(c(485, 565)) +  ylim(c(485, 565)) + 
        #ggtitle('All')+
        theme_bw()
        ggsave('Output//test-40_cor-MB_composite.png',g2, width = 10, height =10, units = 'cm')

    # plot repeatability
          # estimate BASIC 
              lfsim = list()
              lfrpt = list()
              for(i in c('Acrosome','Nucleus','Midpiece','Tail')){
                part_ = i
                # part_ = "Acrosome"
                dd = b[part == part_]
                m = lmer(Pixels ~ 1+(1|id), dd)
                Rf = R_out(part_)
                lfsim[[i]] = Rf[, method_CI:='arm package']

                R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
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
              xy[, part := factor(part, levels=c("Acrosome", "Nucleus", "Midpiece","Tail"))] 
              xy[nchar(CI) == 7 & upr == 10, upr := 100]
          # plot BASIC
              g1 = 
                ggplot(xy, aes(x = part, y = pred, col = method_CI)) +
                  geom_errorbar(aes(ymin = lwr, ymax = upr, col = method_CI), width = 0.1, position = position_dodge(width = 0.25) ) +
                  #ggtitle ("Sim based")+
                  geom_point(position = position_dodge(width = 0.25)) +
                  scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
                  scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
                  labs(x = NULL, y = "Repeatability [%]")+
                  ylim(c(85,100))+
                  coord_flip()+
                  theme_bw() +
                  theme(plot.title = element_text(size=9))
            
              ggsave('Output/test-40_Repeatability-MB_BASIC.png',g1, width = 10, height =7, units = 'cm')
          # estimate COMPOSITE
              lfsim = list()
              lfrpt = list()
              for(i in c("Head", "Flagellum","Total")){
                part_ = i
                # part_ = "Flagellum"
                dd = bc[part == part_]
                m = lmer(Pixels ~ 1+(1|id), dd)
                Rf = R_out(part_)
                lfsim[[i]] = Rf[, method_CI:='arm package']

                R = rpt(Pixels ~ (1 | id), grname = "id", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
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
              xy[, part := factor(part, levels=c("Head","Flagellum","Total"))] 
              xy[nchar(CI) == 7 & upr == 10, upr := 100]
          # plot COMPOSITE
              g1 = 
                ggplot(xy, aes(x = part, y = pred, col = method_CI)) +
                  geom_errorbar(aes(ymin = lwr, ymax = upr, col = method_CI), width = 0.1, position = position_dodge(width = 0.25) ) +
                  #ggtitle ("Sim based")+
                  geom_point(position = position_dodge(width = 0.25)) +
                  scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
                  scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
                  labs(x = NULL, y = "Repeatability [%]")+
                  ylim(c(97,100))+
                  coord_flip()+
                  theme_bw() +
                  theme(plot.title = element_text(size=9))
            
              ggsave('Output/test-40_Repeatability-MB_COMPOSITE.png',g1, width = 10, height =7, units = 'cm')
