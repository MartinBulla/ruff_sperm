#' ---
#' title: "Motility of ruff sperm - preliminary output"
#' author: "Martin Bulla"
#' date: "`r Sys.time()`"
#' output: 
#'     html_document:
#'         toc: true
#'         toc_float: true
#'         toc_depth: 5
#'         code_folding: hide
#' ---

#+ r setup, include=FALSE 
knitr::opts_chunk$set(message = FALSE, warning = FALSE, cache = TRUE)

#' Code to load tools and data
  # TOOLS 
    require(here)
    source(here::here('R/tools.R'))
    
    require(arm) 
    require(effects)
    require(ggpubr)
    require(ggsci) 
    require(gridExtra)
    require(magrittr)
    require(multcomp)
    require(PerformanceAnalytics)
    require(rptR) 
    require(stringi)
    require(viridis)


    # constants
      round_ = 3 # number of decimal places to round model coefficients
      nsim = 5000 # number of simulations to extract estimates and 95%CrI
      ax_lines = "grey60" # defines color of the axis lines
      colors <- c("#999999", "#E69F00", "#56B4E9") #viridis(3)
    
    # functions
      getime = function (x) {ifelse(is.na(x), as.numeric(NA), as.numeric(difftime(x, trunc(x,"day"), units = "hours")))}
      
      getDay = function (x) {as.Date(trunc(x, "day"))}
     
      # for R output based on sim
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
  # DATA  - TO DO add 'manip' column
    # composite measures
      x = fread(here::here('Data/DAT_morpho.csv'))
      setnames(x,old = 'pic', new = 'sperm_ID')
      x[, sample_ID:=as.character(sample_ID)]

      d1 = x[,.(bird_ID, sample_ID, sperm_ID, part, Pixels)]
      
      d2 = d1[, sum(Pixels), by = list(bird_ID, sample_ID, sperm_ID)]
      names(d2)[4] = 'Pixels'
      d2[,part := 'Total']

      d3 = d1[part%in%c('Midpiece', 'Tail'), sum(Pixels), by = list(bird_ID, sample_ID, sperm_ID)]
      names(d3)[4] = 'Pixels'
      d3[,part := 'Flagellum']

      d4 = d1[part%in%c('Acrosome', 'Nucleus'), sum(Pixels), by = list(bird_ID, sample_ID, sperm_ID)]
      names(d4)[4] = 'Pixels'
      d4[,part := 'Head']

      b = merge(d1,d2, all = TRUE)
      b = merge(b,d3, all = TRUE)
      b = merge(b,d4, all = TRUE)
      b[, Length_µm:=Pixels*0.078]

    # add metadata and motility
      b = merge(b, x[,.(bird_ID, sample_ID, sperm_ID, part, manip)], all.x = TRUE, by=c('bird_ID', 'sample_ID', 'sperm_ID', 'part'))
      b[manip%in%"", manip:=NA]
      
      d = data.table(read_excel(here::here('Data/motility.xlsx'), sheet = 1))
      s = data.table(read_excel(here::here('Data/sampling_2021_cleaned.xlsx'), sheet = 1))
      s = s[!is.na(recording)]

      m = data.table(read_excel(here::here('Data/ruff_males_Seewiesen.xlsx'), sheet = 1))#, range = "A1:G161"))
      m[, hatch_year:=as.numeric(substr(Hatch_date,1,4)) ]
      m[, age := 2021-hatch_year]

      m = m[Ind_ID %in%unique(s$DB_ID[!is.na(s$DB_ID)]),.(Ind_ID, Morph, age)]
      s = merge(s,m, by.x = 'DB_ID', by.y = 'Ind_ID', all.x = TRUE)
      b = merge(b, s[,.(sample_ID, bird_ID, Morph, age, datetime, type, sperm, recording, rec_measured, month)], by.x = c('sample_ID', 'bird_ID'), by.y = c('sample_ID', 'bird_ID'), all.x = TRUE)
      b = merge(b, d, by.x = c('bird_ID', 'month'), by.y = c('ID', 'date'), all.x = TRUE)

      #d = (d[!duplicated(d)]) # shouldn't be necessary, but for some reason the above call duplicates the dataset
      
      b[is.na(Morph), Morph := 'Zebra finch']
      b[Morph == 'F', Morph := 'Faeder']
      b[Morph == 'I', Morph := 'Independent']
      b[Morph == 'S', Morph := 'Satellite']
      b[is.na(issues), issues := 'zero']

      b = b[!Morph %in% 'Zebra finch']
      b[, part := factor(part, levels = c('Acrosome', 'Nucleus', 'Midpiece', 'Tail','Head','Flagellum', 'Total'))]
      b[part %in% c('Acrosome', 'Nucleus', 'Midpiece', 'Tail'), measure := 'original']
      b[part %in% c('Head','Flagellum', 'Total'), measure := 'composite']

      #nrow(d) # N 139

    # prepare for correlations and repeatability
      bw = reshape(b[,.(bird_ID,Morph, age, datetime, month, sample_ID, sperm_ID, VAP,VSL,VCL, motileCount, part, Length_µm)], idvar = c('bird_ID','Morph','age','datetime', 'month', 'sample_ID', 'sperm_ID','VAP','VSL','VCL', 'motileCount'), timevar = 'part', direction = "wide")  
      setnames(bw,old = c('Length_µm.Acrosome', 'Length_µm.Nucleus','Length_µm.Flagellum','Length_µm.Head','Length_µm.Midpiece','Length_µm.Tail', 'Length_µm.Total'), new = c('Acrosome', 'Flagellum', 'Head','Nucleus', 'Midpiece', 'Tail','Total'))
    # add relative measures
      bw[, Midpiece_rel := Midpiece/Total]
      bw[, Flagellum_rel := Flagellum/Total]
    
    # mean/male dataset
      a = b[, mean(Length_µm), by = list(bird_ID, Morph, age, part)]
      setnames(a, old = 'V1', new = 'Length_avg')
      a1 = data.table(bird_ID = unique(a$bird_ID), blind = c(rep(c('Independent','Satellite', 'Faeder'), floor(length(unique(a$bird_ID))/3)), 'Independent','Satellite'))
      a = merge(a,a1, all.x = TRUE)

      aw = reshape(a, idvar = c('bird_ID','Morph', 'blind','age'), timevar = "part", direction = "wide")
      names(aw) = c('bird_ID','Morph', 'blind','age',as.character(unique(a$part)))
      aw[, Midpiece_rel := Midpiece/Total]
      aw[, Flagellum_rel := Flagellum/Total]

    # morph as factor
      b[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 
      bw[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 
      a[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 
      aw[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 

#' ## Exploration
 
#+ fig.width=10, fig.height = 14
  b[, order_ := mean(Length_µm), by = bird_ID]
  b_ = b[part =='Total']
  b_[,bird_ID := reorder(bird_ID, Length_µm, mean)]
  b[, bird_ID := factor(bird_ID, levels = levels(b_$bird_ID))]

  g = ggplot(b, aes(x = reorder(as.factor(bird_ID),Length_µm,mean), y = Length_µm)) + facet_wrap(~part, nrow = 7, scales = 'free') + #x = reorder(as.factor(bird_ID),Length_µm,mean)
    geom_boxplot(aes(col = Morph)) + 
    #geom_dotplot(binaxis = 'y', stackdir = 'center',
     #            position = position_dodge(), col = 'red', fill ='red')+
    scale_colour_manual(values = colors) + 
    #scale_color_viridis(discrete=TRUE)+
    xlab('Male ID') +
    theme_bw() +
    theme(axis.text.x = element_blank())
    #$, legend.position = "none")
  g
  #ggsave('Output/morpho_within_male_boxplots_ordered.png', g, width = 20, height = 15, units = 'cm')
  
#+ cor_parts, fig.width=7, fig.height = 7
  chart.Correlation(bw[, c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel')], histogram=TRUE, pch=19)
  mtext("Single sperm", side=3, line=3)
  #dev.copy(png,'Output/VD_corr_single.png')
  #dev.off()
#+ cor_parts_avg, fig.width=7, fig.height = 7  
  chart.Correlation(aw[, c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total'.'Midpiece_rel', 'Flagellum_rel')], histogram=TRUE, pch=19)
  mtext("Male averages", side=3, line=3)
  #dev.copy(png,'Output/VD_corr_avg.png')
  #dev.off()
#+ cor_parts_I, fig.width=7, fig.height = 7  
  chart.Correlation(bw[Morph == 'Independent', c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel')], histogram=TRUE, pch=19)
  mtext("Independent", side=3, line=3)
  #dev.copy(png,'Output/morpho_corr_single_I.png')
  #dev.off()
#+ cor_parts_S, fig.width=7, fig.height = 7    
  chart.Correlation(bw[Morph == 'Satellite', c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel')], histogram=TRUE, pch=19)
  mtext("Satellite", side=3, line=3)
  #dev.copy(png,'Output/morpho_corr_single_S.png')
  #dev.off()
#+ cor_parts_F, fig.width=7, fig.height = 7    
  chart.Correlation(bw[Morph == 'Faeder', c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel')], histogram=TRUE, pch=19)
  mtext("Feader", side=3, line=3)
  #dev.copy(png,'Output/morpho_corr_single_F.png')
  #dev.off()

# repeatability within male
  # estimate
    lfsim = list()
    lfrpt = list()
    for(i in c('Acrosome','Flagellum','Head','Midpiece','Nucleus','Tail','Total')){
      part_ = i
      # part_ = "Acrosome"
      dd = b[part == part_]
      m = lmer(Length_µm ~ 1+(1|bird_ID), dd)
      Rf = R_out(part_)
      lfsim[[i]] = Rf[, method_CI:='arm package']

      R = rpt(Length_µm ~ (1 | bird_ID), grname = "bird_ID", data = dd, datatype = "Gaussian")#, nboot = 0, npermut = 0)
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
      dd = b[part == part_]
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
    facet_wrap(~part, scales = 'free', nrow = 3)+
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
      facet_wrap(~part, scales = 'free_y', nrow = 3)+
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
      facet_wrap(~part, scales = 'free_y', nrow = 3)+
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
   cv_ =  b[, cv(Length_µm), by = list(part, Morph)]
   cv_[ , Morph123 := as.numeric(Morph)]
   names(cv_) [3]='CV' 
   g=
   ggplot(cv_, aes(x =Morph123 , y = CV, col = part, fill = part)) + geom_point() + geom_line() + scale_x_continuous(br eaks =c(1,2,3), labels = c('Independent','Satellite','Faeder')) +
        theme_bw() +
        theme(axis.title.x = element_blank())
   ggsave('Output/VD_CV_overall.png',g, width = 10, height =7, units = 'cm')

   cv_ =  b[, cv(Length_µm), by = list(bird_ID, part, Morph)]
   cv_[ , Morph123 := as.numeric(Morph)]
   names(cv_) [4]='CV' 
   g=
   ggplot(cv_, aes(x =Morph , y = CV, col = part, fill = part)) + 
   #geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(), col = 'grey',  dotsize = 0.5)+
   geom_boxplot(fill = NA) + 
   theme_bw() +
   theme(axis.title.x = element_blank())
   ggsave('Output/VD_CV_male_alternative.png',g, width = 10, height =7, units = 'cm')
   
   g=
   ggplot(cv_, aes(x =Morph , y = CV)) + 
   geom_boxplot(fill = NA) + 
   geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(), col = 'red',  fill = 'grey', dotsize = 1)+
   facet_wrap(~part, nrow = 3, scales = "free_y") +
   guides(x =  guide_axis(angle = -45)) +
   theme_bw() +
   theme(axis.title.x = element_blank())
   ggsave('Output/VD_CV_male.png',g, width = 10, height =10, units = 'cm')
    
  # correlations
    # median with 95%
    bs = b[, quantile(Length_µm, prob = c(0.5)), by = list(part,Morph)]
    names(bs)[3] = 'median'
    bs$lwr = b[, quantile(Length_µm, prob = c(0.025)), by = list(part,Morph)]$V1
    bs$upr = b[, quantile(Length_µm, prob = c(0.975)), by = list(part,Morph)]$V1

    # mean +/-sd
    bs = b[, mean(Length_µm), by = list(part,Morph)]
    names(bs)[3] = 'median'
    bs$lwr = b[, mean(Length_µm)-sd(Length_µm), by = list(part,Morph)]$V1
    bs$upr = b[, mean(Length_µm)+sd(Length_µm), by = list(part,Morph)]$V1

    bs_h = bs[part == 'Head']
    bs_h[,Midpiece_µm:= bs[part == 'Midpiece',.(median)]]

    bs_m = bs[part == 'Midpiece']
    bs_m[,Head_µm:= bs[part == 'Head',.(median)]]
    
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

    bs_h_ = bs[part == 'Head']
    bs_h_[,Flagellum_µm:= bs[part == 'Flagellum',.(median)]]

    bs_f = bs[part == 'Flagellum']
    bs_f[,Head_µm:= bs[part == 'Head',.(median)]]
      
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

    bs_n = bs[part == 'Nucleus']
    bs_n[,Acrosome_µm:= bs[part == 'Acrosome',.(median)]]

    bs_a = bs[part == 'Acrosome']
    bs_a[,Nucleus_µm:= bs[part == 'Nucleus',.(median)]] 
    
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

    as = a[, quantile(Length_avg, prob = c(0.5)), by = list(part,Morph)]
    names(as)[3] = 'median'
    as$lwr = a[, quantile(Length_avg, prob = c(0.025)), by = list(part,Morph)]$V1
    as$upr = a[, quantile(Length_avg, prob = c(0.975)), by = list(part,Morph)]$V1

    # mean +/-sd
    as = a[, mean(Length_avg), by = list(part,Morph)]
    names(as)[3] = 'median'
    as$lwr = a[, mean(Length_avg)-sd(Length_avg), by = list(part,Morph)]$V1
    as$upr = a[, mean(Length_avg)+sd(Length_avg), by = list(part,Morph)]$V1

    as_h = as[part == 'Head']
    as_h[,Midpiece_µm:= as[part == 'Midpiece',.(median)]]

    as_m = as[part == 'Midpiece']
    as_m[,Head_µm:= as[part == 'Head',.(median)]]

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

    as_h_ = as[part == 'Head']
    as_h_[,Flagellum_µm:= as[part == 'Flagellum',.(median)]]

    as_f = as[part == 'Flagellum']
    as_f[,Head_µm:= as[part == 'Head',.(median)]]

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

    as_n = as[part == 'Nucleus']
    as_n[,Acrosome_µm:= as[part == 'Acrosome',.(median)]]

    as_a = as[part == 'Acrosome']
    as_a[,Nucleus_µm:= as[part == 'Nucleus',.(median)]] 
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
        for(i in unique(a$part)){
          #i ='Nucleus'
          m = lm(scale(Length_avg) ~ Morph, a[part == i])
          #summary(m)
          #plot(allEffects(m))
          bsim = sim(m, n.sim=nsim) 
          v = apply(bsim@coef, 2, quantile, prob=c(0.5))
          ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
          l[[i]]=data.frame(response=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

          # get predictions
          m = lm(Length_avg ~ Morph, a[part == i])
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
          newD$part=i
          lp[[i]] = data.table(newD)

          print(i)     
          }          
        
        ll = data.table(do.call(rbind,l) ) 
        ll[, response := factor(response, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
        ll[, effect := factor(effect, levels=rev(c("(Intercept)", "MorphSatellite", "MorphFaeder")))]
        
        llp = data.table(do.call(rbind,lp) ) 
        llp[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
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
          facet_wrap(~part, scales = 'free_y', nrow = 3)+
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
        for(i in unique(a$part)){
          #i ='Nucleus'
          m = lmer(scale(Length_µm) ~ Morph + (1|bird_ID), b[part == i])
          #summary(m)
          #plot(allEffects(m))
          bsim = sim(m, n.sim=5000) 
          v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
          ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
          ls[[i]]=data.frame(response=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])
         
          # get predictions
          m = lmer(Length_µm ~ Morph + (1|bird_ID), b[part == i])
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
          newD$part=i
          lps[[i]] = data.table(newD)

          print(i)         
          }  
        
        lls = data.table(do.call(rbind,ls) ) 
        lls[, response := factor(response, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
        lls[, effect := factor(effect, levels=rev(c("(Intercept)", "MorphSatellite", "MorphFaeder")))]
  
        llps = data.table(do.call(rbind,lps) ) 
        llps[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
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
          facet_wrap(~part, scales = 'free_y', nrow = 3)+
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


# End     