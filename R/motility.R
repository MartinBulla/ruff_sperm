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

# DATA  
  d = data.table(read_excel('Data/motility.xlsx', sheet = 1))
  s = data.table(read_excel('Data/sampling_2021_cleaned.xlsx', sheet = 1))
  s = s[!is.na(recording)]

  # add morph and age
    m = data.table(read_excel('Data/ruff_males_Seewiesen.xlsx', sheet = 1))#, range = "A1:G161"))
    m[, hatch_year:=as.numeric(substr(Hatch_date,1,4)) ]
    m[, age := 2021-hatch_year]

    m = m[Ind_ID %in%unique(s$DB_ID[!is.na(s$DB_ID)]),.(Ind_ID, Morph, age)]
    s = merge(s,m, by.x = 'DB_ID', by.y = 'Ind_ID', all.x = TRUE)
    d = merge(d, s[,.(bird_ID, Morph, age)], by.x = 'ID', by.y = 'bird_ID', all.x = TRUE)
    d = (d[!duplicated(d)]) # shouldn't be necessary, but for some reason the above call duplicates the dataset
    
    d[is.na(Morph), Morph := 'Zebra finch']
    d[Morph == 'F', Morph := 'Faeder']
    d[Morph == 'I', Morph := 'Independent']
    d[Morph == 'S', Morph := 'Satellite']

  # prepare for correlations and repeatability
    dr = d[ID%in%ID[duplicated(ID)]]
    drw = reshape(dr[,.(date,ID,VAP,VSL,VCL, motileCount, Morph, age)], idvar = c('ID','Morph','age'), timevar = 'date', direction = "wide")  

# Sample size
    table(d$ID, d$date)

    j = s[datetime>as.POSIXct('2021-05-15')]
    nrow(j) # 108 recordings
    length(unique(j$bird_ID)) # from 99 sampled individuals

    nrow(d[date == 'June']) # 93 recordings analyzed
    length(unique(d[date == 'June', ID])) # for 93 individuals

    jj = unique(j$bird_ID)
    dd = unique(d[date == 'June', ID])
    jj[!jj%in%dd]
    dd[!dd%in%jj] # all IDs in motility measurements are from the sampled 

    m = s[!datetime>as.POSIXct('2021-05-15')]
    nrow(m) # 53 recordings
    length(unique(m$bird_ID)) # 47 of 
    nrow(d[date == 'May']) # 46 + 1 (A01718) without tracking data

# Exploration
  nrow(d)
  summary(d$motileCount) # N sperm/sample
  densityplot(~motileCount, group = date, data = d, auto.key = TRUE)

  cor(d$VCL,log(d$motileCount))
  
  ggplot(d,aes(x = motileCount, y = VCL)) + 
    geom_point() + 
    stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1)

  # velocity ~ N all
    g1 = 
    ggplot(d,aes(x = motileCount, y = VCL, col = Morph)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1, se = FALSE)+
      ylab('Velocity μm/s') + 
      ylim(c(0,72)) + 
      ggtitle('Curvilinear') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
         plot.title = element_text(size=9, hjust = 0.5))
    
    g2 = 
    ggplot(d,aes(x = motileCount, y = VAP, col = Morph)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1, se = FALSE)+
      ylim(c(0,72)) +
      ggtitle('Average-path') + 
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())
    
    g3 =  
    ggplot(d,aes(x = motileCount, y = VSL, col = Morph)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1, se = FALSE) +
      ylim(c(0,72)) + 
      ggtitle('Straight-line') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())

    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "first"))
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave('Output/Motility-SpermN.png',cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =8*1.5, units = 'cm')  
  # velocity ~ N by sampling date
    g1 = 
    ggplot(d,aes(x = motileCount, y = VCL, col = date)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1, se = FALSE)+
      ylab('Velocity μm/s') + 
      ylim(c(0,72)) + 
      ggtitle('Curvilinear') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_text(size=9, hjust = 0.5))
    
    g2 = 
    ggplot(d,aes(x = motileCount, y = VAP, col = date)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1, se = FALSE)+
      ylim(c(0,72)) +
      ggtitle('Average-path') + 
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())
    
    g3 =  
    ggplot(d,aes(x = motileCount, y = VSL, col = date)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1, se = FALSE) +
      ylim(c(0,72)) + 
      ggtitle('Straight-line') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())

    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "last"))
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave('Output/Motility-SpermN_by_Date.png',cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =8*1.5, units = 'cm')  
  # velocity ~ log(N)       
    g1 = 
    ggplot(d,aes(x = motileCount, y = VCL, col = Morph)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm",size = 1, se = FALSE)+
      scale_x_continuous(trans='log10')+
      ylab('Velocity μm/s') + 
      ylim(c(0,70)) + 
      ggtitle('Curvilinear velocity') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
         plot.title = element_text(size=9, hjust = 0.5))
    
    g2 = 
    ggplot(d,aes(x = motileCount, y = VAP, col = Morph)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm",size = 1, se = FALSE)+
      scale_x_continuous(trans='log10')+
      ylim(c(0,70)) +
      ggtitle('Average-path velocity') + 
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())
    
    g3 =  
    ggplot(d,aes(x = motileCount, y = VSL, col = Morph)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm",size = 1, se = FALSE)+
      scale_x_continuous(trans='log10')+
      ylim(c(0,70)) + 
      ggtitle('Straight-line velocity') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())

    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "first"))
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave('Output/Motility-logSpermN.png',cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =8*1.5, units = 'cm')  
  # N ~ date
    ggplot(d,aes(x = date, y = motileCount, col = Morph)) + 
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(), aes(col = Morph), fill ='grey40', dotsize = 1)+
      geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2,position = position_dodge()) + 
      stat_summary(fun=mean, geom="point", color="red", fill="red") +
      scale_colour_viridis(discrete=TRUE)+
      theme_bw() 
  # velocity ~ date
    g1 = 
    ggplot(d,aes(x = date, y = VCL, col = Morph)) + 
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(), aes(col = Morph), fill ='grey40', dotsize = 2)+
      #geom_boxplot(col = 'grey',alpha = 0.2,position = position_dodge())+
      #geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2,position = position_dodge())+
      stat_summary(fun=mean, geom="point", color="red", fill="red") +
      ylab('Velocity μm/s') + 
      ylim(c(0,72)) + 
      ggtitle('Curvilinear') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
         plot.title = element_text(size=9, hjust = 0.5))
    
    g2 = 
    ggplot(d,aes(x = date, y = VAP)) + 
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(), aes(col = Morph), fill ='grey40', dotsize = 2)+
      #geom_boxplot(col = 'grey',alpha = 0.2,position = position_dodge())+
      stat_summary(fun=mean, geom="point", color="red", fill="red") +
      ylim(c(0,72)) +
      ggtitle('Average-path') + 
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())
    
    g3 =  
    ggplot(d,aes(x = date, y = VSL, col = Morph)) + 
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(), aes(col = Morph), fill ='grey40', dotsize = 2)+
      #geom_boxplot(col = 'grey',alpha = 0.2,position = position_dodge())+
      #geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2,position = position_dodge())+
      stat_summary(fun=mean, geom="point", color="red", fill="red") +
      ylim(c(0,72)) + 
      ggtitle('Straight-line') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())

    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "first"))
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave('Output/Motility-Date.png',cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =8*1.5, units = 'cm')  
  # velocity ~ age
    g1 = 
    ggplot(d[date == 'June'],aes(x = age, y = VCL, col = Morph)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1, se = FALSE)+
      ylab('Velocity μm/s') + 
      ylim(c(0,72)) + 
      ggtitle('Curvilinear') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
         plot.title = element_text(size=9, hjust = 0.5))
    
    g2 = 
    ggplot(d[date == 'June'],aes(x = age, y = VAP, col = Morph)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1, se = FALSE)+
      ylim(c(0,72)) +
      ggtitle('Average-path') + 
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())
    
    g3 =  
    ggplot(d[date == 'June'],aes(x = age, y = VSL, col = Morph)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1, se = FALSE) +
      ylim(c(0,72)) + 
      ggtitle('Straight-line') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())

    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "first"))
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave('Output/Motility-age.png',cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =8*1.5, units = 'cm')  
  # N ~ age
    g1 = 
    ggplot(d[date == 'May' & Morph!='Zebra finch'],aes(x = age, y = motileCount, col = Morph)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm")+
      ylab('N sperm cell tracked') + 
      coord_cartesian(xlim = c(1,13), ylim = c(0,600)) + 
      ggtitle('May') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
         plot.title = element_text(size=9, hjust = 0.5))
    
    g2 = 
    ggplot(d[date == 'June' & Morph!='Zebra finch'],aes(x = age, y = motileCount, col = Morph)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm")+
      coord_cartesian(xlim = c(1,13), ylim = c(0,600)) + 
      ggtitle('June') + 
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())
    
    g3 =  
    ggplot(d[Morph!='Zebra finch'],aes(x = age, y = motileCount, col = Morph)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm", size = 1) +
      coord_cartesian(xlim = c(1,13), ylim = c(0,600)) + 
      ggtitle('Both') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()
          )

    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "first"))
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave('Output/Motility_N-age.png',cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =8*1.5, units = 'cm')  
  
# Correlations
        g1 = 
        ggplot(drw, aes(x = VCL.May, y = VCL.June)) +
        facet_wrap(~Morph, ncol = 1)  +
        stat_smooth(method = 'lm', aes(col = Morph))+geom_point(pch = 21, aes(col = Morph))+
        stat_cor(method="pearson",size = 2) +
        geom_abline(b = 1, col = 'red', lty = 3) + 
        xlim(c(min(c(drw$VCL.May,drw$VCL.June)), max(c(drw$VCL.May,drw$VCL.June)))) +  ylim(c(min(c(drw$VCL.May,drw$VCL.June)), max(c(drw$VCL.May,drw$VCL.June)))) + 
        ggtitle('Curvilinear')+
        ylab('Velocity in June [μm/s[') + 
        theme_bw() + 
        theme(legend.position = "none",
            axis.text = element_text(size=7), 
            axis.title.x = element_text(size = 8, color = 'white'),
            axis.title.y = element_text(size = 8),
            strip.text.x = element_text(size = 7),
            plot.title = element_text(size=8, hjust = 0.5))
        
        g2 = 
        ggplot(drw, aes(x = VAP.May, y = VAP.June)) +
        facet_wrap(~Morph, ncol = 1)  +
        stat_smooth(method = 'lm', aes(col = Morph))+geom_point(pch = 21, aes(col = Morph))+
        stat_cor(method="pearson",size = 2) +
        geom_abline(b = 1, col = 'red', lty = 3) + 
        xlim(c(min(c(drw$VAP.May,drw$VAP.June)), max(c(drw$VAP.May,drw$VAP.June)))) +  ylim(c(min(c(drw$VAP.May,drw$VAP.June)), max(c(drw$VAP.May,drw$VAP.June)))) + 
        ggtitle('Average path')+
        xlab('Velocity in May [μm/s]') + 
        theme_bw() + 
        theme(legend.position = "none",
            axis.text = element_text(size=7), 
            axis.title.y = element_blank(), 
            axis.title.x = element_text(size = 8, hjust = 0.5),
            strip.text.x = element_text(size = 7),
            plot.title = element_text(size=8, hjust = 0.5))

        g3 = 
        ggplot(drw, aes(x = VSL.May, y = VSL.June)) +
        facet_wrap(~Morph, ncol = 1)  +
        stat_smooth(method = 'lm', aes(col = Morph))+geom_point(pch = 21, aes(col = Morph))+
        stat_cor(method="pearson",size = 2) +
        geom_abline(b = 1, col = 'red', lty = 3) + 
        xlim(c(min(c(drw$VSL.May,drw$VSL.June)), max(c(drw$VSL.May,drw$VSL.June)))) +  ylim(c(min(c(drw$VSL.May,drw$VSL.June)), max(c(drw$VSL.May,drw$VSL.June)))) + 
        ggtitle('Straight line')+
        theme_bw() + 
        theme(legend.position = "none",
          axis.text = element_text(size=7), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          strip.text.x = element_text(size = 7),
          plot.title = element_text(size=8, hjust = 0.5))

        grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "first"))
        
        gg1 <- ggplotGrob(g1)
        gg2 <- ggplotGrob(g2) 
        gg3 <- ggplotGrob(g3) 
        ggsave('Output/Motility-cor_MayJune.png',cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =7*1.5, units = 'cm')  
# Repeatability
    # prepare
      velo = 'VAP'
        m = lmer(VAP ~ 1+(1|ID), dr)
        Rf_VAP = R_out(velo)
        Rf_VAP[, method_CI:='arm package']
        names(Rf_VAP)[1] = 'velocity'
        Rf_VAP[, pred:= as.numeric(substr(repeatability,1,nchar(repeatability)-1))]
        Rf_VAP[, lwr:= as.numeric(substr(CI,1,2))]
        Rf_VAP[, upr:= as.numeric(substr(CI,4,5))]

        R = rpt(VAP ~ (1 | ID), grname = "ID", data = dr, datatype = "Gaussian")#, nboot = 0, npermut = 0)
        RR = data.table(merge(data.frame(velocity =velo), paste0(round(R$R*100),'%'))) %>% setnames(new = c('velocity', 'repeatability'))
        RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
        RR[, method_CI := 'rpt package']
        RR[, pred := 100*R$R]
        RR[, lwr := 100*R$CI_emp[1]]
        RR[, upr := 100*R$CI_emp[2]]
        RR_VAP = RR

        VAP = rbind(Rf_VAP, RR_VAP)

      velo = 'VSL'
        m = lmer(VSL ~ 1+(1|ID), dr)
        Rf_VSL = R_out(velo)
        Rf_VSL[, method_CI:='arm package']
        names(Rf_VSL)[1] = 'velocity'
        Rf_VSL[, pred:= as.numeric(substr(repeatability,1,nchar(repeatability)-1))]
        Rf_VSL[, lwr:= as.numeric(substr(CI,1,2))]
        Rf_VSL[, upr:= as.numeric(substr(CI,4,5))]

        R = rpt(VSL ~ (1 | ID), grname = "ID", data = dr, datatype = "Gaussian")#, nboot = 0, npermut = 0)
        RR = data.table(merge(data.frame(velocity =velo), paste0(round(R$R*100),'%'))) %>% setnames(new = c('velocity', 'repeatability'))
        RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
        RR[, method_CI := 'rpt package']
        RR[, pred := 100*R$R]
        RR[, lwr := 100*R$CI_emp[1]]
        RR[, upr := 100*R$CI_emp[2]]
        RR_VSL = RR

        VSL = rbind(Rf_VSL, RR_VSL)

      velo = 'VCL'
        m = lmer(VCL ~ 1+(1|ID), dr)
        Rf_VCL = R_out(velo)
        Rf_VCL[, method_CI:='arm package']
        names(Rf_VCL)[1] = 'velocity'
        Rf_VCL[, pred:= as.numeric(substr(repeatability,1,nchar(repeatability)-1))]
        Rf_VCL[, lwr:= as.numeric(substr(CI,1,2))]
        Rf_VCL[, upr:= as.numeric(substr(CI,4,5))]

        R = rpt(VCL ~ (1 | ID), grname = "ID", data = dr, datatype = "Gaussian")#, nboot = 0, npermut = 0)
        RR = data.table(merge(data.frame(velocity =velo), paste0(round(R$R*100),'%'))) %>% setnames(new = c('velocity', 'repeatability'))
        RR[, CI := paste0(paste(round(R$CI_emp*100)[1], round(R$CI_emp*100)[2], sep = "-"), '%')] 
        RR[, method_CI := 'rpt package']
        RR[, pred := 100*R$R]
        RR[, lwr := 100*R$CI_emp[1]]
        RR[, upr := 100*R$CI_emp[2]]
        RR_VCL = RR

        VCL = rbind(Rf_VCL, RR_VCL)

      r = rbind(VCL,VAP,VSL) 

      r[velocity == 'VCL', velocity := 'Curvilinear']
      r[velocity == 'VAP', velocity := 'Average path']
      r[velocity == 'VSL', velocity := 'Straight line']

    # plot    
      ggplot(r, aes(x = velocity, y = pred, col = method_CI)) +
            geom_errorbar(aes(ymin = lwr, ymax = upr, col = method_CI), width = 0, position = position_dodge(width = 0.4) ) +
            #ggtitle ("Sim based")+
            geom_point(position = position_dodge(width = 0.4)) +
            #scale_fill_brewer(palette = "Set1", guide = guide_legend(reverse = TRUE)) +
            #scale_color_brewer(palette = "Set1", guide = guide_legend(reverse = TRUE))  +
            scale_color_viridis(discrete=TRUE, guide = guide_legend(reverse = TRUE))  +
            scale_fill_viridis(discrete=TRUE, guide = guide_legend(reverse = TRUE)) + 
            #scale_shape(guide = guide_legend(reverse = TRUE)) + 
            #geom_hline(yintercept = 99, col = 'red')+
            labs(x = 'Velocity', y = "Repeatability [%]")+
            scale_y_continuous(expand = c(0, 0), lim = c(0,80)) +
            #scale_y_continuous(breaks = c(65, 70, 80, 85, 90, 95, 100))+
            #geom_text(y =99, x =0.6, label = '99', col = 'red', size = 3)+
            coord_flip()+
            theme_bw() +
            theme(
                #axis.title.x = element_blank(),
                legend.text=element_text(size=6),
                legend.title=element_text(size=7),
                legend.key = element_rect(colour = NA, fill = NA),
                legend.key.height= unit(0.5,"line"),
                legend.key.width = unit(0.25, "cm"),
                legend.margin = margin(0,0,0,0, unit="cm"),
                legend.box.margin = margin(l = -6), #legend.justification = c(-1,0),
                legend.background = element_blank()
                )
      ggsave('Output/Motility-repeatability.png', width = 7*1.5, height =2.75*1.5, units = 'cm')  

# Models - TEST Morph-blind
    densityplot(d$VCL)
    d[VCL<40]
    sample(c('a','b','c'), size = nrow(d), replace = TRUE)
    d[, blind := sample(c('a','b','c'), size = nrow(d), replace = TRUE)]
    ggplot(d,aes(x = blind, y = VCL)) + 
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(), col = 'grey', aes(fill =blind), dotsize = 2)+
      geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
      stat_summary(fun=mean, geom="point", color="red", fill="red") +
      scale_fill_viridis(discrete=TRUE)+
      ylab('Curvilinear velocity') +
      theme_bw() +
      #guides(x =  guide_axis(angle = -45)) +
      theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )

    m = lm(VCL ~ log(motileCount) + blind, d) 
    summary(glht(m)) 
# Models - Real Morph-values 
  # June   
    g1= 
    ggplot(d[date == 'June'],aes(x = Morph, y = VCL)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        ylab('VCL [μm/s]') +
        theme_bw() +
        guides(x =  guide_axis(angle = -45)) +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          plot.title = element_text(size=9)
          )

    g2= 
    ggplot(d[date == 'June'],aes(x = Morph, y = VSL)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        ylab('VSL [μm/s]') +
        theme_bw() +
        guides(x =  guide_axis(angle = -45)) +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          plot.title = element_text(size=9)
          )
    
    g3= 
    ggplot(d[date == 'June'],aes(x = Morph, y = VAP)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        ylab('VAP [μm/s]') +
        theme_bw() +
        guides(x =  guide_axis(angle = -45)) +
        theme(
          legend.position = "none",
          #axis.title.x = element_blank(), 
          #axis.text.x = element_blank(),
          plot.title = element_text(size=9)
          )   

    ggplot(d[Morph!='Zebra finch' & date == 'June'],aes(x = Morph, y = VCL)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        ylab('Curvilinear velocity') +
        theme_bw() +
        guides(x =  guide_axis(angle = -45)) +
        theme(
          legend.position = "none",
          #axis.title.x = element_blank(), 
          #axis.text.x = element_blank(),
          plot.title = element_text(size=9)
          )    
    
    grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "last"))
        #grid.arrange(g1,g2)
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave('Output/Motility_Morp_boxplots.png',rbind(gg1,gg2,gg3, size = "last"), width = 7*1.5, height =15*1.5, units = 'cm')  

    m = lm(VCL ~ log(motileCount) + Morph, d[Morph!='Zebra finch' & date == 'June']) 
    m = lm(VAP ~ log(motileCount) + Morph, d[Morph!='Zebra finch' & date == 'June']) 
    m = lm(VSL ~ log(motileCount) + Morph, d[Morph!='Zebra finch' & date == 'June']) 
    summary(glht(m)) 
    summary(m)
    plot(allEffects(m))
  # May   
    g1= 
    ggplot(d[date == 'June'],aes(x = Morph, y = VCL)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        ylab('VCL [μm/s]') +
        theme_bw() +
        guides(x =  guide_axis(angle = -45)) +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          plot.title = element_text(size=9)
          )

    g2= 
    ggplot(d[date == 'June'],aes(x = Morph, y = VSL)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        ylab('VSL [μm/s]') +
        theme_bw() +
        guides(x =  guide_axis(angle = -45)) +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          plot.title = element_text(size=9)
          )
    
    g3= 
    ggplot(d[date == 'June'],aes(x = Morph, y = VAP)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        ylab('VAP [μm/s]') +
        theme_bw() +
        guides(x =  guide_axis(angle = -45)) +
        theme(
          legend.position = "none",
          #axis.title.x = element_blank(), 
          #axis.text.x = element_blank(),
          plot.title = element_text(size=9)
          )   

    ggplot(d[Morph!='Zebra finch' & date == 'June'],aes(x = Morph, y = VCL)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        ylab('Curvilinear velocity') +
        theme_bw() +
        guides(x =  guide_axis(angle = -45)) +
        theme(
          legend.position = "none",
          #axis.title.x = element_blank(), 
          #axis.text.x = element_blank(),
          plot.title = element_text(size=9)
          )    
    
    grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "last"))
        #grid.arrange(g1,g2)
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave('Output/Motility_Morp_boxplots.png',rbind(gg1,gg2,gg3, size = "last"), width = 7*1.5, height =15*1.5, units = 'cm')  

    m = lm(VCL ~ log(motileCount) + Morph, d[Morph!='Zebra finch' & date == 'June']) 
    m = lm(VAP ~ log(motileCount) + Morph, d[Morph!='Zebra finch' & date == 'June']) 
    m = lm(VSL ~ log(motileCount) + Morph, d[Morph!='Zebra finch' & date == 'June']) 
    summary(glht(m)) 
    summary(m)
    plot(allEffects(m))
  # ALL 
    # plot together 
      g1= 
      ggplot(d,aes(x = Morph, y = VCL)) + 
          geom_dotplot(binaxis = 'y', stackdir = 'center',
                       position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
          geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
          stat_summary(fun=mean, geom="point", color="red", fill="red") +
          scale_fill_viridis(discrete=TRUE)+
          ylab('VCL [μm/s]') +
          theme_bw() +
          guides(x =  guide_axis(angle = -45)) +
          theme(
            legend.position = "none",
            axis.title.x = element_blank(), 
            axis.text.x = element_blank(),
            plot.title = element_text(size=9)
            )

      g2= 
      ggplot(d,aes(x = Morph, y = VSL)) + 
          geom_dotplot(binaxis = 'y', stackdir = 'center',
                       position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
          geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
          stat_summary(fun=mean, geom="point", color="red", fill="red") +
          scale_fill_viridis(discrete=TRUE)+
          ylab('VSL [μm/s]') +
          theme_bw() +
          guides(x =  guide_axis(angle = -45)) +
          theme(
            legend.position = "none",
            axis.title.x = element_blank(), 
            axis.text.x = element_blank(),
            plot.title = element_text(size=9)
            )
      
      g3= 
      ggplot(d,aes(x = Morph, y = VAP)) + 
          geom_dotplot(binaxis = 'y', stackdir = 'center',
                       position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
          geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
          stat_summary(fun=mean, geom="point", color="red", fill="red") +
          scale_fill_viridis(discrete=TRUE)+
          ylab('VAP [μm/s]') +
          theme_bw() +
          guides(x =  guide_axis(angle = -45)) +
          theme(
            legend.position = "none",
            #axis.title.x = element_blank(), 
            #axis.text.x = element_blank(),
            plot.title = element_text(size=9)
            )   

      ggplot(d[Morph!='Zebra finch'],aes(x = Morph, y = VCL)) + 
          geom_dotplot(binaxis = 'y', stackdir = 'center',
                       position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
          geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
          stat_summary(fun=mean, geom="point", color="red", fill="red") +
          scale_fill_viridis(discrete=TRUE)+
          ylab('Curvilinear velocity') +
          theme_bw() +
          guides(x =  guide_axis(angle = -45)) +
          theme(
            legend.position = "none",
            #axis.title.x = element_blank(), 
            #axis.text.x = element_blank(),
            plot.title = element_text(size=9)
            )    
      
      grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "last"))
          #grid.arrange(g1,g2)
      gg1 <- ggplotGrob(g1)
      gg2 <- ggplotGrob(g2) 
      gg3 <- ggplotGrob(g3) 
      ggsave('Output/Motility_Morp_boxplots.png',rbind(gg1,gg2,gg3, size = "last"), width = 7*1.5, height =15*1.5, units = 'cm')  
    # plot by date
      g1= 
      ggplot(d,aes(x = Morph, y = VCL, col = date)) + 
          geom_dotplot(binaxis = 'y', stackdir = 'center',
                       dotsize = 1, aes(col = date, fill = date),
                       position = position_dodge(width = 0.3))+
          geom_boxplot(outlier.color = NA, fill = NA, alpha = 0.2, width = 0.3, position = position_dodge(1.1)) + 
          #stat_summary(fun=mean, geom="point", color="red", fill="red") +
          #scale_color_viridis(discrete=TRUE)+
          scale_fill_npg( alpha = 0.2)+
          scale_color_npg()+
          #scale_color_viridis(discrete=TRUE)+
          #scale_fill_viridis(discrete=TRUE, alpha = 0.2)+
          ylab('VCL [μm/s]') +
          theme_bw() +
          guides(x =  guide_axis(angle = -45)) +
          theme(
            legend.position = "top",
            axis.title.x = element_blank(), 
            axis.text.x = element_blank(),
            plot.title = element_text(size=9)
            )

      g2= 
      ggplot(d,aes(x = Morph, y = VSL, col = date)) + 
          geom_dotplot(binaxis = 'y', stackdir = 'center',
                       dotsize = 1, aes(col = date,fill = date),
                       position = position_dodge(width = 0.3))+
          geom_boxplot(outlier.color = NA, fill = NA, alpha = 0.2, width = 0.3, position = position_dodge(1.1)) + 
          #stat_summary(fun=mean, geom="point", color="red", fill="red") +
          scale_fill_npg( alpha = 0.2)+
          scale_color_npg()+
          #scale_color_viridis(discrete=TRUE)+
          #scale_fill_viridis(discrete=TRUE, alpha = 0.2)+
          ylab('VSL [μm/s]') +
          theme_bw() +
          guides(x =  guide_axis(angle = -45)) +
          theme(
            legend.position = "none",
            axis.title.x = element_blank(), 
            axis.text.x = element_blank(),
            plot.title = element_text(size=9)
            )
      
      g3= 
      ggplot(d,aes(x = Morph, y = VAP, col = date)) + 
          geom_dotplot(binaxis = 'y', stackdir = 'center',
                       dotsize = 1, aes(col = date,fill = date),
                       position = position_dodge(width = 0.3))+
          geom_boxplot(outlier.color = NA, fill = NA, alpha = 0.2, width = 0.3, position = position_dodge(1.1)) + 
          #stat_summary(fun=mean, geom="point", color="red", fill="red") +
          scale_fill_npg( alpha = 0.2)+
          scale_color_npg()+
          #scale_color_viridis(discrete=TRUE)+
          #scale_fill_viridis(discrete=TRUE, alpha = 0.2)+
          ylab('VAP [μm/s]') +
          theme_bw() +
          guides(x =  guide_axis(angle = -45)) +
          theme(
            legend.position = "none",
            #axis.title.x = element_blank(), 
            #axis.text.x = element_blank(),
            plot.title = element_text(size=9)
            )   
  
      grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "last"))
        #grid.arrange(g1,g2)

      gg1 <- ggplotGrob(g1)
      gg2 <- ggplotGrob(g2) 
      gg3 <- ggplotGrob(g3) 
      ggsave('Output/Motility_Morp_boxplots_date.png',rbind(gg1,gg2,gg3, size = "last"), width = 7*1.5, height =15*1.5, units = 'cm') 
    
    d$Morph_F = as.factor(d$Morph)
    m = lmer(VCL ~ log(motileCount) + date + age + Morph + (1|ID), d[Morph!='Zebra finch']) 
    m = lmer(VAP ~ log(motileCount) + date + age + Morph + (1|ID), d[Morph!='Zebra finch']) 
    m = lmer(VSL ~ log(motileCount) + date + age + Morph + (1|ID), d[Morph!='Zebra finch']) 
    #m = lm(VSL ~ log(motileCount) + date + age + Morph, d[Morph!='Zebra finch']) 
    summary(glht(m)) 
    summary(m)
    plot(allEffects(m))

# END    