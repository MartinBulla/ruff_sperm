# TOOLS & DATA
  require(here)
  source(here::here('R/tools.R'))
  
  require(arm) 
  require(effects)
  require(ggpubr)
  require(gridExtra)
  require(magrittr)
  require(multcomp)
  #require(rptR) 
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

  # load data  
    d = data.table(read_excel('Data/ruffs_motility_June2021.xlsx', sheet = 1))
    s = data.table(read_excel('Data/sampling_2021-05-03_05_2021-06-07.xlsx', sheet = 1))
    m = data.table(read_excel('Data/ruff_males_Seewiesen.xlsx', sheet = 1))#, range = "A1:G161"))
    m[, hatch_year:=as.numeric(substr(Hatch_date,1,4)) ]
    m[, age := 2021-hatch_year]

    m = m[Ind_ID %in%unique(s$DB_ID[!is.na(s$DB_ID)]),.(Ind_ID, Morph)]
    s$Morph = m$Morph[match(s$DB_ID, m$Ind_ID)]
    d$Morph = s$Morph[match(d$ID, s$bird_ID)]

    d$ID[is.na(d$Morph)]
    d[is.na(Morph), Morph := 'Zebra finch']
    d[Morph == 'F', Morph := 'Faeder']
    d[Morph == 'I', Morph := 'Independent']
    d[Morph == 'S', Morph := 'Satellite']
   
# EXPLORE
  nrow(d)
  summary(d$motileCount) # N sperm/sample

  ss = s[bird_ID %in% s$bird_ID[duplicated(s$bird_ID)]]
  ss[, date := getDay(datetime)]
  ss[, time := getime(datetime)]
  ss[, date2 := getDay(datetime_2)]
  ss[, time2 := getime(datetime_2)]
  table(ss$bird_ID)

  x = ss[, by = bird_ID, .(min(date), max(date), min(time), max(time))]
  x [V1 == V2 & (V4-V3)>1]
  x [V4-V3>1]

  y = ss[, by = bird_ID, .(min(date2), max(date2), min(time2), max(time2))]
  y [V1 == V2 & (V4-V3)>1]
  y [V4-V3>1]

  cor(d$VCL,log(d$motileCount))
  
  ggplot(d,aes(x = motileCount, y = VCL)) + 
    geom_point() + 
    stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1)

  
  g1 = 
  ggplot(d,aes(x = motileCount, y = VCL, col = Morph)) + 
    geom_point(pch = 21) + 
    stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1, se = FALSE)+
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
    stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1, se = FALSE)+
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
    stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1, se = FALSE) +
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
  ggsave('Output/Motility-SpermN.png',cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =8*1.5, units = 'cm')  
     

# TEST Morph-blind
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


# Real Morph-values    
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


  m = lm(VCL ~ log(motileCount) + Morph, d[Morph!='Zebra finch']) 
  m = lm(VAP ~ log(motileCount) + Morph, d[Morph!='Zebra finch']) 
  m = lm(VSL ~ log(motileCount) + Morph, d[Morph!='Zebra finch']) 
  summary(glht(m)) 
  summary(m)
  plot(allEffects(m))