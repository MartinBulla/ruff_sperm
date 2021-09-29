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
  d = data.table(read_excel(here::here('Data/motility.xlsx'), sheet = 1))
  s = data.table(read_excel(here::here('Data/sampling_2021_cleaned.xlsx'), sheet = 1))
  s = s[!is.na(recording)]

  # add morph and age
    m = data.table(read_excel(here::here('Data/ruff_males_Seewiesen.xlsx'), sheet = 1))#, range = "A1:G161"))
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
    d[is.na(issues), issues := 'zero']

    #nrow(d) # N 139

  # prepare for correlations and repeatability
    dr = d[ID%in%ID[duplicated(ID)]]
    drw = reshape(dr[,.(date,ID,VAP,VSL,VCL, motileCount, Morph, age)], idvar = c('ID','Morph','age'), timevar = 'date', direction = "wide")  

#' ## Sample sizes
#+ results = 'hide'      
    #table(d$ID, d$date)

    # June
    j = s[datetime>as.POSIXct('2021-05-15')]
    nrow(j) # 108 recordings
    length(unique(j$bird_ID)) # from 99 sampled individuals

    nrow(d[date == 'June']) # 93 recordings analyzed
    length(unique(d[date == 'June', ID])) # for 93 individuals

    jj = unique(j$bird_ID)
    dd = unique(d[date == 'June', ID])
    jj[!jj%in%dd]
    dd[!dd%in%jj] # all IDs in motility measurements are from the sampled 

    # May
    m = s[!datetime>as.POSIXct('2021-05-15')]
    nrow(m) # 53 recordings
    length(unique(m$bird_ID)) # 47 of 
    nrow(d[date == 'May']) # 46 + 1 (A01718) without tracking data

#' ## Exploration
#' ### Densities  
#' ####  Velocity
#+ fig.width=6,    
    g1 = 
    ggplot(d,aes(x = VCL, col = Morph)) + 
      geom_density() + 
      xlab('Velocity μm/s') + 
      xlim(c(0,72)) + 
      ggtitle('Curvilinear') +
      theme_bw() +
      theme(
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0.15, 0.65),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.height= unit(0.5,"line"),
        legend.key.width = unit(0.25, "cm"),
        legend.margin = margin(0,0,0,0, unit="cm"),
        legend.box.margin = margin(l = -6), #legend.justification = c(-1,0),
        legend.background = element_blank(),
        plot.title = element_text(size=9, hjust = 0.5))
    g2 = 
    ggplot(d,aes(x = VAP, col = Morph)) + 
      geom_density() + 
      xlim(c(0,72)) +
      ggtitle('Average path') + 
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_text(hjust = 0.5),
          axis.title.x = element_blank(), 
          axis.text.x = element_blank())
    
    g3 =  
    ggplot(d,aes(x =  VSL, col = Morph)) + 
      geom_density() + 
      xlim(c(0,72)) +
      xlab('Velocity μm/s') + 
      ggtitle('Straight line') +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_text(color = 'white', hjust = 0.5)
          )

    grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "last"))
    
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave(here::here('Output/Motility-density.png'),rbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "last"), width = 7*1.5, height =8*1.5, units = 'cm')    
#' ####  Velocity - date
    g1 = 
    ggplot(d,aes(x = VCL, col = Morph)) + 
      geom_density() + 
      xlab('Velocity μm/s') + 
      xlim(c(0,72)) + 
      ggtitle('Curvilinear') +
      facet_wrap(~date, ncol = 2) +
      theme_bw() +
      theme(
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0.15, 0.55),
        legend.text=element_text(size=6),
        legend.title=element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.height= unit(0.5,"line"),
        legend.key.width = unit(0.25, "cm"),
        legend.margin = margin(0,0,0,0, unit="cm"),
        legend.box.margin = margin(l = -6), #legend.justification = c(-1,0),
        legend.background = element_blank(),
        plot.title = element_text(size=9, hjust = 0.5))
    g2 = 
    ggplot(d,aes(x = VAP, col = Morph)) + 
      geom_density() + 
      xlim(c(0,72)) +
      ggtitle('Average-path') + 
       facet_wrap(~date, ncol = 2) +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_text(hjust = 0.5),
          axis.title.x = element_blank(), 
          axis.text.x = element_blank())
    
    g3 =  
    ggplot(d,aes(x =  VSL, col = Morph)) + 
      geom_density() + 
      xlim(c(0,72)) +
      xlab('Velocity μm/s') + 
      ggtitle('Straight-line') +
      facet_wrap(~date, ncol = 2) +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_text(color = 'white', hjust = 0.5)
          )

    grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "last"))
    
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave(here::here('Output/Motility-density-date.png'),rbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "last"), width = 7*1.5, height =8*1.5, units = 'cm')  
#' ####  N tracked sperm
#+ fig.width=5, fig.height=2, dpi = 100    
    g1 =      
    ggplot(d,aes(x = motileCount, col = Morph)) + 
      geom_density() + 
      xlab('N tracked sperm cells') + 
      theme_bw() +
      theme(
        axis.title.y = element_blank(),
        legend.position = c(0.85, 0.75),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.height= unit(0.5,"line"),
        legend.key.width = unit(0.25, "cm"),
        legend.margin = margin(0,0,0,0, unit="cm"),
        legend.box.margin = margin(l = -6), #legend.justification = c(-1,0),
        legend.background = element_blank(),
        plot.title = element_text(size=9, hjust = 0.5))
    
    #densityplot(~motileCount, group = date, data = d, auto.key = TRUE)
    
    g2 = 
    ggplot(d,aes(x = motileCount, col = date)) + 
      geom_density() + 
      xlab('N tracked sperm cells') + 
      scale_colour_discrete(guide = guide_legend(reverse = TRUE)) + 
      theme_bw() +
      theme(
        axis.title.y = element_blank(),
        legend.position = c(0.85, 0.75),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.height= unit(0.5,"line"),
        legend.key.width = unit(0.25, "cm"),
        legend.margin = margin(0,0,0,0, unit="cm"),
        legend.box.margin = margin(l = -6), #legend.justification = c(-1,0),
        legend.background = element_blank(),
        plot.title = element_text(size=9, hjust = 0.5))
    
    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), size = "first"))
    
    ggsave(here::here('Output/Motility-density_N.png'), cbind(ggplotGrob(g1), ggplotGrob(g2), size = "first"), width = 8*1.5, height =4*1.5, units = 'cm')    
    
    summary(d$motileCount) 
#' ####  Age
#+ fig.width=7, fig.height=2, dpi = 100    
    g1 = 
    ggplot(d,aes(x = age, col = Morph)) + 
      geom_density() + 
      xlab('Age') + 
      theme_bw() +
      theme(
        legend.position = c(0.85, 0.85),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.height= unit(0.5,"line"),
        legend.key.width = unit(0.25, "cm"),
        legend.margin = margin(0,0,0,0, unit="cm"),
        legend.box.margin = margin(l = -6), #legend.justification = c(-1,0),
        legend.background = element_blank(),
        plot.title = element_text(size=9, hjust = 0.5))
     
   g2 = 
    ggplot(d,aes(x = age, col = date)) + 
      geom_density() + 
      xlab('Age') + 
      scale_colour_discrete(guide = guide_legend(reverse = TRUE)) + 
      theme_bw() +
      theme(
        axis.title.y = element_blank(),
        legend.position = c(0.85, 0.85),
        legend.text=element_text(size=6),
        legend.title=element_text(size=7),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.height= unit(0.5,"line"),
        legend.key.width = unit(0.25, "cm"),
        legend.margin = margin(0,0,0,0, unit="cm"),
        legend.box.margin = margin(l = -6), #legend.justification = c(-1,0),
        legend.background = element_blank(),
        plot.title = element_text(size=9, hjust = 0.5))
      
    g3 = 
     ggplot(d,aes(x = age, col = Morph)) + 
      geom_density() + 
      xlab('Age (log)') + 
      theme_bw() +
      scale_x_continuous(trans = 'log') + 
      theme(
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size=9, hjust = 0.5))
    
    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g3), ggplotGrob(g2), size = "first"))

    ggsave(here::here('Output/Motility-density_age.png'), cbind(ggplotGrob(g1), ggplotGrob(g3),ggplotGrob(g2), size = "first"), width = 10*1.5, height =4*1.5, units = 'cm') 

#' ### Histograms
#' ####  Velocity
    g1 = 
    ggplot(d,aes(x = VCL, col = Morph)) + 
      geom_histogram() + 
      coord_cartesian(xlim = c(0,72), ylim = c(0,20)) +
      ggtitle('Curvilinear') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
         axis.title.x = element_text(colour= 'white', hjust = 0.75),
         plot.title = element_text(size=9, hjust = 0.5))
    g2 = 
    ggplot(d,aes(x = VAP, col = Morph)) + 
      geom_histogram() + 
       coord_cartesian(xlim = c(0,72), ylim = c(0,20)) +
      xlab('Velocity μm/s') + 
      ggtitle('Average-path') + 
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()
          )
    
    g3 =  
    ggplot(d,aes(x =  VSL, col = Morph)) + 
      geom_histogram() + 
      coord_cartesian(xlim = c(0,72), ylim = c(0,20)) +
      ggtitle('Straight-line') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank()
          )

    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "first"))
    
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave(here::here('Output/Motility-hist.png'),cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =8*1.5, units = 'cm')       
#' ####  N and age
    g1 = ggplot(d,aes(x = motileCount, col = Morph)) + 
      geom_histogram() + 
      #coord_cartesian(xlim = c(0,72), ylim = c(0,20)) +
      xlab('N tracked sperm cells') + 
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none"
         )

    g2 = ggplot(d,aes(x = age, col = Morph)) + 
      geom_histogram() + 
      #coord_cartesian(xlim = c(0,72), ylim = c(0,20)) +
      xlab('Ages') + 
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none"
         )  
    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), size = "first"))
    
    ggsave(here::here('Output/Motility-hist_N&Age.png'), cbind(ggplotGrob(g1), ggplotGrob(g2), size = "first"), width = 7*1.5, height =7*1.5, units = 'cm')  

#' ### Correlations
#' #### Velocity ~ N all
    g1 = 
    ggplot(d,aes(x = motileCount, y = VCL)) + 
      geom_point(pch = 21, aes(col = Morph)) + 
      stat_smooth(method = "lm", aes(col = Morph),formula = y ~ poly(x,2), size = 1, se = FALSE)+
      stat_cor(method="pearson",size = 2) +
      ylab('Velocity μm/s') + 
      ylim(c(0,72)) + 
      ggtitle('Curvilinear') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
         axis.title.x = element_text(color="white"), 
         plot.title = element_text(size=9, hjust = 0.5))
    
    g2 = 
    ggplot(d,aes(x = motileCount, y = VAP)) + 
      geom_point(pch = 21, aes(col = Morph)) + 
      stat_smooth(method = "lm", aes(col = Morph),formula = y ~ poly(x,2), size = 1, se = FALSE)+
      stat_cor(method="pearson",size = 2) +
      ylim(c(0,72)) +
      ggtitle('Average-path') + 
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())
    
    g3 =  
    ggplot(d,aes(x = motileCount, y = VSL)) + 
      geom_point(pch = 21, aes(col = Morph)) + 
      stat_smooth(method = "lm", aes(col = Morph),formula = y ~ poly(x,2), size = 1, se = FALSE)+
      stat_cor(method="pearson",size = 2) +
      ylim(c(0,72)) + 
      ggtitle('Straight-line') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())

    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "first"))
    
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave(here::here('Output/Motility-corSpermN.png'),cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =8*1.5, units = 'cm')  
#' #### Velocity ~ N by sampling date
    g1 = 
    ggplot(d,aes(x = motileCount, y = VCL, col = date)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1, se = FALSE)+
      #stat_cor(method="pearson",size = 2) +
      ylab('Velocity μm/s') + 
      ylim(c(0,72)) + 
      ggtitle('Curvilinear') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
        axis.title.x = element_text(color = "white"),
        plot.title = element_text(size=9, hjust = 0.5))
    
    g2 = 
    ggplot(d,aes(x = motileCount, y = VAP, col = date)) + 
      geom_point(pch = 21) + 
      stat_smooth(method = "lm", formula = y ~ poly(x,2), size = 1, se = FALSE)+
      #stat_cor(method="pearson",size = 2) +
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
      #stat_cor(method="pearson",size = 2) +
      ylim(c(0,72)) + 
      ggtitle('Straight-line') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())

    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "first"))
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave(here::here('Output/Motility-corSpermN_by_Date.png'),cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =8*1.5, units = 'cm')  
#' #### Velocity ~ log(N)       
    g1 = 
    ggplot(d,aes(x = motileCount, y = VCL)) + 
      geom_point(pch = 21, aes(col = Morph)) + 
      stat_smooth(method = "lm", aes(col = Morph),size = 1, se = FALSE)+
      stat_cor(method="pearson",size = 2) +
      scale_x_continuous(trans='log10')+
      ylab('Velocity μm/s') + 
      ylim(c(0,70)) + 
      ggtitle('Curvilinear') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
         plot.title = element_text(size=9, hjust = 0.5))
    
    g2 = 
    ggplot(d,aes(x = motileCount, y = VAP)) + 
      geom_point(pch = 21, aes(col = Morph)) + 
      stat_smooth(method = "lm", aes(col = Morph),size = 1, se = FALSE)+
      stat_cor(method="pearson",size = 2) +
      scale_x_continuous(trans='log10')+
      ylim(c(0,70)) +
      ggtitle('Average path') + 
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())
    
    g3 =  
    ggplot(d,aes(x = motileCount, y = VSL)) + 
      geom_point(pch = 21, aes(col = Morph)) + 
      stat_smooth(method = "lm", aes(col = Morph),size = 1, se = FALSE)+
      stat_cor(method="pearson",size = 2) +
      scale_x_continuous(trans='log10')+
      ylim(c(0,70)) + 
      ggtitle('Straight line') +
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
    ggsave(here::here('Output/Motility-corlogSpermN.png'),cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =8*1.5, units = 'cm')    
#' #### Velocity ~ date
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
        axis.title.x = element_text(color = "white"), 
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
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())

    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "first"))
    
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave(here::here('Output/Motility-corDate.png'),cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =8*1.5, units = 'cm')  
#' #### Velocity ~ age
    g1 = 
    ggplot(d[date == 'June'],aes(x = age, y = VCL)) + 
      geom_point(pch = 21, aes(col = Morph)) + 
      stat_smooth(method = "lm", aes(col = Morph), formula = y ~ poly(x,2), size = 1, se = FALSE)+
      stat_cor(method="pearson",size = 2) +
      ylab('Velocity μm/s') + 
      ylim(c(0,72)) + 
      ggtitle('Curvilinear') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
         axis.title.x = element_text(color = "white"), 
         plot.title = element_text(size=9, hjust = 0.5))
    
    g2 = 
    ggplot(d[date == 'June'],aes(x = age, y = VAP)) + 
      geom_point(pch = 21, aes(col = Morph)) + 
      stat_smooth(method = "lm", aes(col = Morph), formula = y ~ poly(x,2), size = 1, se = FALSE)+
      stat_cor(method="pearson",size = 2) +
      ylim(c(0,72)) +
      ggtitle('Average-path') + 
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())
    
    g3 =  
    ggplot(d[date == 'June'],aes(x = age, y = VSL)) + 
      geom_point(pch = 21, aes(col = Morph)) + 
      stat_smooth(method = "lm", aes(col = Morph), formula = y ~ poly(x,2), size = 1, se = FALSE) +
      stat_cor(method="pearson",size = 2) +
      ylim(c(0,72)) + 
      ggtitle('Straight-line') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())

    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "first"))
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave(here::here('Output/Motility-corage.png'),cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =8*1.5, units = 'cm')  
#' #### Velocity ~ issues - Morph
    g1 = 
    ggplot(d,aes(x = issues, y = VCL, col = Morph)) + 
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(), fill ='grey40', dotsize = 2)+
      #geom_boxplot(col = 'grey',alpha = 0.2,position = position_dodge())+
      #geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2,position = position_dodge())+
      stat_summary(fun=mean, geom="point", color="red", fill="red") +
      ylab('Velocity μm/s') + 
      ylim(c(0,72)) + 
      ggtitle('Curvilinear') +
      facet_wrap(~Morph, ncol = 1)  +
      guides(x =  guide_axis(angle = -45)) +
      theme_bw() +
      theme(legend.position = "none",
        axis.title.x = element_text(color = "white"), 
         plot.title = element_text(size=9, hjust = 0.5))
    
    g2 = 
    ggplot(d,aes(x = issues, y = VAP)) + 
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(), aes(col = Morph), fill ='grey40', dotsize = 2)+
      #geom_boxplot(col = 'grey',alpha = 0.2,position = position_dodge())+
      stat_summary(fun=mean, geom="point", color="red", fill="red") +
      ylim(c(0,72)) +
      ggtitle('Average-path') + 
      facet_wrap(~Morph, ncol = 1)  +
      guides(x =  guide_axis(angle = -45)) +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())
    
    g3 =  
    ggplot(d,aes(x = issues, y = VSL, col = Morph)) + 
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(), aes(col = Morph), fill ='grey40', dotsize = 2)+
      #geom_boxplot(col = 'grey',alpha = 0.2,position = position_dodge())+
      #geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2,position = position_dodge())+
      stat_summary(fun=mean, geom="point", color="red", fill="red") +
      ylim(c(0,72)) + 
      ggtitle('Straight-line') +
      facet_wrap(~Morph, ncol = 1)  +
      guides(x =  guide_axis(angle = -45)) +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())

    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "first"))
    
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave(here::here('Output/Motility-corIssues_Morph.png'),cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =8*1.5, units = 'cm')  
#' #### Velocity ~ issues
#+ fig.width=9, fig.height=2.5    
    g1 = 
    ggplot(d,aes(x = issues, y = VCL, col = issues)) + 
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(), fill ='grey40', dotsize = 1)+
      #geom_boxplot(col = 'grey',alpha = 0.2,position = position_dodge())+
      #geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2,position = position_dodge())+
      ylab('Velocity μm/s') + 
      ylim(c(0,72)) + 
      ggtitle('Curvilinear') +
      guides(x =  guide_axis(angle = -45)) +
      theme_bw() +
      theme(legend.position = "none",
        axis.title.x = element_text(color = "white"), 
         plot.title = element_text(size=9, hjust = 0.5))
    
    g2 = 
    ggplot(d,aes(x = issues, y = VAP, col = issues)) + 
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(), fill ='grey40', dotsize = 1)+
      #geom_boxplot(col = 'grey',alpha = 0.2,position = position_dodge())+
      ylim(c(0,72)) +
      ggtitle('Average-path') + 
      guides(x =  guide_axis(angle = -45)) +
      theme_bw() +
      xlab('Issues during sperm tracking')+
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())
    
    g3 =  
    ggplot(d,aes(x = issues, y = VSL, col = issues)) + 
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(),  fill ='grey40', dotsize = 1)+
      #geom_boxplot(col = 'grey',alpha = 0.2,position = position_dodge())+
      #geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2,position = position_dodge())+
      ylim(c(0,72)) + 
      ggtitle('Straight-line') +
      guides(x =  guide_axis(angle = -45)) +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())

    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "first"))
    
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave(here::here('Output/Motility-corIssues.png'),cbind(gg1,gg2,gg3, size = "first"), width = 18, height =5, units = 'cm')  

#' #### N ~ date
#+ fig.width=5, fig.height=2, dpi = 100       
    ggplot(d,aes(x = date, y = motileCount, col = Morph)) +  
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(), aes(col = Morph), fill ='grey40', dotsize = 1)+
      geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2) + 
      stat_summary(fun=mean, geom="point", color="red", fill="red") +
      scale_colour_viridis(discrete=TRUE)+
      theme_bw() 

    ggsave(here::here('Output/Motility-corSpermN_Date.png'), width = 13, height =6, units = 'cm')  
#' #### N ~ age
    g1 = 
    ggplot(d[date == 'May' & Morph!='Zebra finch'], aes(x = age, y = motileCount)) + 
      geom_point(pch = 21, aes(col = Morph)) + 
      stat_smooth(method = "lm", aes(col = Morph), size = 1) +
      stat_cor(method="pearson",size = 2) +
      ylab('N sperm cell tracked') + 
      coord_cartesian(xlim = c(1,13), ylim = c(0,600)) + 
      ggtitle('May') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
         axis.title.x = element_text(color = "white"),
         plot.title = element_text(size=9, hjust = 0.5))
    
    g2 = 
    ggplot(d[date == 'June' & Morph!='Zebra finch'],aes(x = age, y = motileCount)) + 
      geom_point(pch = 21, aes(col = Morph)) + 
      stat_smooth(method = "lm", aes(col = Morph), size = 1) +
      stat_cor(method="pearson",size = 2) +
      coord_cartesian(xlim = c(1,13), ylim = c(0,600)) + 
      ggtitle('June') + 
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank())
    
    g3 =  
    ggplot(d[Morph!='Zebra finch'],aes(x = age, y = motileCount)) + 
      geom_point(pch = 21, aes(col = Morph)) + 
      stat_smooth(method = "lm", aes(col = Morph), size = 1) +
      stat_cor(method="pearson",size = 2) +
      coord_cartesian(xlim = c(1,13), ylim = c(0,600)) + 
      ggtitle('Both') +
      facet_wrap(~Morph, ncol = 1)  +
      theme_bw() +
      theme(legend.position = "none",
          plot.title = element_text(size=9, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()
          )

    grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "first"))
    
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    ggsave(here::here('Output/Motility_corN-age.png'),cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =8*1.5, units = 'cm')  
  
#'***
#'***
  
#' ## Correlations between May & June
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
    ggsave(here::here('Output/Motility-cor_MayJune.png'),cbind(gg1,gg2,gg3, size = "first"), width = 7*1.5, height =7*1.5, units = 'cm')  
 
#' ## Repeatability
#+ results="hide", warning = FALSE, message = FALSE, cache = TRUE   
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
#+ fig.width=5, fig.height=2, dpi=100, results="markup"     
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
      ggsave(here::here('Output/Motility-repeatability.png'), width = 7*1.5, height =2.75*1.5, units = 'cm')

    # table  
      knitr::kable(r[,1:4])

#'***
#'***

#' ## Analyses of interest 
#' ### Compare May, June, May & June  
#+ fig.width=6, fig.height=6  
  # plot May   
    g1= 
    ggplot(d[date == 'May'],aes(x = Morph, y = VCL)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        coord_cartesian(ylim = c(20,72)) +
        ylab('VCL [μm/s]') +
        theme_bw() +
        guides(x =  guide_axis(angle = -45)) +
        ggtitle("May")+
        theme(
          legend.position = "none",
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          plot.title = element_text(size=9)
          )

    g2= 
    ggplot(d[date == 'May'],aes(x = Morph, y = VSL)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        coord_cartesian(ylim = c(8,47)) +
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
    ggplot(d[date == 'May'],aes(x = Morph, y = VAP)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        coord_cartesian(ylim = c(13,51)) +
        ylab('VAP [μm/s]') +
        theme_bw() +
        guides(x =  guide_axis(angle = -45)) +
        theme(
          legend.position = "none",
          axis.title.x = element_text(color = "white"), 
          #axis.text.x = element_blank(),
          plot.title = element_text(size=9)
          )   

    grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), size = "last"))
        #grid.arrange(g1,g2)
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    #ggsave(here::here('Output/Motility_Morp_boxplots_June.png'),rbind(gg1,gg2,gg3, size = "last"), width = 7*1.5, height =15*1.5, units = 'cm')  
    #summary(glht(m)) 
    #plot(allEffects(m))
  # plot June   
    dj = d[date == 'June']
    # dummy row to include zebra finch in the plot
    df = d[date == 'May' & Morph == 'Zebra finch'][1,]
    df$date = 'June'
    df$VAP = 100
    df$VCL = 100
    df$VSL = 100

    dj = rbind(dj,df)

    g4= 
    ggplot(dj,aes(x = Morph, y = VCL)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 0.85)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        coord_cartesian(ylim = c(20,72)) +
        ylab('VCL [μm/s]') +
        theme_bw() +
        guides(x =  guide_axis(angle = -45)) +
        ggtitle("June")+
        theme(
          legend.position = "none",
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),  
          plot.title = element_text(size=9)
          )

    g5= 
    ggplot(dj,aes(x = Morph, y = VSL)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 0.45)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        coord_cartesian(ylim = c(8,47)) +
        ylab('VSL [μm/s]') +
        theme_bw() +
        guides(x =  guide_axis(angle = -45)) +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),  
          plot.title = element_text(size=9)
          )
    
    g6= 
    ggplot(dj,aes(x = Morph, y = VAP)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 0.45)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        coord_cartesian(ylim = c(13,51)) +
        ylab('VAP [μm/s]') +
        theme_bw() +
        guides(x =  guide_axis(angle = -45)) +
        theme(
          legend.position = "none",
          #axis.title.x = element_blank(), 
          #axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),  
          plot.title = element_text(size=9)
          )   

    
    grid.draw(rbind(ggplotGrob(g4), ggplotGrob(g5), ggplotGrob(g6), size = "last"))
        #grid.arrange(g1,g2)
    gg4 <- ggplotGrob(g4)
    gg5 <- ggplotGrob(g5) 
    gg6 <- ggplotGrob(g6) 
    #ggsave(here::here('Output/Motility_Morp_boxplots_May.png'),rbind(gg4,gg5,gg6, size = "last"), width = 7*1.5, height =15*1.5, units = 'cm')  
  # plot ALL  together 
    g7= 
    ggplot(d,aes(x = Morph, y = VCL)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        coord_cartesian(ylim = c(20,72)) +
        ylab('VCL [μm/s]') +
        theme_bw() +
        guides(x =  guide_axis(angle = -45)) +
        ggtitle("May & June")+
        theme(
          legend.position = "none",
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(), 
          plot.title = element_text(size=9)
          )

    g8= 
    ggplot(d,aes(x = Morph, y = VSL)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        coord_cartesian(ylim = c(8,47)) +
        ylab('VSL [μm/s]') +
        theme_bw() +
        guides(x =  guide_axis(angle = -45)) +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          plot.title = element_text(size=9)
          )
    
    g9= 
    ggplot(d,aes(x = Morph, y = VAP)) + 
        geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
        geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
        stat_summary(fun=mean, geom="point", color="red", fill="red") +
        scale_fill_viridis(discrete=TRUE)+
        coord_cartesian(ylim = c(13,51)) +
        ylab('VAP [μm/s]') +
        theme_bw() +
        guides(x =  guide_axis(angle = -45)) +
        theme(
          legend.position = "none",
          axis.title.x = element_text(color = "white"), 
          #axis.title.y = element_text(color = "white"), 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          #axis.text.x = element_blank(),
          plot.title = element_text(size=9)
          )   

    grid.draw(rbind(ggplotGrob(g7), ggplotGrob(g8), ggplotGrob(g9), size = "last"))
        #grid.arrange(g1,g2)
    gg7 <- ggplotGrob(g7)
    gg8 <- ggplotGrob(g8) 
    gg9 <- ggplotGrob(g9) 
    #ggsave(here::here('Output/Motility_Morp_boxplots.png'),rbind(gg1,gg2,gg3, size = "last"), width = 7*1.5, height =15*1.5, units = 'cm')
  # combine
    may = rbind(gg1,gg2,gg3, size = "last")
    june = rbind(gg4,gg5,gg6, size = "last")
    all = rbind(gg7,gg8,gg9, size = "last")

    grid.draw(cbind(may,june,all, size = "first"))

    ggsave(here::here('Output/Motility_Morp_boxplots_MayJuneAll.png'),cbind(may,june,all, size = "last"), width = 10*1.5, height =10*1.5, units = 'cm')

#' ### Compare May & June
#+ fig.width=4, fig.height=5, dpi=100
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
    ggsave(here::here('Output/Motility_Morp_boxplots_date.png'),rbind(gg1,gg2,gg3, size = "last"), width = 7*1.5, height =15*1.5, units = 'cm') 
    
#' ### Model outputs
#' #### for May only
#' Curvilinear      
      m = lm(VCL ~ log(motileCount) + Morph, d[Morph!='Zebra finch' & date == 'May']) 
      summary(m)
#' Average path      
      m = lm(VAP ~ log(motileCount) + Morph, d[Morph!='Zebra finch' & date == 'May']) 
      summary(m)
#' Straight line     
      m = lm(VSL ~ log(motileCount) + Morph, d[Morph!='Zebra finch' & date == 'May']) 
      summary(m)
      #summary(glht(m)) 
      #plot(allEffects(m))

#' #### for June only
#' Curvilinear 
      m = lm(VCL ~ log(motileCount) + Morph, d[Morph!='Zebra finch' & date == 'June']) 
      summary(m)
#' Average path  
      m = lm(VAP ~ log(motileCount) + Morph, d[Morph!='Zebra finch' & date == 'June']) 
      summary(m)
#' Straight line    
      m = lm(VSL ~ log(motileCount) + Morph, d[Morph!='Zebra finch' & date == 'June']) 
      summary(m)
      
#' #### for May & June 
#' Curvilinear 
      #d$Morph_F = as.factor(d$Morph)
      m = lmer(VCL ~ log(motileCount) + date + age + Morph + (1|ID), d[Morph!='Zebra finch']) 
      summary(m)
#' Average path  
      m = lmer(VAP ~ log(motileCount) + date + age + Morph + (1|ID), d[Morph!='Zebra finch'])
      summary(m) 
#' Straight line  
      m = lmer(VSL ~ log(motileCount) + date + age + Morph + (1|ID), d[Morph!='Zebra finch']) 
      summary(m)
      #m = lm(VSL ~ log(motileCount) + date + age + Morph, d[Morph!='Zebra finch']) 
      #summary(glht(m)) 
      #summary(m)
      #plot(allEffects(m))

#'***
#'***

#' ## Test Morph-blind
#+ eval = FALSE  
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

# END    