#' ---
#' title: "Morphology of ruff sperm"
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
    require(grid)
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
  # DATA 
    # composite measures
      x = fread(here::here('Data/DAT_morpho.csv'))
      setnames(x,old = 'pic', new = 'sperm_ID')
      x[, sample_ID:=as.character(sample_ID)]

        #  mneasurements from minip vs rest are the same - so no need to control for
        #ggplot(x[part == 'Tail'], aes(x = manip, y = Pixels)) + geom_boxplot()
      
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
      a = b[, list(mean(Length_µm), mean(VAP),mean(VSL), mean(VCL), mean(motileCount)), by = list(bird_ID, Morph, age, part)]
      setnames(a, old = c('V1', 'V2','V3','V4','V5'), new = c('Length_avg', 'VAP', 'VSL', 'VCL', 'motileCount'))
      a1 = data.table(bird_ID = unique(a$bird_ID), blind = c(rep(c('Independent','Satellite', 'Faeder'), floor(length(unique(a$bird_ID))/3)), 'Independent','Satellite'))
      a = merge(a,a1, all.x = TRUE)

      aw = reshape(a, idvar = c('bird_ID','Morph', 'blind','age','VAP','VSL','VCL', 'motileCount'), timevar = "part", direction = "wide")
      names(aw) = c('bird_ID','Morph', 'VAP','VSL','VCL', 'motileCount','blind','age',as.character(unique(a$part)))
      aw[, Midpiece_rel := Midpiece/Total]
      aw[, Flagellum_rel := Flagellum/Total]

    # morph as factor
      b[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 
      bw[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 
      a[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 
      aw[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 

    # dataset for relative measurements
      bt = b[part == 'Total',.(Morph, bird_ID, sample_ID, sperm_ID, part, Length_µm, VAP, VSL, VCL, motileCount)]
      
      b1 = b[part%in%c('Midpiece'),.(Morph, bird_ID, sample_ID, sperm_ID, part, Length_µm, VAP, VSL, VCL, motileCount)]
      b1[, Length_rel := Length_µm/bt$Length_µm]
     
      b2 = b[part%in%c('Flagellum'),.(Morph, bird_ID, sample_ID, sperm_ID, part, Length_µm,  VAP, VSL, VCL, motileCount)]
      b2[, Length_rel := Length_µm/bt$Length_µm]
      
      br = rbind(b1,b2)
      br$Length_µm = NULL

      at = a[part == 'Total',.(Morph, bird_ID,  part, Length_avg, VAP, VSL, VCL, motileCount)]
      
      a1 = a[part%in%c('Midpiece'),.(Morph, bird_ID, part, Length_avg, VAP, VSL, VCL, motileCount)]
      a1[, Length_rel := Length_avg/at$Length_avg]
     
      a2 = a[part%in%c('Flagellum'),.(Morph, bird_ID, part, Length_avg, VAP, VSL, VCL, motileCount)]
      a2[, Length_rel := Length_avg/at$Length_avg]
      
      ar = rbind(a1,a2)
      ar$Length_avg = NULL
    
    # dataset for correlations
      h = b[part == 'Head']
      setnames(h, old = "Length_µm", new="Head_µm")
      h[, Acrosome_µm := b[part == 'Acrosome',.(Length_µm)]]
      h[, Nucleus_µm := b[part == 'Nucleus',.(Length_µm)]]
      h[, Midpiece_µm := b[part == 'Midpiece',.(Length_µm)]]
      h[, Tail_µm := b[part == 'Tail',.(Length_µm)]]
      h[, Flagellum_µm := b[part == 'Flagellum',.(Length_µm)]]
      h[, Total_µm := b[part == 'Total',.(Length_µm)]]
      h[, Midpiece_rel := br[part == 'Midpiece',.(Length_rel)]]
      h[, Flagellum_rel := br[part == 'Flagellum',.(Length_rel)]]

      ha = a[part == 'Head']
      setnames(ha, old = "Length_avg", new="Head_µm")
      ha[, Acrosome_µm := a[part == 'Acrosome',.(Length_avg)]]
      ha[, Nucleus_µm := a[part == 'Nucleus',.(Length_avg)]]
      ha[, Midpiece_µm := a[part == 'Midpiece',.(Length_avg)]]
      ha[, Tail_µm := a[part == 'Tail',.(Length_avg)]]
      ha[, Flagellum_µm := a[part == 'Flagellum',.(Length_avg)]]
      ha[, Total_µm := a[part == 'Total',.(Length_avg)]]
      ha[, Midpiece_rel := ar[part == 'Midpiece',.(Length_rel)]]
      ha[, Flagellum_rel := ar[part == 'Flagellum',.(Length_rel)]]
  
#' ## Exploration
 
#' ## Correlations
#+ cor, fig.width=10, fig.height = 6
  ga1 = 
  ggplot(aw, aes(x = Acrosome, y = VAP, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      annotate(geom="text", x=3, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      annotate(geom="text", x=3, y=13, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      annotate(geom="text", x=3, y=11, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      ggtitle('Male means') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  ga2 = 
  ggplot(aw, aes(x = Acrosome, y = VSL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )    
  ga3 = 
  ggplot(aw, aes(x = Acrosome, y = VCL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 

  gn1 = 
  ggplot(aw, aes(x = Nucleus, y = VAP, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
             #scale_colour_viridis(discrete=TRUE)+
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  gn2 = 
  ggplot(aw, aes(x = Nucleus, y = VSL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )    
  gn3 = 
  ggplot(aw, aes(x = Nucleus, y = VCL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 

  gm1 = 
  ggplot(aw, aes(x = Midpiece, y = VAP, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
            #scale_colour_viridis(discrete=TRUE)+
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  gm2 = 
  ggplot(aw, aes(x = Midpiece, y = VSL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )    
  gm3 = 
  ggplot(aw, aes(x = Midpiece, y = VCL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  
  gt1 = 
  ggplot(aw, aes(x = Tail, y = VAP, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
            #scale_colour_viridis(discrete=TRUE)+
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  gt2 = 
  ggplot(aw, aes(x = Tail, y = VSL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )    
  gt3 = 
  ggplot(aw, aes(x = Tail, y = VCL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 

  gh1 = 
  ggplot(aw, aes(x = Head, y = VAP, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  gh2 = 
  ggplot(aw, aes(x = Head, y = VSL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        plot.title = element_text(size=9)
        )    
  gh3 = 
  ggplot(aw, aes(x = Head, y = VCL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 

  gf1 = 
  ggplot(aw, aes(x = Flagellum, y = VAP, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  gf2 = 
  ggplot(aw, aes(x = Flagellum, y = VSL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )    
  gf3 = 
  ggplot(aw, aes(x = Flagellum, y = VCL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 

  go1 = 
  ggplot(aw, aes(x = Total, y = VAP, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
            #scale_colour_viridis(discrete=TRUE)+
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  go2 = 
  ggplot(aw, aes(x = Total, y = VSL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )    
  go3 = 
  ggplot(aw, aes(x = Total, y = VCL, col = Morph, fill = Morph)) +
      geom_point(pch = 21)+
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 

  grid.draw(cbind(
    rbind(ggplotGrob(ga1), ggplotGrob(ga2), ggplotGrob(ga3), size = "first"),
    rbind(ggplotGrob(gn1), ggplotGrob(gn2), ggplotGrob(gn3), size = "first"),
    rbind(ggplotGrob(gm1), ggplotGrob(gm2), ggplotGrob(gm3), size = "first"),
    rbind(ggplotGrob(gt1), ggplotGrob(gt2), ggplotGrob(gt3), size = "first"),
    rbind(ggplotGrob(gh1), ggplotGrob(gh2), ggplotGrob(gh3), size = "first"),
    rbind(ggplotGrob(gf1), ggplotGrob(gf2), ggplotGrob(gf3), size = "first"),
    rbind(ggplotGrob(go1), ggplotGrob(go2), ggplotGrob(go3), size = "first"),
    size = "first")
    )
  
  ggsave(here::here('Output/morpho_motil.png'),cbind(
          rbind(ggplotGrob(ga1), ggplotGrob(ga2), ggplotGrob(ga3), size = "first"),
          rbind(ggplotGrob(gn1), ggplotGrob(gn2), ggplotGrob(gn3), size = "first"),
          rbind(ggplotGrob(gm1), ggplotGrob(gm2), ggplotGrob(gm3), size = "first"),
          rbind(ggplotGrob(gt1), ggplotGrob(gt2), ggplotGrob(gt3), size = "first"),
          rbind(ggplotGrob(gh1), ggplotGrob(gh2), ggplotGrob(gh3), size = "first"),
          rbind(ggplotGrob(gf1), ggplotGrob(gf2), ggplotGrob(gf3), size = "first"),
          rbind(ggplotGrob(go1), ggplotGrob(go2), ggplotGrob(go3), size = "first"),
          size = "first"), 
          width = 18, height =8, units = 'cm') 
#+ cor_l, fig.width=10, fig.height = 6
  ga1 = 
  ggplot(aw, aes(x = Acrosome, y = VAP, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      annotate(geom="text", x=3, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      annotate(geom="text", x=3, y=13, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      annotate(geom="text", x=3, y=11, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      ggtitle('Male means') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  ga2 = 
  ggplot(aw, aes(x = Acrosome, y = VSL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )    
  ga3 = 
  ggplot(aw, aes(x = Acrosome, y = VCL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 

  gn1 = 
  ggplot(aw, aes(x = Nucleus, y = VAP, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
             #scale_colour_viridis(discrete=TRUE)+
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  gn2 = 
  ggplot(aw, aes(x = Nucleus, y = VSL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )    
  gn3 = 
  ggplot(aw, aes(x = Nucleus, y = VCL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 

  gm1 = 
  ggplot(aw, aes(x = Midpiece, y = VAP, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
            #scale_colour_viridis(discrete=TRUE)+
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  gm2 = 
  ggplot(aw, aes(x = Midpiece, y = VSL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )    
  gm3 = 
  ggplot(aw, aes(x = Midpiece, y = VCL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  
  gt1 = 
  ggplot(aw, aes(x = Tail, y = VAP, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
            #scale_colour_viridis(discrete=TRUE)+
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  gt2 = 
  ggplot(aw, aes(x = Tail, y = VSL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )    
  gt3 = 
  ggplot(aw, aes(x = Tail, y = VCL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 

  gh1 = 
  ggplot(aw, aes(x = Head, y = VAP, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  gh2 = 
  ggplot(aw, aes(x = Head, y = VSL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        plot.title = element_text(size=9)
        )    
  gh3 = 
  ggplot(aw, aes(x = Head, y = VCL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 

  gf1 = 
  ggplot(aw, aes(x = Flagellum, y = VAP, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  gf2 = 
  ggplot(aw, aes(x = Flagellum, y = VSL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )    
  gf3 = 
  ggplot(aw, aes(x = Flagellum, y = VCL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 

  go1 = 
  ggplot(aw, aes(x = Total, y = VAP, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
            #scale_colour_viridis(discrete=TRUE)+
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  go2 = 
  ggplot(aw, aes(x = Total, y = VSL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )    
  go3 = 
  ggplot(aw, aes(x = Total, y = VCL, col = Morph, fill = Morph)) +
      geom_smooth(method = 'rlm') +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=28, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=28, y=14, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=28, y=13, label="Faeder", size = 3, col = colors[3], hjust = 0) +
            #scale_colour_viridis(discrete=TRUE)+
      #ggtitle('Single measurements') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 

  grid.draw(cbind(
    rbind(ggplotGrob(ga1), ggplotGrob(ga2), ggplotGrob(ga3), size = "first"),
    rbind(ggplotGrob(gn1), ggplotGrob(gn2), ggplotGrob(gn3), size = "first"),
    rbind(ggplotGrob(gm1), ggplotGrob(gm2), ggplotGrob(gm3), size = "first"),
    rbind(ggplotGrob(gt1), ggplotGrob(gt2), ggplotGrob(gt3), size = "first"),
    rbind(ggplotGrob(gh1), ggplotGrob(gh2), ggplotGrob(gh3), size = "first"),
    rbind(ggplotGrob(gf1), ggplotGrob(gf2), ggplotGrob(gf3), size = "first"),
    rbind(ggplotGrob(go1), ggplotGrob(go2), ggplotGrob(go3), size = "first"),
    size = "first")
    )
  
     ggsave(here::here('Output/morpho_motil_l.png'),cbind(
          rbind(ggplotGrob(ga1), ggplotGrob(ga2), ggplotGrob(ga3), size = "first"),
          rbind(ggplotGrob(gn1), ggplotGrob(gn2), ggplotGrob(gn3), size = "first"),
          rbind(ggplotGrob(gm1), ggplotGrob(gm2), ggplotGrob(gm3), size = "first"),
          rbind(ggplotGrob(gt1), ggplotGrob(gt2), ggplotGrob(gt3), size = "first"),
          rbind(ggplotGrob(gh1), ggplotGrob(gh2), ggplotGrob(gh3), size = "first"),
          rbind(ggplotGrob(gf1), ggplotGrob(gf2), ggplotGrob(gf3), size = "first"),
          rbind(ggplotGrob(go1), ggplotGrob(go2), ggplotGrob(go3), size = "first"),
          size = "first"), 
          width = 18, height =8, units = 'cm') 

#' ## Model outcomes
#+ pred, results = "hide" 
      l = list()
      lp =list()
      # VAP
         for(i in unique(a$part)){
            #i ='Nucleus'
            ai = a[part == i]
            m = lm(scale(VAP) ~ Morph*scale(Length_avg), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            l[[paste('VAP',i)]]=data.frame(mot = 'VAP', response=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VAP ~ Morph*Length_avg, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
              }
            newD=do.call(rbind,nd)

            # values to predict for
            X <- model.matrix(~ Morph*Length_avg,data=newD) # exactly the model which was used has to be specified here 
            newD$VAP <-(X%*%v) 
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                            predmatrix[predmatrix < 0] <- 0
                            newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                            newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                            #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
            newD$part=i
            newD$mot = 'VAP'
            setnames(newD, old = 'VAP', new = 'motility')
            lp[[paste(i,newD$mot[1])]] = data.table(newD)

            print(paste(i,newD$mot[1]))     
            }          
      # VSL
         for(i in unique(a$part)){
            #i ='Nucleus'
            ai = a[part == i]
            m = lm(scale(VSL) ~ Morph*scale(Length_avg), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            l[[paste('VSL',i)]]=data.frame(mot = 'VSL', response=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VSL ~ Morph*Length_avg, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
              }
            newD=do.call(rbind,nd)
            X <- model.matrix(~ Morph*Length_avg,data=newD) # exactly the model which was used has to be specified here 
            newD$VSL <-(X%*%v) 
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                            predmatrix[predmatrix < 0] <- 0
                            newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                            newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                            #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
            newD$part=i
            newD$mot = 'VSL'
            setnames(newD, old = 'VSL', new = 'motility')
            lp[[paste(i,newD$mot[1])]] = data.table(newD)

            print(paste(i,newD$mot[1]))     
            }          
      # VCL
         for(i in unique(a$part)){
            #i ='Nucleus'
            ai = a[part == i]
            m = lm(scale(VCL) ~ Morph*scale(Length_avg), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            l[[paste('VCL',i)]]=data.frame(mot = 'VCL', response=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VCL ~ Morph*Length_avg, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
              }
            newD=do.call(rbind,nd)
            X <- model.matrix(~ Morph*Length_avg,data=newD) # exactly the model which was used has to be specified here 
            newD$VCL <-(X%*%v) 
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                            predmatrix[predmatrix < 0] <- 0
                            newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                            newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                            #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
            newD$part=i
            newD$mot = 'VCL'
            setnames(newD, old = 'VCL', new = 'motility')
            lp[[paste(i,newD$mot[1])]] = data.table(newD)

            print(paste(i,newD$mot[1]))     
            }          
                 
      ll = data.table(do.call(rbind,l) ) 
      ll[, response := factor(response, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
      ll[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      ll[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      ll[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      ll[effect == 'scale(Length_avg)', effect := 'Length_µm for Independent']
      ll[effect == 'MorphSatellite:scale(Length_avg)', effect := 'Length_µm Satellite\n(relative to Independent)']
      ll[effect == 'MorphFaeder:scale(Length_avg)', effect := 'Length_µm Faeder\n(relative to Independent)']

      ll[, effect := factor(effect, levels=c('Length_µm Faeder\n(relative to Independent)','Length_µm Satellite\n(relative to Independent)','Length_µm for Independent',"Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llp = data.table(do.call(rbind,lp) ) 
      llp[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
      llp[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
#+ effects, fig.width =5, fig.height = 7           
        g = 
        ggplot(ll, aes(y = effect, x = estimate, col = response, shape = mot)) +
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = response), width = 0.1, position = position_dodge(width = 0.9) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.9)) +
          geom_vline(xintercept = 0, col = "red", lty =1)+
          #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          scale_shape(guide = guide_legend(reverse = TRUE)) + 
          #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
          labs(y = NULL ,x = "Standardized effect size") +
          #coord_fixed(ratio = 1.5)+
          #ylim(c(0,100))+
          #coord_flip()+
          theme_bw() +
          theme( #legend.position ="right",
                plot.title = element_text(size=7),
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
                axis.ticks.length = unit(1, "pt"),
                axis.text.x = element_text(colour="black", size = 7),
                axis.text.y=element_text(colour="black", size = 7),
                axis.title=element_text(size=9)
                )
        g
        ggsave(here::here('Output/morpho_motil_effectSizes_virid.png'),g, width = 10, height =14, units = 'cm')

#end        