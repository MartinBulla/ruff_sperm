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
    require(readxl)


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
      
      # for each bird 10 sperm measured from one of the two sampling occasions
        #bb = b[part == 'Tail']
        #bbx = data.table(table(bb$bird_ID,bb$month))
        #bbx[!N%in%c(0,10), unique(V1)]

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

    # add metadata 
      b = merge(b, x[,.(bird_ID, sample_ID, sperm_ID, part, manip)], all.x = TRUE, by=c('bird_ID', 'sample_ID', 'sperm_ID', 'part'))
      b[manip%in%"", manip:=NA]
      # adjust multiple samples from the same time to the one recorded for motility
        b[sample_ID==51, sample_ID:=52]
        b[sample_ID==3, sample_ID:=4]
      s = data.table(read_excel(here::here('Data/sampling_2021_cleaned.xlsx'), sheet = 1))
      s = s[!is.na(recording)]
      m = data.table(read_excel(here::here('Data/ruff_males_Seewiesen.xlsx'), sheet = 1))#, range = "A1:G161"))
      m[, hatch_year:=as.numeric(substr(Hatch_date,1,4)) ]
      m[, age := 2021-hatch_year]

      s = merge(s,m[,.(Ind_ID, Morph, age)], by.x = 'DB_ID', by.y = 'Ind_ID', all.x = TRUE)
      # adjust IDs (where missed typed) for smooth merging
        b[bird_ID=='A027121649', bird_ID:='AO27121649']
        b[bird_ID=='A03999183', bird_ID:='AO3999183']
        b[bird_ID=='A0414173NL', bird_ID:='AO414-17-3NL']
        b[bird_ID=='A0414175', bird_ID:='AO414-17-5']
        b[bird_ID=='A059381830', bird_ID:='AO59381830']
        b[bird_ID=='A07422181', bird_ID:='AO7422181NL']
        b[bird_ID=='A079561654', bird_ID:='AO79561654']
        b[bird_ID=='A079561656', bird_ID:='AO79561656']
        b[bird_ID=='A079561718', bird_ID:='AO7956-17-18']
        b[bird_ID=='AIFA-016507', bird_ID:='AIFAO16507']
        b[bird_ID=='AIFA01511', bird_ID:='AIFAO15-11']
        b[bird_ID=='A079561718', bird_ID:='AO79561718']
        b[bird_ID=='G20005', bird_ID:='G200055']
        b[bird_ID=='cz005', bird_ID:='CZ005']
        b[bird_ID=='AO414-17-3NL', bird_ID:='AO414-17-3-NL']

      b = merge(b, s[,.(sample_ID, bird_ID, Morph, age, datetime, type, sperm, recording, rec_measured, month)], by.x = c('sample_ID', 'bird_ID'), by.y = c('sample_ID', 'bird_ID'), all.x = TRUE)
      
      # fwrite(ss[!duplicated(bird_ID), .(DB_ID, bird_ID)], file = 'Data/used_males_for_Clemens.csv')
      #b = (b[!duplicated(b)]) # shouldn't be necessary,
    
    # add motility on to individual morpho measurements (only one motility value per individual)
      d = data.table(read_excel(here::here('Data/motility.xlsx'), sheet = 1))
      d[, motileCount_ln := log(motileCount)]
      setnames(d, old=c('date','ID'), new = c('month','bird_ID'))

      # for bird 1339 motility measured from May sample, but sperm from June, so meta data changed to May sample 
         b[sample_ID==183, month:='May']
         b[sample_ID==183, datetime:=as.POSIXct('2021-05-04 12:33')]
         b[sample_ID==183, sample_ID:=95]
      
      b = merge(b, d, by = c('bird_ID', 'month'), all.x = TRUE)

      #b[is.na(Morph), unique(bird_ID)]
     
      b[Morph == 'F', Morph := 'Faeder']
      b[Morph == 'I', Morph := 'Independent']
      b[Morph == 'S', Morph := 'Satellite']
      b[is.na(issues), issues := 'zero']

      b[, part := factor(part, levels = c('Acrosome', 'Nucleus', 'Midpiece', 'Tail','Head','Flagellum', 'Total'))]
      b[part %in% c('Acrosome', 'Nucleus', 'Midpiece', 'Tail'), measure := 'original']
      b[part %in% c('Head','Flagellum', 'Total'), measure := 'composite']

      #nrow(d) # N 139
   
    # add metadata to motility dataset
      # check for multiple samples per individuals at a single sampling period & keep only one for merging
        ajm = data.table(table(s$bird_ID,s$month))
        #s[bird_ID%in%ajm[N>1, unique(V1)]]
        ss = s[rec_measured %in% 'yes']    
      d = merge(d, ss[,.(sample_ID, bird_ID, Morph, age, datetime, type, sperm, recording, month)], by = c('month','bird_ID'), all.x = TRUE)
     
      d[is.na(Morph), Morph := 'Zebra finch']
      d[Morph == 'F', Morph := 'Faeder']
      d[Morph == 'I', Morph := 'Independent']
      d[Morph == 'S', Morph := 'Satellite']
      d[is.na(issues), issues := 'zero']

    # prepare for correlations and repeatability
      bw = reshape(b[,.(bird_ID,Morph, age, datetime, month, sample_ID, sperm_ID, VAP,VSL,VCL, motileCount, part, Length_µm)], idvar = c('bird_ID','Morph','age','datetime', 'month', 'sample_ID', 'sperm_ID','VAP','VSL','VCL', 'motileCount'), timevar = 'part', direction = "wide")  
      setnames(bw,old = c('Length_µm.Acrosome', 'Length_µm.Nucleus','Length_µm.Head','Length_µm.Midpiece','Length_µm.Tail','Length_µm.Flagellum', 'Length_µm.Total'), new = c('Acrosome', 'Nucleus', 'Head','Midpiece', 'Tail','Flagellum','Total'))
      # add relative measures
        bw[, Midpiece_rel := Midpiece/Total]
        bw[, Flagellum_rel := Flagellum/Total]

     dw = reshape(d[bird_ID%in%d[duplicated(bird_ID), bird_ID],.(bird_ID,species, Morph, age, month, VAP,VSL,VCL, motileCount, motileCount_ln)], idvar = c('bird_ID','species','Morph','age'), timevar = 'month', direction = "wide")  
   
    # mean/male dataset
      a = b[, list(mean(Length_µm), mean(VAP),mean(VSL), mean(VCL), mean(motileCount)), by = list(month, bird_ID, Morph, age, part)]
       setnames(a, old = c('V1', 'V2','V3','V4','V5'), new = c('Length_avg', 'VAP', 'VSL', 'VCL', 'motileCount'))
      a1 = data.table(bird_ID = unique(a$bird_ID), blind = c(rep(c('Independent','Satellite', 'Faeder'), floor(length(unique(a$bird_ID))/3)), 'Independent','Satellite'))
      a = merge(a,a1, all.x = TRUE)

      aw = reshape(a, idvar = c('month','bird_ID','Morph', 'blind','age','VAP','VSL','VCL', 'motileCount'), timevar = "part", direction = "wide")
      names(aw) = c('bird_ID','month','Morph', 'age','VAP','VSL','VCL', 'motileCount', 'blind', as.character(unique(a$part)))
      # add relative measures
      aw[, Midpiece_rel := Midpiece/Total]
      aw[, Flagellum_rel := Flagellum/Total]
    # add inbreeding & HL    
      # prepare
        n = data.table(read_excel(here::here("Data/Ruffs2020_inbreeding.xlsx"), sheet = "allBins", na = "NA"))
        n = n[Sex == 1]
        n = n[!is.na(OriginalRing)]
        n[, c('Ruff8a', 'Ruff8b') := NULL]
        n[ Cme9a==Cme9b, Cme9 := 1]
        n[ Cme9a!=Cme9b, Cme9 := 0]
        n[ Ppu47a==Ppu47b, Ppu47 := 1]
        n[ Ppu47a!=Ppu47b, Ppu47 := 0]
        n[ Ruff1a==Ruff1b, Ruff1 := 1]
        n[ Ruff1a!=Ruff1b, Ruff1 := 0]
        n[ Ruff12a==Ruff12b, Ruff12 := 1]
        n[ Ruff12a!=Ruff12b, Ruff12 := 0]
        n[ Ruff50a==Ruff50b, Ruff50 := 1]
        n[ Ruff50a!=Ruff50b, Ruff50 := 0]
        n[ Ppu20a==Ppu20b, Ppu20 := 1]
        n[ Ppu20a!=Ppu20b, Ppu20 := 0]
        n[ Ppu22a==Ppu22b, Ppu22 := 1]
        n[ Ppu22a!=Ppu22b, Ppu22 := 0]
        n[ Ppu24a==Ppu24b, Ppu24 := 1]
        n[ Ppu24a!=Ppu24b, Ppu24 := 0]
        n[ Ppu25a==Ppu25b, Ppu25 := 1]
        n[ Ppu25a!=Ppu25b, Ppu25 := 0]
        n[ Ppu28a==Ppu28b, Ppu28 := 1]
        n[ Ppu28a!=Ppu28b, Ppu28 := 0]
        n[ Ppu48a==Ppu24b, Ppu48 := 1]
        n[ Ppu48a!=Ppu24b, Ppu48 := 0]
        n[ Cme2a==Cme2b, Cme2 := 1]
        n[ Cme2a!=Cme2b, Cme2 := 0]   
        n[ Cme6a==Cme6b, Cme6 := 1]
        n[ Cme6a!=Cme6b, Cme6 := 0]
        n[ Ruff6a==Ruff6b, Ruff6 := 1]
        n[ Ruff6a!=Ruff6b, Ruff6 := 0]
        n[ Ppu9a==Ppu9b, Ppu9 := 1]
        n[ Ppu9a!=Ppu9b, Ppu9 := 0]
        n[ Ppu19a==Ppu19b, Ppu19 := 1]
        n[ Ppu19a!=Ppu19b, Ppu19 := 0]
        n[ Ppu31a==Ppu31b, Ppu31 := 1]
        n[ Ppu31a!=Ppu31b, Ppu31 := 0]
        n[ Ppu56a==Ppu56b, Ppu56 := 1]
        n[ Ppu56a!=Ppu56b, Ppu56 := 0]
        n[ Tgu06_Slatea==Tgu06_Slateb, Tgu06_Slate := 1]
        n[ Tgu06_Slatea!=Tgu06_Slateb, Tgu06_Slate := 0]
        n[ Phil2a==Phil2b, Phil2 := 1]
        n[ Phil2a!=Phil2b, Phil2 := 0]  
        
        nn = melt(n, id.vars = c("OriginalRing"),
                measure.vars = c("Cme9", "Ppu47", "Ruff1","Ruff12","Ruff50","Ppu20","Ppu22","Ppu24","Ppu25","Ppu28","Ppu48","Cme2","Cme6","Ruff6","Ppu9","Ppu19","Ppu31",'Ppu56',"Tgu06_Slate","Phil2"),
                variable.name = "marker", value.name = "homo")
        nn = nn[!(OriginalRing %in% '7 -04 - 105' & marker %in% 'Ppu9') ] # remove uncertain assignment
        h = nn[, .(sum(homo)/.N, .N), by = OriginalRing]
        #h = nn[, .(sum(homo)/.N), by = OriginalRing]
        setnames(h, old = 'V1', new = 'inbreeding') 
           
        n = data.table(read_excel(here::here("Data/RuffsCK_2010-2019_RunData6.xlsx"), sheet = "allAlleles", na = "NA", col_type = 'text'))
        n = n[sex == 1]
        n = n[!is.na(OriginalRing)]
        n[, c('Ruff8a', 'Ruff8b') := NULL]
        n[ Cme9a==Cme9b, Cme9 := 1]
        n[ Cme9a!=Cme9b, Cme9 := 0]
        n[ Ppu47a==Ppu47b, Ppu47 := 1]
        n[ Ppu47a!=Ppu47b, Ppu47 := 0]
        n[ Ruff1a==Ruff1b, Ruff1 := 1]
        n[ Ruff1a!=Ruff1b, Ruff1 := 0]
        n[ Ruff12a==Ruff12b, Ruff12 := 1]
        n[ Ruff12a!=Ruff12b, Ruff12 := 0]
        n[ Ruff50a==Ruff50b, Ruff50 := 1]
        n[ Ruff50a!=Ruff50b, Ruff50 := 0]
        n[ Ppu20a==Ppu20b, Ppu20 := 1]
        n[ Ppu20a!=Ppu20b, Ppu20 := 0]
        n[ Ppu22a==Ppu22b, Ppu22 := 1]
        n[ Ppu22a!=Ppu22b, Ppu22 := 0]
        n[ Ppu24a==Ppu24b, Ppu24 := 1]
        n[ Ppu24a!=Ppu24b, Ppu24 := 0]
        n[ Ppu25a==Ppu25b, Ppu25 := 1]
        n[ Ppu25a!=Ppu25b, Ppu25 := 0]
        n[ Ppu28a==Ppu28b, Ppu28 := 1]
        n[ Ppu28a!=Ppu28b, Ppu28 := 0]
        n[ Ppu48a==Ppu24b, Ppu48 := 1]
        n[ Ppu48a!=Ppu24b, Ppu48 := 0]
        n[ Cme2a==Cme2b, Cme2 := 1]
        n[ Cme2a!=Cme2b, Cme2 := 0]   
        n[ Cme6a==Cme6b, Cme6 := 1]
        n[ Cme6a!=Cme6b, Cme6 := 0]
        n[ Ruff6a==Ruff6b, Ruff6 := 1]
        n[ Ruff6a!=Ruff6b, Ruff6 := 0]
        n[ Ppu9a==Ppu9b, Ppu9 := 1]
        n[ Ppu9a!=Ppu9b, Ppu9 := 0]
        n[ Ppu19a==Ppu19b, Ppu19 := 1]
        n[ Ppu19a!=Ppu19b, Ppu19 := 0]
        n[ Ppu31a==Ppu31b, Ppu31 := 1]
        n[ Ppu31a!=Ppu31b, Ppu31 := 0]
        n[ Ppu56a==Ppu56b, Ppu56 := 1]
        n[ Ppu56a!=Ppu56b, Ppu56 := 0]
        n[ Tgu06_Slatea==Tgu06_Slateb, Tgu06_Slate := 1]
        n[ Tgu06_Slatea!=Tgu06_Slateb, Tgu06_Slate := 0]
        n[ Phil2a==Phil2b, Phil2 := 1]
        n[ Phil2a!=Phil2b, Phil2 := 0]  
        
        nn = melt(n, id.vars = c("OriginalRing"),
                measure.vars = c("Cme9", "Ppu47", "Ruff1","Ruff12","Ruff50","Ppu20","Ppu22","Ppu24","Ppu25","Ppu28","Ppu48","Cme2","Cme6","Ruff6","Ppu9","Ppu19","Ppu31",'Ppu56',"Tgu06_Slate","Phil2"),
                variable.name = "marker", value.name = "homo")
        h2 = nn[, .(sum(homo)/.N, .N), by = OriginalRing]
        #h2 = nn[, .(sum(homo)/.N), by = OriginalRing]
        setnames(h2, old = 'V1', new = 'inb2014_19') 

        hh = merge(h,h2, by = 'OriginalRing', all.x = TRUE)
        hh[N.x!=N.y, inbreeding := inb2014_19]
        h = hh[,.(OriginalRing,inbreeding)]
        h2 = h2[!OriginalRing%in%h$OriginalRing]
        h2[, N := NULL]
        setnames(h2, old = 'inb2014_19', new = 'inbreeding')
        h = rbind(h,h2)

        h[ OriginalRing== '7 -04 - 105', OriginalRing := 704105]
        h[ OriginalRing== 'AIF AO - 15 - 11', OriginalRing := 'AIFAO15-11']
        h[ OriginalRing== 'A 8209 AIF AO 16 507', OriginalRing := 'AIFAO16507']
        h[ OriginalRing== 'AO 2712-16-49-B 6.5', OriginalRing := 'AO27121649']
        h[ OriginalRing== 'AO 3999-18-3-NL 6.5', OriginalRing := 'AO3999183']
        h[ OriginalRing== 'AO 414-17-3-NL 6.5', OriginalRing := 'AO414-17-3-NL']
        h[ OriginalRing== 'AO 414-17-5-NL 6.5', OriginalRing := 'AO414-17-5']
        h[ OriginalRing== 'AO 5938-18-30-NL 6.5', OriginalRing := 'AO59381830']
        h[ OriginalRing== 'AO 7422-18-1-NL 6.5', OriginalRing := 'AO7422181NL']
        h[ OriginalRing== 'AO 7422-18-1-NL 6.5', OriginalRing := 'AO7422181NL']
        h[ OriginalRing== 'AO 7956-17-18-NL 6.5', OriginalRing := 'AO7956-17-18']
        h[ OriginalRing== 'AO 7956-16-54-NL 6.5', OriginalRing := 'AO79561654']
        h[ OriginalRing== 'AO 7956-16-56-NL 6.5', OriginalRing := 'AO79561656']
        h[ OriginalRing== 'CZ 6.5 005', OriginalRing := 'CZ005']
        h[ OriginalRing== 'G 5.0 - 14 - 223', OriginalRing := 'G14223']
        h[ OriginalRing== 'G 5.0 - 17 - 267', OriginalRing := 'G17267']
        h[substr(OriginalRing, 1,3) == 'G20', OriginalRing := paste0('G20', substring(OriginalRing, 5,8))]
      
      # merge
        aw = merge(aw,h, by.x = 'bird_ID', by.y = 'OriginalRing')
        a = merge(a,h, by.x = 'bird_ID', by.y = 'OriginalRing')
        z = fread("Data/Dat_HL.txt")
        aw = merge(aw,z[,.(sampleid, HL)], by.x = 'bird_ID', by.y = 'sampleid')
        a = merge(a,z[,.(sampleid, HL)], by.x = 'bird_ID', by.y = 'sampleid')

    # morph as factor
      b[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 
      bw[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 
      a[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 
      aw[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder"))] 
      
      d[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder", "Zebra finch"))] 
      dw[, Morph := factor(Morph, levels=c("Independent", "Satellite", "Faeder", "Zebra finch"))] 
      #ggplot(dw, aes(x = VAP.May, y = VAP.June))+geom_point() + stat_smooth(method = 'lm') + geom_abline(intercept = 0, slope = 1, lty = 3, col = 'red')
    
    # CV dataset
        cv_ =  b[, cv(Length_µm), by = list(bird_ID, part, Morph)]
        cv_[ , Morph123 := as.numeric(Morph)]
        names(cv_) [4]='CV'
        cv_ = merge(cv_,h, by.x = 'bird_ID', by.y = 'OriginalRing')
        cv_ = merge(cv_,z[,.(sampleid, HL)], by.x = 'bird_ID', by.y = 'sampleid')

    # dataset for relative measurements
      bt = b[part == 'Total',.(Morph, bird_ID, sample_ID, sperm_ID, part, Length_µm, VAP, VSL, VCL, motileCount, motileCount_ln)]
      
      b1 = b[part%in%c('Midpiece'),.(Morph, bird_ID, sample_ID, sperm_ID, part, Length_µm, VAP, VSL, VCL, motileCount, motileCount_ln)]
      b1[, Length_rel := Length_µm/bt$Length_µm]
     
      b2 = b[part%in%c('Flagellum'),.(Morph, bird_ID, sample_ID, sperm_ID, part, Length_µm, VAP, VSL, VCL, motileCount, motileCount_ln)]
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
      ar = merge(ar,h, by.x = 'bird_ID', by.y = 'OriginalRing')
      ar = merge(ar,z[,.(sampleid, HL)], by.x = 'bird_ID', by.y = 'sampleid')
    
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

    # export for Mihai
      a$CV = cv_$CV[match(a$bird_ID, cv_$bird_ID)]
      #save(b,br,a,ar, file = 'Data/ruff_sperm_for_Mihai.RData')

#' ## Basic info
#' - for each male 10 sperm measured from one of the two sampling occasions
#' - motility for most individuals and for many in May and June
#' - for all morpho data we have a motility measurement that corresponds with the sampling date, except for one male, where motility measured in May, but morpho from June sample
#' - initial exploration of the motility is available from [here](https://rawcdn.githack
#'.com/MartinBulla/ruff_sperm/60918bbc69163e2820df7c712c33d0ad11ea6e61/Output/motility.html).
#' - for motility ~ morph analyses I, use
#'      1. June motility values and for 4 males without June, May values. Is this ok or better to use male averages? 
#'      2. bird_ID as random intercept - 50 males measured once 42 twice (June & May) - TO DO
#' - to decide: whether to control for pedigree, relatedness matrix based on microsatellite genotypes or nothing
#' - to do: Run models taking the recordings with issues out and check whether the model outcomes differ

#' ***
#' ## Exploration
#+ inbr_HL,fig.width=4, fig.height = 4
  ggplot(aw, aes(x = HL, y = inbreeding)) + geom_point() + geom_abline(slope = 1, col = 'red')
#+ inbr_morph, fig.width=4, fig.height = 4
  ggplot(aw, aes(x = Morph, y = inbreeding)) +
          geom_dotplot(binaxis = 'y', stackdir = 'center',
                        position = position_dodge(), col = 'grey', fill = 'lightgrey', dotsize = 0.75)+
          geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2) + 
          stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
          #stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
          scale_color_viridis(discrete=TRUE)+
          theme_bw() +
          theme(legend.position = "none",
            plot.title = element_text(size=9)
            )  
#+ inbr_part, fig.width=8, fig.height = 4
  ggplot(aw, aes(x = Morph, y = inbreeding)) +
          geom_dotplot(binaxis = 'y', stackdir = 'center',
                        position = position_dodge(), col = 'grey', fill = 'lightgrey', dotsize = 0.75)+
          geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2) + 
          stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
          #stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
          scale_color_viridis(discrete=TRUE)+
          theme_bw() +
          theme(legend.position = "none",
            plot.title = element_text(size=9)
            )  
#+ inbr_CV, fig.width=8, fig.height = 4 
 ggplot(cv_, aes(x = inbreeding, y = CV)) +
          stat_smooth(method = 'rlm') +
          geom_point(aes(col = Morph), alpha = 0.5) + 
          scale_color_viridis(discrete=TRUE)+
          facet_wrap(~part, scales = 'free_y', nrow = 2)+
          #labs(title = 'Predictions from sperm part specific models\non single sperm') +
          xlab('Inbreeding') +
          ylab('CV') +
          theme_bw() +
          guides(x =  guide_axis(angle = -45)) +
          theme(legend.position = "none",
            #axis.title.x = element_blank(), 
            #axis.text.x = element_blank(),
            plot.title = element_text(size=9)
            ) #panel.spacing.y = unit(0, "mm")) #,     
#+ inbr_velo fig.width=8, fig.height = 4     
  ga1 = 
    ggplot(aw, aes(x = inbreeding, y = VAP)) +
      geom_smooth(method = 'rlm') +
      geom_point(aes(col = Morph), alpha = 0.5) + 
      scale_color_viridis(discrete=TRUE)+
     # xlab('Inbreeding') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  
  ga2 = 
  ggplot(aw, aes(x = inbreeding, y = VSL)) +
      geom_smooth(method = 'rlm') +
      geom_point(aes(col = Morph), alpha = 0.5) + 
      scale_color_viridis(discrete=TRUE)+
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )    
  ga3 = 
  ggplot(aw, aes(x = inbreeding, y = VCL)) +
      geom_smooth(method = 'rlm') +
      geom_point(aes(col = Morph), alpha = 0.5) + 
      scale_color_viridis(discrete=TRUE)+
      xlab("Inbreeding") + 
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  grid.draw(
    rbind(ggplotGrob(ga1), ggplotGrob(ga2), ggplotGrob(ga3), size = "first")
    )

#+ HL_morph, fig.width=4, fig.height = 4
  ggplot(aw, aes(x = Morp, y = HL)) +
          geom_dotplot(binaxis = 'y', stackdir = 'center',
                        position = position_dodge(), col = 'grey', fill = 'lightgrey', dotsize = 0.75)+
          geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2) + 
          stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
          #stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
          scale_color_viridis(discrete=TRUE)+
          theme_bw() +
          theme(legend.position = "none",
            plot.title = element_text(size=9)
            )  
#+ HL_part, fig.width=8, fig.height = 4
  ggplot(a, aes(x = HL, y = Length_avg)) +
    stat_smooth(method = 'rlm') +
    geom_point(aes(col = Morph), alpha = 0.5) + 
    stat_cor(aes(label = ..r.label..),  label.x = 0.3, size = 2)+
    scale_color_viridis(discrete=TRUE)+
    facet_wrap(~part, scales = 'free_y', nrow = 2)+
    #stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
    scale_color_viridis(discrete=TRUE)+
    xlab('Homozygousity by locus') +
    ylab('Length [µm]') +
    theme_bw() +
    theme(legend.position = "none",
      plot.title = element_text(size=9)
      )  
#+ HL_CV, fig.width=8, fig.height = 4 
 ggplot(cv_, aes(x = HL, y = CV)) +
          stat_smooth(method = 'rlm') +
          geom_point(aes(col = Morph), alpha = 0.5) + 
          stat_cor(aes(label = ..r.label..),  label.x = 0.3, size = 2)+
          scale_color_viridis(discrete=TRUE)+
          facet_wrap(~part, scales = 'free_y', nrow = 2)+
          #labs(title = 'Predictions from sperm part specific models\non single sperm') +
          xlab('Homozygousity by locus') +
          ylab('CV') +
          theme_bw() +
          #guides(x =  guide_axis(angle = -45)) +
          theme(legend.position = "none",
            #axis.title.x = element_blank(), 
            #axis.text.x = element_blank(),
            plot.title = element_text(size=9)
            ) #panel.spacing.y = unit(0, "mm")) #,     
#+ HL_velo fig.width=8, fig.height = 4     
  ga1 = 
    ggplot(aw, aes(x = HL, y = VAP)) +
      geom_smooth(method = 'rlm') +
      geom_point(aes(col = Morph), alpha = 0.5) + 
      stat_cor(aes(label = ..r.label..),  label.x = 0.3, size = 2)+
      scale_color_viridis(discrete=TRUE)+
     # xlab('Inbreeding') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  
  ga2 = 
  ggplot(aw, aes(x = HL, y = VSL)) +
      geom_smooth(method = 'rlm') +
      geom_point(aes(col = Morph), alpha = 0.5) + 
      stat_cor(aes(label = ..r.label..),  label.x = 0.3, size = 2)+
      scale_color_viridis(discrete=TRUE)+
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )    
  ga3 = 
  ggplot(aw, aes(x = HL, y = VCL)) +
      geom_smooth(method = 'rlm') +
      geom_point(aes(col = Morph), alpha = 0.5) + 
      stat_cor(aes(label = ..r.label..),  label.x = 0.3, size = 2)+
      scale_color_viridis(discrete=TRUE)+
      xlab("Homozygousity by locus") + 
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 
  grid.draw(
    rbind(ggplotGrob(ga1), ggplotGrob(ga2), ggplotGrob(ga3), size = "first")
    )

#+ fig.width=8, fig.height = 10
  b[, order_ := mean(Length_µm), by = bird_ID]
  b_ = b[part =='Total']
  b_[,bird_ID := reorder(bird_ID, Length_µm, mean)]
  b[, bird_ID := factor(bird_ID, levels = levels(b_$bird_ID))]

  g = ggplot(b, aes(x = reorder(as.factor(bird_ID),Length_µm,mean), y = Length_µm)) + facet_wrap(~part, nrow = 7, scales = 'free') + #x = reorder(as.factor(bird_ID),Length_µm,mean)
    geom_boxplot(aes(col = Morph)) + 
    labs(subtitle = 'Distribution within & across males (ordered by total length)')+
    #geom_dotplot(binaxis = 'y', stackdir = 'center',
     #            position = position_dodge(), col = 'red', fill ='red')+
    scale_colour_manual(values = colors) + 
    #scale_color_viridis(discrete=TRUE)+
    xlab('Male ID') +
    theme_bw() +
    theme(axis.text.x = element_blank())
    #$, legend.position = "none")
  g
  ggsave(here::here('Output/morpho_within_male_boxplots_ordered.png'), g, width = 20, height = 15, units = 'cm')
  
#+ cor_parts, figures-side, fig.show="hold", out.width="50%", fig.dim = c(8, 8)
  #fig.width=5, fig.height = 5
  chart.Correlation(bw[, c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel')], histogram=TRUE, pch=19)
  mtext("Single sperm", side=3, line=3)
  #dev.copy(png,'Output/VD_corr_single.png')
  #dev.off()

  # cor_parts_avg, fig.width=5, fig.height=5  
  chart.Correlation(aw[, c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total','Midpiece_rel', 'Flagellum_rel')], histogram=TRUE, pch=19)
  mtext("Male averages", side=3, line=3)
  #dev.copy(png,'Output/VD_corr_avg.png')
  #dev.off()

#' 
#' <span style="color: red;">!!! Nucleus strongly predict head and tail total sperm length. !!!  </span>  
#' <br>
#' 
#' ***
#+ cor_parts_I, figures-side, fig.show="hold", out.width="33%", fig.dim = c(6, 6)
  #fig.width=5, fig.height = 5  
  chart.Correlation(bw[Morph == 'Independent', c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel')], histogram=TRUE, pch=19)
  mtext("Independent", side=3, line=3)
  #dev.copy(png,'Output/morpho_corr_single_I.png')
  #dev.off()  
  #cor_parts_S, fig.width=5, fig.height = 5    
  chart.Correlation(bw[Morph == 'Satellite', c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel')], histogram=TRUE, pch=19)
  mtext("Satellite", side=3, line=3)
  #dev.copy(png,'Output/morpho_corr_single_S.png')
  #dev.off()
  #cor_parts_F, fig.width=5, fig.height = 5    
  chart.Correlation(bw[Morph == 'Faeder', c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel')], histogram=TRUE, pch=19)
  mtext("Feader", side=3, line=3)
  #dev.copy(png,'Output/morpho_corr_single_F.png')
  #dev.off()
#' 
#'   
#' <span style="color: red;">!!! Correlations look similar across morphs, except for midpiece in faeders where it more strongly correlates with tail and hence also total sperm length  !!!</span>    
#' 
#' ***
#' ## Repeatability 
#+ estimate, results = "hide"  
  # male 
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
  # morph
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
    xy2 = rbind(x,y)
    xy2[, part := factor(part, levels=c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total"))] 
#+ R_plot, fig.width=6, fig.height = 3
  g1 = 
      ggplot(xy, aes(x = part, y = pred, col = method_CI)) +
        geom_errorbar(aes(ymin = lwr, ymax = upr, col = method_CI), width = 0.1, position = position_dodge(width = 0.25) ) +
        #ggtitle ("Sim based")+
        geom_point(position = position_dodge(width = 0.25)) +
        scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
        scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
        labs(x = NULL, y = "Repeatability [%]", subtitle = "within male")+
        ylim(c(0,100))+
        coord_flip()+
        theme_bw() +
        theme(plot.title = element_text(size=9),
          legend.position = "none")

  #ggsave(here::here('Output/morpho_Repeatability_within-male.png'),g1, width = 10, height =7, units = 'cm')

  g2 = 
      ggplot(xy2, aes(x = part, y = pred, col = method_CI)) +
        geom_errorbar(aes(ymin = lwr, ymax = upr, col = method_CI), width = 0.1, position = position_dodge(width = 0.25) ) +
        #ggtitle ("Sim based")+
        geom_point(position = position_dodge(width = 0.25)) +
        scale_color_viridis(discrete=TRUE, begin=0, end = 0.5, guide = guide_legend(reverse = TRUE)) +
        #scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
        labs(x = NULL, y = "Repeatability [%]", subtitle = "within morph")+
        ylim(c(0,100))+
        coord_flip()+
        theme_bw() +
        theme(plot.title = element_text(size=9),
          axis.title.y = element_blank(), 
        axis.text.y = element_blank())
    
    #ggsave(here::here('Output/morpho_Repeatability_within-morph.png'),g, width = 10, height =7, units = 'cm')

  grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), size = "last"))
  #grid.arrange(g1,g2)
  ggsave(here::here('Output/morpho_Repeatability.png'),cbind(ggplotGrob(g1),ggplotGrob(g2) , size = "last"), width = 6*2, height =3*1.5, units = 'cm')  

#' ## Differences - "raw"
#' Red indicates mean
#+ boxplot, fig.width=7, fig.height = 7
    g1 =  # dummy to extract variables for median calculation
      ggplot(b, aes(x = Morph, y = Length_µm)) +
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(), aes(col = Morph, fill =Morph), dotsize = 0.5)+
      geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
      stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
      scale_fill_viridis(discrete=TRUE, alpha = 0.4)+
      scale_color_viridis(discrete=TRUE, alpha = 0.8)+
      facet_wrap(~part, scales = 'free_y', nrow = 2)+
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
                   position = position_dodge(), aes(col = Morph, fill =Morph), dotsize = 1)+
      #              position = position_dodge(), col = 'grey', aes(fill =Morph), dotsize = 1)+
      scale_fill_viridis(discrete=TRUE, alpha = 0.4)+
      scale_color_viridis(discrete=TRUE, alpha = 0.8)+
      stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
      facet_wrap(~part, scales = 'free_y', nrow = 2)+
      ggtitle('Male means') +
      ylab('Length [µm]') +
      theme_bw() +
      guides(x =  guide_axis(angle = -45)) +
      theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) #panel.spacing.y = unit(0, "mm")) #axis.title.x = element_blank(), axis.text.x = element_blank(), 
    
    cv_ =  b[, cv(Length_µm), by = list(bird_ID, part, Morph)]
    cv_[ , Morph123 := as.numeric(Morph)]
    names(cv_) [4]='CV' 
    
    g3 = 
     ggplot(cv_, aes(x =Morph , y = CV)) + 
     geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
     geom_dotplot(binaxis = 'y', stackdir = 'center',
                     position = position_dodge(), aes(col = Morph, fill =Morph), dotsize = 1)+
     stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
     scale_fill_viridis(discrete=TRUE, alpha = 0.4)+
     scale_color_viridis(discrete=TRUE, alpha = 0.8)+
     facet_wrap(~part, nrow = 2, scales = "free_y") +
     guides(x =  guide_axis(angle = -45)) +
     ylab('Coefficient of variation') +
     ggtitle('Male values') +
     theme_bw() +
     theme(legend.position = "none",
          plot.title = element_text(size=9),
          axis.title.x = element_blank()
          )
    
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    
    grid.draw(rbind(gg1, gg2, gg3, size = "last"))
      #grid.arrange(g1,g2)
 
    ggsave(here::here('Output/morpho_boxplots_.png'), rbind(gg1,gg2,gg3, size = "last"), width = 7*3, height =10*1.5, units = 'cm')  

    #ggsave(here::here('Output/morpho_per_male.png'),g, width = 10, height =10, units = 'cm')
#+ CV_box_alt, fig.width=6, fig.height = 3   
   g=
   ggplot(cv_, aes(x =Morph , y = CV, col = part, fill = part)) + 
   #geom_dotplot(binaxis = 'y', stackdir = 'center', position = position_dodge(), col = 'grey',  dotsize = 0.5)+
   geom_boxplot(fill = NA) + 
   theme_bw() +
   theme(axis.title.x = element_blank())
   g
   ggsave(here::here('Output/CV_per_male_alternative.png'),g, width = 10, height =7, units = 'cm')
#+ relative, fig.width=7/2, fig.height = 3.5
    g1 =  # dummy to extract variables for median calculation
      ggplot(br, aes(x = Morph, y = Length_rel)) +
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(), aes(col = Morph, fill =Morph), dotsize = 0.5)+
      geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
      stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
      scale_fill_viridis(discrete=TRUE, alpha = 0.4)+
      scale_color_viridis(discrete=TRUE, alpha = 0.8)+
      facet_wrap(~part, scales = 'free_y', nrow = 2)+
      ggtitle('Single measurements') +
      ylab('Relative to total sperm length') +
      theme_bw() +
      guides(x =  guide_axis(angle = -45)) +
      theme(legend.position = "none",
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) #panel.spacing.y = unit(0, "mm")) #, 
    g2 =  # dummy to extract variables for median calculation
      ggplot(ar, aes(x = Morph, y = Length_rel)) +
      geom_dotplot(binaxis = 'y', stackdir = 'center',
                   position = position_dodge(), aes(col = Morph, fill =Morph), dotsize = 0.5)+
      geom_boxplot(col = 'grey40', fill = NA, alpha = 0.2) + 
      stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
      scale_fill_viridis(discrete=TRUE, alpha = 0.4)+
      scale_color_viridis(discrete=TRUE, alpha = 0.8)+
      facet_wrap(~part, scales = 'free_y', nrow = 2)+
      ggtitle('Male means') +
      ylab('Relative to total sperm length') +
      theme_bw() +
      guides(x =  guide_axis(angle = -45)) +
      theme(legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        plot.title = element_text(size=9)
        ) 
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    grid.draw(cbind(gg1, gg2, size = "last"))
      #grid.arrange(g1,g2)
      
    ggsave(here::here('Output/morpho_boxplots_rel.png'),cbind(gg1,gg2, size = "last"), width = 7*2, height =10*1.5, units = 'cm')    

#' ## Morpho-space
#+ cor_with_95, fig.width=4, fig.height = 8
  # single sperm  
    bs = b[, quantile(Length_µm, prob = c(0.5)), by = list(part,Morph)]
    names(bs)[3] = 'median'
    bs$lwr = b[, quantile(Length_µm, prob = c(0.025)), by = list(part,Morph)]$V1
    bs$upr = b[, quantile(Length_µm, prob = c(0.975)), by = list(part,Morph)]$V1

    # mean +/-sd
    bs = b[, mean(Length_µm), by = list(part,Morph)]
    names(bs)[3] = 'median'
    bs$lwr = b[, mean(Length_µm)-sd(Length_µm), by = list(part,Morph)]$V1
    bs$upr = b[, mean(Length_µm)+sd(Length_µm), by = list(part,Morph)]$V1

    bs_f = bs[part == 'Flagellum']
    bs_f[,Midpiece_µm:= bs[part == 'Midpiece',.(median)]]

    bs_fm = bs[part == 'Midpiece']
    bs_fm[,Flagellum_µm:= bs[part == 'Flagellum',.(median)]]
    
    g0 =
    ggplot(h, aes(x = Flagellum_µm, y = Midpiece_µm, col = Morph, fill = Morph)) +
      geom_point(pch = 21)+
      geom_point(data = bs_f, aes(x = median, y = Midpiece_µm), col = 'black', size = 2) +
      geom_segment(data = bs_f, aes(x = lwr, xend = upr, y =Midpiece_µm, yend =Midpiece_µm), col = "black" ) +
      geom_segment(data = bs_fm, aes(x = Flagellum_µm, xend = Flagellum_µm, y =lwr, yend =upr), col = "black" )+
      geom_point(data = bs_f, aes(x = median, y = Midpiece_µm, col = Morph), size = 1) +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_f, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      annotate(geom="text", x=80, y=29.5, label="Independent", size = 2, col = colors[1], hjust = 0) +
      annotate(geom="text", x=80, y=29, label="Satellite", size = 2, col = colors[2], hjust = 0) +
      annotate(geom="text", x=80, y=28.5, label="Faeder", size = 2, col = colors[3], hjust = 0) +
      annotate(geom="text", x=80, y=28, label="Mean +/-SD", size = 2, hjust = 0) +
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

    bs_h = bs[part == 'Head']
    bs_h[,Midpiece_µm:= bs[part == 'Midpiece',.(median)]]

    bs_m = bs[part == 'Midpiece']
    bs_m[,Head_µm:= bs[part == 'Head',.(median)]]
    
    g1 =
    ggplot(h, aes(x = Head_µm, y = Midpiece_µm, col = Morph, fill = Morph)) +
      geom_point(pch = 21)+
      geom_point(data = bs_h, aes(x = median, y = Midpiece_µm), col = 'black', size = 2) +
      geom_segment(data = bs_h, aes(x = lwr, xend = upr, y =Midpiece_µm, yend =Midpiece_µm), col = "black" ) +
      geom_segment(data = bs_m, aes(x = Head_µm, xend = Head_µm, y =lwr, yend =upr), col = "black" )+
      geom_point(data = bs_h, aes(x = median, y = Midpiece_µm, col = Morph), size = 1) +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=25, y=29.5, label="Independent", size = 2, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=25, y=29, label="Satellite", size = 2, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=25, y=28.5, label="Faeder", size = 2, col = colors[3], hjust = 0) +
      #annotate(geom="text", x=25, y=28, label="Mean +/-SD", size = 2, hjust = 0) +
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
    ggplot(h, aes(x = Head_µm, y = Flagellum_µm, col = Morph, fill = Morph)) +
      geom_point(pch = 21)+
      geom_point(pch = 21)+
      geom_point(data = bs_h_, aes(x = median, y = Flagellum_µm), col = 'black', size = 2) +
      geom_segment(data = bs_h_, aes(x = lwr, xend = upr, y =Flagellum_µm, yend =Flagellum_µm), col = "black" ) +
      geom_segment(data = bs_f, aes(x = Head_µm, xend = Head_µm, y =lwr, yend =upr), col = "black" )+
      geom_point(data = bs_h_, aes(x = median, y = Flagellum_µm, col = Morph), size = 1) +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h_, aes(x = median, y = Flagellum_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
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
    ggplot(h, aes(x = Nucleus_µm, y = Acrosome_µm, col = Morph, fill = Morph)) +
     geom_point(pch = 21)+
      geom_point(data = bs_n, aes(x = median, y = Acrosome_µm), col = 'black', size = 2) +
      geom_segment(data = bs_n, aes(x = lwr, xend = upr, y =Acrosome_µm, yend =Acrosome_µm), col = "black" ) +
      geom_segment(data = bs_a, aes(x = Nucleus_µm, xend = Nucleus_µm, y =lwr, yend =upr), col = "black" )+
      geom_point(data = bs_n, aes(x = median, y = Acrosome_µm, col = Morph), size = 1) +
      # geom_point(data = bs_n, aes(x = median, y = Acrosome_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #scale_colour_viridis(discrete=TRUE)+
      theme_bw() +
      theme(legend.position = "none",
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        )  
  # male means
    as = a[, quantile(Length_avg, prob = c(0.5)), by = list(part,Morph)]
    names(as)[3] = 'median'
    as$lwr = a[, quantile(Length_avg, prob = c(0.025)), by = list(part,Morph)]$V1
    as$upr = a[, quantile(Length_avg, prob = c(0.975)), by = list(part,Morph)]$V1

    # mean +/-sd
    as = a[, mean(Length_avg), by = list(part,Morph)]
    names(as)[3] = 'median'
    as$lwr = a[, mean(Length_avg)-sd(Length_avg), by = list(part,Morph)]$V1
    as$upr = a[, mean(Length_avg)+sd(Length_avg), by = list(part,Morph)]$V1

    as_f = as[part == 'Flagellum']
    as_f[,Midpiece_µm:= as[part == 'Midpiece',.(median)]]

    as_fm = as[part == 'Midpiece']
    as_fm[,Flagellum_µm:= as[part == 'Flagellum',.(median)]]
    
    g0a =
    ggplot(ha, aes(x = Flagellum_µm, y = Midpiece_µm, col = Morph, fill = Morph)) +
      geom_point(pch = 21)+
      geom_point(data = as_f, aes(x = median, y = Midpiece_µm), col = 'black', size = 2) +
      geom_segment(data = as_f, aes(x = lwr, xend = upr, y =Midpiece_µm, yend =Midpiece_µm), col = "black" ) +
      geom_segment(data = as_fm, aes(x = Flagellum_µm, xend = Flagellum_µm, y =lwr, yend =upr), col = "black" )+
      geom_point(data = as_f, aes(x = median, y = Midpiece_µm, col = Morph), size = 1) +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = aas_f, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
     #scale_colour_viridis(discrete=TRUE)+
      ggtitle('Male means') +
      theme_bw() +
      theme(
        #legend.position=c(.9,.9),  
        legend.position = "none",
        #axis.title.x = element_blank(), 
        #axis.text.x = element_blank(),
        plot.title = element_text(size=9)
        ) 

    as_h = as[part == 'Head']
    as_h[,Midpiece_µm:= as[part == 'Midpiece',.(median)]]

    as_m = as[part == 'Midpiece']
    as_m[,Head_µm:= as[part == 'Head',.(median)]]

    g1a =
    ggplot(ha, aes(x = Head_µm, y = Midpiece_µm, col = Morph, fill = Morph)) +
      geom_point(pch = 21)+
      geom_point(data = as_h, aes(x = median, y = Midpiece_µm), col = 'black', size = 2) +
      geom_segment(data = as_h, aes(x = lwr, xend = upr, y =Midpiece_µm, yend =Midpiece_µm), col = "black" ) +
      geom_segment(data = as_m, aes(x = Head_µm, xend = Head_µm, y =lwr, yend =upr), col = "black" )+
      geom_point(data = as_h, aes(x = median, y = Midpiece_µm, col = Morph), size = 1) +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = as_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ggtitle('Male means') +
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
    ggplot(ha, aes(x = Head_µm, y = Flagellum_µm, col = Morph, fill = Morph)) +
      geom_point(pch = 21)+
      geom_point(data = as_h_, aes(x = median, y = Flagellum_µm), col = 'black', size = 2) +
      geom_segment(data = as_h_, aes(x = lwr, xend = upr, y =Flagellum_µm, yend =Flagellum_µm), col = "black" ) +
      geom_segment(data = as_f, aes(x = Head_µm, xend = Head_µm, y =lwr, yend =upr), col = "black" )+
      geom_point(data = as_h_, aes(x = median, y = Flagellum_µm, col = Morph), size = 1) +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
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
    ggplot(ha, aes(x = Nucleus_µm, y = Acrosome_µm, col = Morph, fill = Morph)) +
      geom_point(pch = 21)+
      geom_point(data = as_n, aes(x = median, y = Acrosome_µm), col = 'black', size = 2) +
      geom_segment(data = as_n, aes(x = lwr, xend = upr, y =Acrosome_µm, yend =Acrosome_µm), col = "black" ) +
      geom_segment(data = as_a, aes(x = Nucleus_µm, xend = Nucleus_µm, y =lwr, yend =upr), col = "black" )+
      geom_point(data = as_n, aes(x = median, y = Acrosome_µm, col = Morph), size = 1) +
      scale_colour_manual(values = colors)+
      scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = as_n, aes(x = median, y = Acrosome_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      theme_bw() +
      theme(legend.position = "none",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        plot.title = element_text(size=9)
        ) 
  # export
    gg0 <- ggplotGrob(g0)
    gg1 <- ggplotGrob(g1)
    gg2 <- ggplotGrob(g2) 
    gg3 <- ggplotGrob(g3) 
    gg0a <- ggplotGrob(g0a)
    gg1a <- ggplotGrob(g1a)
    gg2a <- ggplotGrob(g2a) 
    gg3a <- ggplotGrob(g3a) 

    grid.draw(cbind(rbind(gg0,gg1,gg2, gg3, size = "first"), rbind(gg0a,gg1a,gg2a, gg3a, size = "first"), size = "first"))
      #grid.arrange(g1,g2)
    
    ggsave(here::here('Output/morpho_corWithMeans.png'),cbind(rbind(gg0,gg1,gg2, gg3, size = "first"), rbind(gg0a,gg1a,gg2a, gg3a, size = "first"), size = "first"), width = 10, height =20, units = 'cm') 
#+ cor_3D, fig.width=10 
    bw[,morph123 :=ifelse(Morph == 'Independent', 1, ifelse(Morph == 'Satellite', 2,3))]             

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
    update(p[rep(1, prod(npanel))], layout = npanel,
           panel = function(..., screen) {
              panel.cloud(..., screen = list(z = rotz[current.column()],
                                                 x = rotx[current.row()]))
          })
#+ export, results = "hide"    
    png(here::here('Output/morpho_cor_3D.png'),width = 20, height = 10, units = "cm", res = 600)
      update(p[rep(1, prod(npanel))], layout = npanel,
       panel = function(..., screen) {
        panel.cloud(..., screen = list(z = rotz[current.column()],
                x = rotx[current.row()]))
            })
       dev.off()

#' ## Models
#' ### Effect sizes
#+ est_pred, results = "hide"
    # morpho  
      # for averages
        l = list()
        lp =list()
        lpr = list()
        lcv = list()
        lpcv = list()
        
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
        for(i in unique(ar$part)){
          if(i == 'Midpiece'){ii = 'Midpiece_rel'}
          if(i == 'Flagellum'){ii = 'Flagellum_rel'}
          m = lm(scale(Length_rel) ~ Morph, ar[part == i])
          #summary(m)
          #plot(allEffects(m))
          bsim = sim(m, n.sim=nsim) 
          v = apply(bsim@coef, 2, quantile, prob=c(0.5))
          ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
          l[[ii]]=data.frame(response=ii,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

          # get predictions
          m = lm(Length_rel ~ Morph, ar[part == i])
          bsim = sim(m, n.sim=nsim) 
          v = apply(bsim@coef, 2, quantile, prob=c(0.5))
          newD=data.frame(Morph = unique(a$Morph)) # values to predict for
          X <- model.matrix(~ Morph,data=newD) # exactly the model which was used has to be specified here 
          newD$Length_rel <-(X%*%v) 
          predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
          for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                          predmatrix[predmatrix < 0] <- 0
                          newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                          newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                          #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
          newD$part=i
          lpr[[i]] = data.table(newD)

          print(i)     
          }          
        for(i in unique(cv_$part)){
          #i ='Nucleus'
          m = lm(scale(CV) ~ Morph, cv_[part == i])
          #summary(m)
          #plot(allEffects(m))
          bsim = sim(m, n.sim=nsim) 
          v = apply(bsim@coef, 2, quantile, prob=c(0.5))
          ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
          lcv[[i]]=data.frame(response=paste('CV',i),effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

          # get predictions
          m = lm(CV ~ Morph, cv_[part == i])
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
          lpcv[[i]] = data.table(newD)

          print(i)     
          }          
        
        ll = data.table(do.call(rbind,l) ) 
        ll[, response := factor(response, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
        ll[, effect := factor(effect, levels=rev(c("(Intercept)", "MorphSatellite", "MorphFaeder")))]
        ll[, unit := 'male average']

        llcv = data.table(do.call(rbind,lcv) ) 
        llcv[, response := factor(response, levels=rev(c("CV Acrosome", "CV Nucleus", "CV Head", "CV Midpiece","CV Tail","CV Flagellum","CV Total")))] 
        llcv[, effect := factor(effect, levels=rev(c("(Intercept)", "MorphSatellite", "MorphFaeder")))]
        
        llp = data.table(do.call(rbind,lp) ) 
        llp[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total")))] 
        llp[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]

        llpr = data.table(do.call(rbind,lpr) ) 
        llpr[, part := factor(part, levels=rev(c("Midpiece","Flagellum")))] 
        llpr[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]

      # for single values 
        ls = list()
        lps = list()
        lpsr = list()
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
        for(i in unique(ar$part)){
          if(i == 'Midpiece'){ii = 'Midpiece_rel'}
          if(i == 'Flagellum'){ii = 'Flagellum_rel'}
          m = lmer(scale(Length_rel) ~ Morph + (1|bird_ID), br[part == i])
          #summary(m)
          #plot(allEffects(m))
          bsim = sim(m, n.sim=5000) 
          v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
          ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
          ls[[ii]]=data.frame(response=ii,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])
         
          # get predictions
          m = lmer(Length_rel ~ Morph + (1|bird_ID), br[part == i])
          bsim = sim(m, n.sim=nsim) 
          v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
          newD=data.frame(Morph = unique(ar$Morph)) # values to predict for
          X <- model.matrix(~ Morph,data=newD) # exactly the model which was used has to be specified here 
          newD$Length_rel <-(X%*%v) 
          predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
          for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                          predmatrix[predmatrix < 0] <- 0
                          newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                          newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                          #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
          newD$part=i
          lpsr[[i]] = data.table(newD)

          print(i)     
          }          
        
        lls = data.table(do.call(rbind,ls) ) 
        lls[, response := factor(response, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
        lls[, effect := factor(effect, levels=rev(c("(Intercept)", "MorphSatellite", "MorphFaeder")))]
        lls[, unit := 'single sperm']

        llps = data.table(do.call(rbind,lps) ) 
        llps[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
        llps[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]

        llpsr = data.table(do.call(rbind,lpsr) ) 
        llpsr[, part := factor(part, levels=rev(c("Midpiece","Flagellum")))] 
        llpsr[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]

    # motility
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lv = list()
      lvp =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      #dd = a[part =='Acrosome']
      #dd[, motileCount_ln:=scale(log(motileCount))]
      # VAP
        m = lm(scale(VAP) ~ scale(log(motileCount)) + Morph, dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lv[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VAP ~ motileCount_ln + Morph, dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvp[['vap']] = data.table(newD)   
      # VSL
        m = lm(scale(VSL) ~ scale(log(motileCount))  + Morph, dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lv[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VSL ~ motileCount_ln + Morph, dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln),Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvp[['VSL']] = data.table(newD)
      # VCL
        m = lm(scale(VCL) ~scale(log(motileCount)) + Morph, dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lv[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VCL ~ motileCount_ln + Morph, dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvp[['VCL']] = data.table(newD) 
                 
      llv = data.table(do.call(rbind,lv) ) 
      llv = llv[effect != 'scale(log(motileCount))' ]
      llv[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llv[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llv[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llv[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      
      llv[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvp = data.table(do.call(rbind,lvp) ) 
      llvp[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvp[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
    # motility 2 - only one measure per bird
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lvx = list()
      lvpx =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      dd1 = dd[month == 'June']
      dd2 = dd[month == 'May']
      ddx = rbind(dd1,dd2[!bird_ID%in%dd1$bird_ID])
      # VAP
        m = lm(scale(VAP) ~ scale(log(motileCount)) + Morph, ddx)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lvx[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VAP ~ motileCount_ln + Morph, ddx)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(ddx$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvpx[['vap']] = data.table(newD)   
      # VSL
        m = lm(scale(VSL) ~ scale(log(motileCount))  + Morph, ddx)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lvx[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VSL ~ motileCount_ln + Morph, ddx)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(ddx$motileCount_ln),Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvpx[['VSL']] = data.table(newD)
      # VCL
        m = lm(scale(VCL) ~scale(log(motileCount)) + Morph, ddx)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
        lvx[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lm(VCL ~ motileCount_ln + Morph, ddx)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@coef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(ddx$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvpx[['VCL']] = data.table(newD) 
                 
      llvx = data.table(do.call(rbind,lvx) ) 
      llvx = llvx[effect != 'scale(log(motileCount))' ]
      llvx[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvx[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llvx[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llvx[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llvx[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvpx = data.table(do.call(rbind,lvpx) ) 
      llvpx[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvpx[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
    # motility 3 - mixed model
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lvq = list()
      lvpq =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      # VAP
        m = lmer(scale(VAP) ~ scale(log(motileCount)) + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvq[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VAP ~ motileCount_ln + Morph+ (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvpq[['vap']] = data.table(newD)   
      # VSL
        m = lmer(scale(VSL) ~ scale(log(motileCount))  + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvq[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VSL ~ motileCount_ln + Morph + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln),Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvpq[['VSL']] = data.table(newD)
      # VCL
        m = lmer(scale(VCL) ~scale(log(motileCount)) + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvq[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VCL ~ motileCount_ln + Morph + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvpq[['VCL']] = data.table(newD) 
                 
      llvq = data.table(do.call(rbind,lvq) ) 
      llvq = llvq[effect != 'scale(log(motileCount))' ]
      llvq[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvq[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llvq[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llvq[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llvq[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvpq = data.table(do.call(rbind,lvpq) ) 
      llvpq[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvpq[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
    # motility 4 - mixed model with month
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lvm = list()
      lvpm =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      # VAP
        m = lmer(scale(VAP) ~ scale(log(motileCount)) + month + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvm[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VAP ~ motileCount_ln+ month + Morph+ (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5, Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvpm[['vap']] = data.table(newD)   
      # VSL
        m = lmer(scale(VSL) ~ scale(log(motileCount))+ month  + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvm[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VSL ~ motileCount_ln + Morph+ month + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5,Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvpm[['VSL']] = data.table(newD)
      # VCL
        m = lmer(scale(VCL) ~scale(log(motileCount))+ month + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvm[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VCL ~ motileCount_ln+ month + Morph + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5, Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvpm[['VCL']] = data.table(newD) 
                 
      llvm = data.table(do.call(rbind,lvm) ) 
      llvm = llvm[!effect %in% c('scale(log(motileCount))','monthMay') ]
      llvm[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvm[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llvm[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llvm[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llvm[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvpm = data.table(do.call(rbind,lvpm) ) 
      llvpm[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvpm[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
    # motility 5 - mixed model, only data without issues
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lvi = list()
      lvpi =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      ddi = dd[issues == 'zero']
      # VAP
        m = lmer(scale(VAP) ~ scale(log(motileCount)) + Morph + (1|bird_ID), ddi)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvi[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VAP ~ motileCount_ln + Morph+ (1|bird_ID), ddi)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(ddi$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvpi[['vap']] = data.table(newD)   
      # VSL
        m = lmer(scale(VSL) ~ scale(log(motileCount))  + Morph + (1|bird_ID), ddi)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvi[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VSL ~ motileCount_ln + Morph + (1|bird_ID), ddi)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(ddi$motileCount_ln),Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvpi[['VSL']] = data.table(newD)
      # VCL
        m = lmer(scale(VCL) ~scale(log(motileCount)) + Morph + (1|bird_ID), ddi)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvi[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VCL ~ motileCount_ln + Morph + (1|bird_ID), ddi)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(ddi$motileCount_ln), Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvpi[['VCL']] = data.table(newD)              
      llvi = data.table(do.call(rbind,lvi) ) 
      llvi = llvi[effect != 'scale(log(motileCount))' ]
      llvi[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvi[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llvi[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llvi[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llvi[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvpi = data.table(do.call(rbind,lvpi) ) 
      llvpi[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvpi[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
    # motility 6- mixed model with month and issuee
      #chart.Correlation(d[, c('VAP', 'VSL','VCL', 'motileCount','motileCount_ln', 'NumberFields')], histogram=TRUE, pch=19)
      #mtext("Single sperm", side=3, line=3)
      lvs = list()
      lvps =list()
      d[, motileCount_ln:=scale(log(motileCount))]
      dd = d[!Morph%in%'Zebra finch']
      dd[issues == 'zero', issue:='no']
      dd[issues != 'zero', issue:='yes']
      # VAP
        m = lmer(scale(VAP) ~ scale(log(motileCount)) + month + issue +  Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvs[['VAP']]=data.frame(response='Average-path (VAP)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VAP ~ motileCount_ln+ month + issue + Morph+ (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5, issue = 0.5,  Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + issue + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Average-path (VAP)'
        lvps[['vap']] = data.table(newD)   
      # VSL
        m = lmer(scale(VSL) ~ scale(log(motileCount))+ month + issue  + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvs[['VSL']]=data.frame(response='Straight-line (VSL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VSL ~ motileCount_ln + Morph+ month + issue + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5, issue = 0.5, Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + issue + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Straight-line (VSL)'
        lvps[['VSL']] = data.table(newD)
      # VCL
        m = lmer(scale(VCL) ~scale(log(motileCount))+ month + issue + Morph + (1|bird_ID), dd)
        #summary(m)
        #plot(allEffects(m))
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
        lvs[['VCL']]=data.frame(response='Curvilinear (VCL)',effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

        # get predictions
        m = lmer(VCL ~ motileCount_ln+ month + Morph + issue + (1|bird_ID), dd)
        bsim = sim(m, n.sim=nsim) 
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
        newD=data.frame(motileCount_ln = mean(dd$motileCount_ln), month = 0.5, issue = 0.5, Morph = unique(b$Morph)) # values to predict for
        X <- model.matrix(~ motileCount_ln+ month + issue + Morph,data=newD) # exactly the model which was used has to be specified here 
        newD$pred <-(X%*%v) 
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@fixef[j,])}
                    predmatrix[predmatrix < 0] <- 0
                    newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                    newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                    #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
        newD$motility = 'Curvilinear (VCL)'
        lvps[['VCL']] = data.table(newD) 
                 
      llvs = data.table(do.call(rbind,lvs) ) 
      llvs = llvs[!effect %in% c('scale(log(motileCount))','monthMay','issueyes') ]
      llvs[, response := factor(response, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvs[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llvs[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llvs[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llvs[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llvps = data.table(do.call(rbind,lvps) ) 
      llvps[, motility := factor(motility, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
      llvps[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]
    
#+ effect_sizes, fig.width =4, fig.height = 4.5        
        llll = rbind(ll,lls)  
        llll[, unit := factor(unit, levels=rev(c("single sperm", "male average")))] 
        llll[effect == '(Intercept)', effect:='Independent\n(Intercept)']
        llll[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
        llll[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
        llll[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 
        
        g = 
        ggplot(llll, aes(y = effect, x = estimate, col = response, shape = unit)) +
          geom_vline(xintercept = 0, col = "grey30", lty =3)+
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = response), width = 0.1, position = position_dodge(width = 0.4) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.4)) +
          #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          scale_shape(guide = guide_legend(reverse = TRUE)) + 
          #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
          labs(y = NULL ,x = "Standardized effect size",title = 'Sperm-part specific models') +
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
        ggsave(here::here('Output/morpho_effectSizes_virid.png'),g, width = 12, height =10, units = 'cm')
        #ggsave('Output/morpho_effectSizes_pair.png',g, width = 12, height =10, units = 'cm')
#'   
#' <span style="color: red;">!!! Satellites seem to have the longest tails (and flagellums), faeders the longest midpieces, which translates into the differences in relative measurements.  !!!</span>    
#' 
#' ***
#'
#+ effect_sizes_CV, fig.width =4, fig.height = 3  
        llcv[effect == '(Intercept)', effect:='Independent\n(Intercept)']
        llcv[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
        llcv[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
        llcv[, effect := factor(effect, levels=c("Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 
        
        g = 
        ggplot(llcv, aes(y = effect, x = estimate, col = response)) +
          geom_vline(xintercept = 0, col = "grey30", lty =3)+
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = response), width = 0.1, position = position_dodge(width = 0.4) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.4)) +
          #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          #scale_shape(guide = guide_legend(reverse = TRUE)) + 
          #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
          labs(y = NULL ,x = "Standardized effect size",title = 'Sperm-part specific models on male CV') +
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
        ggsave(here::here('Output/morpho_effectSizes_virid_CV.png'),g, width = 12, height =10, units = 'cm')
#+ effect_sizes_motil_ALL, fig.width =4, fig.height = 2 
        llvx[,model := 'linear on June recordings']
        llvq[,model := 'mixed, all recordings']
        llvm[,model := 'mixed, all recordings & control for month']
        llvs[,model := 'mixed, all recordings & control for month & issues']
        llvi[,model := 'mixed, recordings without issues']
        

        xqm = rbind(llvx,llvq, llvm, llvi, llvs)
        xqm[, model := factor(model, levels=c('mixed, recordings without issues', 'mixed, all recordings & control for month & issues','mixed, all recordings & control for month', 'mixed, all recordings',"linear on June recordings"))] 
        #table(dd$issues,dd$Morph) 
        g2 = 
        ggplot(xqm, aes(y = effect, x = estimate, col = response, fill = response, shape = model)) +
          geom_vline(xintercept = 0, col = "grey30", lty =3)+
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = response), width = 0.1, position = position_dodge(width = 0.6) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.6)) +
          #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          scale_shape_manual(guide = guide_legend(reverse = TRUE, override.aes = list(fill = c('black'))), values =c(25, 24,23,22,21))+
          #scale_shape_manual(guide = guide_legend(reverse = TRUE), values =c(5,2,0,1)) + 
          #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
          labs(y = NULL ,x = "Standardized effect size",title = 'Velocity-type specific models') +
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
        g2
        ggsave(here::here('Output/motility_effectSizes_virid_all.png'),g2, width = 14, height =10, units = 'cm')
# effect_sizes_motility, fig.width =4, fig.height = 2 
        g2 = 
        ggplot(llvx, aes(y = effect, x = estimate, col = response)) +
          geom_vline(xintercept = 0, col = "grey30", lty =3)+
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = response), width = 0.1, position = position_dodge(width = 0.4) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.4)) +
          #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          #scale_shape(guide = guide_legend(reverse = TRUE)) + 
          #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
          labs(y = NULL ,x = "Standardized effect size",title = 'Velocity-type specific models on male velocity\nusing single recording per male (for all but three males from June)') +
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
        #g2
        ggsave(here::here('Output/motility_effectSizes_virid.png'),g2, width = 10, height =8, units = 'cm')
# effect_sizes_motil_mix, fig.width =4, fig.height = 2 
        g2 = 
        ggplot(llvq, aes(y = effect, x = estimate, col = response)) +
          geom_vline(xintercept = 0, col = "grey30", lty =3)+
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = response), width = 0.1, position = position_dodge(width = 0.4) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.4)) +
          #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          #scale_shape(guide = guide_legend(reverse = TRUE)) + 
          #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
          labs(y = NULL ,x = "Standardized effect size",title = 'Velocity-type specific models on male velocity\nusing all recording and bird_ID as random intercept') +
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
        #g2
        ggsave(here::here('Output/motility_effectSizes_virid_mixed.png'),g2, width = 10, height =8, units = 'cm')
# effect_sizes_motil_mix_m, fig.width =4, fig.height = 2 
        g2 = 
        ggplot(llvm, aes(y = effect, x = estimate, col = response)) +
          geom_vline(xintercept = 0, col = "grey30", lty =3)+
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = response), width = 0.1, position = position_dodge(width = 0.4) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.4)) +
          #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          #scale_shape(guide = guide_legend(reverse = TRUE)) + 
          #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
          labs(y = NULL ,x = "Standardized effect size",title = 'Velocity-type specific models on male velocity\nusing all recordings, controlled for month and bird_ID (random intercept)') +
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
        #g2
        ggsave(here::here('Output/motility_effectSizes_virid_mixed_m.png'),g2, width = 10, height =8, units = 'cm')
#'
#' <span style="color: red;">!!! Against prediction about selection for the fastest sperm, faeders might have most variable & slowest sperm !!!</span>    
#' 
#' ***
#' ### Predictions with raw data
#+ effectSizes&dist, fig.width =7 
      g1 = 
          ggplot(a, aes(x = Morph, y = Length_avg)) +
          geom_dotplot(binaxis = 'y', stackdir = 'center',
                        position = position_dodge(), col = 'grey', fill = 'lightgrey', dotsize = 1.5)+
          geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2) + 
          geom_errorbar(data = llp, aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(width = 0.25), col = 'red' ) +
          geom_point(data = llp, aes(x = Morph, y =Length_avg), position = position_dodge(width = 0.25), col = 'red', size = 0.75) +
          #stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
          scale_color_viridis(discrete=TRUE)+
          facet_wrap(~part, scales = 'free_y', nrow = 2)+
          labs(title = 'Predictions from sperm part specific models\non male means') +
          ylab('Length [µm]') +
          theme_bw() +
          guides(x =  guide_axis(angle = -45)) +
          theme(legend.position = "none",
            axis.title.x = element_blank(), 
            #axis.text.x = element_blank(),
            plot.title = element_text(size=9)
            ) #panel.spacing.y = unit(0, "mm")) #,    
 
      g2 = 
          ggplot(b, aes(x = Morph, y = Length_µm)) +
          geom_dotplot(binaxis = 'y', stackdir = 'center',
                        position = position_dodge(), col = 'grey', fill = 'lightgrey', dotsize = 1)+
          geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2) + 
          geom_errorbar(data = llps, aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(width = 0.25), col = 'red' ) +
          geom_point(data = llps, aes(x = Morph, y =Length_µm), position = position_dodge(width = 0.25), col = 'red', size = 0.75) +
          #stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
          scale_color_viridis(discrete=TRUE)+
          facet_wrap(~part, scales = 'free_y', nrow = 2)+
          labs(title = 'Predictions from sperm part specific models\non single sperm') +
          ylab('Length [µm]') +
          theme_bw() +
          guides(x =  guide_axis(angle = -45)) +
          theme(legend.position = "none",
            axis.title.x = element_blank(), 
            axis.text.x = element_blank(),
            plot.title = element_text(size=9)
            ) #panel.spacing.y = unit(0, "mm")) #,     

      grid.draw(rbind(ggplotGrob(g2), ggplotGrob(g1), size = "last"))
      #grid.arrange(g1,g2)
      gg1 <- ggplotGrob(g1)
      gg2 <- ggplotGrob(g2) 
      ggsave(here::here('Output/morpho_predictions+boxPlots.png'),rbind(gg2,gg1, size = "last"), width = 15, height =15, units = 'cm')    

#+ effectSizesRel&dist, fig.width =7/2, fig.height =7/2 
      g1 = 
          ggplot(br, aes(x = Morph, y = Length_rel)) +
          geom_dotplot(binaxis = 'y', stackdir = 'center',
                        position = position_dodge(), col = 'grey', fill = 'lightgrey', dotsize = 1)+
          geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2) + 
          geom_errorbar(data = llpsr, aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(width = 0.25), col = 'red' ) +
          geom_point(data = llpsr, aes(x = Morph, y =Length_rel), position = position_dodge(width = 0.25), col = 'red', size = 0.75) +
          #stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
          scale_color_viridis(discrete=TRUE)+
          facet_wrap(~part, scales = 'free_y', nrow = 2)+
          ylab('Relative length to total sperm length') +
          labs(title = 'Predictions from sperm part specific models\non single sperm') +
          theme_bw() +
          guides(x =  guide_axis(angle = -45)) +
          theme(legend.position = "none",
            plot.title = element_text(size=9)
            ) #panel.spacing.y = unit(0, "mm")) #,     

      g2 = 
          ggplot(ar, aes(x = Morph, y = Length_rel)) +
          geom_dotplot(binaxis = 'y', stackdir = 'center',
                        position = position_dodge(), col = 'grey', fill = 'lightgrey', dotsize = 1.5)+
          geom_boxplot(aes(col = Morph), fill = NA, alpha = 0.2) + 
          geom_errorbar(data = llpr, aes(ymin = lwr, ymax = upr), width = 0, position = position_dodge(width = 0.25), col = 'red' ) +
          geom_point(data = llpr, aes(x = Morph, y =Length_rel), position = position_dodge(width = 0.25), col = 'red', size = 0.75) +
          #stat_summary(fun.y=mean, geom="point", color="red", fill="red") +
          scale_color_viridis(discrete=TRUE)+
          facet_wrap(~part, scales = 'free_y', nrow = 2)+
          labs(title = '\non male means') +
          theme_bw() +
          guides(x =  guide_axis(angle = -45)) +
          theme(legend.position = "none",
            axis.title.y = element_blank(), 
            axis.text.y = element_blank(),
            plot.title = element_text(size=9)
            ) #panel.spacing.y = unit(0, "mm")) #,    
 
      gg1 <- ggplotGrob(g1)
      gg2 <- ggplotGrob(g2) 

      grid.draw(cbind(gg1, gg2, size = "last"))
      #grid.arrange(g1,g2)
      
      ggsave(here::here('Output/morpho_predictions+boxPlots_rel.png'),rbind(gg1,gg2, size = "last"), width = 7.5, height =15/4, units = 'cm')   
#'
#' ***
#' ### Does morpho predict motility?
#+ pred, results = "hide"   
    # not controlled for motile count  
      l = list()
      lp =list()
      # VAP
         for(i in unique(a$part)){
            #i ='Nucleus'
            ai = a[part == i]
            m = lm(scale(VAP) ~ Morph+scale(Length_avg), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            l[[paste('VAP',i)]]=data.frame(mot = 'VAP', part=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VAP ~ Morph+Length_avg, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
              }
            newD=do.call(rbind,nd)

            # values to predict for
            X <- model.matrix(~ Morph+Length_avg,data=newD) # exactly the model which was used has to be specified here 
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
            m = lm(scale(VSL) ~ Morph+scale(Length_avg), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            l[[paste('VSL',i)]]=data.frame(mot = 'VSL', part=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VSL ~ Morph+Length_avg, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
              }
            newD=do.call(rbind,nd)
            X <- model.matrix(~ Morph+Length_avg,data=newD) # exactly the model which was used has to be specified here 
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
            m = lm(scale(VCL) ~ Morph+scale(Length_avg), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            l[[paste('VCL',i)]]=data.frame(mot = 'VCL', part=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VCL ~ Morph+Length_avg, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
              }
            newD=do.call(rbind,nd)
            X <- model.matrix(~ Morph+Length_avg,data=newD) # exactly the model which was used has to be specified here 
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
      ll[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
      ll[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      ll[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      ll[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      ll[effect == 'scale(Length_avg)', effect := 'Length_µm']
      
      ll[, effect := factor(effect, levels=c('Length_µm',"Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      ll[mot=='VCL', mot:='Curvilinear (VCL)']
      ll[mot=='VSL', mot:='Straight-line (VSL)']
      ll[mot=='VAP', mot:='Average-path (VAP)']
      ll[, mot := factor(mot, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 

      llp = data.table(do.call(rbind,lp) ) 
      llp[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
      llp[, Morph := factor(Morph, levels=rev(c("Independent", "Satellite", "Faeder")))]  
    # controlled for motile count  
      lz = list()
      lpz =list()
      lprz =list()
      a[,motileCount_ln_z := scale(log(motileCount))]
      ar[,motileCount_ln_z := scale(log(motileCount))]
      # VAP
         for(i in unique(a$part)){
            #i ='Nucleus'
            ai = a[part == i]
            m = lm(scale(VAP) ~ motileCount_ln_z+ Morph+scale(Length_avg), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lz[[paste('VAP',i)]]=data.frame(mot = 'VAP', part=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VAP ~ motileCount_ln_z + Morph+Length_avg, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z),Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
              }
            newD=do.call(rbind,nd)

            # values to predict for
            X <- model.matrix(~ motileCount_ln_z + Morph+Length_avg,data=newD) # exactly the model which was used has to be specified here 
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
            lpz[[paste(i,newD$mot[1])]] = data.table(newD)

            print(paste(i,newD$mot[1]))     
            }          
         for(i in unique(ar$part)){
            #i ='Nucleus'
            if(i == 'Midpiece'){ii = 'Midpiece_rel'}
            if(i == 'Flagellum'){ii = 'Flagellum_rel'}
            ai = ar[part == i]
            m = lm(scale(VAP) ~ motileCount_ln_z+ Morph+scale(Length_rel), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lz[[paste('VAP',ii)]]=data.frame(mot = 'VAP', part=ii,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VAP ~ motileCount_ln_z + Morph+Length_rel, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z),Morph = j, Length_rel = seq(min(ai$Length_rel), max(ai$Length_rel), length.out = 200)) 
              }
            newD=do.call(rbind,nd)

            # values to predict for
            X <- model.matrix(~ motileCount_ln_z + Morph+Length_rel,data=newD) # exactly the model which was used has to be specified here 
            newD$VAP <-(X%*%v) 
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                            predmatrix[predmatrix < 0] <- 0
                            newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                            newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                            #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
            newD$part=ii
            newD$mot = 'VAP'
            setnames(newD, old = 'VAP', new = 'motility')
            lpz[[paste(ii,newD$mot[1])]] = data.table(newD)

            print(paste(ii,newD$mot[1]))     
            }            
      # VSL
         for(i in unique(a$part)){
            #i ='Nucleus'
            ai = a[part == i]
            m = lm(scale(VSL) ~ motileCount_ln_z + Morph+scale(Length_avg), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lz[[paste('VSL',i)]]=data.frame(mot = 'VSL', part=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VSL ~ motileCount_ln_z + Morph+Length_avg, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z), Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
              }
            newD=do.call(rbind,nd)
            X <- model.matrix(~ motileCount_ln_z + Morph+Length_avg,data=newD) # exactly the model which was used has to be specified here 
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
            lpz[[paste(i,newD$mot[1])]] = data.table(newD)

            print(paste(i,newD$mot[1]))     
            }          
         for(i in unique(ar$part)){
            #i ='Nucleus'
            if(i == 'Midpiece'){ii = 'Midpiece_rel'}
            if(i == 'Flagellum'){ii = 'Flagellum_rel'}
            ai = ar[part == i]
            m = lm(scale(VSL) ~ motileCount_ln_z+ Morph+scale(Length_rel), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lz[[paste('VSL',ii)]]=data.frame(mot = 'VSL', part=ii,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VSL ~ motileCount_ln_z + Morph+Length_rel, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z),Morph = j, Length_rel = seq(min(ai$Length_rel), max(ai$Length_rel), length.out = 200)) 
              }
            newD=do.call(rbind,nd)

            # values to predict for
            X <- model.matrix(~ motileCount_ln_z + Morph+Length_rel,data=newD) # exactly the model which was used has to be specified here 
            newD$VSL <-(X%*%v) 
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                            predmatrix[predmatrix < 0] <- 0
                            newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                            newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                            #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
            newD$part=ii
            newD$mot = 'VSL'
            setnames(newD, old = 'VSL', new = 'motility')
            lpz[[paste(ii,newD$mot[1])]] = data.table(newD)

            print(paste(ii,newD$mot[1]))     
            }            
      # VCL
         for(i in unique(a$part)){
            #i ='Nucleus'
            ai = a[part == i]
            m = lm(scale(VCL) ~ motileCount_ln_z + Morph+scale(Length_avg), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lz[[paste('VCL',i)]]=data.frame(mot = 'VCL', part=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VCL ~ motileCount_ln_z + Morph+Length_avg, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z), Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
              }
            newD=do.call(rbind,nd)
            X <- model.matrix(~ motileCount_ln_z + Morph+Length_avg,data=newD) # exactly the model which was used has to be specified here 
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
            lpz[[paste(i,newD$mot[1])]] = data.table(newD)

            print(paste(i,newD$mot[1]))     
            }          
         for(i in unique(ar$part)){
            #i ='Nucleus'
            if(i == 'Midpiece'){ii = 'Midpiece_rel'}
            if(i == 'Flagellum'){ii = 'Flagellum_rel'}
            ai = ar[part == i]
            m = lm(scale(VCL) ~ motileCount_ln_z+ Morph+scale(Length_rel), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lz[[paste('VCL',ii)]]=data.frame(mot = 'VCL', part=ii,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VCL ~ motileCount_ln_z + Morph+Length_rel, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z),Morph = j, Length_rel = seq(min(ai$Length_rel), max(ai$Length_rel), length.out = 200)) 
              }
            newD=do.call(rbind,nd)

            # values to predict for
            X <- model.matrix(~ motileCount_ln_z + Morph+Length_rel,data=newD) # exactly the model which was used has to be specified here 
            newD$VCL <-(X%*%v) 
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                            predmatrix[predmatrix < 0] <- 0
                            newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                            newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                            #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
            newD$part=ii
            newD$mot = 'VCL'
            setnames(newD, old = 'VCL', new = 'motility')
            lpz[[paste(ii,newD$mot[1])]] = data.table(newD)

            print(paste(ii,newD$mot[1]))     
            }          
               
      llz = data.table(do.call(rbind,lz) ) 
      llz[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
      llz[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llz[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llz[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llz[effect == 'scale(Length_avg)', effect := 'Length_µm']
      llz[effect == 'scale(Length_rel)', effect := 'Length_µm'] # dummy variable
      
      llz[, effect := factor(effect, levels=c('Length_µm','motileCount_ln_z',"Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llz[mot=='VCL', mot:='Curvilinear (VCL)']
      llz[mot=='VSL', mot:='Straight-line (VSL)']
      llz[mot=='VAP', mot:='Average-path (VAP)']
      llz[, mot := factor(mot, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
    # controlled for motile count & inbreeding  
      lzi = list()
      lpzi =list()
      lprzi =list()
      a[,motileCount_ln_z := scale(log(motileCount))]
      a[,inbreeding_z := scale(inbreeding)]
      ar[,motileCount_ln_z := scale(log(motileCount))]
      ar[,inbreeding_z := scale(inbreeding)]
      # VAP
         for(i in unique(a$part)){
            #i ='Nucleus'
            ai = a[part == i]
            m = lm(scale(VAP) ~ motileCount_ln_z+ inbreeding_z + Morph+scale(Length_avg), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lzi[[paste('VAP',i)]]=data.frame(mot = 'VAP', part=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VAP ~ motileCount_ln_z + inbreeding_z + Morph+Length_avg, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z),inbreeding_z = mean(ai$inbreeding_z), Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
              }
            newD=do.call(rbind,nd)

            # values to predict for
            X <- model.matrix(~ motileCount_ln_z + inbreeding_z + Morph+Length_avg,data=newD) # exactly the model which was used has to be specified here 
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
            lpzi[[paste(i,newD$mot[1])]] = data.table(newD)

            print(paste(i,newD$mot[1]))     
            }          
         for(i in unique(ar$part)){
            #i ='Nucleus'
            if(i == 'Midpiece'){ii = 'Midpiece_rel'}
            if(i == 'Flagellum'){ii = 'Flagellum_rel'}
            ai = ar[part == i]
            m = lm(scale(VAP) ~ motileCount_ln_z + inbreeding_z+ Morph+scale(Length_rel), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lzi[[paste('VAP',ii)]]=data.frame(mot = 'VAP', part=ii,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VAP ~ motileCount_ln_z +inbreeding_z+ Morph+Length_rel, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z),inbreeding_z = mean(ai$inbreeding_z), Morph = j, Length_rel = seq(min(ai$Length_rel), max(ai$Length_rel), length.out = 200)) 
              }
            newD=do.call(rbind,nd)

            # values to predict for
            X <- model.matrix(~ motileCount_ln_z +inbreeding_z + Morph+Length_rel,data=newD) # exactly the model which was used has to be specified here 
            newD$VAP <-(X%*%v) 
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                            predmatrix[predmatrix < 0] <- 0
                            newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                            newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                            #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
            newD$part=ii
            newD$mot = 'VAP'
            setnames(newD, old = 'VAP', new = 'motility')
            lpzi[[paste(ii,newD$mot[1])]] = data.table(newD)

            print(paste(ii,newD$mot[1]))     
            }            
      # VSL
         for(i in unique(a$part)){
            #i ='Nucleus'
            ai = a[part == i]
            m = lm(scale(VSL) ~ motileCount_ln_z + inbreeding_z + Morph+scale(Length_avg), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lzi[[paste('VSL',i)]]=data.frame(mot = 'VSL', part=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VSL ~ motileCount_ln_z + inbreeding_z + Morph+Length_avg, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z), inbreeding_z = mean(ai$inbreeding_z), Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
              }
            newD=do.call(rbind,nd)
            X <- model.matrix(~ motileCount_ln_z +inbreeding_z + Morph+Length_avg,data=newD) # exactly the model which was used has to be specified here 
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
            lpzi[[paste(i,newD$mot[1])]] = data.table(newD)

            print(paste(i,newD$mot[1]))     
            }          
         for(i in unique(ar$part)){
            #i ='Nucleus'
            if(i == 'Midpiece'){ii = 'Midpiece_rel'}
            if(i == 'Flagellum'){ii = 'Flagellum_rel'}
            ai = ar[part == i]
            m = lm(scale(VSL) ~ motileCount_ln_z+ Morph+scale(Length_rel), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lzi[[paste('VSL',ii)]]=data.frame(mot = 'VSL', part=ii,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VSL ~ motileCount_ln_z + Morph+Length_rel, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z), inbreeding_z = mean(ai$inbreeding_z),Morph = j, Length_rel = seq(min(ai$Length_rel), max(ai$Length_rel), length.out = 200)) 
              }
            newD=do.call(rbind,nd)

            # values to predict for
            X <- model.matrix(~ motileCount_ln_z + Morph+Length_rel,data=newD) # exactly the model which was used has to be specified here 
            newD$VSL <-(X%*%v) 
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                            predmatrix[predmatrix < 0] <- 0
                            newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                            newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                            #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
            newD$part=ii
            newD$mot = 'VSL'
            setnames(newD, old = 'VSL', new = 'motility')
            lpzi[[paste(ii,newD$mot[1])]] = data.table(newD)

            print(paste(ii,newD$mot[1]))     
            }            
      # VCL
         for(i in unique(a$part)){
            #i ='Nucleus'
            ai = a[part == i]
            m = lm(scale(VCL) ~ motileCount_ln_z + inbreeding_z + Morph+scale(Length_avg), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lzi[[paste('VCL',i)]]=data.frame(mot = 'VCL', part=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VCL ~ motileCount_ln_z + inbreeding_z + Morph+Length_avg, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z), inbreeding_z = mean(ai$inbreeding_z), Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
              }
            newD=do.call(rbind,nd)
            X <- model.matrix(~ motileCount_ln_z + inbreeding_z + Morph+Length_avg,data=newD) # exactly the model which was used has to be specified here 
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
            lpzi[[paste(i,newD$mot[1])]] = data.table(newD)

            print(paste(i,newD$mot[1]))     
            }          
         for(i in unique(ar$part)){
            #i ='Nucleus'
            if(i == 'Midpiece'){ii = 'Midpiece_rel'}
            if(i == 'Flagellum'){ii = 'Flagellum_rel'}
            ai = ar[part == i]
            m = lm(scale(VCL) ~ motileCount_ln_z+ inbreeding_z + Morph+scale(Length_rel), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lzi[[paste('VCL',ii)]]=data.frame(mot = 'VCL', part=ii,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VCL ~ motileCount_ln_z +inbreeding_z + Morph+Length_rel, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z), inbreeding_z = mean(ai$inbreeding_z),Morph = j, Length_rel = seq(min(ai$Length_rel), max(ai$Length_rel), length.out = 200)) 
              }
            newD=do.call(rbind,nd)

            # values to predict for
            X <- model.matrix(~ motileCount_ln_z +inbreeding_z + Morph+Length_rel,data=newD) # exactly the model which was used has to be specified here 
            newD$VCL <-(X%*%v) 
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                            predmatrix[predmatrix < 0] <- 0
                            newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                            newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                            #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
            newD$part=ii
            newD$mot = 'VCL'
            setnames(newD, old = 'VCL', new = 'motility')
            lpzi[[paste(ii,newD$mot[1])]] = data.table(newD)

            print(paste(ii,newD$mot[1]))     
            }          
               
      llzi = data.table(do.call(rbind,lzi) ) 
      llzi[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
      llzi[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llzi[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llzi[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llzi[effect == 'scale(Length_avg)', effect := 'Length_µm']
      llzi[effect == 'scale(Length_rel)', effect := 'Length_µm'] # dummy variable
      
      llzi[, effect := factor(effect, levels=c('Length_µm','motileCount_ln_z','inbreeding_z',"Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llzi[mot=='VCL', mot:='Curvilinear (VCL)']
      llzi[mot=='VSL', mot:='Straight-line (VSL)']
      llzi[mot=='VAP', mot:='Average-path (VAP)']
      llzi[, mot := factor(mot, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 
    # controlled for motile count & HL  
      lzhl = list()
      lpzhl =list()
      lprzhl =list()
      a[,motileCount_ln_z := scale(log(motileCount))]
      a[,HL_z := scale(HL)]
      ar[,motileCount_ln_z := scale(log(motileCount))]
      ar[,HL_z := scale(HL)]
      # VAP
         for(i in unique(a$part)){
            #i ='Nucleus'
            ai = a[part == i]
            m = lm(scale(VAP) ~ motileCount_ln_z+ HL_z + Morph+scale(Length_avg), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lzhl[[paste('VAP',i)]]=data.frame(mot = 'VAP', part=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VAP ~ motileCount_ln_z + HL_z + Morph+Length_avg, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z),HL_z = mean(ai$HL_z), Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
              }
            newD=do.call(rbind,nd)

            # values to predict for
            X <- model.matrix(~ motileCount_ln_z + HL_z + Morph+Length_avg,data=newD) # exactly the model which was used has to be specified here 
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
            lpzhl[[paste(i,newD$mot[1])]] = data.table(newD)

            print(paste(i,newD$mot[1]))     
            }          
         for(i in unique(ar$part)){
            #i ='Nucleus'
            if(i == 'Midpiece'){ii = 'Midpiece_rel'}
            if(i == 'Flagellum'){ii = 'Flagellum_rel'}
            ai = ar[part == i]
            m = lm(scale(VAP) ~ motileCount_ln_z + HL_z+ Morph+scale(Length_rel), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lzhl[[paste('VAP',ii)]]=data.frame(mot = 'VAP', part=ii,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VAP ~ motileCount_ln_z +HL_z+ Morph+Length_rel, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z),HL_z = mean(ai$HL_z), Morph = j, Length_rel = seq(min(ai$Length_rel), max(ai$Length_rel), length.out = 200)) 
              }
            newD=do.call(rbind,nd)

            # values to predict for
            X <- model.matrix(~ motileCount_ln_z +HL_z + Morph+Length_rel,data=newD) # exactly the model which was used has to be specified here 
            newD$VAP <-(X%*%v) 
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                            predmatrix[predmatrix < 0] <- 0
                            newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                            newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                            #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
            newD$part=ii
            newD$mot = 'VAP'
            setnames(newD, old = 'VAP', new = 'motility')
            lpzhl[[paste(ii,newD$mot[1])]] = data.table(newD)

            print(paste(ii,newD$mot[1]))     
            }            
      # VSL
         for(i in unique(a$part)){
            #i ='Nucleus'
            ai = a[part == i]
            m = lm(scale(VSL) ~ motileCount_ln_z + HL_z + Morph+scale(Length_avg), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lzhl[[paste('VSL',i)]]=data.frame(mot = 'VSL', part=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VSL ~ motileCount_ln_z + HL_z + Morph+Length_avg, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z), HL_z = mean(ai$HL_z), Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
              }
            newD=do.call(rbind,nd)
            X <- model.matrix(~ motileCount_ln_z +HL_z + Morph+Length_avg,data=newD) # exactly the model which was used has to be specified here 
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
            lpzhl[[paste(i,newD$mot[1])]] = data.table(newD)

            print(paste(i,newD$mot[1]))     
            }          
         for(i in unique(ar$part)){
            #i ='Nucleus'
            if(i == 'Midpiece'){ii = 'Midpiece_rel'}
            if(i == 'Flagellum'){ii = 'Flagellum_rel'}
            ai = ar[part == i]
            m = lm(scale(VSL) ~ motileCount_ln_z+ Morph+scale(Length_rel), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lzhl[[paste('VSL',ii)]]=data.frame(mot = 'VSL', part=ii,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VSL ~ motileCount_ln_z + Morph+Length_rel, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z), HL_z = mean(ai$HL_z),Morph = j, Length_rel = seq(min(ai$Length_rel), max(ai$Length_rel), length.out = 200)) 
              }
            newD=do.call(rbind,nd)

            # values to predict for
            X <- model.matrix(~ motileCount_ln_z + Morph+Length_rel,data=newD) # exactly the model which was used has to be specified here 
            newD$VSL <-(X%*%v) 
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                            predmatrix[predmatrix < 0] <- 0
                            newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                            newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                            #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
            newD$part=ii
            newD$mot = 'VSL'
            setnames(newD, old = 'VSL', new = 'motility')
            lpzhl[[paste(ii,newD$mot[1])]] = data.table(newD)

            print(paste(ii,newD$mot[1]))     
            }            
      # VCL
         for(i in unique(a$part)){
            #i ='Nucleus'
            ai = a[part == i]
            m = lm(scale(VCL) ~ motileCount_ln_z + HL_z + Morph+scale(Length_avg), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lzhl[[paste('VCL',i)]]=data.frame(mot = 'VCL', part=i,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VCL ~ motileCount_ln_z + HL_z + Morph+Length_avg, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z), HL_z = mean(ai$HL_z), Morph = j, Length_avg = seq(min(ai$Length_avg), max(ai$Length_avg), length.out = 200)) 
              }
            newD=do.call(rbind,nd)
            X <- model.matrix(~ motileCount_ln_z + HL_z + Morph+Length_avg,data=newD) # exactly the model which was used has to be specified here 
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
            lpzhl[[paste(i,newD$mot[1])]] = data.table(newD)

            print(paste(i,newD$mot[1]))     
            }          
         for(i in unique(ar$part)){
            #i ='Nucleus'
            if(i == 'Midpiece'){ii = 'Midpiece_rel'}
            if(i == 'Flagellum'){ii = 'Flagellum_rel'}
            ai = ar[part == i]
            m = lm(scale(VCL) ~ motileCount_ln_z+ HL_z + Morph+scale(Length_rel), ai)
            #summary(m)
            #plot(allEffects(m))
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 
            lzhl[[paste('VCL',ii)]]=data.frame(mot = 'VCL', part=ii,effect=rownames(coef(summary(m))),estimate=v, lwr=ci[1,], upr=ci[2,])

            # get predictions
            m = lm(VCL ~ motileCount_ln_z +HL_z + Morph+Length_rel, ai)
            bsim = sim(m, n.sim=nsim) 
            v = apply(bsim@coef, 2, quantile, prob=c(0.5))
            nd = list()
            for(j in unique(ai$Morph)){
              #j=1
              nd[[j]] = data.frame(motileCount_ln_z = mean(ai$motileCount_ln_z), HL_z = mean(ai$HL_z),Morph = j, Length_rel = seq(min(ai$Length_rel), max(ai$Length_rel), length.out = 200)) 
              }
            newD=do.call(rbind,nd)

            # values to predict for
            X <- model.matrix(~ motileCount_ln_z +HL_z + Morph+Length_rel,data=newD) # exactly the model which was used has to be specified here 
            newD$VCL <-(X%*%v) 
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(j in 1:nsim) {predmatrix[,j] <- (X%*%bsim@coef[j,])}
                            predmatrix[predmatrix < 0] <- 0
                            newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                            newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
                            #newD$pred <- apply(predmatrix, 1, quantile, prob=0.5)
            newD$part=ii
            newD$mot = 'VCL'
            setnames(newD, old = 'VCL', new = 'motility')
            lpzhl[[paste(ii,newD$mot[1])]] = data.table(newD)

            print(paste(ii,newD$mot[1]))     
            }          
               
      llzhl = data.table(do.call(rbind,lzhl) ) 
      llzhl[, part := factor(part, levels=rev(c("Acrosome", "Nucleus", "Head", "Midpiece","Tail","Flagellum","Total","Midpiece_rel","Flagellum_rel")))] 
      llzhl[effect == '(Intercept)', effect:='Independent\n(Intercept)']
      llzhl[effect == 'MorphSatellite', effect := 'Satellite\n(relative to Independent)']
      llzhl[effect == 'MorphFaeder', effect := 'Faeder\n(relative to Independent)']
      llzhl[effect == 'scale(Length_avg)', effect := 'Length_µm']
      llzhl[effect == 'scale(Length_rel)', effect := 'Length_µm'] # dummy variable
      
      llzhl[, effect := factor(effect, levels=c('Length_µm','motileCount_ln_z','HL_z',"Faeder\n(relative to Independent)","Satellite\n(relative to Independent)","Independent\n(Intercept)"))] 

      llzhl[mot=='VCL', mot:='Curvilinear (VCL)']
      llzhl[mot=='VSL', mot:='Straight-line (VSL)']
      llzhl[mot=='VAP', mot:='Average-path (VAP)']
      llzhl[, mot := factor(mot, levels=rev(c("Curvilinear (VCL)", "Straight-line (VSL)", "Average-path (VAP)")))] 

#+ effects_mot_morp, fig.width =7/2, fig.height = 2.5                
        g10 = 
        ggplot(llz[effect%in%'Length_µm'], aes(y = part, x = estimate, col = part, shape = mot)) +
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = part), width = 0.1, position = position_dodge(width = 0.5) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.5)) +
          geom_vline(xintercept = 0, col = "red", lty =1)+
          #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          scale_shape(guide = guide_legend(reverse = TRUE, name = 'Motility')) + 
          #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
          labs(y = 'Predictor' , x = "Standardized effect size") +
          guides(color=FALSE, fill = FALSE) +
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
          g10
          ggsave(here::here('Output/morpho_motil_effectSizes_virid_simple_controlled_MotileCount.png'),g10, width = 10, height =7, units = 'cm')
#+ cor, fig.width=10, fig.height = 5
  ga1 = 
  ggplot(aw, aes(x = Acrosome, y = VAP)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
      #scale_colour_manual(values = colors)+
      #scale_fill_manual(values = alpha(colors, 0.4))+
      #geom_point(data = bs_h, aes(x = median, y = Midpiece_µm),position = position_dodge(width = 0.25), col = 'black', size = 2) +
      #ylim(c(20,28)) +
      #xlim(c(24,36)) +
      #annotate(geom="text", x=3, y=15, label="Independent", size = 3, col = colors[1], hjust = 0) +
      #annotate(geom="text", x=3, y=13, label="Satellite", size = 3, col = colors[2], hjust = 0) +
      #annotate(geom="text", x=3, y=11, label="Faeder", size = 3, col = colors[3], hjust = 0) +
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
  ggplot(aw, aes(x = Acrosome, y = VSL)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
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
  ggplot(aw, aes(x = Acrosome, y = VCL)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
     
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
  ggplot(aw, aes(x = Nucleus, y = VAP)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
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
  ggplot(aw, aes(x = Nucleus, y = VSL)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
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
  ggplot(aw, aes(x = Nucleus, y = VCL)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
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
  ggplot(aw, aes(x = Midpiece, y = VAP)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
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
  ggplot(aw, aes(x = Midpiece, y = VSL)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
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
  ggplot(aw, aes(x = Midpiece, y = VCL)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
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
  ggplot(aw, aes(x = Tail, y = VAP)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
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
  ggplot(aw, aes(x = Tail, y = VSL)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
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
  ggplot(aw, aes(x = Tail, y = VCL)) +
      geom_smooth(method = 'rlm') +
      geom_point(pch = 21)+
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
  ggplot(aw, aes(x = Head, y = VAP)) +
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
  ggplot(aw, aes(x = Head, y = VSL)) +
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
  ggplot(aw, aes(x = Head, y = VCL)) +
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
  ggplot(aw, aes(x = Flagellum, y = VAP)) +
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
  ggplot(aw, aes(x = Flagellum, y = VSL)) +
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
  ggplot(aw, aes(x = Flagellum, y = VCL)) +
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
  ggplot(aw, aes(x = Total, y = VAP)) +
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
  ggplot(aw, aes(x = Total, y = VSL)) +
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
  ggplot(aw, aes(x = Total, y = VCL)) +
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
  
  ggsave(here::here('Output/morpho_motil_simple.png'),cbind(
          rbind(ggplotGrob(ga1), ggplotGrob(ga2), ggplotGrob(ga3), size = "first"),
          rbind(ggplotGrob(gn1), ggplotGrob(gn2), ggplotGrob(gn3), size = "first"),
          rbind(ggplotGrob(gm1), ggplotGrob(gm2), ggplotGrob(gm3), size = "first"),
          rbind(ggplotGrob(gt1), ggplotGrob(gt2), ggplotGrob(gt3), size = "first"),
          rbind(ggplotGrob(gh1), ggplotGrob(gh2), ggplotGrob(gh3), size = "first"),
          rbind(ggplotGrob(gf1), ggplotGrob(gf2), ggplotGrob(gf3), size = "first"),
          rbind(ggplotGrob(go1), ggplotGrob(go2), ggplotGrob(go3), size = "first"),
          size = "first"), 
          width = 18, height =8, units = 'cm') 

#' <span style="color: red;">!!! Perhaps a tendency for longer acrosomes making slower sperm and longer other traits making faster sperm !!!</span>  

#' Clustering
  require('apcluster') 
  require('Rcpp')                                                                                                                       
  bw[, ID := .I]
  aw[, ID := .I]

  get_clusters <- function(apc) {
    x <- lapply(apc@clusters, function(x) data.table(ID = names(x) %>% as.numeric()))
    for (i in 1:length(x)) x[[i]][, clust_ID := i]
    rbindlist(x)
  }

  #  all data
    apc = apclusterK
    apc <- apcluster(corSimMat(r = 1), bw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)], q = 0 ,  details = TRUE)
    apc = apclusterK(corSimMat(r = 1), bw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)],  K = 3,  details = TRUE)
    
    apc <- apcluster(corSimMat(r = 1), bw[, .(Head, Midpiece, Tail, Total)], q = 0 ,  details = TRUE)
    apc = apclusterK(corSimMat(r = 1), bw[, .(Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total)],  K = 3,  details = TRUE)


    apc <- apcluster(negDistMat(r=2) , bw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)], q = 0 ,  details = TRUE)
    apc = apclusterK(negDistMat(r = 2), bw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)],  K = 3,  details = TRUE)
    apc <- apcluster(expSimMat(r=2) , bw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)], q = 0 ,  details = TRUE)
    #apc = apclusterK(expSimMat(r = 2), bw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)],  K = 3,  details = TRUE)
    apc <- apcluster(linKernel() , bw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)], q = 0 ,  details = TRUE)

    cc <- get_clusters(apc)

    O <- merge(cc, bw, by = "ID")

    #O[, .N, .(clust_ID, Morph)]

    #dev.new()
    ggplot(O, aes(x = clust_ID, y = Morph)) +
    geom_count()

    plot(apc, bw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)])
    heatmap(apc)

#  averages   
    apc <- apcluster(corSimMat(r = 1), aw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)], q = 0 ,  details = TRUE)
     apc = apclusterK(corSimMat(r = 1), aw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)],  K = 3,  details = TRUE)
    
    apc <- apcluster(corSimMat(r = 1), aw[, .(Head, Midpiece, Tail, Total)], q = 0 ,  details = TRUE)
    apc <- apclusterK(corSimMat(r = 1), aw[, .(Head, Midpiece, Tail, Total)], K = 3 ,  details = TRUE)

    apc <- apcluster(corSimMat(r = 1), aw[, .(VAP, VSL, VCL)], q = 0 ,  details = TRUE)
    apc <- apclusterK(corSimMat(r = 1), aw[, .(VAP, VSL, VCL)], K = 3 ,  details = TRUE)

    apc <- apcluster(negDistMat(r=2) , aw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)], q = 0.5 ,  details = TRUE)
     apc = apclusterK(negDistMat(r = 2), aw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)],  K = 3,  details = TRUE)
    
    apc <- apcluster(expSimMat(r=2) , aw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)], q = 0 ,  details = TRUE)
     
     apc = apclusterK(expSimMat(r = 2), aw[, .(VAP, VSL, VCL, Acrosome, Flagellum, Head, Midpiece, Nucleus, Tail, Total, Midpiece_rel, Flagellum_rel)],  K = 3,  details = TRUE)

    cc <- get_clusters(apc)

    O <- merge(cc, aw, by = "ID")

    O[, .N, .(clust_ID, Morph)]

    #dev.new()
    ggplot(O, aes(x = clust_ID, y = Morph)) +
    geom_count()

    plot(apc, aw[, .(Head, Midpiece, Tail, Total)])
    heatmap(apc)
# End     