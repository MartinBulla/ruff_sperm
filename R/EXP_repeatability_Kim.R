# TOOLS & DATA
  require(here)
  source(here::here('R/tools.R'))
  library(ggpubr)
  require(arm)
  require(rptR) 

rpt.corr <- function(y, groups, CI=0.95, nboot=1000, npermut=1000) {
    # initial checks
    if(length(y)!= length(groups)) 
        stop("y and group have to be of equal length")
    if(nboot < 0)   nboot <- 0
    if(npermut < 1) npermut <- 1
    if(any(is.na(y))) 
        stop("missing values in y ")
    if(any(is.na(groups)))
        stop("missing values in groups ")
    if(!all(table(groups)==2))     
        stop("not exactly two data points per group")
    # preparation
    sortix <- sort.int(as.numeric(groups),index.return=TRUE)$ix
    y1     <- y[sortix][seq(1,length(y),by=2)]
    y2     <- y[sortix][seq(2,length(y),by=2)]
    k      <- length(unique(groups))
    # functions: point estimates of R
    R.pe  <- function(y1, y2, k) {
        y <- c(y1, y2)
        R <- 1/(k-1)*sum((y1-mean(y))*(y2-mean(y))) / var(y)
        return(R) 
    }   
    # point estimation according to equations 4 and 5
    R      <- R.pe(y1, y2, k)
    # confidence interval estimation according to equation 6 and 7
    bootstr <- function(y1, y2, k) {
        samp<- sample(1:k, k, replace=TRUE)
        R.pe(y1[samp], y2[samp], k)
    }
    if(nboot > 0)
        R.boot <- replicate(nboot, bootstr(y1, y2, k), simplify=TRUE) 
    else
        R.boot <- NA
    CI.R   <- quantile(R.boot, c((1-CI)/2, 1-(1-CI)/2), na.rm=TRUE)
    se     <- sd(R.boot, na.rm=TRUE)
    # significance test by permutation
    permut   <- function(y1, y2, k) {
        samp <- sample(1:k, k)
        R.pe(y1[samp], y2, k)
    }
    if(npermut > 1) {
        R.permut <- c(R, replicate(npermut-1, permut(y1, y2, k), simplify=TRUE))
        P.permut <- sum(R.permut >= R)/npermut
    }
    else {
        R.permut <- R
        P.permut <- NA
    }
    # return of results
    res <- list(call=match.call(), datatype="Gaussian", method="corr", CI=CI, 
                R=R, se=se, CI.R=CI.R, P=P.permut, R.permut=R.permut) 
    class(res) <- "rpt"
    return(res) 
}

    # functions: point estimates of R
    R.pe  <- function(y1, y2, k) {
        y <- c(y1, y2)
        R <- 1/(k-1)*sum((y1-mean(y))*(y2-mean(y))) / var(y)
        return(R) 
    }   
    # point estimation according to equations 4 and 5
   

  d = fread('Data/Compare_MH_vs_TL.csv')
  d1 = d[,.(sperm, part, length_A)]
  setnames(d1, c('sperm_ID','part','length'))
  d1$measur = 'A'
  d2 = d[,.(sperm, part, length_B)]
  setnames(d2, c('sperm_ID','part','length'))
  d2$measur = 'B'
  dd = rbind(d1,d2)


  ggplot(d, aes(x = length_A, y = length_B)) + 
    stat_smooth(method = 'lm') + geom_point() + 
    stat_cor(method="pearson") + 
    facet_wrap(~part, scales = 'free') + 
    theme_light() + theme(axis.text.x=element_text(angle = 45, hjust=1))

  

  ggplot(dd, aes(x = measur, y = length, fill = measur)) + geom_boxplot() +facet_wrap(~part, scales = 'free') + theme_light()
  
  ggplot(dd, aes(x = measur, y = length, fill = measur)) + 
    facet_wrap(~part, scales = 'free') + 
    theme_light()+
    geom_dotplot(binaxis = 'y', stackdir = 'center',position = position_dodge())

  ggplot(dd[part == 'head' & measur == 'A'], aes(x = measur, y = length)) + 
    theme_light()+
    geom_dotplot(binaxis = 'y', stackdir = 'center',position = position_dodge())

  m = lmer(length ~ 1+(1|measur), dd[part == 'head'])
    l=data.frame(summary(m)$varcor)
    l=l[is.na(l$var2),]
    ri1=data.frame(model='HEAD', type='% of var',effect=l$grp, estimate_r=round(100*l$vcov/sum(l$vcov)))


  m = lmer(length ~ 1+(1|sperm_ID), dd[part == 'head'])
    l=data.frame(summary(m)$varcor)
    l=l[is.na(l$var2),]
    ri1=data.frame(model='HEAD', type='% of var',effect=l$grp, estimate_r=round(100*l$vcov/sum(l$vcov)))

  h = d[part == 'head']  
  R = R.pe(y1 = h$length_A, y2 = h$length_B, k = length(unique(h$sperm)))
    
  rpt(length ~ (1 | sperm_ID), grname = "sperm_ID", data = dd[part == 'head'], datatype = "Gaussian", nboot = 0, npermut = 0)

  rpt(length ~ (1 | measur), grname = "measur", data = dd[part == 'head'], datatype = "Gaussian", nboot = 0, npermut = 0)

  m = lmer(length ~ 1+(1|measur), dd[part == 'midpiece'])
    l=data.frame(summary(m)$varcor)
    l=l[is.na(l$var2),]
    ri1=data.frame(model='midpiece', type='% of var',effect=l$grp, estimate_r=round(100*l$vcov/sum(l$vcov)))

      
  m = lmer(length ~ 1+(1|sperm_ID), dd[part == 'midpiece'])
    l=data.frame(summary(m)$varcor)
    l=l[is.na(l$var2),]
    ri2=data.frame(model='midpiece', type='% of var',effect=l$grp, estimate_r=round(100*l$vcov/sum(l$vcov)))

  h = d[part == 'midpiece']  
  R.pe(y1 = h$length_A, y2 = h$length_B, k = length(unique(h$sperm)))

  rpt(length ~ (1 | measur), grname = "measur", data = dd[part == 'midpiece'], datatype = "Gaussian", nboot = 0, npermut = 0)
   
    m = lmer(length ~ 1+(1|measur), dd[part == 'tail'])
    l=data.frame(summary(m)$varcor)
    l=l[is.na(l$var2),]
    ri1=data.frame(model='midpiece', type='% of var',effect=l$grp, estimate_r=round(100*l$vcov/sum(l$vcov)))
     
  m = lmer(length ~ 1+(1|sperm_ID), dd[part == 'tail'])
    l=data.frame(summary(m)$varcor)
    l=l[is.na(l$var2),]
    ri3=data.frame(model='tail', type='% of var',effect=l$grp, estimate_r=round(100*l$vcov/sum(l$vcov)))    
    
   h = d[part == 'tail']  
  R.pe(y1 = h$length_A, y2 = h$length_B, k = length(unique(h$sperm)))  

  rpt(length ~ (1 | measur), grname = "measur", data = dd[part == 'tail'], datatype = "Gaussian", nboot = 0, npermut = 0)
  rbind(ri1,ri2,ri3)    