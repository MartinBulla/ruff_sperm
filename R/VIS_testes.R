# TOOLS & DATA
  require(here)
  source(here::here('R/tools.R'))
  require(gridExtra)
  require(viridis)

# load data  
  d = data.table(read_excel('Data/GSI_for_Liam.xlsx', sheet = 1))#, range = "A1:G161"))
  d[, Morph := factor(Morph, levels=c("Res", "Sat", "Faed"))] 
  d[Morph == 'Res', Morph := 'Independent']
  d[Morph == 'Sat', Morph := 'Satellite']
  d[Morph == 'Faed', Morph := 'Faeder']

# plot
  g1 = ggplot(d, aes(x = Morph, y = Gonadmass)) + 
    geom_boxplot() + 
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(), col = 'darkgrey', aes(fill =Morph))+
    scale_fill_viridis(discrete=TRUE)+
    ylab('Gonad Mass (g)') +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none")

  g2 = ggplot(d, aes(x = Morph, y = GSI)) + 
    geom_boxplot() +
    geom_dotplot(binaxis = 'y', stackdir = 'center',
                 position = position_dodge(), col = 'darkgrey', aes(fill =Morph))+
    scale_fill_viridis(discrete=TRUE)+
    ylab('Gonadosomatic Index') +
    theme_bw()+theme(legend.position = "none")  

  grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), size = "last"))
  #grid.arrange(g1,g2)
  gg1 <- ggplotGrob(g1)
  gg2 <- ggplotGrob(g2) 

  ggsave('Output/testes.png',rbind(gg1,gg2, size = "last"), width = 7, height =10, units = 'cm')  

  # End