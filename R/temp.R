require(data.table)
require(ggplot2)

d = fread('Data/results_inv_KT_test.csv')
d2 = fread('Data/results_gre_KT_test.csv')

names(d) = c('pic','inv')
names(d2) = c('pic','gre')
merge(d,d2)