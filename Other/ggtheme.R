### A clean theme that gets rid of unnecessary gridlines and colored backgrounds.

require(ggplot2)
ggtheme = theme(legend.background = element_blank(), 
                legend.key = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.grid.major = element_blank(), 
                panel.background = element_blank(), 
                panel.border = element_blank(), 
                strip.background = element_blank(), 
                plot.background = element_blank())
detach(package:ggplot2)