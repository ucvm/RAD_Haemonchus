library(png)
library(grid)
library(gridExtra)
library(cowplot)

###############################
# IMPORTANT!!!!!
# Download RAD_Haemonchus folder and as set working directory
setwd("~/Desktop/RAD_Haemonchus/")
###############################

LEFT <- rasterGrob(readPNG("data/left_fig_1.png", native = FALSE),
           interpolate = FALSE)
RIGHT <- rasterGrob(readPNG("data/pie.png", native = FALSE),
                    interpolate = FALSE)

LIST <- list()
LIST[[1]] <- LEFT
LIST[[2]] <- RIGHT


#Combining left and right

total <- plot_grid(LIST[[1]],LIST[[2]], align = "h", nrow = 1,  rel_widths = c(2.5,1))#, rel_heights = c(1/4, 1/4, 1/2))
ggsave("Wit_ea_Fig2.pdf", total, width = 7.5, height =3.5 )
