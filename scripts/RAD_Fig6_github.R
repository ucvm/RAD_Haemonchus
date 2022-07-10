library("tidyr")  
library("ggplot2")
library("patchwork")
library("dplyr")
library("plyr")
library("cowplot")
library("grid")


###############################
# IMPORTANT!!!!!
# Download RAD_Haemonchus folder and as set working directory
setwd("~/Desktop/RAD_Haemonchus/")
###############################


## Data
assign("LongData",get(load("data/PooledFstData.rda")))

LongData %>%
  dplyr::group_by(Pops)

#Aligning X axis across comparisons
chrom_length <- data.frame("Chr" = c("Chr1","Chr2","Chr3","Chr4","Chr5","ChrX"), "length" = c(45800,47400,43600,51800,48900,46000),
                           "start"= c(0,45800,93200,136800,188600,237500))

ToPlot <- LongData %>%
  dplyr::left_join(chrom_length, by = "Chr")  %>%
  dplyr::group_by(Pops,Chr) %>%
  dplyr::arrange(Chr) %>%
  dplyr::mutate(newX = row_number()) %>%
  dplyr::mutate(normXC = ((newX/max(newX)*length)+start)/7) %>%
  dplyr::mutate(Colour = as.factor(Colour))

#Make beta-tub isotype-1 data frame
Isotypes = LongData[c(1:4,4471:4474,8942:8945,13413:13416,17884:17887,22355:22358),c(1,2,4)]
Positions = c(7027095) #,5611669,1273069 13433296,
Chromosomes = c("Chr1") #,"ChrX","ChrX""Chr2",
Populations = c(1,5,9,13,17,21)
PopNames = c("S1.S2","S1.R1","S1.R2","S2.R1","S2.R2","R1.R2")
for(a in (c(1:4))){
  add = a-1
  Isotypes[(Populations + add),"Pos"] = Positions[a]
  Isotypes[(Populations + add),"Chr"] = Chromosomes[a]
}
for(b in c(1:6)){
  Isotypes[c((Populations[b]):(Populations[b] + 3)),"Pops"] = PopNames[b]
}
Isotypes <- Isotypes %>% na.omit(Pos)

PlotIsotypes <- Isotypes %>%
  dplyr::left_join(chrom_length, by = "Chr")  %>%
  dplyr::group_by(Pops,Chr) %>%
  dplyr::arrange(Chr) %>%
  dplyr::mutate(newX = row_number()) %>%
  dplyr::mutate(normXC = (Pos/7074951*924))

#Plot theme
extra_RAD <- ggplot2::theme(panel.background = element_rect(fill = NA, colour = "black", size = 0.5),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            legend.position="none", 
                            axis.text.y = element_text(size = 6),
                            axis.title.x = element_text( size=10, face = "bold"), 
                            axis.title.y = element_text( size=8, face = "bold", vjust = 1.5), 
                            title = element_text(size=8, face ="bold"),
                            strip.text = element_text(colour = 'black'),
                            axis.ticks = element_line(colour = "black", size = 0.2), 
                            axis.ticks.length = unit(0.05, "cm"))
pointsize = 0.05
cols <- c("1"=0.1, "2" = 1)
sizes = c("1"=0.5, "2" = .5)

# Making plots
plots <- 1
plot_list = list()

for(i in c(1:4)){
  if(i == 1){
    comp <- "S1.R1"
    lab = "T1-U1"
  }
  if(i == 2){
    comp <- "S2.R1"
    lab = "T1-U2"
  }
  if(i == 3){
    comp <- "S1.R2"
    lab = "T2-U1"
  }
  if(i == 4){
    comp <- "S2.R2"
    lab = "T2-U2"
  }
  plot_list[[plots]] <- ToPlot[ToPlot$Pops == comp,] %>%
  ggplot(aes(x = normXC, y = Fst, axes = FALSE, colour = Chr, alpha = Colour, size = Colour)) +
  geom_vline(data = PlotIsotypes, aes(xintercept = normXC), color = "blue", linetype="dashed", size=0.1) +
  geom_point(size = pointsize) + #shape = "."
  scale_x_continuous(name = NULL, breaks=NULL, labels = NULL, expand =c(0.02,0.02)) +
  scale_y_continuous(name = NULL, breaks=c(0,.25,.5,.75,1), labels = c("0.00","0.25","0.50","0.75","1.00"),limits = c(0,1), expand =c(0.02,0.02)) + #name = expression("U1 - F"[ST]), 
  scale_colour_manual(values = c("red","maroon","mediumpurple4","steelblue","turquoise4","mediumseagreen"),  
                      aesthetics = c("colour", "fill")) +
  scale_alpha_manual(values = cols) +
  extra_RAD + 
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.4), "cm"))
    plots= plots+1
}

## Adding labelling
T1U1 <- plot_grid(NULL,plot_list[[1]], rel_widths = c(0.02,1)) + 
  draw_label(expression(bold(paste("F"[ST],sep = ""))), 
             color = "black", size = 10, angle = 90, x = 0.015, y = 0.5)
T1U2 <- plot_grid(NULL,plot_list[[2]], rel_widths = c(0.02,1)) + 
  draw_label(expression(bold(paste("F"[ST],sep = ""))), 
             color = "black", size = 10, angle = 90, x = 0.015, y = 0.5)
T2U1 <- plot_grid(NULL,plot_list[[3]], rel_widths = c(0.02,1)) + 
  draw_label(expression(bold(paste("F"[ST],sep = ""))), 
             color = "black", size = 10, angle = 90, x = 0.015, y = 0.5)
T2U2 <- plot_grid(NULL,plot_list[[4]], rel_widths = c(0.02,1)) + 
  draw_label(expression(bold(paste("F"[ST],sep = ""))), 
             color = "black", size = 10, angle = 90, x = 0.015, y = 0.5)

## Combining all pairwise comparisons
all_plots <- plot_grid(NULL,T1U1,NULL,NULL,T1U2,NULL,NULL,T2U1,NULL,NULL,T2U2,NULL,ncol=3, rel_widths = c(0.02,1,0,0.02,1,0,0.02,1,0,0.02,1,0),align = "hv", 
             labels = c("","A","","","B","","","C","","","D",""), label_size =10, label_x = -0.01) +
  theme(plot.background = element_rect(fill="white", color = NA))

## Adding more r0om below X axis and labels
labroom <- plot_grid(all_plots,NULL,ncol =1, rel_heights = c(1,.025))+
  theme(plot.background = element_rect(fill="white", color = NA))

withlabel <- plot_grid(labroom) +
  draw_label("Genomic position", fontface ="bold",
            color = "black", size = 10, angle = 0, x = 0.52, y = 0.015) 

rect <- rectGrob(
  x = .85,
  y = c(.235,0.722,0.965,.478),
  width = unit(0.6, "in"),
  height = unit(0.3, "in"),
  gp = gpar(fill = "white", alpha=1)
)

labels <- ggdraw(withlabel) +
  draw_grob(rect) +
  draw_text(c("T1-U1","T1-U2","T2-U1","T2-U2"), x = 0.85, y= c(.965,0.722,.478,0.235),fontface = "bold", color = "black", size =10)

##Save figure 6
#ggsave2("Wit_ea_Fig6.png", labels , height = 7.5, width = 6)
ggsave2("Wit_ea_Fig6.pdf", labels , height = 7.5, width = 5)


  