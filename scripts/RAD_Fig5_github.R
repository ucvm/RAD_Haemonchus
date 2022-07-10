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
assign("all_data", get(load("data/AllFstData.rda")))


all_data <- all_data[all_data$Chr =="Chr1",]
all_data$Chr <- factor(all_data$Chr, levels = c("Chr1","Chr2","Chr3","Chr4","Chr5","ChrX"))
all_data$Pop <- factor(all_data$Pop , levels = c("R1S1", "R1S2","R2S1","R2S2","R1R2","S1S2"))
all_data$Analysis <- factor(all_data$Analysis , levels = c("Read_2M","Read_1.5M","Read_1M","Read_.5M", "Read_.25M", "Read_0.1M",
                                                      "Indiv_Stringent", "Indiv_All","Indiv_15","Indiv_10","Indiv_5"))
cols <- c("1"=0.1, "2" = 1)
sizes = c("1"=0.5, "2" = .5)
all_data$Colour <- as.factor(all_data$Colour)

#Plot theme

Fig5 <- ggplot2::theme(panel.background = element_rect(fill = NA, colour = "black", size = 0.5),
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

##Plot all data

allplots <- 1
allplot_list = list()

for(k in c(1:11)){
  if(k == 1){
    sub <- "Read_2M"
  }
  if(k == 2){
    sub <- "Read_1.5M"
  }
  if(k == 3){
    sub <- "Read_1M"
  }
  if(k == 4){
    sub <- "Read_.5M"
  }
  if(k == 5){
    sub <- "Read_.25M"
  }
  if(k == 6){
    sub <- "Read_0.1M"
  }
  if(k == 7){
    sub <-  "Indiv_Stringent"
  }
  if(k == 8){
    sub <- "Indiv_All"
  }
  if(k == 9){
    sub <- "Indiv_15"
  }
  if(k == 10){
    sub <- "Indiv_10"
  }
  if(k == 11){
    sub <- "Indiv_5"
  }
  allplot_list[[allplots]] <- all_data %>%
  dplyr::filter(Pop == "R1S1") %>%
  dplyr::filter(Analysis == sub) %>%
  ggplot( aes(x = BP, y = Fst, axes = FALSE, colour = Chr, alpha = Colour, size = Colour)) +
  geom_point(size = pointsize) + 
  Fig5 +
  ylim(0,1) +
  scale_x_continuous(NULL , breaks= NULL) + 
  labs(y = NULL) +
  theme(legend.position="none", axis.title.x = element_text( size=10), axis.title.y = element_text( size=10), title = element_text(size=8, face ="bold")) + #face="bold",
  scale_colour_manual(values = c("red","maroon","mediumpurple4","steelblue","turquoise4","mediumseagreen"),  aesthetics = c("colour", "fill")) +
  scale_alpha_manual(values = cols) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.4), "cm"))
  allplots = allplots+1
}

##Figure extra's
triangle <- cowplot::ggdraw() + cowplot::draw_image("~/Desktop/DESKTOP_Cleanup/Bureaublad_Aug21/20210809_Triangle3.png",scale = .98) # scale = x = 1, y = 1, hjust = 1, vjust = 1, width = 0.13, height = 0.2
rectLab <- rectGrob(
  x = c(.437,.935),
  y = .96,
  width = unit(.95, "in"),
  height = unit(0.3, "in")
)

#Make figure A
MakeA <- plot_grid(allplot_list[[7]],allplot_list[[8]],allplot_list[[9]],allplot_list[[10]],allplot_list[[11]],align = "hv", 
                   axis = "b",rel_heights = c(1,1,1), ncol =1, greedy = FALSE) +
  theme(plot.background = element_rect(fill="white", color = NA))

MakeAfull <- plot_grid(NULL,MakeA,triangle, rel_widths = c(0.02666,1.333,0.45), ncol = 3) + 
  draw_label(expression(bold(paste("Genetic differentiation (F"[ST],")",sep = ""))), 
             color = "black", size = 8, angle = 90, x = 0.015, y = 0.5)+
  theme(plot.background = element_rect(fill="white", color = NA))

#Make figure B
MakeB <- plot_grid(allplot_list[[1]],allplot_list[[2]],allplot_list[[3]],allplot_list[[4]],allplot_list[[5]],allplot_list[[6]],align = "hv", 
                      axis = "b",rel_heights = c(1,1,1), ncol =1, greedy = FALSE) +
  theme(plot.background = element_rect(fill="white", color = NA))

MakeBfull <- plot_grid(MakeB,triangle, rel_widths = c(1.333,0.45)) +
  theme(plot.background = element_rect(fill="white", color = NA))

#Make full figure and add labels

plot_FIVE <- plot_grid(MakeAfull,MakeBfull, rel_heights = c(1,1), ncol=2, labels = c("A","B"), label_size =10, label_x = -0.01) + 
  draw_grob(rectLab) +
  draw_label("No. of individuals", fontface = "bold", color = "black", size = 8, x = 0.437, y = 0.97)+
  draw_label("Read depth", fontface = "italic", color = "black", size = 8, x = 0.437, y = 0.945)+
  draw_label("20*", color = "black", size = 7.5, x = 0.437, y = 0.9, fontface = "bold") +
  draw_label("9.1", color = "black", size = 7.5, x = 0.437, y = 0.87, fontface = "italic") +
  draw_label("20", color = "black", size = 7.5, x = 0.437, y = .7, fontface = "bold")+
  draw_label("9.1", color = "black", size = 7.5, x = 0.437, y = .67, fontface = "italic")+
  draw_label("15", color = "black", size = 7.5, x = 0.437, y = .5, fontface = "bold")+
  draw_label("9.0", color = "black", size = 7.5, x = 0.437, y = .47, fontface = "italic")+
  draw_label("10", color = "black", size = 7.5, x = 0.437, y = .3, fontface = "bold")+
  draw_label("9.4", color = "black", size = 7.5, x = 0.437, y = .27, fontface = "italic")+
  draw_label("5", color = "black", size = 7.5, x = 0.437, y = .1, fontface = "bold") +
  draw_label("9.0", color = "black", size = 7.5, x = 0.437, y = .07, fontface = "italic") +
  draw_label("No. of reads (M)", fontface = "bold", color = "black", size = 8, x =.935, y = 0.97)+
  draw_label("Read depth", fontface = "italic", color = "black", size = 8, x = 0.935, y = 0.945)+
  draw_label("2.00", color = "black", size = 7.5, x = .935, y = 0.91, fontface = "bold") +
  draw_label("7.4", color = "black", size = 7.5, x = .935, y = 0.88, fontface= "italic") +
  draw_label("1.50", color = "black", size = 7.5, x = .935, y = .75, fontface = "bold")+
  draw_label("6.3", color = "black", size = 7.5, x = .935, y = .72, fontface= "italic")+
  draw_label("1.00", color = "black", size = 7.5, x = .935, y = .585, fontface = "bold")+
  draw_label("5.0", color = "black", size = 7.5, x = .935, y = .555, fontface= "italic")+
  draw_label("0.50", color = "black", size = 7.5, x = .935, y = .42, fontface = "bold")+
  draw_label("3.3", color = "black", size = 7.5, x = .935, y = .39, fontface= "italic")+
  draw_label("0.25", color = "black", size = 7.5, x = .935, y = .255, fontface = "bold")+
  draw_label("2.3", color = "black", size = 7.5, x = .935, y = .225, fontface= "italic")+
  draw_label("0.10", color = "black", size = 7.5, x = .935, y = .09, fontface = "bold")+
  draw_label("1.6", color = "black", size = 7.5, x = .935, y = .06, fontface= "italic")+
  theme(plot.background = element_rect(fill="white", color = NA))

five_LABEL <- plot_grid(plot_FIVE,NULL,ncol=1, rel_heights = c(1,0.03)) + 
  draw_label("Genomic position (Chr I)", fontface ="bold",
             color = "black", size = 8, angle = 0, x = 0.21, y = 0.02) +
  draw_label("Genomic position (Chr I)", fontface ="bold",
             color = "black", size = 8, angle = 0, x = 0.71, y = 0.02) +
  theme(plot.background = element_rect(fill="white", color = NA))

#Save figure 5
#ggsave2("Wit_ea_Fig5.png", five_LABEL, height = 5, width = 7.5)
ggsave2("Wit_ea_Fig5.pdf", five_LABEL, height = 5, width = 7.5)
