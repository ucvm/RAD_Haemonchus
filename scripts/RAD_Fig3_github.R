library("tidyr")  
library("ggplot2")
library("patchwork")
library("dplyr")
library("RcppRoll")
library("cowplot")
library("grid")

###############################
# IMPORTANT!!!!!
# Download RAD_Haemonchus folder and as set working directory
setwd("~/Desktop/RAD_Haemonchus/")
###############################


## Data
AllData <- read.table(file= "data/FstData.tsv", header = TRUE) #Fst
Divs <- read.table(file= "data/ExpHetData.tsv", header = TRUE) #Exp Het


##Fst
#Location of beta tub isotypes for plotting
Length = c(1:24)
Isotypes = AllData[Length,c(3,4)]
Positions = c(13433296,7027095,5611669,1273069)
Chromosomes = c("Chr2","Chr1","ChrX","ChrX")
Populations = c(1,5,9,13,17,21)
PopNames = c("R1S1","R2S1","R1S2","R2S2","S1S2","R1R2")
for(a in (c(1:4))){
  add = a-1
  Isotypes[(Populations + add),"Pos"] = Positions[a]
  Isotypes[(Populations + add),"Chr"] = Chromosomes[a]
}
for(b in c(1:6)){
  Isotypes[c((Populations[b]):(Populations[b] + 3)),"Pop"] = PopNames[b]
}
Isotypes$Pop <- factor(Isotypes$Pop, levels = c("R1S1", "R1S2","R2S1","R2S2","R1R2","S1S2"))

AllData$Chr <- factor(AllData$Chr, levels = c("Chr1","Chr2","Chr3","Chr4","Chr5","ChrX"))
AllData$Pop <- factor(AllData$Pop , levels = c("R1S1", "R1S2","R2S1","R2S2","R1R2","S1S2"))
labelsPop = c("R1S1" = "T1-U1","R2S1" = "T2-U1", "R1S2" = "T1-U2","R2S2" = "T2-U2","S1S2" = "U1-U2", "R1R2" = "T1-T2")
labelsChrom = c("Chr1"= "I","Chr2" = "II","Chr3" = "III","Chr4" = "IV","Chr5" ="V","ChrX" = "X")

cols <- c("1"=0.1, "2" = 1, "iso" = 1)
sizes = c("1"=0.8, "2" = .8, "iso" = 1)
AllData$Colour <- as.factor(AllData$Colour)

#Exp Het

ChromCol <- c("Chr1" = "red","Chr2" = "maroon","Chr3" = "mediumpurple4","Chr4" = "steelblue","Chr5" = "turquoise4","ChrX" = "mediumseagreen")
Divs_cols <- c("1"=0.02, "2" = 1) 
Divs_sizes = c("1"=0.8, "2" = 1) 
Divs$Chr <- factor(Divs$Chr , levels = c("Chr1","Chr2","Chr3","Chr4","Chr5","ChrX"))
Divs$PopID <- factor(Divs$PopID , levels = c("Res1", "Res2","Susc1", "Susc2"))
Divs$Colour <- factor(Divs$Colour , levels = c("1", "2"))
labelsPopID = c("Res1" = "T1","Res2" = "T2", "Susc1" = "U1","Susc2" = "U2")


## Ensuring equal length of chromosomes across populations
chrom_length <- data.frame("Chr" = c("Chr1","Chr2","Chr3","Chr4","Chr5","ChrX"), "length" = c(45800,47400,43600,51800,48900,46000),
                          "start"= c(0,45800,93200,136800,188600,237500))

AllData %>%
  dplyr::group_by(Pop) %>%
  dplyr::summarise(markers = n())

Divs %>%
  dplyr::group_by(PopID) %>%
  dplyr::summarise(markers = n())

lengths_div <- Divs %>%
  dplyr::group_by(PopID) %>%
  dplyr::arrange(Chr) %>%
  dplyr::mutate(newX = row_number()) %>%
  dplyr::summarise(max = max(newX))

lengths_all <- AllData %>%
  dplyr::group_by(Pop) %>%
  dplyr::arrange(Chr) %>%
  dplyr::mutate(newX = row_number()) %>%
  dplyr::summarise(max = max(newX))

plot_Divs <- Divs %>%
  dplyr::left_join(chrom_length, by = "Chr")  %>%
  dplyr::group_by(PopID, Chr) %>%
  dplyr::arrange(Chr) %>%
  dplyr::mutate(newX = row_number()) %>%
  dplyr::mutate(normX = ((newX/max(newX)*length)+start)/7)

plot_All <- AllData %>%
  dplyr::left_join(chrom_length, by = "Chr")  %>%
  dplyr::group_by(Pop,Chr) %>%
  dplyr::arrange(Chr) %>%
  dplyr::mutate(newX = row_number()) %>%
  dplyr::mutate(normX = ((newX/max(newX)*length)+start)/7)

isotypes_all <- plot_All %>%
  dplyr::filter(Chr %in% c("Chr1","Chr2","ChrX"),
                BP %in% c(7021619,13430265,5611464,1275896)) %>%
  dplyr::ungroup() %>%
  dplyr::select(Pop,Chr,BP,normX) %>%#,1275896 
  dplyr::left_join(Isotypes, by = c("Chr","Pop")) %>%
  dplyr::mutate(correct = normX) %>%
  dplyr::mutate(normX = BP/Pos*correct) %>%
  dplyr::filter(!(Chr == "ChrX" & normX <34000)) %>%
  dplyr::filter(!(Chr == "ChrX" & normX >36000)) %>%
  dplyr::mutate(Colour = as.factor("iso"))

Isotype2 <- Isotypes %>%
  dplyr::mutate(PopID = case_when(Pop == "R1S1" ~ "Susc1",
                                Pop == "R2S1" ~ "Susc2",
                                Pop == "R1S2" ~ "Res1",
                                Pop == "R2S2" ~ "Res2",
                                Pop %in% c("S1S2","R1R2") ~ "None")) %>%
  dplyr::filter(!PopID == "None")

isotypes_divs <- plot_Divs  %>%
  dplyr::filter(Chr %in% c("Chr1","Chr2","ChrX"),
                BP %in% c(7028677,13434326,5611637,1281530)) %>%
  dplyr::ungroup() %>%
  dplyr::select(PopID,Chr,BP,normX) %>%#,1275896 
  dplyr::left_join(Isotype2, by = c("Chr","PopID")) %>%
  dplyr::mutate(correct = normX) %>%
  dplyr::mutate(normX = BP/Pos*correct) %>%
  dplyr::filter(!(Chr == "ChrX" & normX <34000)) %>%
  dplyr::filter(!(Chr == "ChrX" & normX >36000)) %>%
  dplyr::mutate(Colour = as.factor("iso"))

##Making figures
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


##Fst loop##
plots <- 1
plot_list = list()


for(i in c(1:6)){
  if(i == 1){
    comp <- "R1S1"
    lab = "T1-U1"
  }
  if(i == 2){
    comp <- "R1S2"
    lab = "T1-U2"
  }
  if(i == 3){
    comp <- "R2S1"
    lab = "T2-U1"
  }
  if(i == 4){
    comp <- "R2S2"
    lab = "T2-U2"
  }
  if(i == 5){
    comp <- "R1R2"
    lab = "T1-T2"
  }
  if(i == 6){
    comp <- "S1S2"
    lab = "U1-U2"
  }
  plot_list[[plots]] <- plot_All %>%
    dplyr::filter(Pop == comp) %>%
    ggplot(aes(x = normX, y = Fst, axes = FALSE, colour = Chr, alpha = Colour, size = Colour)) +
    geom_point(shape = ".") +
    scale_x_continuous(name = NULL, breaks=NULL, labels = NULL, expand =c(0.02,0.02)) +
    scale_y_continuous(name = NULL, breaks=c(0,.25,.5,.75,1), labels = c("0.00","0.25","0.50","0.75","1.00"), expand =c(0.02,0.02)) + #name = expression("U1 - F"[ST]), 
    scale_colour_manual(values = c("red","maroon","mediumpurple4","steelblue","turquoise4","mediumseagreen"),  aesthetics = c("colour", "fill")) +
    scale_alpha_manual(values = cols) +
    scale_size_manual(values = sizes) +
    #geom_vline(data = isotypes_all[isotypes_all$Pop == comp,], aes(xintercept = normX), colour = "blue", linetype = 2, size =0.2) +
    geom_segment(data = isotypes_all[isotypes_all$Pop == comp,], aes(x = normX, xend = normX, y = 0, yend=1), colour = "blue", linetype = 2, size =0.2) +
    extra_RAD +
   # geom_text(aes(x=30000, y=0.85), label = lab, color = "black", size = 8 ) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.4), "cm"))
  plots= plots+1
}

##Exp het loop##
divplots <- 1
divplot_list = list()


for(j in c(1:6)){
  if(j == 1){
    pop <- "Res1"
  }
  if(j == 2){
    pop <- "Res2"
  }
  if(j == 3){
    pop <- "Susc1"
  }
  if(j == 4){
    pop <- "Susc2"
  }
  ROLL <- plot_Divs[plot_Divs$PopID==pop,] 
  ROLL$rolling <- roll_mean(ROLL$Exp_Het,n=500, align = "center",fill=NA, by=20)
  ROLL <- ROLL %>%
    dplyr::filter(!is.na(rolling))
  divplot_list[[divplots]] <- plot_Divs %>%
    dplyr::filter(PopID == c(pop)) %>%
    ggplot(aes(x = normX, y = Exp_Het, axes = FALSE, colour = Chr, alpha = Colour, size = Colour)) +
    geom_point(shape = ".") + 
    geom_line(data = ROLL,aes(x = normX, y = rolling,colour = Chr), size =0.4, alpha =1) +
    scale_x_continuous(name = NULL, breaks=NULL, labels = NULL, expand =c(0.02,0.02)) +
    scale_y_continuous(name = NULL, breaks=c(0,0.25,0.5), labels = c("0.00","0.25","0.50"), limits = c(0,0.52632), expand =c(0.02,0.02)) + #expression(paste("Diversity (", pi, ")", sep =""))
    scale_colour_manual(values = c("red","maroon","mediumpurple4","steelblue","turquoise4","mediumseagreen"),  aesthetics = c("colour", "fill")) +
    scale_alpha_manual(values = cols) +
    geom_segment(data = isotypes_divs[isotypes_divs$PopID == pop,], aes(x = normX, xend = normX, y = 0.15, yend=0.52), colour = "blue", linetype = 2, size =0.2) +
    extra_RAD +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.4), "cm"))
  divplots= divplots+1
}

complete_one <- plot_grid(plot_list[[1]],plot_list[[3]],plot_list[[2]],
                          plot_list[[4]],plot_list[[5]],plot_list[[6]],
                          NULL, NULL,
                          divplot_list[[1]],divplot_list[[3]], 
                          divplot_list[[2]],divplot_list[[4]], ncol = 2, rel_heights = c(1,1,1,0.1,1,1), align = "hv") +
  theme(plot.background = element_rect(fill="white", color = NA))
y_space <- plot_grid(NULL, complete_one, ncol =2, rel_widths = c(0.03,1))+
  theme(plot.background = element_rect(fill="white", color = NA))
x_space <- plot_grid(y_space, NULL, ncol =1, rel_heights = c(1, 0.03))+
  theme(plot.background = element_rect(fill="white", color = NA))

rect_right <- rectGrob(
  x = .92,
  y = c(.96,0.77,0.58,.37,.18),
  width = unit(0.35, "in"),
  height = unit(0.2, "in"),
  gp = gpar(fill = "white", alpha=1)
)

rect_left <- rectGrob(
  x = .44,
  y = c(.96,0.77,0.58,.37,.18),
  width = unit(0.35, "in"),
  height = unit(0.2, "in"),
  gp = gpar(fill = "white", alpha=1)
)

labelled <- x_space +
  draw_label(expression(bold(paste("Genetic differentiation (F"[ST],")",sep = ""))), 
             color = "black", size = 8, angle = 90, x = 0.02, y = 0.72) + 
  draw_label("Expected Heterozygosity", fontface = "bold",
             color = "black", size = 8, angle = 90, x = 0.02, y = 0.23) +
  draw_label("Genomic position", color = "black", size = 8, angle = 0, x = .28, y = 0.02,fontface = "bold") +
  draw_label("Genomic position", color = "black", size = 8, angle = 0, x = .77, y = 0.02,fontface = "bold") +
  draw_label("A", color = "black", size = 10, angle = 0, x = 0.02, y = 0.98,fontface = "bold") +
  draw_label("B", color = "black", size = 10, angle = 0, x = 0.02, y = 0.40,fontface = "bold") +
  draw_grob(rect_right ) +
  draw_grob(rect_left ) +
  draw_text(c("T2-U1","T2-U2","U1-U2","U1","U2"), x = 0.92, y= c(.96,0.77,0.58,.37,.18),fontface = "bold",color = "black", size = 8) +
  draw_text(c("T1-U1","T1-U2","T1-T2","T1","T2"), x = 0.44, y= c(.96,0.77,0.58,.37,.18),fontface = "bold",color = "black", size = 8)

##Save as png or pdf
#ggsave2("~/Desktop/Wit_ea_Fig3.png", labelled, height = 5, width = 7.5)
ggsave2("Wit_ea_Fig3.pdf", labelled, height = 5, width = 7.5)
