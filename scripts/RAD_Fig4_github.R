library("tidyr")  
library("ggplot2")
library("patchwork")
library("dplyr")
library("RcppRoll")
library("cowplot")

# Data loading exp het
###############################
# IMPORTANT!!!!!
# Download RAD_Haemonchus folder and as set working directory
setwd("~/Desktop/RAD_Haemonchus/")
###############################


DIV =  read.table(file= "data/DataxpHetData.tsv", header = TRUE) #Exp Het
DIV$Type = "DIV"

##Figure setup
extra_RAD_x <- ggplot2::theme(panel.background = element_rect(fill = NA, colour = "black", size = 0.5),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            legend.position="none", 
                            axis.text.y = element_text(size = 6),
                            axis.title.x = element_text( size=10, face = "bold"), 
                            axis.title.y = element_text( size=10, face = "bold", vjust = 1.5), 
                            title = element_text(size=8, face ="bold"),
                            strip.text = element_text(colour = 'black'),
                            axis.ticks = element_line(colour = "black", size = 0.2), 
                            axis.ticks.length = unit(0.05, "cm"))

extra_RAD_nox <- ggplot2::theme(panel.background = element_rect(fill = NA, colour = "black", size = 0.5),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            legend.position="none", 
                            axis.text.x = element_blank(), 
                            axis.text.y = element_text(size = 6),
                            axis.title.x = element_blank(), 
                            axis.title.y = element_text( size=10, face = "bold", vjust = 1.5), 
                            title = element_text(size=8, face ="bold"),
                            strip.text = element_text(colour = 'black'),
                            axis.ticks = element_line(colour = "black", size = 0.2), 
                            axis.ticks.length = unit(0.05, "cm"))

extra_RAD_xnoy <- ggplot2::theme(panel.background = element_rect(fill = NA, colour = "black", size = 0.5),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              legend.position="none", 
                              axis.text.y = element_text(size = 6),
                              axis.title.x = element_text( size=10, face = "bold"), 
                              axis.title.y = element_blank(), 
                              title = element_text(size=8, face ="bold"),
                              strip.text = element_text(colour = 'black'),
                              axis.ticks = element_line(colour = "black", size = 0.2), 
                              axis.ticks.length = unit(0.05, "cm"))

extra_RAD_noxnoy <- ggplot2::theme(panel.background = element_rect(fill = NA, colour = "black", size = 0.5),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                legend.position="none", 
                                axis.text.x = element_blank(), 
                                axis.text.y = element_text(size = 6),
                                axis.title.x = element_blank(), 
                                axis.title.y = element_blank(),  
                                title = element_text(size=8, face ="bold"),
                                strip.text = element_text(colour = 'black'),
                                axis.ticks = element_line(colour = "black", size = 0.2), 
                                axis.ticks.length = unit(0.05, "cm"))


chr_col <- c("1"= "red","2" ="maroon","3" ="mediumpurple4","4" ="steelblue","5" = "turquoise4","X" = "mediumseagreen")

##Loop through all populations and both expected heterozygosity and LD. Can loop through multiple chromosomes at a time

for(i in c(1:4)){
  if(i == 1 | i == 2){
    pop1 <- paste("R",i,"_", sep="")
    pop2 <- paste("Res",i,sep= "")
  }
  if(i == 3 | i == 4){
    pop1 <- paste("S",i-2,"_", sep="")
    pop2 <- paste("Susc",i-2,sep= "")
  }
  for(j in 1){ #Chromosome of interest, can be more than one.
    chr <- paste("chr",j,sep = "")
    Chr <- paste("Chr", j, sep ="")
    ld_file <- paste("data/LDLine_",pop1,Chr,".tsv", sep="")
    col_sel <- chr_col[j] 

    ##LD
    LD1 = read.table(file= ld_file, header = TRUE)
    LD1$Type = "LD"
    LD1$BP = LD1$Bin*150000-75000
    LD1$Value = LD1$RealHits/LD1$MaxHits
    LDSelect = LD1[,c(13,11,12,14)]
    
    ##DIV
    DivCh1 = DIV[DIV$Chr == Chr & DIV$PopID == pop2,]
    DIVSelect = DivCh1[,c(6,2,9,8)]#[,c(5,25,26,16)]
    ROLL <- DIVSelect
    
    ROLL$rolling <- roll_mean(ROLL$Exp_Het,n=500, align = "center",fill=NA, by=20)
    ROLL <- ROLL %>%
      dplyr::filter(!is.na(rolling))

    ##File Headigns and merging
    colnames(DIVSelect) = c("BP","Sign","Type","Value")
    colnames(LDSelect) = c("BP","Sign","Type","Value")
    ALL = rbind(DIVSelect,LDSelect) 

    ALL$BP = as.integer(ALL$BP)
    ALL$Value = as.numeric(ALL$Value)
    ALL$Sign  <- as.factor(ALL$Sign)

    ##Re-labelling
    ALL$Type <- factor(ALL$Type, levels = c("DIV","LD"))
    cols <- c("1"=0.05, "2" = 1)
    sizes = c("1"=0.8, "2" = .8)

    #Looping through individual populations and measures
    if(i == 1){
      labelsType = c("FST1" = "T1-U1","FST2" = "T1-U2", "FST3" = "T1-T2","DIV" = "pi" ,"LD" = "LD (R2)")
      RES1LD <- ggplot(ALL[ALL$Type == "LD", ], aes(x= BP,y= Value, axes = FALSE, alpha = Sign, size = Sign)) + 
        theme_light() +
        geom_point(data = subset(ALL, Type == "LD" & Sign == 2), colour = col_sel) +
        geom_line(data = subset(ALL, Type == "LD" & Sign == 1), color = col_sel) +
        geom_vline(xintercept=c(7027095,26596797), color = "blue", linetype="dashed", size=0.2) +
        scale_alpha_manual(values = cols) +
        scale_size_manual(values = sizes) +
        labs(y = expression(bold(paste("LD (R"^2,")")))) +
        scale_x_continuous(name = "Nucleotide position", breaks=c(0,20000000,40000000), labels = c("0","20M","40M")) +
        scale_y_continuous(breaks=c(0,0.2,0.4,0.6), labels = c("0.00","0.20","0.40","0.60"), limits = c(0,.67)) +
        extra_RAD_nox+
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
     RES1DIV <- ggplot(ALL[ALL$Type == "DIV", ], aes(x= BP,y= Value, axes = FALSE, alpha = Sign, size = Sign)) + 
        theme_light() +
        geom_point(colour = col_sel) +
        geom_line(data = ROLL,aes(x = BP, y = rolling), colour = "red", size =0.6, alpha =1) +
        scale_alpha_manual(values = cols) +
        scale_size_manual(values = sizes) +
        geom_vline(xintercept=c(7027095,26596797) , color = "blue", linetype="dashed", size=0.2) +
        scale_x_continuous(name = "Nucleotide position", breaks=c(0,20000000,40000000), labels = c("0","20M","40M")) +
        scale_y_continuous(name = "Exp het", breaks=c(0,0.25,0.5), labels = c("0.00","0.25","0.50"), limits = c(0,0.52632), expand =c(0.02,0.02)) + #,expression(bold(paste("Diversity (", pi, ")", sep ="")))
        extra_RAD_nox+
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
      }
    if(i == 2){
      labelsType = c("FST1" = "T2-U1","FST2" = "T2-U2", "FST3" = "T1-T2","DIV" = "pi" ,"LD" = "LD (R2)", "TAJD" = "Tajima's D")
      RES2LD <- ggplot(ALL[ALL$Type == "LD", ], aes(x= BP,y= Value, axes = FALSE, alpha = Sign, size = Sign)) + 
        theme_light() +
        geom_point(data = subset(ALL, Type == "LD" & Sign == 2), colour = col_sel) +
        geom_line(data = subset(ALL, Type == "LD" & Sign == 1), color = col_sel) +
        geom_vline(xintercept=c(7027095,26596797) , color = "blue", linetype="dashed", size=0.2) +
        scale_alpha_manual(values = cols) +
        scale_size_manual(values = sizes) +
        scale_x_continuous(name = "Nucleotide position", breaks=c(0,20000000,40000000), labels = c("0","20M","40M")) +
        scale_y_continuous(breaks=c(0,0.2,0.4,0.6), labels = c("0.00","0.20","0.40","0.60"), limits = c(0,.67)) +
        extra_RAD_noxnoy+ #extra_RAD_xnoy + #
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
     RES2DIV <- ggplot(ALL[ALL$Type == "DIV", ], aes(x= BP,y= Value, axes = FALSE, alpha = Sign, size = Sign)) + 
        theme_light() +
        geom_point(colour = col_sel) +
        geom_line(data = ROLL,aes(x = BP, y = rolling), colour = "red", size =0.6, alpha =1) +
        #geom_vline(xintercept=7027095 , color = "blue", linetype="dashed", size=0.2) +
        scale_alpha_manual(values = cols) +
        scale_size_manual(values = sizes) +
        geom_vline(xintercept=c(7027095,26596797) , color = "blue", linetype="dashed", size=0.2) +
        scale_x_continuous(name = "Nucleotide position", breaks=c(0,20000000,40000000), labels = c("0","20M","40M")) +
        scale_y_continuous(name = NULL, breaks=c(0,0.25,0.5), labels = c("0.00","0.25","0.50"), limits = c(0,0.52632), expand =c(0.02,0.02)) + 
        extra_RAD_noxnoy +
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
    }
    if(i == 3){
      labelsType = c("FST3" = "U1-U2","DIV" = "pi" ,"LD" = "LD (R2)", "TAJD" = "Tajima's D")
      ALL <- ALL %>%
        dplyr::filter(!Type %in% c("FST1","FST2"))
      Susc1LD <- ggplot(ALL[ALL$Type == "LD", ], aes(x= BP,y= Value, axes = FALSE, alpha = Sign, size = Sign)) + 
        theme_light() +
        geom_point(data = subset(ALL, Type == "LD" & Sign == 2), colour = col_sel) +
        geom_line(data = subset(ALL, Type == "LD" & Sign == 1), color = col_sel) +
        geom_vline(xintercept=c(7027095,26596797) , color = "blue", linetype="dashed", size=0.2) +
        scale_alpha_manual(values = cols) +
        scale_size_manual(values = sizes) +
        labs(y = expression(bold(paste("LD (R"^2,")")))) +
        scale_x_continuous(name = "Nucleotide position", breaks=c(0,20000000,40000000), labels = c("0","20M","40M")) +
        scale_y_continuous(breaks=c(0,0.2,0.4,0.6), labels = c("0.00","0.20","0.40","0.60"), limits = c(0,.67)) +
        extra_RAD_x+ #extra_RAD_nox+
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
      Susc1DIV <- ggplot(ALL[ALL$Type == "DIV", ], aes(x= BP,y= Value, axes = FALSE, alpha = Sign, size = Sign)) + 
        theme_light() +
        geom_point(colour = col_sel) +
        geom_line(data = ROLL,aes(x = BP, y = rolling), colour = "red", size =0.6, alpha =1) +
        scale_alpha_manual(values = cols) +
        scale_size_manual(values = sizes) +
        geom_vline(xintercept=c(7027095,26596797) , color = "blue", linetype="dashed", size=0.2) +
        scale_x_continuous(name = "Nucleotide position", breaks=c(0,20000000,40000000), labels = c("0","20M","40M")) +
        scale_y_continuous(name = "Exp het", breaks=c(0,0.25,0.5), labels = c("0.00","0.25","0.50"), limits = c(0,0.52632), expand =c(0.02,0.02)) + #,expression(bold(paste("Diversity (", pi, ")", sep ="")))
        labs(y= "Tajima's D")+ 
        extra_RAD_nox+
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
    }
    if(i == 4){
      labelsType = c("FST3" = "U1-U2","DIV" = "pi" ,"LD" = "LD (R2)", "TAJD" = "Tajima's D")
      ALL <- ALL %>%
        dplyr::filter(!Type %in% c("FST1","FST2"))
      Susc2LD <- ggplot(ALL[ALL$Type == "LD", ], aes(x= BP,y= Value, axes = FALSE, alpha = Sign, size = Sign)) + 
        theme_light() +
        geom_point(data = subset(ALL, Type == "LD" & Sign == 2), colour = col_sel) +
        geom_line(data = subset(ALL, Type == "LD" & Sign == 1), color = col_sel) +
        geom_vline(xintercept=c(7027095,26596797) , color = "blue", linetype="dashed", size=0.2) +
        scale_alpha_manual(values = cols) +
        scale_size_manual(values = sizes) +
        labs(y = "Untreated 2") +
        scale_x_continuous(name = "Nucleotide position", breaks=c(0,20000000,40000000), labels = c("0","20M","40M")) +
        scale_y_continuous(breaks=c(0,0.2,0.4,0.6), labels = c("0.00","0.20","0.40","0.60"), limits = c(0,.67)) +
        extra_RAD_xnoy +
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
      Susc2DIV <- ggplot(ALL[ALL$Type == "DIV", ], aes(x= BP,y= Value, axes = FALSE, alpha = Sign, size = Sign)) + 
        theme_light() +
        geom_point(colour = col_sel) +
        geom_line(data = ROLL,aes(x = BP, y = rolling), colour = "red", size =0.6, alpha =1) +
        scale_alpha_manual(values = cols) +
        scale_size_manual(values = sizes) +
        geom_vline(xintercept=c(7027095,26596797) , color = "blue", linetype="dashed", size=0.2) +
        scale_x_continuous(name = "Nucleotide position", breaks=c(0,20000000,40000000), labels = c("0","20M","40M")) +
        scale_y_continuous(name = "Exp het", breaks=c(0,0.25,0.5), labels = c("0.00","0.25","0.50"), limits = c(0,0.52632), expand =c(0.02,0.02)) + #expression(paste("Diversity (", pi, ")", sep =""))
        extra_RAD_noxnoy +
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
    }

  }
}


##Generating figure for chrmosome 1

Res1 <- plot_grid(RES1DIV,RES1LD, ncol =1, align = "v", rel_heights = c(1,1) )  
Res2 <- plot_grid(RES2DIV,RES2LD, ncol =1, align = "v", rel_heights = c(1,1) )  
Susc1 <- plot_grid(Susc1DIV,Susc1LD, ncol =1, align = "v", rel_heights = c(1,1.25) )
Susc2 <- plot_grid(Susc2DIV,Susc2LD, ncol =1, align = "v", rel_heights = c(1,1.25) ) 


All <- plot_grid(Res1,NULL,Res2,NULL,NULL,NULL,Susc1,NULL,Susc2, ncol= 3, rel_widths = c(1,0.05,1,1,0.05,1,1,0.05,1), 
                 labels = c("A","B","","","","","C","D",""), label_size =10, rel_heights = c(0.86,0.05,1))+
  theme(plot.background = element_rect(fill="white", color = NA))

#Save as png or pdf
#ggsave2("Wit_ea_Fig4.png", All, height = 5, width = 7.5) 
ggsave2("Wit_ea_Fig4.pdf", All, height = 5, width = 7.5) 
