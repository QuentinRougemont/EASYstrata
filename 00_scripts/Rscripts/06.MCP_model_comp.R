#!/usr/bin/env Rscript

################################################################################
######  Script to plot Ds values and perform changepoint analyses     ##########
#Author: QR
#Date: 26-01-24
################################################################################
argv <- commandArgs(T)
if (argv[1]=="-h" || length(argv)==0){
    cat("run script as:\n \tRscript ./06.MCP_model_comp.R YES/NO \n\n")
    cat(paste0("\033[0;41m","compulsory parameter:","\033[0m","\n")) 
    cat("\ta 'YES' or 'NO' argument stating wether or not an ancestral species was used in previous steps or not\n\n")
    cat("other input includes the dS values ordered along the ancestral proxy (generated from the previous steps)\n")
    cat("\tif not provided this input will be read automatically from file 02_results/dS.values.forchangepoint.txt\n")
} else {
    #take the single exepected argument:
    is_anc <- argv[1] #a simple YES/NO to state wether ancestral species is used or not
    ds_arg <- argv[2]

    if (exists("is_anc")) {
        print(paste0("was an ancestral genome used?", is_anc))
    } else {
        print("error no argument related to ancestral species provided")
        quit("no")
    }
    #---- load data ---- # 
    dsfile <- "02_results/dS.values.forchangepoint.txt"
    if (file.exists(dsfile)){
        print(paste0("reading ", dsfile)) 
        df <- read.table(dsfile, h = T) #a table with two columns : Ds and order 
        df <- na.omit(df)
    } else if (exists("ds_arg")) {
        print(paste0("reading ", ds_arg)) 
        df <- read.table(ds_arg, h = T) 
    } else {
        print("error no ds file provided")
        print("please provide a file of ds")
        quit("no")
    } 

    #--------------- check if library are installed -------------------------------#
    libs <- c('mcp','dplyr','ggplot2','magrittr','cowplot','ggstatsplot' )
    install.packages(setdiff(libs, rownames(installed.packages())), repos="https://cloud.r-project.org" )
    
    #---------------- load libraries ---------------------------------------------#
    invisible(lapply(libs, suppressWarnings(suppressMessages(suppressPackageStartupMessages(library))), character.only = TRUE))
        
    #create dir if not present:
    if (!dir.exists("02_results/modelcomp")){
      dir.create("02_results/modelcomp")
    }
    
    if (!dir.exists("02_results/modelcomp/noprior")){
      dir.create("02_results/modelcomp/noprior")
    }
    
    ## use my usual theme:
    ## ------------------ GGPLOT  CUSTOMISATION ------------------------------------------------##
    th_plot <-  theme(axis.title.x=element_text(size=14, family="Helvetica",face="bold"),
        axis.text.x=element_text(size=14,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
        axis.title.y=element_text(size=18, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
        axis.text.y=element_text(size=14,family="Helvetica",face="bold"),
        strip.text.x = element_text(size=18),panel.grid.major = element_blank())
    
    ## same but reduced police size:
    th_plot2 <-  theme(axis.title.x=element_text(size=8, family="Helvetica",face="bold"),
        axis.text.x=element_text(size=8,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
        axis.title.y=element_text(size=8, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
        axis.text.y=element_text(size=8,family="Helvetica",face="bold"),
        strip.text.x = element_text(size=7),panel.grid.major = element_blank())
    
    th_plot3<-  theme(axis.title.x=element_text(size=10, family="Helvetica",face="bold"),
        axis.text.x=element_text(size=10,family="Helvetica",face="bold", angle=0, hjust=0, vjust=0.5),
        axis.title.y=element_text(size=12, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
        axis.text.y=element_text(size=10,family="Helvetica",face="bold"),
        strip.text.x = element_text(size=12),panel.grid.major = element_blank())
    
    ## color:
    mycolor <-c("#E69F00",  "#0072B2" ,"#5B1A79",  "#CC79A7", "#D55E00", "red", "black","yellow")
    
    ## ------------------ DECLARE USEFULL FUNCTION----------------------------------------------##
    
    # simple color plot along the genomic coordinates of ancestral species: 
    dplot <- function(df_of_ds, nstrata, columnstrata) {
        columnstrata=sym(columnstrata)
        ggplot(df_of_ds, aes(x = start, y = Ds, colour = !!columnstrata)) + 
        geom_point( size = .5) + 
        facet_wrap(~scaff, scale="free_x") +
        theme_classic() +
        ylim(c(0,0.3)) +
        xlab("position along chr (bp)") +
        ylab(expression(italic(d[s]))) +
        th_plot2 + 
        theme(legend.position = "none") + 
        scale_colour_manual(values=mycolor[1:nstrata])  
    
    }
    
    # simple color plot along the ancestral order:
    dplot2 <- function(df_of_ds, nstrata, columnstrata) {
        columnstrata=sym(columnstrata)
        ggplot(df_of_ds, aes(x = orderchp, y = Ds, colour = !!columnstrata)) + 
        geom_point( size = .5) + 
        theme_classic() +
        ylim(c(0,0.4)) +
        xlab("position along gene order") +
        ylab(expression(italic(d[s]))) +
        th_plot2 + 
        theme(legend.position = "none") + 
        scale_colour_manual(values=mycolor[1:nstrata])  
    
    }
    
    # a simple function to plot violin plot multiple times:
    vplot <- function(data, nstrata, columnstrata) {
        columnstrata=sym(columnstrata)
        ggplot(data, aes(x = !!columnstrata,  y = Ds, fill = !!columnstrata)) + 
        geom_violin(trim = FALSE) + 
        geom_jitter(shape=16, position=position_jitter(0.2)) +
        theme_classic() + 
        th_plot  + 
        ylab(expression(italic(d[s]))) +
        scale_fill_manual(values=mycolor[1:nstrata])  + 
        theme(legend.position="none")
    }
    
    
    xl <- expression(paste("order along ", italic("ancestral species"), " mating chromsome"))
    plotcp <- function(cpmodel, title) {
     plot(cpmodel, q_fit = TRUE) + ggplot2::ggtitle(title) + 
         theme_classic() + 
         th_plot +
         geom_point(color = "darkblue", size = 0.1) +
         xlab(xl) +
         ylab(expression(italic(d[s])))
    }
    
    ################################################################################
    #                   perform the changepoint analyis here: 
    ################################################################################
    #define the model we want to test:
    #model10cp = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1, 1 ~ 1 , 1 ~ 1 , 1 ~ 1, 1  ~ 1, 1 ~ 1, 1 ~ 1,1 ~ 1) 
    #we probably don't need 10 !
    model9cp = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1, 1 ~ 1 , 1 ~ 1 , 1 ~ 1, 1  ~ 1 , 1 ~ 1, 1 ~ 1) 
    model8cp = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1, 1 ~ 1 , 1 ~ 1 , 1 ~ 1, 1  ~ 1 , 1 ~ 1) 
    model7cp = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1, 1 ~ 1 , 1 ~ 1 , 1 ~ 1, 1  ~ 1 ) 
    model6cp = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1, 1 ~ 1 , 1 ~ 1 , 1 ~ 1)
    model5cp = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1, 1 ~ 1 , 1 ~ 1)
    model4cp = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1, 1 ~ 1)
    model3cp = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1)
    model2cp = list(Ds ~ 1, 1~ 1, 1 ~ 1)
    model1cp = list(Ds ~ 1, 1~ 1)
    
    ## NOTE :
    #we don't use any prior for now
    path <- "02_results/modelcomp/noprior/"
    
    modelcp <- list(model1cp, model2cp, model3cp, model4cp, 
                    model5cp, model6cp, model7cp, model8cp)
    
    #note : the whole good below will be changed
    maxchgp = 8
    figcp <- vector('list', max(maxchgp))
    fitcp <- vector('list', max(maxchgp))
    m <- list()
    
    
    for(i in 1:maxchgp){
        message(i)
        fitcp[[i]] <- mcp(modelcp[[i]], 
                     data = df, 
                     par_x = "orderchp", 
                     iter = 5e5, 
                     adapt = 1.5e3,  
                     chains = 5, 
                     cores = 5 )
    
        figcp[[i]] <- plotcp(fitcp[[i]], paste0("Posterior fit ", i ,"  changepoint")) 
        
        pdf(file = paste0(path,"/Strata_comparison_", i , "chpt.pdf"), 10,6)
        plot_grid(print(figcp[[i]]), labels = "AUTO", ncol = 1)
        dev.off()
    
       m[[i]] <- summary(fitcp[[i]])
       write.table(m[[i]], paste0(path, "/model",i,"chpt.txt"), quote =F )
    
      if (i == 1){
        pdf(file = paste(path,"pars_1cp.pdf"))
        print(plot_pars(fitcp[[i]], pars = c("cp_1")))
        dev.off()
        df$two_strata <- ifelse(df$orderchp < m[[i]]$mean[1], "strata1", "strata2")
    
        fitcp[[i]]$loo <- loo(fitcp[[i]])
        
      } else if (i == 2){
        pdf(file =paste0(path, "pars_2cp.pdf"))
        print(plot_pars(fitcp[[i]], pars = c("cp_1" ,"cp_2")))
        dev.off()
       
         df$three_strata <- ifelse(df$orderchp < m[[2]]$mean[1], "strata1",
              ifelse(df$orderchp > m[[2]]$mean[2], "strata3", "strata2"))
        
        fitcp[[i]]$loo <- loo(fitcp[[i]])
        vp3 <- vplot(df,3,"three_strata")
        pdf(file = paste0(path, "violin_plot3strata.pdf"), 10,6)
        print(vp3)
        dev.off()
        
      } else if (i == 3){
          pdf(file = paste0(path, "pars_3cp.pdf"))
          print(plot_pars(fitcp[[i]], pars = c("cp_1" ,"cp_2","cp_3")))
          dev.off()
      
          df$four_strata <- ifelse(df$orderchp < m[[3]]$mean[1], "strata1", 
                            ifelse(df$orderchp > m[[3]]$mean[3], "strata4",
                            ifelse(df$orderchp > m[[3]]$mean[1] & df$orderchp < 
                                     m[[3]]$mean[2], 
                                   "strata2","strata3"))) 
    
          fitcp[[i]]$loo <- loo(fitcp[[i]])
          
          vp4 <- vplot(df,4, "four_strata")
          pdf(file = paste0(path, "violin_plot4Strata.pdf"),10,6)
          print(vp4)
          dev.off()
          
      } else if (i == 4) {
          pdf(file = paste0(path, "pars_4cp.pdf"))
          print(plot_pars(fitcp[[i]], 
                          pars = c("cp_1" ,"cp_2","cp_3","cp_4")))
          dev.off()
         df$five_strata <- ifelse(df$orderchp < m[[4]]$mean[1] , "strata1", 
                           ifelse(df$orderchp > m[[4]]$mean[4], "strata5",
                           ifelse(df$orderchp > m[[4]]$mean[1] & df$orderchp < 
                                    m[[4]]$mean[2], "strata2",
                           ifelse(df$orderchp > m[[4]]$mean[2] & df$orderchp < 
                                    m[[4]]$mean[3], "strata3",
                                  "strata4"))))  
          
        fitcp[[i]]$loo <- loo(fitcp[[i]])
        
        vp5 <- vplot(df,5,"five_strata")
        pdf(file = paste0(path, "violin_plot5strata.pdf"),10,6)
        print(vp5)
        dev.off()
        
      } else if (i == 5){
          pdf(file = paste0(path, "pars_5cp.pdf"))
          print(plot_pars(fitcp[[i]], 
                          pars = c("cp_1" ,"cp_2","cp_3","cp_4","cp_5")))
          dev.off()
            
          fitcp[[i]]$loo <- loo(fitcp[[i]])
        
          df$six_strata <- ifelse(df$orderchp < m[[5]]$mean[1], "strata1", 
                           ifelse(df$orderchp > m[[5]]$mean[5], "strata6",
                           ifelse(df$orderchp > m[[5]]$mean[1] & df$orderchp < m[[5]]$mean[2], "strata2",
                           ifelse(df$orderchp > m[[5]]$mean[2] & df$orderchp < m[[5]]$mean[3], "strata3",
                           ifelse(df$orderchp > m[[5]]$mean[3] & df$orderchp < m[[5]]$mean[4], "strata4", 
                                  "strata5" )))))
            
          vp6 <- vplot(df,6,"six_strata")
          pdf(file = paste0(path, "violin_plot6strata.pdf"), 10,6)
          print(vp6)
          dev.off()
          
      } else if (i == 6) {
            pdf(file = paste0(path, "pars_6cp.pdf"))
            print(plot_pars(fitcp[[i]], 
                            pars = c("cp_1" ,"cp_2","cp_3",
                                    "cp_4","cp_5","cp_6")))
            dev.off()
            
            fitcp[[i]]$loo <- loo(fitcp[[i]])
         
            df$seven_strata <- ifelse(df$orderchp < m[[6]]$mean[1], "strata1", 
                               ifelse(df$orderchp > m[[6]]$mean[6], "strata7",
                               ifelse(df$orderchp > m[[6]]$mean[1] & df$orderchp < m[[6]]$mean[2], "strata2",
                               ifelse(df$orderchp > m[[6]]$mean[2] & df$orderchp < m[[6]]$mean[3], "strata3",
                               ifelse(df$orderchp > m[[6]]$mean[3] & df$orderchp < m[[6]]$mean[4], "strata4",
                               ifelse(df$orderchp > m[[6]]$mean[4] & df$orderchp < m[[6]]$mean[5], "strata5",
                                      "strata6"))))))  
            
            vp7 <- vplot(df,7,"seven_strata")
            pdf(file = paste0(path, "violin_plot7strata.pdf"), 10,6)
            print(vp7)
            dev.off()
            
      } else if (i == 7){
            pdf(file = paste0(path, "pars_7cp.pdf"))
            print(plot_pars(fitcp[[i]], pars = c("cp_1" ,"cp_2","cp_3","cp_4",
                                           "cp_5","cp_6", "cp_7")))
            dev.off()
            
            fitcp[[i]]$loo <- loo(fitcp[[i]])
           
            df$eight_strata <- ifelse(df$orderchp < m[[7]]$mean[1], "strata1", 
                               ifelse(df$orderchp > m[[7]]$mean[7], "strata8",
                               ifelse(df$orderchp > m[[7]]$mean[1] & df$orderchp < m[[7]]$mean[2], "strata2",
                               ifelse(df$orderchp > m[[7]]$mean[2] & df$orderchp < m[[7]]$mean[3], "strata3",
                               ifelse(df$orderchp > m[[7]]$mean[3] & df$orderchp < m[[7]]$mean[4], "strata4",
                               ifelse(df$orderchp > m[[7]]$mean[4] & df$orderchp < m[[7]]$mean[5], "strata5",
                               ifelse(df$orderchp > m[[7]]$mean[5] & df$orderchp < m[[7]]$mean[6], "strata6",
                                        "strata7" )))))))
             
            
            vp8 <- vplot(df,8,"eight_strata")
            pdf(file = paste0(path, "violin_plot8strata.pdf"), 10,6)
            print(vp8)
            dev.off()
            
      } else { #assuming (i = 8)
            pdf(file = paste0(path, "pars_8cp.pdf"))
            print(plot_pars(fitcp[[i]], 
                            pars = c("cp_1" ,"cp_2","cp_3","cp_4",
                                     "cp_5","cp_6", "cp_7","cp_8")))
            dev.off()
            
            fitcp[[i]]$loo <- loo(fitcp[[i]])
     
            df$nine_strata <- ifelse(df$orderchp <m[[8]]$mean[1], "strata1", 
                              ifelse(df$orderchp > m[[8]]$mean[8], "strata9",
                              ifelse(df$orderchp > m[[8]]$mean[1] & df$orderchp < m[[8]]$mean[2], "strata2",
                              ifelse(df$orderchp > m[[8]]$mean[2] & df$orderchp < m[[8]]$mean[3], "strata3",
                              ifelse(df$orderchp > m[[8]]$mean[3] & df$orderchp < m[[8]]$mean[4], "strata4",
                              ifelse(df$orderchp > m[[8]]$mean[4] & df$orderchp < m[[8]]$mean[5], "strata5",
                              ifelse(df$orderchp > m[[8]]$mean[5] & df$orderchp < m[[8]]$mean[6], "strata6",
                              ifelse(df$orderchp > m[[8]]$mean[6] & df$orderchp < m[[8]]$mean[7], "strata7",
                                              "strata8" ))))))))
            vp9 <- vplot(df,9,"nine_strata")
            pdf(file = paste0(path, "violin_plot9strata.pdf"), 10,6)
            print(vp9)
            dev.off()
            
      }
    }
    
    
    #save.image( file = "02_results/modelcomp/changepoint_analysis.RData")
    
    # part 2: 
    writeLines("\nn~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    writeLines("\n\nperforming model choice \n\n")
    writeLines("\nn~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    
    #perform model choice and extract weights:
    loo_list = list()
    
    for (i in 1:maxchgp) {
    loo_list[i] <- list(fitcp[[i]]$loo)
    }
    
    #
    m.choice <- loo::loo_compare(loo_list)
    weights <- loo::loo_model_weights(loo_list, method="pseudobma")
    
    write.table(m.choice,"02_results/modelcomp/noprior/model_choice.txt",quote=F)
    write.table(weights,"02_results/modelcomp/noprior/model_weights.txt",quote=F, col.names=("weights"))
    write.table(df,"02_results/modelcomp/noprior/df.txt",quote=F,row.names=F,col.names=T,sep="\t")
    
    # ----- some hypothesis testing regarding differences among intervals -------- #
    #testing hypothesis
    #see more here: https://lindeloev.github.io/mcp/articles/comparison.html
    #below 3 changepoints, there is little relevance, so we test this directly
    
    writeLines("\nn~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    writeLines("\n\n compare adjacent strata value through BF \n\n")
    writeLines("\nn~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    
    
    for (i in 3:maxchgp){
      if (i == 3){
      hyp3 <- data.frame(matrix(ncol = 6, nrow = 4)) %>%
        set_colnames(.,c('hypothesis','mean','lower','upper','p','BF'))
      hyp3[1,] <-  hypothesis(fitcp[[i]],  c("int_1 < int_2"))
      hyp3[2,] <-  hypothesis(fitcp[[i]],  c("int_2 < int_3"))
      hyp3[3,] <-  hypothesis(fitcp[[i]],  c("int_1 > int_2"))
      hyp3[4,] <-  hypothesis(fitcp[[i]],  c("int_2 > int_3"))
        
      } else if (i == 4){
        hyp4 <- data.frame(matrix(ncol = 6, nrow = 6))%>%
          set_colnames(.,c('hypothesis','mean','lower','upper','p','BF'))
      
        hyp4[1,] <-  hypothesis(fitcp[[i]],  c("int_1 < int_2"))
        hyp4[2,] <-  hypothesis(fitcp[[i]],  c("int_2 < int_3"))
        hyp4[3,] <-  hypothesis(fitcp[[i]],  c("int_3 < int_4"))
        hyp4[4,] <-  hypothesis(fitcp[[i]],  c("int_1 > int_2"))
        hyp4[5,] <-  hypothesis(fitcp[[i]],  c("int_2 > int_3"))
        hyp4[6,] <-  hypothesis(fitcp[[i]],  c("int_3 > int_4"))
        
      } else if (i == 5){
        hyp5 <- data.frame(matrix(ncol = 6, nrow = 8))%>%
          set_colnames(.,c('hypothesis','mean','lower','upper','p','BF'))
        
        hyp5[1,] <-  hypothesis(fitcp[[i]],  c("int_1 < int_2"))
        hyp5[2,] <-  hypothesis(fitcp[[i]],  c("int_2 < int_3"))
        hyp5[3,] <-  hypothesis(fitcp[[i]],  c("int_3 < int_4"))
        hyp5[4,] <-  hypothesis(fitcp[[i]],  c("int_4 < int_5"))
        hyp5[5,] <-  hypothesis(fitcp[[i]],  c("int_1 > int_2"))
        hyp5[6,] <-  hypothesis(fitcp[[i]],  c("int_2 > int_3"))
        hyp5[7,] <-  hypothesis(fitcp[[i]],  c("int_3 > int_4"))
        hyp5[8,] <-  hypothesis(fitcp[[i]],  c("int_4 > int_5"))
    
     } else if (i == 6){
        hyp6 <- data.frame(matrix(ncol = 6, nrow = 10))%>%
          set_colnames(.,c('hypothesis','mean','lower','upper','p','BF'))
        
        hyp6[1,] <-  hypothesis(fitcp[[i]],  c("int_1 < int_2"))
        hyp6[2,] <-  hypothesis(fitcp[[i]],  c("int_2 < int_3"))
        hyp6[3,] <-  hypothesis(fitcp[[i]],  c("int_3 < int_4"))
        hyp6[4,] <-  hypothesis(fitcp[[i]],  c("int_4 < int_5"))
        hyp6[5,] <-  hypothesis(fitcp[[i]],  c("int_5 < int_6"))
        hyp6[6,] <-  hypothesis(fitcp[[i]],  c("int_1 > int_2"))
        hyp6[7,] <-  hypothesis(fitcp[[i]],  c("int_2 > int_3"))
        hyp6[8,] <-  hypothesis(fitcp[[i]],  c("int_3 > int_4"))
        hyp6[9,] <-  hypothesis(fitcp[[i]],  c("int_4 > int_5"))
        hyp6[10,] <-  hypothesis(fitcp[[i]],  c("int_5 > int_6"))
        
     } else if (i == 7){
      
       hyp7 <- data.frame(matrix(ncol = 6, nrow = 12))%>%
         set_colnames(.,c('hypothesis','mean','lower','upper','p','BF'))
       
       hyp7[1,] <-  hypothesis(fitcp[[i]],  c("int_1 < int_2"))
       hyp7[2,] <-  hypothesis(fitcp[[i]],  c("int_2 < int_3"))
       hyp7[3,] <-  hypothesis(fitcp[[i]],  c("int_3 < int_4"))
       hyp7[4,] <-  hypothesis(fitcp[[i]],  c("int_4 < int_5"))
       hyp7[5,] <-  hypothesis(fitcp[[i]],  c("int_5 < int_6"))
       hyp7[6,] <-  hypothesis(fitcp[[i]],  c("int_6 < int_7"))
       hyp7[7,] <-  hypothesis(fitcp[[i]],  c("int_1 > int_2"))
       hyp7[8,] <-  hypothesis(fitcp[[i]],  c("int_2 > int_3"))
       hyp7[9,] <-  hypothesis(fitcp[[i]],  c("int_3 > int_4"))
       hyp7[10,] <-  hypothesis(fitcp[[i]],  c("int_4 > int_5"))
       hyp7[11,] <-  hypothesis(fitcp[[i]],  c("int_5 > int_6"))
       hyp7[12,] <-  hypothesis(fitcp[[i]],  c("int_6 > int_7"))
    
     } else if (i == 8){
       
       hyp8 <- data.frame(matrix(ncol = 6, nrow = 14))%>%
         set_colnames(.,c('hypothesis','mean','lower','upper','p','BF'))
       
       hyp8[1,] <-  hypothesis(fitcp[[i]],  c("int_1 < int_2"))
       hyp8[2,] <-  hypothesis(fitcp[[i]],  c("int_2 < int_3"))
       hyp8[3,] <-  hypothesis(fitcp[[i]],  c("int_3 < int_4"))
       hyp8[4,] <-  hypothesis(fitcp[[i]],  c("int_4 < int_5"))
       hyp8[5,] <-  hypothesis(fitcp[[i]],  c("int_5 < int_6"))
       hyp8[6,] <-  hypothesis(fitcp[[i]],  c("int_6 < int_7"))
       hyp8[7,] <-  hypothesis(fitcp[[i]],  c("int_7 < int_8"))
       hyp8[8,] <-  hypothesis(fitcp[[i]],  c("int_1 > int_2"))
       hyp8[9,] <-  hypothesis(fitcp[[i]],  c("int_2 > int_3"))
       hyp8[10,] <-  hypothesis(fitcp[[i]],  c("int_3 > int_4"))
       hyp8[11,] <-  hypothesis(fitcp[[i]],  c("int_4 > int_5"))
       hyp8[12,] <-  hypothesis(fitcp[[i]],  c("int_5 > int_6"))
       hyp8[13,] <-  hypothesis(fitcp[[i]],  c("int_6 > int_7"))
       hyp8[14,] <-  hypothesis(fitcp[[i]],  c("int_7 > int_8"))
       
       
     }
    }
    
    #hypothesis testing attempt: 
    write.table(hyp3, paste0(path, "hypothesis3strata"), quote= F)
    write.table(hyp4, paste0(path,  "hypothesis4strata.txt"), quote= F)
    write.table(hyp5, paste0(path, "hypothesis5strata.txt"), quote= F)
    write.table(hyp6, paste0(path, "hypothesis6strata.txt"), quote= F)
    write.table(hyp7, paste0(path, "hypothesis7strata.txt"), quote= F)
    write.table(hyp8, paste0(path, "hypothesis8strata.txt"), quote= F)
    
    writeLines("\nn~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    writeLines("\n\n exporting some more plots \n\n")
    writeLines("\nn~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    
    
    #other plots: 
    ds3 <- dplot(df, nstrata=3, "three_strata")
    ds4 <- dplot(df, nstrata=4, "four_strata")
    ds5 <- dplot(df, nstrata=5, "five_strata")
    ds6 <- dplot(df, nstrata=6, "six_strata")
    ds7 <- dplot(df, nstrata=7, "seven_strata")
    ds8 <- dplot(df, nstrata=8, "eight_strata")
    
    pdf(file=paste0(path,"plot_dS_all_position.pdf"),8,12)
    plot_grid(ds3, ds4, ds5, ds6, ds7, ds8, 
              labels = c("A - three changepoint",
                         "B - four changepoint" ,
                         "C - five changepoint" ,
                         "D - six changepoint"  ,
                         "E - seven changepoint" ,
                         "F - eight changepoint") , 
             label_size = 7,
             hjust = -0.5, vjust = -0.5,
             ncol = 1)
    dev.off()
    
    ds3 <- dplot2(df, nstrata=3, "three_strata")
    ds4 <- dplot2(df, nstrata=4, "four_strata")
    ds5 <- dplot2(df, nstrata=5, "five_strata")
    ds6 <- dplot2(df, nstrata=6, "six_strata")
    ds7 <- dplot2(df, nstrata=7, "seven_strata")
    ds8 <- dplot2(df, nstrata=8, "eight_strata")
    
    pdf(file=paste0(path,"plot_Ds_along_order.pdf"),8,12)
    plot_grid(ds3, ds4, ds5, ds6, ds7, ds8, 
              labels = c("A - three changepoint",
                         "B - four changepoint" ,
                         "C - five changepoint" ,
                         "D - six changepoint"  ,
                         "E - seven changepoint" ,
                         "F - eight changepoint") , 
             label_size = 7,
             hjust = -0.5, vjust = -0.5,
             ncol = 1)
    dev.off()
    
    
    #finally: 
    if(is_anc=="YES"){ 
    s3.anc.h1 <- select(df, gene, geneX, three_strata)
    s4.anc.h1 <- select(df, gene, geneX, four_strata)
    s5.anc.h1 <- select(df, gene, geneX, five_strata)
    s6.anc.h1 <- select(df, gene, geneX, five_strata)
    s6.anc.h1 <- select(df, gene, geneX, six_strata)
    s7.anc.h1 <- select(df, gene, geneX, seven_strata)
    s8.anc.h1 <- select(df, gene, geneX, eight_strata)
    
    s3.h1.h2 <- select(df, geneX, geneY.x, three_strata)
    s4.h1.h2 <- select(df, geneX, geneY.x, four_strata)
    s5.h1.h2 <- select(df, geneX, geneY.x, five_strata)
    s6.h1.h2 <- select(df, geneX, geneY.x, five_strata)
    s6.h1.h2 <- select(df, geneX, geneY.x, six_strata)
    s7.h1.h2 <- select(df, geneX, geneY.x, seven_strata)
    s8.h1.h2 <- select(df, geneX, geneY.x, eight_strata)
    
    write.table(s3.anc.h1,paste0(path,"classif.s3.ancestral.haplo1"),
        quote=F,row.names=F,col.names=F,sep="\t")
    write.table(s4.anc.h1,paste0(path,"classif.s4.ancestral.haplo1"),
        quote=F,row.names=F,col.names=F,sep="\t")
    write.table(s5.anc.h1,paste0(path,"classif.s5.ancestral.haplo1"),
        quote=F,row.names=F,col.names=F,sep="\t")
    write.table(s6.anc.h1,paste0(path,"classif.s6.ancestral.haplo1"),
        quote=F,row.names=F,col.names=F,sep="\t")
    write.table(s7.anc.h1,paste0(path,"classif.s7.ancestral.haplo1"),
        quote=F,row.names=F,col.names=F,sep="\t")
    write.table(s8.anc.h1,paste0(path,"classif.s8.ancestral.haplo1"),
        quote=F,row.names=F,col.names=F,sep="\t")
    
    }else{
            #geneY.x ortho   geneY.y
    s3.h1.h2 <- select(df, gene, geneY.y, three_strata)
    s4.h1.h2 <- select(df, gene, geneY.y, four_strata)
    s5.h1.h2 <- select(df, gene, geneY.y, five_strata)
    s6.h1.h2 <- select(df, gene, geneY.y, five_strata)
    s6.h1.h2 <- select(df, gene, geneY.y, six_strata)
    s7.h1.h2 <- select(df, gene, geneY.y, seven_strata)
    s8.h1.h2 <- select(df, gene, geneY.y, eight_strata)
    
    }
    
    write.table(s3.h1.h2,paste0(path,"classif.s3.haplo1.haplo2"),
        quote=F,row.names=F,col.names=F,sep="\t")
    write.table(s4.h1.h2,paste0(path,"classif.s4.haplo1.haplo2"),
        quote=F,row.names=F,col.names=F,sep="\t")
    write.table(s5.h1.h2,paste0(path,"classif.s5.haplo1.haplo2"),
        quote=F,row.names=F,col.names=F,sep="\t")
    write.table(s6.h1.h2,paste0(path,"classif.s6.haplo1.haplo2"),
        quote=F,row.names=F,col.names=F,sep="\t")
    write.table(s7.h1.h2,paste0(path,"classif.s7.haplo1.haplo2"),
        quote=F,row.names=F,col.names=F,sep="\t")
    write.table(s8.h1.h2,paste0(path,"classif.s8.haplo1.haplo2"),
        quote=F,row.names=F,col.names=F,sep="\t")
    
    ################################################################################
    #finally construct some combined plot:
    plot <- list()
    for (i in 1:maxchgp) {
    plot[[i]] <- plot_grid(print(figcp[[i]]) + th_plot3, labels = "AUTO", ncol = 1)
    }
    
    writeLines("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    writeLines("\n\n comparing models with gstatsplot\n\n")
    writeLines("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    
    vp3s <- ggbetweenstats(df, three_strata, Ds)  + th_plot3
    vp4s <- ggbetweenstats(df, four_strata, Ds)   + th_plot3
    vp5s <- ggbetweenstats(df, five_strata, Ds)   + th_plot3
    vp6s <- ggbetweenstats(df, six_strata, Ds)    + th_plot3
    vp7s <- ggbetweenstats(df, seven_strata, Ds)  + th_plot3
    #this is probably too much: 
    vp8s <- ggbetweenstats(df, eight_strata, Ds, palette = "Paired")  + th_plot3
    vp9s <- ggbetweenstats(df, nine_strata, Ds, palette = "Paired")  + th_plot3
    
    #comment if you don't want some:
    #TO DO: modify to plot each separately and to remove boxplot when n is low 
    
    pdf(file = paste0(path, "viobox_ds_strata_distribution_priors.pdf"), 10,28)
    plot_grid(
        #vp3s, 
        vp4s,
        vp5s,
        vp6s,
        vp7s,
        vp8s, 
        vp9s,
        labels = "AUTO", ncol = 1)
    dev.off()
    
    pdf(file = paste0(path, "all_comp.pdf"), 10,20)
    plot_grid(
      #plot[[3]],
      plot[[4]],
      plot[[5]],
      plot[[6]],
      plot[[7]],
      plot[[8]],
      labels = "AUTO", ncol = 1)
    dev.off()
    
    
    #comment those that are not wanted:
    pdf(file = paste0(path, "strata_and_viobox_ds_strata_distribution_priors.pdf"), 20,20)
    plot_grid(
      #plot[[2]],vp3s,
      #plot[[3]],vp4s,
      plot[[4]],vp5s,
      plot[[5]],vp6s,
      plot[[6]],vp6s,
      plot[[7]],vp8s,
      plot[[8]],vp9s,
      labels = "AUTO", ncol = 2, rel_heights = c(.4,1,.4,1,.4,1))
    dev.off()
    
    writeLines("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    writeLines("\n analyses finished !! \n\n")
    writeLines("\nn~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    
}
