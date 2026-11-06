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
        xlab("Position along chromosome (bp)") +
        ylab(expression(italic(d[S]))) +
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
        xlab("Position along inferred gene rank") +
        ylab(expression(italic(d[S]))) +
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
        ylab(expression(italic(d[S]))) +
        scale_fill_manual(values=mycolor[1:nstrata])  + 
        theme(legend.position="none")
    }
    
    
    xl <- expression(paste("Order along  mating chromsome"))
    plotcp <- function(cpmodel, title) {
     plot(cpmodel, interval = "fitted") + ggplot2::ggtitle(title) + 
         theme_classic() + 
         th_plot +
         geom_point(color = "darkblue", size = 0.1) +
         xlab(xl) +
         ylab(expression(italic(d[S])))
    }
   
    ################################################################################
    #                   perform the changepoint analyis here: 
    ################################################################################
    #define the model we want to test:
    #model10cp = list(Ds ~ 1, ~ 1, ~ 1, ~ 1, ~ 1 , ~ 1 , ~ 1, ~ 1, ~ 1, ~ 1, ~ 1) 
    #we probably don't need 10 !
    model9cp = list(Ds ~ 1, ~ 1, ~ 1, ~ 1, ~ 1 , ~ 1 , ~ 1, ~ 1 , ~ 1, ~ 1) #10 strata (8 strata if PAR on each side) 
    model8cp = list(Ds ~ 1, ~ 1, ~ 1, ~ 1, ~ 1 , ~ 1 , ~ 1, ~ 1 , ~ 1)        #9 strata  (7 strata if PAR on each side)  
    model7cp = list(Ds ~ 1, ~ 1, ~ 1, ~ 1, ~ 1 , ~ 1 , ~ 1, ~ 1 )               #8 strata  
    model6cp = list(Ds ~ 1, ~ 1, ~ 1, ~ 1, ~ 1 , ~ 1 , ~ 1)                        #7 strata
    model5cp = list(Ds ~ 1, ~ 1, ~ 1, ~ 1, ~ 1 , ~ 1)                                #....
    model4cp = list(Ds ~ 1, ~ 1, ~ 1, ~ 1, ~ 1)
    model3cp = list(Ds ~ 1, ~ 1, ~ 1, ~ 1)
    model2cp = list(Ds ~ 1, ~ 1, ~ 1)
    model1cp = list(Ds ~ 1, ~ 1)
    
    ## NOTE :
    #we don't use any prior for now
    path <- "02_results/modelcomp/noprior/"
    
    modelcp <- list(model1cp, model2cp, model3cp, model4cp, 
                    model5cp, model6cp, model7cp, model8cp, model9cp)
    
    #note : the whole good below will be changed
    maxchgp <- 9
    figcp <- vector('list', max(maxchgp))
    fitcp <- vector('list', max(maxchgp))
    m <- list()
    

    future::plan(future::multisession, workers = 6)  # run chains in parallel
    for(i in 1:maxchgp){
        message(i)
        fitcp[[i]] <- mcp(modelcp[[i]], 
                     data = df, 
                     par_x = "orderchp", 
                     iter = 8e4, 
                     warmup = 1e3,  
                     chains = 6 )
    
        figcp[[i]] <- plotcp(fitcp[[i]], paste0("Posterior fit ", i ,"  changepoint")) 
        
        pdf(file = paste0(path,"/Strata_comparison_", i , "chpt.pdf"), 10,5)
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
            
      } else if (i == 8){ 
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
      } else { #assuming (i == 9) - 10 strata
            pdf(file = paste0(path, "pars_9cp.pdf"))
            print(plot_pars(fitcp[[i]], 
                            pars = c("cp_1" ,"cp_2","cp_3","cp_4",
                                     "cp_5","cp_6", "cp_7","cp_8","cp_9")))
            dev.off()
            
            fitcp[[i]]$loo <- loo(fitcp[[i]])
     
            df$ten_strata <- ifelse(df$orderchp <m[[9]]$mean[1], "strata1", 
                              ifelse(df$orderchp > m[[9]]$mean[9], "strata10",
                              ifelse(df$orderchp > m[[9]]$mean[1] & df$orderchp < m[[9]]$mean[2], "strata2",
                              ifelse(df$orderchp > m[[9]]$mean[2] & df$orderchp < m[[9]]$mean[3], "strata3",
                              ifelse(df$orderchp > m[[9]]$mean[3] & df$orderchp < m[[9]]$mean[4], "strata4",
                              ifelse(df$orderchp > m[[9]]$mean[4] & df$orderchp < m[[9]]$mean[5], "strata5",
                              ifelse(df$orderchp > m[[9]]$mean[5] & df$orderchp < m[[9]]$mean[6], "strata6",
                              ifelse(df$orderchp > m[[9]]$mean[6] & df$orderchp < m[[9]]$mean[7], "strata7",
                              ifelse(df$orderchp > m[[9]]$mean[7] & df$orderchp < m[[9]]$mean[8], "strata8",
                                              "strata9" )))))))))
            vp10 <- vplot(df,10,"ten_strata")
            pdf(file = paste0(path, "violin_plot9strata.pdf"), 10,6)
            print(vp10)
            dev.off()

      }
    }
    
    # part 2: 
    writeLines("\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    writeLines("\n\nperforming model choice \n\n")
    writeLines("\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    
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
    
    writeLines("\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    writeLines("\n\n compare adjacent strata value through BF \n\n")
    writeLines("\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    
    for (i in 1:maxchgp){
      if (i == 1){
        hyp2 <- data.frame(matrix(ncol = 7, nrow = 2)) %>%
            set_colnames(.,c('hypothesis','mean','sd','lower','upper','prob','BF'))
        hyp2[1,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 < Intercept_2"))
        hyp2[2,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 < Intercept_1"))
        
        write.table(hyp2, paste0(path, "hypothesis2strata"), quote= F)
        #plot_pos_list[[1]] <- dplot(df, nstrata=2, "two_strata") + ggtitle("A - two strata")
        #plot_order_list[[1]] <- dplot2(df, nstrata=2, "two_strata") + ggtitle("A - two strata")
        #viobox_list[[1]] <- ggbetweenstats(df, two_strata, Ds)  + ylab(expression(italic(d[s]))) + th_plot3
    
        ds2.1 <- dplot(df, nstrata=2, "two_strata") + ggtitle("A - one changepoint")
        ds2.2 <- dplot2(df, nstrata=2, "two_strata") + ggtitle("A - one changepoint")
        vp2s <- ggbetweenstats(df, two_strata, Ds)  + ylab(expression(italic(d[s]))) + th_plot3

     } else if (i == 2){
        hyp3 <- data.frame(matrix(ncol = 7, nrow = 4)) %>%
            set_colnames(.,c('hypothesis','mean','sd','lower','upper','prob','BF'))
        hyp3[1,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 < Intercept_2"))
        hyp3[2,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 < Intercept_3"))
        hyp3[3,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 > Intercept_2"))
        hyp3[4,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 > Intercept_3"))
        
        write.table(hyp3, paste0(path, "hypothesis3strata"), quote= F)
        #plot_pos_list[[2]] <- dplot(df, nstrata=3, "three_strata") + ggtitle("B - three changepoint")
        #plot_order_list[[2]] <- dplot2(df, nstrata=3, "three_strata") +
        #viobox_list[[2]] <- ggbetweenstats(df, two_strata, Ds)  + ylab(expression(italic(d[s]))) + th_plot3
        ds3.1 <- dplot(df, nstrata=3, "three_strata")  + ggtitle("B - two changepoint")
        ds3.2 <- dplot2(df, nstrata=3, "three_strata")  + ggtitle("B - two changepoint")
        vp3s <- ggbetweenstats(df, three_strata, Ds)  + ylab(expression(italic(d[s]))) + th_plot3

      } else if (i == 3){
        hyp4 <- data.frame(matrix(ncol = 7, nrow = 6))%>%
          set_colnames(.,c('hypothesis','mean','sd','lower','upper','prob','BF'))
      
        hyp4[1,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 < Intercept_2"))
        hyp4[2,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 < Intercept_3"))
        hyp4[3,] <-  hypothesis(fitcp[[i]],  c("Intercept_3 < Intercept_4"))
        hyp4[4,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 > Intercept_2"))
        hyp4[5,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 > Intercept_3"))
        hyp4[6,] <-  hypothesis(fitcp[[i]],  c("Intercept_3 > Intercept_4"))
        
        write.table(hyp4, paste0(path,  "hypothesis4strata.txt"), quote= F)
        ds4.1 <- dplot(df, nstrata=4, "four_strata")  + ggtitle("C - three changepoint")
        ds4.2 <- dplot2(df, nstrata=4, "four_strata") + ggtitle("C - three changepoint")
        vp4s <- ggbetweenstats(df, four_strata, Ds)   + ylab(expression(italic(d[s]))) + th_plot3
    
      } else if (i == 4){
        hyp5 <- data.frame(matrix(ncol = 7, nrow = 8))%>%
          set_colnames(.,c('hypothesis','mean','sd','lower','upper','prob','BF'))
        
        hyp5[1,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 < Intercept_2"))
        hyp5[2,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 < Intercept_3"))
        hyp5[3,] <-  hypothesis(fitcp[[i]],  c("Intercept_3 < Intercept_4"))
        hyp5[4,] <-  hypothesis(fitcp[[i]],  c("Intercept_4 < Intercept_5"))
        hyp5[5,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 > Intercept_2"))
        hyp5[6,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 > Intercept_3"))
        hyp5[7,] <-  hypothesis(fitcp[[i]],  c("Intercept_3 > Intercept_4"))
        hyp5[8,] <-  hypothesis(fitcp[[i]],  c("Intercept_4 > Intercept_5"))
    
        write.table(hyp5, paste0(path, "hypothesis5strata.txt"), quote= F)
        ds5.1 <- dplot(df, nstrata=5, "five_strata")  + ggtitle("D - four changepoint")
        ds5.2 <- dplot2(df, nstrata=5, "five_strata") + ggtitle("D - four changepoint")
        vp5s <- ggbetweenstats(df, five_strata, Ds)   + ylab(expression(italic(d[s]))) + th_plot3

     } else if (i == 5){
        hyp6 <- data.frame(matrix(ncol = 7, nrow = 10))%>%
          set_colnames(.,c('hypothesis','mean','sd','lower','upper','prob','BF'))
        
        hyp6[1,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 < Intercept_2"))
        hyp6[2,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 < Intercept_3"))
        hyp6[3,] <-  hypothesis(fitcp[[i]],  c("Intercept_3 < Intercept_4"))
        hyp6[4,] <-  hypothesis(fitcp[[i]],  c("Intercept_4 < Intercept_5"))
        hyp6[5,] <-  hypothesis(fitcp[[i]],  c("Intercept_5 < Intercept_6"))
        hyp6[6,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 > Intercept_2"))
        hyp6[7,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 > Intercept_3"))
        hyp6[8,] <-  hypothesis(fitcp[[i]],  c("Intercept_3 > Intercept_4"))
        hyp6[9,] <-  hypothesis(fitcp[[i]],  c("Intercept_4 > Intercept_5"))
        hyp6[10,] <-  hypothesis(fitcp[[i]],  c("Intercept_5 > Intercept_6"))

        write.table(hyp6, paste0(path, "hypothesis6strata.txt"), quote= F)
        ds6.1 <- dplot(df, nstrata=6, "six_strata")  + ggtitle("E - five changepoint")
        ds6.2 <- dplot2(df, nstrata=6, "six_strata") + ggtitle("E - five changepoint")
        vp6s <- ggbetweenstats(df, six_strata, Ds)    + ylab(expression(italic(d[s]))) +th_plot3

     } else if (i == 6){
      
       hyp7 <- data.frame(matrix(ncol = 7, nrow = 12))%>%
       set_colnames(.,c('hypothesis','mean','sd','lower','upper','prob','BF'))
       
       hyp7[1,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 < Intercept_2"))
       hyp7[2,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 < Intercept_3"))
       hyp7[3,] <-  hypothesis(fitcp[[i]],  c("Intercept_3 < Intercept_4"))
       hyp7[4,] <-  hypothesis(fitcp[[i]],  c("Intercept_4 < Intercept_5"))
       hyp7[5,] <-  hypothesis(fitcp[[i]],  c("Intercept_5 < Intercept_6"))
       hyp7[6,] <-  hypothesis(fitcp[[i]],  c("Intercept_6 < Intercept_7"))
       hyp7[7,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 > Intercept_2"))
       hyp7[8,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 > Intercept_3"))
       hyp7[9,] <-  hypothesis(fitcp[[i]],  c("Intercept_3 > Intercept_4"))
       hyp7[10,] <-  hypothesis(fitcp[[i]],  c("Intercept_4 > Intercept_5"))
       hyp7[11,] <-  hypothesis(fitcp[[i]],  c("Intercept_5 > Intercept_6"))
       hyp7[12,] <-  hypothesis(fitcp[[i]],  c("Intercept_6 > Intercept_7"))

        write.table(hyp7, paste0(path, "hypothesis7strata.txt"), quote= F)
        ds7.1 <- dplot(df, nstrata=7, "seven_strata")  + ggtitle("F - six changepoint")
        ds7.2 <- dplot2(df, nstrata=7, "seven_strata") + ggtitle("F - six changepoint")
        vp7s <- ggbetweenstats(df, seven_strata, Ds)  + ylab(expression(italic(d[s]))) + th_plot3

     } else if (i == 7){
        
        hyp8 <- data.frame(matrix(ncol = 7, nrow = 14))%>%
         set_colnames(.,c('hypothesis','mean','sd','lower','upper','prob','BF'))
       
        hyp8[1,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 < Intercept_2"))
        hyp8[2,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 < Intercept_3"))
        hyp8[3,] <-  hypothesis(fitcp[[i]],  c("Intercept_3 < Intercept_4"))
        hyp8[4,] <-  hypothesis(fitcp[[i]],  c("Intercept_4 < Intercept_5"))
        hyp8[5,] <-  hypothesis(fitcp[[i]],  c("Intercept_5 < Intercept_6"))
        hyp8[6,] <-  hypothesis(fitcp[[i]],  c("Intercept_6 < Intercept_7"))
        hyp8[7,] <-  hypothesis(fitcp[[i]],  c("Intercept_7 < Intercept_8"))
        hyp8[8,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 > Intercept_2"))
        hyp8[9,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 > Intercept_3"))
        hyp8[10,] <-  hypothesis(fitcp[[i]],  c("Intercept_3 > Intercept_4"))
        hyp8[11,] <-  hypothesis(fitcp[[i]],  c("Intercept_4 > Intercept_5"))
        hyp8[12,] <-  hypothesis(fitcp[[i]],  c("Intercept_5 > Intercept_6"))
        hyp8[13,] <-  hypothesis(fitcp[[i]],  c("Intercept_6 > Intercept_7"))
        hyp8[14,] <-  hypothesis(fitcp[[i]],  c("Intercept_7 > Intercept_8"))
    
        write.table(hyp8, paste0(path, "hypothesis8strata.txt"), quote= F)
        ds8.1 <- dplot(df, nstrata=8, "eight_strata")  + ggtitle("H - seven changepoint")
        ds8.2 <- dplot2(df, nstrata=8, "eight_strata") + ggtitle("H - seven changepoint")
        vp8s <- ggbetweenstats(df, eight_strata, Ds, palette = "Paired")  + th_plot3

     } else if (i == 8){
       
        hyp9 <- data.frame(matrix(ncol = 7, nrow = 16))%>%
          set_colnames(.,c('hypothesis','mean','sd','lower','upper','prob','BF'))
        
        hyp9[1,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 < Intercept_2"))
        hyp9[2,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 < Intercept_3"))
        hyp9[3,] <-  hypothesis(fitcp[[i]],  c("Intercept_3 < Intercept_4"))
        hyp9[4,] <-  hypothesis(fitcp[[i]],  c("Intercept_4 < Intercept_5"))
        hyp9[5,] <-  hypothesis(fitcp[[i]],  c("Intercept_5 < Intercept_6"))
        hyp9[6,] <-  hypothesis(fitcp[[i]],  c("Intercept_6 < Intercept_7"))
        hyp9[7,] <-  hypothesis(fitcp[[i]],  c("Intercept_7 < Intercept_8"))
        hyp9[8,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 > Intercept_2"))
        hyp9[9,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 > Intercept_3"))
        hyp9[10,] <-  hypothesis(fitcp[[i]],  c("Intercept_3 > Intercept_4"))
        hyp9[11,] <-  hypothesis(fitcp[[i]],  c("Intercept_4 > Intercept_5"))
        hyp9[12,] <-  hypothesis(fitcp[[i]],  c("Intercept_5 > Intercept_6"))
        hyp9[13,] <-  hypothesis(fitcp[[i]],  c("Intercept_6 > Intercept_7"))
        hyp9[14,] <-  hypothesis(fitcp[[i]],  c("Intercept_7 > Intercept_8"))
        hyp9[15,] <-  hypothesis(fitcp[[i]],  c("Intercept_8 < Intercept_9"))
        hyp9[16,] <-  hypothesis(fitcp[[i]],  c("Intercept_8 > Intercept_9"))
    
        write.table(hyp9, paste0(path, "hypothesis9strata.txt"), quote= F)
        ds8.1 <- dplot(df, nstrata=9, "nine_strata")  + ggtitle("I - eight changepoint")
        ds8.2 <- dplot2(df, nstrata=9, "nine_strata") + ggtitle("I - eight changepoint")
        vp9s <- ggbetweenstats(df, nine_strata, Ds, palette = "Paired")  + th_plot3

     } else if (i == 9){
       
        hyp10 <- data.frame(matrix(ncol = 7, nrow = 18))%>%
          set_colnames(.,c('hypothesis','mean','sd','lower','upper','prob','BF'))
        
        hyp9[1,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 < Intercept_2"))
        hyp9[2,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 < Intercept_3"))
        hyp9[3,] <-  hypothesis(fitcp[[i]],  c("Intercept_3 < Intercept_4"))
        hyp9[4,] <-  hypothesis(fitcp[[i]],  c("Intercept_4 < Intercept_5"))
        hyp9[5,] <-  hypothesis(fitcp[[i]],  c("Intercept_5 < Intercept_6"))
        hyp9[6,] <-  hypothesis(fitcp[[i]],  c("Intercept_6 < Intercept_7"))
        hyp9[7,] <-  hypothesis(fitcp[[i]],  c("Intercept_7 < Intercept_8"))
        hyp9[8,] <-  hypothesis(fitcp[[i]],  c("Intercept_1 > Intercept_2"))
        hyp9[9,] <-  hypothesis(fitcp[[i]],  c("Intercept_2 > Intercept_3"))
        hyp9[10,] <-  hypothesis(fitcp[[i]],  c("Intercept_3 > Intercept_4"))
        hyp9[11,] <-  hypothesis(fitcp[[i]],  c("Intercept_4 > Intercept_5"))
        hyp9[12,] <-  hypothesis(fitcp[[i]],  c("Intercept_5 > Intercept_6"))
        hyp9[13,] <-  hypothesis(fitcp[[i]],  c("Intercept_6 > Intercept_7"))
        hyp9[14,] <-  hypothesis(fitcp[[i]],  c("Intercept_7 > Intercept_8"))
        hyp9[15,] <-  hypothesis(fitcp[[i]],  c("Intercept_8 < Intercept_9"))
        hyp9[16,] <-  hypothesis(fitcp[[i]],  c("Intercept_8 > Intercept_9"))
        hyp9[17,] <-  hypothesis(fitcp[[i]],  c("Intercept_9 < Intercept_10"))
        hyp9[18,] <-  hypothesis(fitcp[[i]],  c("Intercept_9 > Intercept_10"))

        write.table(hyp10, paste0(path, "hypothesis10strata.txt"), quote= F)
        ds9.1 <- dplot(df, nstrata=10, "ten_strata") + ggtitle("J - nine changepoint")
        ds9.2 <- dplot2(df, nstrata=10, "ten_strata") + ggtitle("J - nine changepoint")
        vp10s <- ggbetweenstats(df, ten_strata, Ds, palette = "Paired")  + th_plot3

     }
    }
    
    writeLines("\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    writeLines("\n\n exporting some more plots \n\n")
    writeLines("\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    ds_pos_files <- mget(ls(pattern="ds[0-9].1"))
    ds_order_files <- mget(ls(pattern="ds[0-9].2"))
    viobox_files <-  mget(ls(pattern="vp[0-9]*.s"))

    for(i in 1:length(ds_pos_files)){
        pdf(file=paste0(path,"plot_ds_position_",i,"changepoint.pdf"),8,6)
        print(plot_grid(ds_pos_files[[i]]))
 
        pdf(file=paste0(path,"plot_ds_order_",i,"changepoint.pdf"),8,6)
        print(plot_grid(ds_order_files[[i]]))
    }

    for(i in 1:length(viobox_files)){
        pdf(file=paste0(path,"viobox_ds__",i,"strata.pdf"),8,6)
        print(plot_grid(viobox_files[[i]]))
    }

    
    ############################################################################
    #finally construct some combined plot:
    pdf(file=paste0(path,"plot_dS_all_position.pdf"),8,24)
    print(plot_grid(plotlist=ds_pos_files, ncol=1))
    dev.off()
    
    pdf(file=paste0(path,"plot_ds_along_order.pdf"),8,22)
    print(plot_grid(plotlist=ds_order_files, ncol=1))
    dev.off()

    pdf(file=paste0(path,"plot_viobox_all_strata.pdf"),8,24)
    print(plot_grid(plotlist=viobox_files, ncol=1, 
             label_size = 7,
             hjust = -0.5, vjust = -0.5,
            labels = "AUTO"))
    dev.off()

    plot <- list()
    for (i in 1:maxchgp) {
    plot[[i]] <- plot_grid(print(figcp[[i]]) + th_plot3, ncol = 1)
    }
    
    pdf(file=paste0(path,"plot_all_changepoint.pdf"), 10, 20)
    print(plot_grid(plotlist=plot, ncol=1, labels = "AUTO"))
    dev.off()

    #finally: 
    if(is_anc=="YES"){ 
    s2.anc.h1 <- select(df, gene, geneX, two_strata)
    s3.anc.h1 <- select(df, gene, geneX, three_strata)
    s4.anc.h1 <- select(df, gene, geneX, four_strata)
    s5.anc.h1 <- select(df, gene, geneX, five_strata)
    s6.anc.h1 <- select(df, gene, geneX, five_strata)
    s6.anc.h1 <- select(df, gene, geneX, six_strata)
    s7.anc.h1 <- select(df, gene, geneX, seven_strata)
    s8.anc.h1 <- select(df, gene, geneX, eight_strata)
    
    s2.h1.h2 <- select(df, geneX, geneY.x, two_strata)
    s3.h1.h2 <- select(df, geneX, geneY.x, three_strata)
    s4.h1.h2 <- select(df, geneX, geneY.x, four_strata)
    s5.h1.h2 <- select(df, geneX, geneY.x, five_strata)
    s6.h1.h2 <- select(df, geneX, geneY.x, five_strata)
    s6.h1.h2 <- select(df, geneX, geneY.x, six_strata)
    s7.h1.h2 <- select(df, geneX, geneY.x, seven_strata)
    s8.h1.h2 <- select(df, geneX, geneY.x, eight_strata)
   write.table(s2.anc.h1,paste0(path,"classif.s2.ancestral.haplo1"),
        quote=F,row.names=F,col.names=F,sep="\t")
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
    s2.h1.h2 <- select(df, gene, geneY.y, two_strata)
    s3.h1.h2 <- select(df, gene, geneY.y, three_strata)
    s4.h1.h2 <- select(df, gene, geneY.y, four_strata)
    s5.h1.h2 <- select(df, gene, geneY.y, five_strata)
    s6.h1.h2 <- select(df, gene, geneY.y, five_strata)
    s6.h1.h2 <- select(df, gene, geneY.y, six_strata)
    s7.h1.h2 <- select(df, gene, geneY.y, seven_strata)
    s8.h1.h2 <- select(df, gene, geneY.y, eight_strata)
    
    }
    
    write.table(s2.h1.h2,paste0(path,"classif.s2.haplo1.haplo2"),
        quote=F,row.names=F,col.names=F,sep="\t")
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
    
     
    writeLines("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    writeLines("\n analyses finished !! \n\n")
    writeLines("\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    writeLines("\n exporting Rsession !! \n\n")

#save.image( file = "02_results/modelcomp/noprior/changepoint_analysis.RData")

}
