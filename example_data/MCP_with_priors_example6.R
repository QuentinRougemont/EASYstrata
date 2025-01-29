#!/usr/bin/env Rscript

#Author QR
#Purpose: perform changepoint analysis with prior 

#--------------- check if library are installed -------------------------------#
libs <- c('mcp','dplyr','ggplot2','magrittr','cowplot','ggstatsplot' )
install.packages(setdiff(libs, rownames(installed.packages())), repos="https://cloud.r-project.org" )
    
#---------------- load libraries ---------------------------------------------#
invisible(lapply(libs, suppressWarnings(suppressMessages(suppressPackageStartupMessages(library))), character.only = TRUE))

#---- load data ---- # 
dsfile <- "02_results/dS.values.forchangepoint.txt"
df <- na.omit(df)

#create dir if not present:
if (!dir.exists("02_results/modelcomp")){
dir.create("02_results/modelcomp")
}
 
 
## use my usual theme:
## ------------------ GGPLOTCUSTOMISATION ------------------------------------------------##
th_plot <-theme(axis.title.x=element_text(size=14, family="Helvetica",face="bold"),
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
#   perform the changepoint analyis here: 
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

###########################################################################
#### changepoint with prior example
###########################################################################

##using prior: 

if (!dir.exists("02_results/modelcomp/withprior")){
  dir.create("02_results/modelcomp/withprior")
}


path="02_results/modelcomp/withprior"
#unfortunately it seems that this needs to be run manually given the dependance on user assumed priors

prior2=list(cp_1="dnorm(100, 10)",cp_2="dnorm(330, 10)")
priorfit_2cp = mcp(model2cp,
       data = df,
       par_x = "orderchp",
       iter = 8e3,
       adapt = 1.5e3,
       prior = prior2,
       chains = 5,
       cores = 5 )

prior3=list(cp_1="dnorm(100, 10)",cp_2="dnorm(300, 10)",cp_3="dnorm(331, 10)")
priorfit_3cp = mcp(model3cp,
       data = df,
       par_x = "orderchp",
       iter = 8e3,
       adapt = 1.5e3,
       prior = prior3,
       chains = 5,
       cores = 5 )

prior4=list(cp_1="dnorm(85, 10)",cp_2="dnorm(130, 10)",cp_3="dnorm(280, 10)",cp_4="dnorm(331,10)")
priorfit_4cp = mcp(model4cp,
       data = df,
       par_x = "orderchp",
       iter = 8e3,
       adapt = 1.5e3,
       prior=prior4,
       chains = 5,
       cores = 5 )

prior5=list(cp_1="dnorm(85,10)", cp_2="dnorm(130,10)", cp_3="dnorm(250, 10)", 
            cp_4="dnorm(285,5)", cp_5="dnorm(327,5)")
priorfit_5cp = mcp(model5cp,
       data = df,
       par_x = "orderchp",
       iter = 8e3,
       adapt = 1.5e3,
       prior=prior5,
       chains = 5,
       cores = 5 )

prior6=list(cp_1="dnorm(90,10)", cp_2="dnorm(125,10)", cp_3="dnorm(230,5)", 
            cp_4="dnorm(265,5)", cp_5="dnorm(285,5)", cp_6="dnorm(335,5)")

bck.prior6 <- priorfit_6cp


prior6=list(cp_1="dnorm(90,10)", cp_2="dnorm(125,10)", cp_3="dnorm(220,5)", 
            cp_4="dnorm(255,5)", cp_5="dnorm(285,5)", cp_6="dnorm(335,5)")

priorfit_6cp = mcp(model6cp,
       data = df,
       par_x = "orderchp",
       iter = 8e3,
       adapt = 1.5e3,
       prior = prior6,
       chains = 5,
       cores = 5 )

prior7=list(cp_1="dnorm(90, 10)",cp_2="dnorm(125, 10)",cp_3="dnorm(170, 10)",cp_4="dnorm(200, 10)",
            cp_5="dnorm(250,5)",cp_6="dnorm(285,5)",cp_7="dnorm(335,5)")

bck.prior7 <- priorfit_7cp

prior7=list(cp_1="dnorm(90, 10)",cp_2="dnorm(125, 10)",cp_3="dnorm(170, 10)",cp_4="dnorm(210, 10)",
            cp_5="dnorm(250,5)",cp_6="dnorm(285,5)",cp_7="dnorm(335,5)")

priorfit_7cp = mcp(model7cp,
       data = df,
       par_x = "orderchp",
       iter = 8e3,
       adapt = 1.5e3,
       prior = prior7,
       chains = 5,
       cores = 5 )

prior8=list(cp_1="dnorm(90, 10)",cp_2="dnorm(125, 10)",cp_3="dnorm(170, 10)",cp_4="dnorm(200, 10)",
            cp_5="dnorm(250,5)", cp_6="dnorm(280,5)",cp_7="dnorm(300,5)",cp_8="dnorm(335,5)")
priorfit_8cp = mcp(model8cp,
       data = df,
       par_x = "orderchp",
       iter = 8e3,
       adapt = 1.5e3,
       prior = prior8,
       chains = 5,
       cores = 5 )

prior9=list(cp_1="dnorm(90, 10)",cp_2="dnorm(125, 10)",cp_3="dnorm(150, 10)",cp_4="dnorm(170, 10)",
            cp_5="dnrom(200,5)", cp_6="dnorm(250,5)",cp_7="dnorm(280,5)",cp_8="dnorm(305,5)",cp_9="dnorm(335,5)"))
priorfit_9cp = mcp(model8cp, data = df, par_x = "orderchp",   iter = 8e3, adapt = 1.5e3, prior = prior9, chains = 5, cores = 5 )
#
##summarise them all:
Fig2cpt <- plotcp(priorfit_2cp, "Posterior fit - 2 changepoint") 
Fig3cpt <- plotcp(priorfit_3cp, "Posterior fit - 3 changepoint") 
Fig4cpt <- plotcp(priorfit_4cp, "Posterior fit - 4 changepoint") 
Fig5cpt <- plotcp(priorfit_5cp, "Posterior fit - 5 changepoint") 
Fig6cpt <- plotcp(priorfit_6cp, "Posterior fit - 6 changepoint") 
Fig7cpt <- plotcp(priorfit_7cp, "Posterior fit - 7 changepoint") 
Fig8cpt <- plotcp(priorfit_8cp, "Posterior fit - 8 changepoint") 
Fig9cpt <- plotcp(priorfit_9cp, "Posterior fit - 9 changepoint") 

pdf(file = paste0(path, "changepoint_comparison_withpriors.pdf"), 12,18)
plot_grid(
   Fig2cpt, 
   Fig3cpt, 
   Fig4cpt, 
   Fig5cpt, 
   Fig6cpt, 
   Fig7cpt, 
   Fig8cpt,
   labels = "AUTO", ncol = 1)
dev.off()


#----------------  a bit of model choice --------------------------------------
#to be considered precautionously!
priorfit_1cp$loo <- loo(priorfit_1cp)
priorfit_2cp$loo <- loo(priorfit_2cp)
priorfit_3cp$loo <- loo(priorfit_3cp)
priorfit_4cp$loo <- loo(priorfit_4cp)
priorfit_5cp$loo <- loo(priorfit_5cp)
priorfit_6cp$loo <- loo(priorfit_6cp)
priorfit_7cp$loo <- loo(priorfit_7cp)
priorfit_8cp$loo <- loo(priorfit_8cp)

m.choice <- loo::loo_compare(priorfit_1cp$loo, priorfit_2cp$loo, priorfit_3cp$loo, 
                             priorfit_4cp$loo, priorfit_5cp$loo, priorfit_6cp$loo,
                             priorfit_7cp$loo, priorfit_8cp$loo)


# ----- rule of thumb - to be used as a helper for model choice ---- : 
# you should read the manual (https://lindeloev.github.io/mcp/articles/comparison.html)
# and associated paper to help in choosing a model.
#finally compute the weight:
loo_list = list(priorfit_1cp$loo, priorfit_2cp$loo, priorfit_3cp$loo, 
                priorfit_4cp$loo, priorfit_5cp$loo, priorfit_6cp$loo,
                priorfit_7cp$loo, priorfit_8cp$loo)


weights <- loo::loo_model_weights(loo_list, method="pseudobma")

##export:
write.table(weights,"02_results/modelcomp/model_weights_with_priors.txt",quote=F, col.names=("weights"))
#
##model choice attempt: 
write.table(m.choice ,"02_results/modelcomp/model.choice.priors.txt", quote = F)
#
m2 <- summary(priorfit_2cp)
m3 <- summary(priorfit_3cp)
m4 <- summary(priorfit_4cp)
m5 <- summary(priorfit_5cp)
m6 <- summary(priorfit_6cp)
m7 <- summary(priorfit_7cp)
m8 <- summary(priorfit_8cp)
#
### make some plot based on mean changes:
##getting mean-value along the X of the first changepoint 
##2-changepoint model:
cp1m2 <- m2$mean[1]
cp2m2 <- m2$mean[2]
#
##3-changepoint :
cp1m3 <- m3$mean[1]
cp2m3 <- m3$mean[2]
cp3m3 <- m3$mean[3]
#
##4 - changepoint model :
cp1m4 <- m4$mean[1]
cp2m4 <- m4$mean[2]
cp3m4 <- m4$mean[3]
cp4m4 <- m4$mean[4]
#
##5 - changepoint model :
cp1m5 <- m5$mean[1] ##alternatively: summary(fit_5st)[2][[1]][1]
cp2m5 <- m5$mean[2]
cp3m5 <- m5$mean[3]
cp4m5 <- m5$mean[4]
cp5m5 <- m5$mean[5]
#
##6 - changepoint model :
cp1m6 <- m6$mean[1] ##alternatively: summary(fit_5st)[2][[1]][1]
cp2m6 <- m6$mean[2]
cp3m6 <- m6$mean[3]
cp4m6 <- m6$mean[4]
cp5m6 <- m6$mean[5]
cp6m6 <- m6$mean[6]
#
##7 - changepoint model :
cp1m7 <- m7$mean[1] ##alternatively: summary(fit_5st)[2][[1]][1]
cp2m7 <- m7$mean[2]
cp3m7 <- m7$mean[3]
cp4m7 <- m7$mean[4]
cp5m7 <- m7$mean[5]
cp6m7 <- m7$mean[6]
cp7m7 <- m7$mean[7]
#
##8 - changepoint model :
cp1m8 <- m8$mean[1] ##alternatively: summary(fit_5st)[2][[1]][1]
cp2m8 <- m8$mean[2]
cp3m8 <- m8$mean[3]
cp4m8 <- m8$mean[4]
cp5m8 <- m8$mean[5]
cp6m8 <- m8$mean[6]
cp7m8 <- m8$mean[7]
cp8m8 <- m8$mean[8]
#
##9 - changepoint model :
cp1m9 <- m9$mean[1] ##alternatively: summary(fit_5st)[2][[1]][1]
cp2m9 <- m9$mean[2]
cp3m9 <- m9$mean[3]
cp4m9 <- m9$mean[4]
cp5m9 <- m9$mean[5]
cp6m9 <- m9$mean[6]
cp7m9 <- m9$mean[7]
cp8m9 <- m9$mean[8]
cp9m9 <- m9$mean[9]
#
#
df$three_strata <- ifelse(df$orderchp < cp1m2, "strata1", 
                        ifelse(df$orderchp > cp1m2 & df$orderchp < cp2m2, "strata2", "strata3"))

df$four_strata <- ifelse(df$orderchp < cp1m3, "strata1", 
                  ifelse(df$orderchp > cp3m3, "strata4", 
                  ifelse(df$orderchp > cp2m3 & df$orderchp < cp3m3, "strata3", "strata2")))

df$five_strata <- ifelse(df$orderchp < cp1m4, "strata1", 
                  ifelse(df$orderchp > cp4m4, "strata5",
                  ifelse(df$orderchp > cp1m4 & df$orderchp < cp2m4, "strata2",
                  ifelse(df$orderchp > cp2m4 & df$orderchp < cp3m4, "strata3","strata4"))))

df$six_strata <- ifelse(df$orderchp < cp1m5, "strata1", 
                 ifelse(df$orderchp > cp5m5, "strata6",
                 ifelse(df$orderchp > cp1m5 & df$orderchp < cp2m5, "strata2",
                 ifelse(df$orderchp > cp2m5 & df$orderchp < cp3m5, "strata3",
                 ifelse(df$orderchp > cp3m5 & df$orderchp < cp4m5, "strata4", "strata5" )))))

df$seven_strata <- ifelse(df$orderchp < cp1m6, "strata1", 
                 ifelse(df$orderchp > cp6m6, "strata7",
                 ifelse(df$orderchp > cp1m6 & df$orderchp < cp2m6, "strata2",
                 ifelse(df$orderchp > cp2m6 & df$orderchp < cp3m6, "strata3",
                 ifelse(df$orderchp > cp3m6 & df$orderchp < cp4m6, "strata4",
                 ifelse(df$orderchp > cp4m6 & df$orderchp < cp5m6, "strata5","strata6"))))))

df$eight_strata <- ifelse(df$orderchp < cp1m7, "strata1", 
                  ifelse(df$orderchp > cp7m7, "strata8",
                  ifelse(df$orderchp > cp1m7 & df$orderchp < cp2m7, "strata2",
                  ifelse(df$orderchp > cp2m7 & df$orderchp < cp3m7, "strata3",
                  ifelse(df$orderchp > cp3m7 & df$orderchp < cp4m7, "strata4",
                  ifelse(df$orderchp > cp4m7 & df$orderchp < cp5m7, "strata5",
                  ifelse(df$orderchp > cp5m7 & df$orderchp < cp6m7, "strata6",
                           "strata7" )))))))

df$nine_strata <- ifelse(df$orderchp < cp1m8, "strata1", 
                  ifelse(df$orderchp > cp8m8, "strata9",
                  ifelse(df$orderchp > cp1m8 & df$orderchp < cp2m8, "strata2",
                  ifelse(df$orderchp > cp2m8 & df$orderchp < cp3m8, "strata3",
                  ifelse(df$orderchp > cp3m8 & df$orderchp < cp4m8, "strata4",
                  ifelse(df$orderchp > cp4m8 & df$orderchp < cp5m8, "strata5",
                  ifelse(df$orderchp > cp5m8 & df$orderchp < cp6m8, "strata6",
                  ifelse(df$orderchp > cp6m8 & df$orderchp < cp7m8, "strata7",
                           "strata8" ))))))))

write.table(df,paste0(path, "df_with_prior.txt"),quote=F,row.names=F,col.names=T,sep="\t")
#
vp2 <- vplot(df, 2, "two_strata")
vp3 <- vplot(df, 3, "three_strata")
vp4 <- vplot(df, 4, "four_strata")
vp5 <- vplot(df, 5, "five_strata")
vp6 <- vplot(df, 6, "six_strata")
vp7 <- vplot(df, 7, "seven_strata")
vp8 <- vplot(df, 8, "eight_strata")


##plot them all
pdf(file = paste0(path, "vioplot_ds_strata_distribution_priors.pdf", 12,18)
plot_grid(
    vp3,
    vp4,
    vp5,
    vp6,
    vp7,
    vp8, 
    labels = "AUTO", ncol = 1)
dev.off()
#
##same data but with a bit more stats:
vp3s <- ggbetweenstats(df, three_strata, Ds)  + th_plot3
vp4s <- ggbetweenstats(df, four_strata, Ds)   + th_plot3
vp5s <- ggbetweenstats(df, five_strata, Ds)   + th_plot3
vp6s <- ggbetweenstats(df, six_strata, Ds)    + th_plot3
vp7s <- ggbetweenstats(df, seven_strata, Ds)  + th_plot3
vp8s <- ggbetweenstats(df, eight_strata, Ds, palette = "paired")  + th_plot3
#
pdf(file = paste0(path, "viobox_ds_strata_distribution_priors.pdf"), 12,28)
plot_grid(
   vp3s,
   vp4s,
   vp5s,
   vp6s,
   vp7s,
   vp8s, 
   labels = "AUTO", 
   ncol = 1)
dev.off()

pdf(file = paste0(path, "Figure.pdf"), 24,14)
plot_grid(
   Fig5cpt, vp5s, 
   Fig6cpt, vp6s, 
   Fig7cpt, vp7s,
   labels = "AUTO", 
   ncol = 2)
dev.off()

write.table(m1, paste0(path, "summary1strata.priors.txt"), quote = F)
write.table(m2, paste0(path, "summary2strata.priors.txt"), quote = F)
write.table(m3, paste0(path, "summary3strata.priors.txt"), quote = F)
write.table(m4, paste0(path, "summary4strata.priors.txt"), quote = F)
write.table(m5, paste0(path, "summary5strata.priors.txt"), quote = F)
write.table(m6, paste0(path, "summary6strata.priors.txt"), quote = F)
write.table(m7, paste0(path, "summary7strata.priors.txt"), quote = F)
write.table(m8, paste0(path, "summary8strata.priors.txt"), quote = F)

ds3 <- dplot(df, nstrata=3, "three_strata")
ds4 <- dplot(df, nstrata=4, "four_strata")
ds5 <- dplot(df, nstrata=5, "five_strata")
ds6 <- dplot(df, nstrata=6, "six_strata")
ds7 <- dplot(df, nstrata=7, "seven_strata")
ds8 <- dplot(df, nstrata=8, "eight_strata")

pdf(file=paste0(path, "plot_dS_all_position_priors.pdf" ),8,12)
plot_grid(
    ds3, 
    ds4, 
    ds5, 
    ds6, 
    ds7, 
    ds8, 
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


pdf(file=paste0(path, "plot_Ds_along_order_priors.pdf" ),8,12)
plot_grid(
    ds3, 
    ds4, 
    ds5, 
    ds6, 
    ds7, 
    ds8, 
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

write.table(s3.anc.h1,"02_results/modelcomp/classif.s3.ancestral.haplo1.priors",
    quote=F,row.names=F,col.names=F,sep="\t")
write.table(s4.anc.h1,"02_results/modelcomp/classif.s4.ancestral.haplo1.priors",
    quote=F,row.names=F,col.names=F,sep="\t")
write.table(s5.anc.h1,"02_results/modelcomp/classif.s5.ancestral.haplo1.priors",
    quote=F,row.names=F,col.names=F,sep="\t")
write.table(s6.anc.h1,"02_results/modelcomp/classif.s6.ancestral.haplo1.priors",
    quote=F,row.names=F,col.names=F,sep="\t")
write.table(s7.anc.h1,"02_results/modelcomp/classif.s7.ancestral.haplo1.priors",
    quote=F,row.names=F,col.names=F,sep="\t")
write.table(s8.anc.h1,"02_results/modelcomp/classif.s8.ancestral.haplo1",
    quote=F,row.names=F,col.names=F,sep="\t")

write.table(s3.h1.h2,"02_results/modelcomp/classif.s3.haplo1.haplo2.priors",
    quote=F,row.names=F,col.names=F,sep="\t")
write.table(s4.h1.h2,"02_results/modelcomp/classif.s4.haplo1.haplo2.priors",
    quote=F,row.names=F,col.names=F,sep="\t")
write.table(s5.h1.h2,"02_results/modelcomp/classif.s5.haplo1.haplo2.priors",
    quote=F,row.names=F,col.names=F,sep="\t")
write.table(s6.h1.h2,"02_results/modelcomp/classif.s6.haplo1.haplo2.priors",
    quote=F,row.names=F,col.names=F,sep="\t")
write.table(s7.h1.h2,"02_results/modelcomp/classif.s7.haplo1.haplo2.priors",
    quote=F,row.names=F,col.names=F,sep="\t")
write.table(s8.h1.h2,"02_results/modelcomp/classif.s8.haplo1.haplo2",
    quote=F,row.names=F,col.names=F,sep="\t")
