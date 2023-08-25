library("ggplot2")
library("ggpubr")
library("dunn.test")

df = read.csv( "galleria_infection.csv" )
df_nc = df[df$Strain != "Ctrl",]
# df = subset(df,duplicated(df$Strain) | duplicated(df$Strain, fromLast=TRUE))

df = df[df$"remove"=="n",] # Data frame will be cleared of all the strains to be removed

new_col = c()
setty = c()

for (row_num in 1:nrow(df)) {
    # Creating the line
    days = c(0,1,2,3,4,5,6,7,8)
    live = c(10, df[row_num,6], df[row_num,7], df[row_num,8], df[row_num,9], df[row_num,10], df[row_num,11], df[row_num,12], df[row_num,13])

    intercept = 10
    pazuzu = lm(I(live - intercept) ~ 0 + days )

    new_thing = summary(pazuzu)$coefficients["days","Estimate"]

    new_col = c(new_col, new_thing)

}

df$slopes = new_col
write.csv(df, file = "galleria_infection_corr.csv")

new_df = data.frame(matrix(ncol = 13, nrow = 0))
new_colnames = colnames(df)
colnames(new_df) = c("strain","D0","D1","D2","D3","D4","D5","D6","D7","D8","auc", "slope", "color")

for (strains in unique(df$Strain)){ 

    # Make step graphs for each strain from df_nc (I will also probably have to make control graphs too.
    strain_df = df[df$Strain == strains,]
    days = c(0,1,2,3,4,5,6,7,8)
    colors = c("red", "blue", "green", "purple")

    # Piggy backing my analysis stuff in here
    werms_per_day = c(30,sum(strain_df[,6]), sum(strain_df[,7]), sum(strain_df[,8]), sum(strain_df[,9]), sum(strain_df[,10]), sum(strain_df[,11]), sum(strain_df[,12]), sum(strain_df[,13]))
    auc_val = sum(werms_per_day)

    intercept = 30
    pazuzu = lm(I(werms_per_day - intercept) ~ 0 + days )

    slope = summary(pazuzu)$coefficients["days","Estimate"]

    color = "1"

    if ( strains == "73" ){ 

        color = "4"

    } else if ( strains == "85" ) {

        color = "3"

    } else if ( strains == "Ctrl" ) { 

        color = "2"

    }

    new_row = c(strains, 30, as.numeric(sum(strain_df[,6])), as.numeric(sum(strain_df[,7])), as.numeric(sum(strain_df[,8])), as.numeric(sum(strain_df[,9])), as.numeric(sum(strain_df[,10])), as.numeric(sum(strain_df[,11])), as.numeric(sum(strain_df[,12])), as.numeric(sum(strain_df[,13])), auc_val, slope, color )
    new_df[nrow(new_df)+1,] = new_row
 
    live_1 = c(10, strain_df[1,6], strain_df[1,7], strain_df[1,8], strain_df[1,9], strain_df[1,10], strain_df[1,11], strain_df[1,12], strain_df[1,13])
    live_2 = c(10, strain_df[2,6], strain_df[2,7], strain_df[2,8], strain_df[2,9], strain_df[2,10], strain_df[2,11], strain_df[2,12], strain_df[2,13])
    live_3 = c(10, strain_df[3,6], strain_df[3,7], strain_df[3,8], strain_df[3,9], strain_df[3,10], strain_df[3,11], strain_df[3,12], strain_df[3,13])

    x = ggplot( ) +
        geom_step( aes( x = days , y = as.numeric(new_row[2:10])) ) + 
        lims(x = c(0,8), y =c(0,30) ) +
        geom_abline(intercept = 30, slope = slope, colour = "red" ) + 
        theme_bw() +
        theme(  
                    # axis.text = element_blank(), 
                    # axis.title = element_blank(), 
                    # legend.text = element_blank(),
                    # axis.text.x = element_blank(),
                    # panel.grid.major.x = element_blank(),
                    panel.grid.minor.x = element_blank(), 
                    axis.ticks.x = element_blank(), 
                    axis.ticks.y = element_blank(), 
                    # legend.position = "none",
                    # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                    panel.border = element_blank() ) + 
        labs( x = "Days", y = "Living Worms" )

    ggsave(paste0("strain_kms/",strains,"_TOTAL_km_curve.png"), height = 7, width = 7)

    y = ggplot( ) + 
        geom_step( aes(x = days, y = live_1, color = colors[1]) ) + 
        geom_text(aes(label = strain_df[1,2], x = 8, y = live_1[9], color = colors[1] )) + 
        geom_step( aes(x = days, y = live_2, color = colors[2]) ) +
        geom_text(aes(label = strain_df[2,2], x = 8, y = live_2[9], color = colors[2] )) + 
        geom_step( aes(x = days, y = live_3, color = colors[3]) ) + 
        geom_text(aes(label = strain_df[3,2], x = 8, y = live_3[9], color = colors[3] )) + 
        lims(x = c(0,8), y =c(0,10) ) +

        theme(legend.position = "none")

    # if ( (strains == "Ctrl") |
    #      (strains == "73") |
    #      (strains == "85") |
    #      (strains == "120") |
    #      (strains == "123") |
    #      (strains == "125") |
    #      (strains == "126") |
    #      (strains == "127") |
    #      (strains == "128") |
    #      (strains == "129") |
    #      (strains == "130") |
    #      (strains == "137") |
    #      (strains == "138") |
    #      (strains == "139") |
    #      (strains == "140") |
    #      (strains == "142") |
    #      (strains == "143") |
    #      (strains == "144") |
    #      (strains == "147") |
    #      (strains == "148") |
    #      (strains == "152") |
    #      (strains == "889") |
    #      (strains == "890") |
    #      (strains == "891") |
    #      (strains == "892") |
    #      (strains == "893") |
    #      (strains == "894") |
    #      (strains == "899") |
    #      (strains == "902") |
    #      (strains == "903") |
    #      (strains == "912") |
    #      (strains == "914") |
    #      (strains == "915") |
    #      (strains == "923") |
    #      (strains == "931") |
    #      (strains == "936") |
    #      (strains == "953") ){

    #     live_4 = c(10, strain_df[4,6], strain_df[4,7], strain_df[4,8], strain_df[4,9], strain_df[4,10], strain_df[4,11], strain_df[4,12], strain_df[4,13])
    #     y + geom_step( aes(x = days, y = live_4, color = colors[4] )) + 
    #         geom_text(aes(label = "4", x = 8, y = live_4[9], color = colors[4] ))
        
    # }

    ggsave(paste0("strain_kms/",strains,"_km_curve.png"), height = 7, width = 7)

}

new_df$auc = as.numeric(as.character(new_df$auc))
new_df$slope = as.numeric(as.character(new_df$slope))

write.csv( new_df, "strains_auc_slope.csv" )

new_df = new_df[order(new_df$"color"),]
new_df$"color" = factor(new_df$"color", levels = c("1","2","3","4"))

xx_la_bigboyo = ggplot( data = new_df, aes(x = `auc`, y = `slope`)) + 
    geom_point( shape = 21, size = 5, aes(fill = `color`)) + 
    scale_fill_manual( values = c("black", "red", "#825CA6", "#E7872B") ) +  
    theme_bw() + 
    labs( x = "Area Under Curve (Worm x Days)", y = "Slope of Linear Best-Fit") + 
    theme(
        # panel.grid.major.x = element_blank(), 
        # panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(), 
        legend.position = "None"
        )

    ggsave("auc_v_slope.png", width = 7, height = 3.31)



