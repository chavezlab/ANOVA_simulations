library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# Set up simulation ------------------------------------------------------------


# Noise level vector
noise <- c(8.00, 4.00, 2.00, 1.00)  

# Sample size vector
ssize <- c(20,40,60,80,100)

# Define standard error function------
se_func <- function (x, na.rm = FALSE) {
  if (na.rm) 
    x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

# Set random seed for reproducibility
set.seed(6)

# Simulation for a within-subject ANOVA  ----------------------------------------

# Create empty data frame
full_df <- data.frame()

#(one within-subject factor "Time", one between-subjects factor "Condition")
for(i in noise){
  for(j in ssize){
    n <- 0
    repeat{
      
      # Update index for each loop
      n <- n + 1
      
      # Create data
      Subject <- as.factor(1:j)
      Con1_Lev1 <- rnorm(j, mean = 0 , sd = i) 
      Con2_Lev1 <- rnorm(j, mean = 0 , sd = i) 
      Con1_Lev2 <- rnorm(j, mean = 1 , sd = i) 
      Con2_Lev2 <- rnorm(j, mean = 0 , sd = i) 
      
      
      # Create data frame and put into long format for aov()
      df <- as.data.frame(cbind(Subject, Con1_Lev1, Con2_Lev1, Con1_Lev2, Con2_Lev2))
      df_melt <- melt(df,id.vars = "Subject", variable.name = "Cell")
      df_melt$Condition <- ifelse(df_melt$Cell == "Con1_Lev1" | df_melt$Cell == "Con1_Lev2", "Con1","Con2")
      df_melt$Level <- ifelse(df_melt$Cell == "Con1_Lev1" | df_melt$Cell == "Con2_Lev1", "Lev1","Lev2")
      
      # Calculate means and SEMs for plotting later
      gg_df <- df_melt %>% 
        select(Condition, Level, value) %>% 
        group_by(Condition,Level) %>% 
        dplyr::summarise(Value = mean(value), se = se_func(value))
      
      # Assign means/SEMs 
      con1lev1_mean <- gg_df[[3]][1]
      con1lev1_se <- gg_df[[4]][1]
      
      con1lev2_mean <- gg_df[[3]][2]
      con1lev2_se <- gg_df[[4]][2]
      
      con2lev1_mean <- gg_df[[3]][3]
      con2lev1_se <- gg_df[[4]][3]
      
      con2lev2_mean <- gg_df[[3]][4]
      con2lev2_se <- gg_df[[4]][4]
      

      # Run ANOVAs and pull F-values and p-values for interaction and main effects.
      F_int <- summary(aov(value ~ Condition*Level + Error(as.factor(Subject)),data = df_melt))[[2]][[1]][[4]][3]
      p_int <- summary(aov(value ~ Condition*Level + Error(as.factor(Subject)),data = df_melt))[[2]][[1]][[5]][3]
      
      F_condition <- summary(aov(value ~ Condition*Level + Error(as.factor(Subject)),data = df_melt))[[2]][[1]][[4]][1]
      p_condition <- summary(aov(value ~ Condition*Level + Error(as.factor(Subject)),data = df_melt))[[2]][[1]][[5]][1]
      
      F_lev <- summary(aov(value ~ Condition*Level + Error(as.factor(Subject)),data = df_melt))[[2]][[1]][[4]][2]
      p_lev <- summary(aov(value ~ Condition*Level + Error(as.factor(Subject)),data = df_melt))[[2]][[1]][[5]][2]
      
    
    # Assign variables to a one row data frame
    sim_df <- data.frame(n, 
                          i,
                          j,
                          F_int,
                          F_condition,
                          F_lev,
                          p_int,
                          p_condition,
                          p_lev,
                          con1lev1_mean,
                          con1lev1_se,
                          con1lev2_mean,
                          con1lev2_se,
                          con2lev1_mean,
                          con2lev1_se,
                          con2lev2_mean,
                          con2lev2_se)  
   
      
      # Concatenate results from each simulation in a dataframe 
      full_df <- rbind(full_df,sim_df)  
      
      # Break after 1000 simulations 
      if(n==1000){break
      }
    }
  }
}

# Post-simulation analyses --------------------------------------------------------------

# Rename variables in data frame
full_df <- full_df %>% dplyr::rename(SNR = i,N = j) 
full_df$SNR <- 1/full_df$SNR   #Transforms simulation standard deviation into SNR value.

# Sort data by interaction significance
full_df_sort <- full_df %>% arrange(SNR, N, p_int)
full_df_sort$sort <- rep(1:1000,20)

# Take 50 most significant interactions
top_df_sort_all <- full_df_sort %>% filter(sort <= 50)


# Plot top twenty most significant interactions
# Make data frame for ggplot means
topfifty_df <- top_df_sort_all %>% 
  select(c1_l1 =con1lev1_mean, 
         c1_l2 = con1lev2_mean, 
         c2_l1 = con2lev1_mean, 
         c2_l2 = con2lev2_mean, 
         sort,
         SNR,
         N) %>%
  melt(id.vars= c("SNR","N","sort"),measure.vars= c("c1_l1","c2_l1","c1_l2","c2_l2"))

# Make data frame for ggplot SEMs
topfifty_df_se <- top_df_sort_all %>% 
  select(c1_l1_se = con1lev1_se, 
         c1_l2_se = con1lev2_se,
         c2_l1_se = con2lev1_se,
         c2_l2_se = con2lev2_se,
         sort,
         SNR,
         N) %>%
  melt(id.vars= c("SNR","N","sort"),measure.vars= c("c1_l1_se","c1_l2_se","c2_l1_se","c2_l2_se"))

# Add SEM vector to the main ggplot data frame 
topfifty_df$se <- topfifty_df_se$value

# Lable conditions for ggplot
topfifty_df$Condition <- ifelse(topfifty_df$variable == "c1_l1" | 
                                   topfifty_df$variable == "c1_l2","Condition 1","Condition 2")
topfifty_df$Level <- ifelse(topfifty_df$variable == "c1_l1" | 
                              topfifty_df$variable == "c2_l1","Level 1","Level 2")


#########################################

# Aggregate top 50 results
topfifty_summary <- topfifty_df %>%
  group_by(Condition,Level,SNR,N) %>%
  summarise(se = mean(se), value = mean(value)) %>%
  filter(SNR == .125 |SNR == .25 | SNR == .5 | SNR == 1)


# Create ggplot--------------------------------------
# Creat main plot
p2 <- ggplot(topfifty_summary, aes(Level,value, group=Condition, color=Condition)) + 
  geom_line() + 
  geom_point() + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1) +
  scale_color_brewer(palette="Set1")  + theme_bw() + 
  facet_grid(N ~ SNR ) + 
  theme(legend.position="bottom", legend.title=element_blank(),panel.grid.minor=element_blank(), panel.grid.major.x = element_blank()) + 
  xlab(NULL) + ylab(NULL) + ggtitle("Aggregated Top 50 Most Significant Interactions \n with True-Effect Present") 

# Add general facet labels
library(gtable)

# get gtable object
z <- ggplot_gtable(ggplot_build(p2))

# add label for right strip
z <- gtable_add_cols(z, z$widths[[8]])
z <- gtable_add_cols(z, unit(1, "line"), 13)
z <- gtable_add_cols(z, unit(1, "line"), 14)
z <- gtable_add_rows(z, unit(.6, "line"), 2)
z <- gtable_add_grob(z, list(rectGrob(gp = gpar(col = NA, fill = gray(0.8))),
                          textGrob("Sample Size", rot = -90, gp = gpar(col = gray(0)))),
                     8, 13, 16, 14, name = paste(runif(2)))

# add label for top strip
z <- gtable_add_rows(z, z$heights[[3]], 2)
z <- gtable_add_grob(z, 
                     list(rectGrob(gp = gpar(col = NA, fill = gray(0.8))),
                          textGrob("SNR", gp = gpar(col = gray(0)))),
                     3, 4, 4, 10, name = paste(runif(2)))

# add margins
z <- gtable_add_cols(z, unit(1/8, "line"), 11)
z <- gtable_add_rows(z, unit(1/8, "line"), 4)

# Draw plot
grid.newpage()
grid.draw(z)



### summary of significant interactions--------------------------------
summary_df <- full_df_sort %>% 
  filter(p_int < .05) %>% arrange(sort)

## Calculate number of disordinal interactions of p < .05
disord_df <- full_df_sort %>% 
  filter(p_int < .05) %>% arrange(sort)


disord_df$disord <-  ifelse(disord_df$con1lev1_mean <= disord_df$con2lev1_mean + disord_df$con2lev1_se & disord_df$con1lev1_mean >= disord_df$con2lev1_mean - disord_df$con2lev1_se | 
                              disord_df$con1lev2_mean <= disord_df$con2lev2_mean + disord_df$con2lev2_se & disord_df$con1lev2_mean >= disord_df$con2lev2_mean - disord_df$con2lev2_se |
                              disord_df$con2lev1_mean <= disord_df$con1lev1_mean + disord_df$con1lev1_se & disord_df$con1lev2_mean >= disord_df$con1lev1_mean - disord_df$con1lev1_se | 
                              disord_df$con2lev2_mean <= disord_df$con1lev2_mean + disord_df$con1lev2_se & disord_df$con1lev1_mean >= disord_df$con1lev2_mean - disord_df$con1lev2_se,"Other","Disordinal")

prop.table(table(disord_df$disord))

disord_df %>% group_by(SNR, N) %>% 
  summarise(percent_disordinal = prop.table(table(disord))[1])

disord_df %>% group_by(SNR, N) %>% 
  summarise(percent_disordinal = prop.table(table(disord))[1]) %>% 
  ggplot(aes(N,percent_disordinal, color=as.factor(SNR))) + 
  geom_point(size =2) + 
  geom_line() + 
  coord_cartesian(ylim = c(0,1)) +
  ggtitle("Percentage of disordinal interactions at p < .05")


## Calculate disordinal interactions of p < .001
disord_df <- full_df_sort %>% 
  filter(p_int < .001) %>% arrange(sort)


disord_df$disord <-  ifelse(disord_df$con1lev1_mean <= disord_df$con2lev1_mean + disord_df$con2lev1_se & disord_df$con1lev1_mean >= disord_df$con2lev1_mean - disord_df$con2lev1_se | 
                              disord_df$con1lev2_mean <= disord_df$con2lev2_mean + disord_df$con2lev2_se & disord_df$con1lev2_mean >= disord_df$con2lev2_mean - disord_df$con2lev2_se |
                              disord_df$con2lev1_mean <= disord_df$con1lev1_mean + disord_df$con1lev1_se & disord_df$con1lev2_mean >= disord_df$con1lev1_mean - disord_df$con1lev1_se | 
                              disord_df$con2lev2_mean <= disord_df$con1lev2_mean + disord_df$con1lev2_se & disord_df$con1lev1_mean >= disord_df$con1lev2_mean - disord_df$con1lev2_se,"Other","Disordinal")

prop.table(table(disord_df$disord))

disord_df %>% group_by(SNR, N) %>% 
  summarise(percent_disordinal = prop.table(table(disord))[1])


disord_df %>% group_by(SNR, N) %>% 
  summarise(percent_disordinal = prop.table(table(disord))[1]) %>% 
  ggplot(aes(N,percent_disordinal, color=as.factor(SNR))) + 
  geom_point(size =2) + 
  geom_line() +
  coord_cartesian(ylim = c(0,1)) +
  ggtitle("Percentage of disordinal interactions at p < .001")


# Model plot ----------------------------------
# A plot of the basic effect we simulating. 

Time <- as.factor(c("Level 1","Level 2","Level 1","Level 2"))
Group <- as.factor(c("Condition 1","Condition 1","Condition 2","Condition 2"))
Value <- c(0,0,0,1)

df <- data.frame(Time, Group, Value)

ggplot(df, aes(Time,Value, group=Group, color=Group)) + 
  geom_line( size=1) + 
  geom_point( size = 2) + 
  scale_color_brewer(palette="Set1") + theme_bw() +
  theme(legend.position="bottom",legend.title=element_blank(),panel.grid.minor=element_blank(), panel.grid.major.x = element_blank()) + 
  xlab("") + ylab("") + coord_cartesian(ylim = c(-.5,1.2)) +
  scale_y_continuous(breaks=seq(-.5, 1.5, .5)) +
  ggtitle("Simulated True-Effect")
