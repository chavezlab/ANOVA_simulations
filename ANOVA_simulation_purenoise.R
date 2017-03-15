library(dplyr)
library(ggplot2)
library(robR)
library(reshape2)
library(RColorBrewer)
#####################

# Create empty data frame
full_df <- data.frame()

# Create simulation number index
n <- 0

set.seed(5)
# Simulation for a mixed effects ANOVA (one within-subject factor "Time", one between-subjects factor "Condition")
repeat{
  
  # Update index for each loop
  n <- n + 1
  
  # Create random data from uniform distribution. 
  Subject <- as.factor(1:40)
  Condition <- as.factor(c(rep("Con1",40),rep("Con2",40)))
  Time1 <- runif(80,0,10)
  Time2 <- runif(80,0,10)
  
  # Create data frame and put into long format for aov()
  df <- as.data.frame(cbind(Subject,Condition,Time1,Time2))
  df_melt <- melt(df,id.vars = c("Subject","Condition"),variable.name = "Time")
  
  # Calculate means and SEMs for plotting later
  gg_df <- df_melt %>% 
    select(Condition, Time, value) %>% 
    group_by(Condition,Time) %>% 
    dplyr::summarise(Value = mean(value), se = se.rob(value))
  
  # Assign means/SEMs 
  con1time1_mean <- gg_df[[3]][1]
  con1time1_se <- gg_df[[4]][1]
  
  con1time2_mean <- gg_df[[3]][2]
  con1time2_se <- gg_df[[4]][2]
  
  con2time1_mean <- gg_df[[3]][3]
  con2time1_se <- gg_df[[4]][3]
  
  con2time2_mean <- gg_df[[3]][4]
  con2time2_se <- gg_df[[4]][4]
  
  # Run ANOVAs and pull F-values and p-values for interaction and main effects.
  F_int <- summary(aov(value ~ Condition*Time + Error(as.factor(Subject)),data = df_melt))[[2]][[1]][[4]][3]
  p_int <- summary(aov(value ~ Condition*Time + Error(as.factor(Subject)),data = df_melt))[[2]][[1]][[5]][3]
  
  F_condition <- summary(aov(value ~ Condition*Time + Error(as.factor(Subject)),data = df_melt))[[2]][[1]][[4]][1]
  p_condition <- summary(aov(value ~ Condition*Time + Error(as.factor(Subject)),data = df_melt))[[2]][[1]][[5]][1]
  
  F_time <- summary(aov(value ~ Condition*Time + Error(as.factor(Subject)),data = df_melt))[[2]][[1]][[4]][2]
  p_time <- summary(aov(value ~ Condition*Time + Error(as.factor(Subject)),data = df_melt))[[2]][[1]][[5]][2]
  
  # Assign variables to a one row data frame
  sim_df <- sim <- data.frame(n, F_int,
                              F_condition,
                              F_time,
                              p_int,
                              p_condition,
                              p_time,
                              con1time1_mean,
                              con1time1_se,
                              con1time2_mean,
                              con1time2_se,
                              con2time1_mean,
                              con2time1_se,
                              con2time2_mean,
                              con2time2_se)  
  
  # Concatenate results from each simulation in a dataframe 
  full_df <- rbind(full_df,sim_df)  
  
  # Run 10000 simulations 
  if(n==10000){break
  }
}

# Write data file for graphing later
#write.csv(full_df,"C:/Users/chavez.95/Google Drive/R/Studies/ANOVA_sims/full_df_aov_noise_withinsub.csv")


# Plotting --------------------------------------
full_df <- read.csv("C:/Users/chavez.95/Google Drive/R/Studies/ANOVA_sims/full_df_aov_noise_withinsub.csv")

# Sort simulation data frame and index it
full_df <- full_df[order(full_df$p_int), ]
full_df$sort <- 1:10000

# Plot top fifty most significant interactions
topfifty_df <- full_df %>% 
  filter(sort <= 50) %>% 
  select(c1t1 =con1time1_mean, 
         c1t2 = con1time2_mean, 
         c2t1 = con2time1_mean, 
         c2t2 = con2time2_mean,
         sort) %>%
  melt(id.vars= "sort",measure.vars= c("c1t1","c2t1","c1t2","c2t2"))

topfifty_df_se <- full_df %>% 
  filter(sort <= 50) %>% 
  select(c1t1_se = con1time1_se, 
         c1t2_se = con1time2_se,
         c2t1_se = con2time1_se,
         c2t2_se = con2time2_se,
         sort) %>%
  melt(id.vars= "sort",measure.vars= c("c1t1_se","c1t2_se","c2t1_se","c2t2_se"))

topfifty_df$se <- topfifty_df_se$value

topfifty_df_pval <- full_df %>% 
  filter(sort <= 50) %>% 
  select(p_value = p_int,
         sort) %>%
  melt(id.vars= "sort",measure.vars= c("p_value"))

topfifty_df$p_value <- topfifty_df_pval$value

# Define aesthetics 
topfifty_df$Group <- ifelse(topfifty_df$variable == "c1t1" | topfifty_df$variable == "c1t2","Condition 1","Condition 2")
topfifty_df$Time <- ifelse(topfifty_df$variable == "c1t1" | topfifty_df$variable == "c2t1","Level 1","Level 2")
topfifty_df$Group_col <- ifelse(topfifty_df$variable == "c1t1" | topfifty_df$variable == "c1t2","black","gray50")
topfifty_df$p_col <- ifelse(topfifty_df$p_value < .001, "p < .001","p < .005")

# Plot
ggplot(topfifty_df, aes(Time,value, group=Group, color=p_col, linetype=Group)) + 
  geom_line() + 
  geom_point() + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, size=.5) +
  scale_color_brewer(palette="Dark2") +
  facet_wrap( ~ sort , ncol = 10) + theme_bw() +
  theme(legend.position="bottom", legend.title=element_blank(),panel.grid.minor=element_blank(), panel.grid.major.x = element_blank()) + 
  scale_y_continuous(breaks=seq(2, 8, 1)) +
  xlab(NULL) + ylab(NULL) + ggtitle("Top 50 Most Significant Interactions from Pure Noise")



### summary of significant interactions --------------------
topfifty_df <- full_df %>% 
  filter(p_int < .05) %>% arrange(sort)

attach(topfifty_df)
topfifty_df$disord <-  ifelse(con1time1_mean <= con2time1_mean + con2time1_se & con1time1_mean >= con2time1_mean - con2time1_se | 
                            con1time2_mean <= con2time2_mean + con2time2_se & con1time2_mean >= con2time2_mean - con2time2_se |
                            con2time1_mean <= con1time1_mean + con1time1_se & con1time2_mean >= con1time1_mean - con1time1_se | 
                            con2time2_mean <= con1time2_mean + con1time2_se & con1time1_mean >= con1time2_mean - con1time2_se,"Other","Disordinal")

prop.table(table(topfifty_df$disord))


# Model plot ----------------------------------
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

