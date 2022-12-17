# This script calculates the importance of each feature (modes of variability and local variable) in proivding
# useful information for predicting drought probability

rm(list = ls())
library(ranger)
library(reshape2)
library(ggplot2)


data_dir <- "E:/Data/drought_paper/data/"


# load the database in a dataframe
db <- read.csv(paste0(data_dir, "ML_Database_All_AWRA_MOf_and_3MPrecip.csv"))


# column names of predictors
predictors <- c("ENSO", "IOD", "SAM", "P_acc_3M", "Qtot",  "Soil_M_root_zone", "Deep_Drainage",
                 "PET_Actual", "E_Actual", "Rainfall" )

# column name of target variable
target <- "Drought...No.Drought"

# keep only the predictors and the target variable
training_testing_all <- db[,c(predictors, target)]

# remove 'no data'
#training_testing_all <- na.omit(training_testing_all)

# Give simpler column names 
predictors_new <- c("ENSO", "IOD", "SAM", "P.Acc.3M", "Q",  "Soil.Moisture", "Deep.Drainage",
                    "PET", "ET", "P")

colnames(training_testing_all) <- c(predictors_new, "Target")

#------ Calculate features importance using Random Forest. Include all years of the drought database

# initialisation
permutation_df <- data.frame()

# Repeat across 50 seed values
for(seed_id in 1:50){
  
  set.seed(seed_id)
  
  # variable importance with RF
  rf.modeli_all <- ranger(Target ~ ., data = training_testing_all, write.forest = TRUE, importance = "permutation", probability=TRUE, oob.error=TRUE)
  
  # extract variable importance
  imp <- rf.modeli_all$variable.importance
  
  # summarise importance for each predictor in a single row 
  row <- data.frame(P.Acc.3M = imp["P.Acc.3M"], Q = imp["Q"], ENSO = imp["ENSO"],
                    IOD=imp["IOD"], SAM = imp["SAM"], Soil.Moisture=imp["Soil.Moisture"], 
                    Deep.Drainage = imp["Deep.Drainage"], PET = imp["PET"],
                    ET = imp["ET"], P = imp["P"])
  
  permutation_df <- rbind(permutation_df, row)
}

# normalise variable importance by dividing each importance by the sum of all variable importance
permutation_df_norm <- permutation_df
permutation_df_norm$total <- rowSums(permutation_df_norm)
for(col_ in predictors_new){
  permutation_df_norm[, col_] <- permutation_df_norm[, col_]/permutation_df_norm$total
}
permutation_df_norm$total <- NULL

# average importance across all seed values
permutation_df_all <- data.frame(Predictor = row.names(data.frame(colMeans(permutation_df_norm))),Score = colMeans(permutation_df_norm))


# generate a plot for normalised variable importance (based on all years) using ggplot
permutation_df_all$Predictor <- factor(permutation_df_all$Predictor, 
                                       levels = rev(c("SAM","ENSO", "IOD",  "P.Acc.3M", "P",  "Soil.Moisture", "Deep.Drainage","Q",
                                                      "ET", "PET" ) ))

varImp_gp <- ggplot(permutation_df_all) +
  aes(x = Predictor, weight = Score) +
  geom_bar(fill = "white", col="black") +
  labs(x = " ", y = "", title = "") +
  coord_flip() +
  theme_minimal()+
  theme(axis.text.x=element_blank())

varImp_gp


#------ Calculate features importance using Random Forest. Include only year 2019

training_testing_only2019 <- db[(db$Year %in% c(2019)), c(predictors, target)]
colnames(training_testing_only2019) <- c(predictors_new, "Target")

permutation_df <- data.frame()

# Repeat across 50 seed values
for(seed_id in 1:50){
  
  set.seed(seed_id)
  
  # variable importance with RF
  rf.modeli_only2019 <- ranger(Target ~ ., data = training_testing_only2019, write.forest = TRUE, importance = "permutation", probability=TRUE, oob.error=TRUE)
  
  imp <- rf.modeli_only2019$variable.importance
  
  # extract importance for each predictor - create a new row 
  row <- data.frame(P.Acc.3M = imp["P.Acc.3M"], Q = imp["Q"], ENSO = imp["ENSO"],
                    IOD=imp["IOD"], SAM = imp["SAM"], Soil.Moisture=imp["Soil.Moisture"], 
                    Deep.Drainage = imp["Deep.Drainage"], PET = imp["PET"],
                    ET = imp["ET"], P = imp["P"])
  
  permutation_df <- rbind(permutation_df, row)
}

# normalise variable importance by dividing each importance by the sum of all variable importance
permutation_df_norm <- permutation_df
permutation_df_norm$total <- rowSums(permutation_df_norm)

for(col_ in predictors_new){
  permutation_df_norm[, col_] <- permutation_df_norm[, col_]/permutation_df_norm$total
}
permutation_df_norm$total <- NULL

# average importance across all seed values
permutation_df_only2019 <- data.frame(Predictor = row.names(data.frame(colMeans(permutation_df_norm))),Score = colMeans(permutation_df_norm))

# generate a plot for normalised variable importance (2019 only) using ggplot
permutation_df_only2019$Predictor <- factor(permutation_df_only2019$Predictor, 
                                            levels = rev(c("SAM","ENSO", "IOD",  "P.Acc.3M", "P",  "Soil.Moisture", "Deep.Drainage","Q",
                                                         "ET", "PET" ) ))
varImp_only_2019_gp <- ggplot(permutation_df_only2019) +
  aes(x = Predictor, weight = Score) +
  geom_bar(fill = "white", col="black") +
  labs(x = " ", y = "", title = "") +
  coord_flip() +
  theme_minimal()+
  theme(axis.text.x=element_blank())

varImp_only_2019_gp

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#--------- How did variable importance change in 2019 compared to all years

changes_df <- permutation_df_only2019
colnames(changes_df) <- c("Predictor", "Score_2019")

# merge importance in all years will change in importance
changes_df <- merge(changes_df, permutation_df_all, by = "Predictor")
changes_df$Score_2019_change <- changes_df$Score_2019 - changes_df$Score 

# convert change into percentage
changes_df[,c("Score_2019","Score", "Score_2019_change")] <- 100* changes_df[,c("Score_2019","Score", "Score_2019_change"  )] 


# In one plot, show variable importance based on all years along with change in importance in 2019

to_plot_df <- changes_df[,c("Predictor", "Score", "Score_2019_change")]
to_plot_df <- melt(to_plot_df, id.vars = "Predictor")
to_plot_df$variable <- as.character(to_plot_df$variable )
to_plot_df[to_plot_df$value < 0, "variable"] <- "Score_neg"
to_plot_df$variable <- factor(to_plot_df$variable, levels = c("Score_2019_change", "Score", "Score_neg") )

varImp_gp <- ggplot(data=to_plot_df, aes(x=Predictor, y=value, fill=variable)) +
  geom_bar(stat="identity")+
  #scale_fill_brewer(palette="Paired")+
  labs(x = " ", y = "(%)", title = "", fill = "") +
  coord_flip() +
  theme_minimal(base_size = 8)+
  theme(legend.position = "top", legend.text =  element_text(size = 4),
        legend.key.size = unit(0.25, 'cm'),
        legend.spacing.x = unit(0.02, 'cm'))+
  scale_fill_manual(labels=c('Imp. increase in 2019', 'Relative Imp.', 'Imp. decrease in 2019'), 
                    values = c("#a6cee3", "#1f78b4", "#b2df8a" ))

varImp_gp
