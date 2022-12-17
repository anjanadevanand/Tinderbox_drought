## this script evaluates the performance of Random Forest drought indicator across 6 metrics of performance

rm(list = ls())

library(ranger) # Random Forest algorithm
library(reshape2)
library(ggplot2)
library(caret) # to create an confusion matrix

data_dir <- "E:/Data/drought_paper/data/"

# calculates the performance of Random Forest across a range of metrics
Get_performance_scores <- function(obs, m.pred){
  
  # create a table for true values and predicted values
  predtable_m_allvar <- table(pred = m.pred, true = obs)
  
  # calculate a confusion matrix
  result_m <- confusionMatrix(predtable_m_allvar, positive="1")
  
  # initialisation
  A <- hits <- predtable_m_allvar[1,1]
  B <- false_alarm <- predtable_m_allvar[1,2]
  C <- misses <- predtable_m_allvar[2,1]
  D <- correct_neg <- predtable_m_allvar[2,2]
  
  accuracy <- (A+D)/(A+B+C+D)
  
  # Metrics of performance 
  
  #kappa <- result_m$overall[2]
  #bias_score <- (A + B)/(A+C)
  true_positive <- A/(A+C) ## Also known as probability of detection rate, Recall,  and sensitivity and hit rate
  false_negative <- C/(A+C)
  false_alarm <- B/(A+B) ## false alarm ratio , (perfect score is 0)
  success_ratio <- A/(A+B) ## also known as precision
  false_positive <- B/(B+D) #Probability of false detection, false alarm rate , (perfect score is 0)
  F1.score <- (true_positive)/(true_positive + 0.5*(false_positive + false_negative))
  #threat_score <- A/(A+B+C)
  #hits_random <- (A+C)*(A+B)/(A+B+C+D)
  #equitable_threat_score <- (A-hits_random)/(A+B+C-hits_random) ## also known as Gillbert Skill score
  #hanssen_Kuipers_disc <- A/(A+C) - B/(B+D) #Hanssen and Kuipers discriminant 
  #true_negative <- D/(D+B) ## Also known as sensitivity
  #balanced_accuracy <- 0.5*(true_positive + true_negative)
  
  # summarise the scores in a single row
  row <- data.frame(Accuracy= accuracy, False.neg.rate = false_negative, False.alarm.rate = false_alarm,
                    F1.score = F1.score, Precision = success_ratio, Recall = true_positive)
  
  return(row)
  
}

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

# Split the data 80% for training and 20% for testing. Repeat for 50 random sampling of training and testing samples 
perf_df <- data.frame()

for(seed in 51:100){
  
  set.seed(seed)
  
  ids_count <- nrow(training_testing_all)
  ## split data into a train and test set
  testids <- sample(ids_count, trunc(ids_count/5))
  
  # samples used for testing
  testseti <- training_testing_all[testids,]
  
  # samples used for training
  trainseti <- training_testing_all[-testids,]
  
  # build a random forest classification model with 'with ranger'
  rf.modeli <- ranger(Target ~ ., data = trainseti, write.forest = TRUE,importance = "impurity",classification = TRUE,   oob.error=TRUE) #classification =TRUE,
  
  # predict drought/no drought for samples in the testing set 
  m.predi <- predict(rf.modeli, testseti[, predictors_new])
  m.pred=m.predi$predictions
  
  # evaluate the performance of Random Forest in the current iteration
  row_df <- Get_performance_scores(obs=testseti[,"Target"], m.pred=m.pred)
  row_df$seed <- seed
  perf_df <- rbind(perf_df, row_df)
  
} 


perf_df <- perf_df[, c("Accuracy", "Precision", "Recall", "F1.score", "False.alarm.rate", "False.neg.rate")]

# reshape the data for plotting with ggplot
perf_df_m <- melt(perf_df)
colnames(perf_df_m) <- c("Metric", "Score")
perf_df_m$Score <- round(perf_df_m$Score,2)
perf_df_m$Metric <- factor(perf_df_m$Metric)


# generate a plot for performance results show 6 metrics
gp_perf <- ggplot(perf_df_m) +
  aes(x = Metric, y = Score) +
  geom_boxplot(shape = "circle") +
  labs(x = " ",  y = " ", title ="") +
  coord_flip() +
  theme_minimal(base_size = 8)+
  scale_y_continuous(breaks = seq(0, 1, by = 0.1))+
  theme(axis.text.y = element_text(size = 7))

gp_perf
