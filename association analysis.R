# Data Processing

# Load data
data <- read.table("data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Inspect data structure
str(data)
summary(data)

# Filter data
data <- subset(data, 
               Func.refGene %in% c("exonic", "splicing") & 
                 ExonicFunc.refGene %in% c("nonsynonymous SNV") &
                 as.numeric(ExAC_ALL) < 0.05)

# Rename columns related to "Otherinfo"
colnames(data)[grepl("Otherinfo", colnames(data))] <- c(
  "Chromosome", "Start", "End", "Ref", "Alt", 
  "Quality", "Filter", "INFO", 
  "GT", "Coverage depth of reference and alternative alleles", "Sequencing depth", "Genotype quality", "Likelihood value of genotype"
)

# Ensure unique column names
data <- data[, !duplicated(colnames(data))]
anyDuplicated(colnames(data))  # Should return 0

# Handle missing values
data <- data %>%
  mutate(across(everything(), ~ifelse(is.na(.), 0, .)))

data$CADD_phred[data$CADD_phred == "."] <- NA
data$CADD_phred <- as.numeric(data$CADD_phred)
data$`GERP.._RS`[data$`GERP.._RS` == "."] <- NA
data$`GERP.._RS` <- as.numeric(data$`GERP.._RS`)
data$ExAC_ALL[data$ExAC_ALL == "."] <- NA
data$ExAC_ALL <- as.numeric(data$ExAC_ALL)

# Gene-based grouping and summarization
gene_summary <- data %>%
  group_by(Gene.refGene) %>%
  summarise(
    SNP_count = n(),
    Avg_CADD = mean(CADD_phred, na.rm = TRUE),
    Avg_GERP = mean(`GERP.._RS`, na.rm = TRUE),
    Avg_AF = mean(ExAC_ALL, na.rm = TRUE)
  ) %>%
  arrange(desc(SNP_count))

# Correlation analysis
cor_test <- cor.test(data$ExAC_ALL, data$CADD_phred, use = "complete.obs")
print(cor_test)

# Visualization of allele frequency and CADD score
library(ggplot2)
ggplot(data, aes(x = ExAC_ALL, y = CADD_phred)) +
  geom_point(color = "black", size = 2, alpha = 0.4) +
  geom_smooth(method = "lm", color = "darkred", fill = "pink", size = 1.2, alpha = 0.3) +
  labs(title = "Allele Frequency vs CADD Score",
       x = "Allele Frequency (ExAC_ALL)",
       y = "CADD Score (phred)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90", linetype = "dashed")
  )

# Consistency analysis for prediction tools
tools <- data %>%
  select(SIFT_pred, Polyphen2_HDIV_pred, Polyphen2_HVAR_pred, MutationTaster_pred) %>%
  mutate(across(everything(), ~ifelse(. %in% c("Damaging", "deleterious"), 1, 0)))

consistency <- rowMeans(tools, na.rm = TRUE)
data$Prediction_Consistency <- consistency

ggplot(data, aes(x = Prediction_Consistency)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "black") +
  labs(title = "Prediction Consistency Distribution",
       x = "Consistency (%)",
       y = "Count") +
  theme_minimal()

# Linear regression model
model <- lm(CADD_phred ~ ExAC_ALL + `GERP.._RS` + phyloP46way_placental, data = data)
data$Predicted <- predict(model, newdata = data)
summary(model)

ggplot(data, aes(x = Predicted, y = CADD_phred)) +
  geom_point(alpha = 0.4, color = "black", size = 2) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Predicted vs Actual",
       x = "Predicted Values",
       y = "Actual Values") +
  theme_minimal()

# Decision tree model
library(rpart)
tree_model <- rpart(CADD_phred ~ ExAC_ALL + `GERP.._RS` + phyloP46way_placental, data = data, method = "anova")
library(rpart.plot)
rpart.plot(tree_model, main = "Decision Tree for CADD Score Prediction")
data$tree_pred <- predict(tree_model)

ggplot(data, aes(x = tree_pred, y = CADD_phred)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "blue") +
  labs(title = "Decision Tree: Predicted vs Observed",
       x = "Predicted CADD Score",
       y = "Observed CADD Score") +
  theme_minimal()

# SVM model
library(e1071)
svm_model <- svm(CADD_phred ~ ExAC_ALL + `GERP.._RS` + phyloP46way_placental, data = data)
data$svm_pred <- predict(svm_model)

ggplot(data, aes(x = svm_pred, y = CADD_phred)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "purple") +
  labs(title = "SVM: Predicted vs Observed",
       x = "Predicted CADD Score",
       y = "Observed CADD Score") +
  theme_minimal()

# Model performance evaluation
evaluate_model <- function(actual, predicted) {
  mse <- mean((actual - predicted)^2)
  r2 <- 1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
  return(list(MSE = mse, R2 = r2))
}

lm_perf <- evaluate_model(data$CADD_phred, data$Predicted)
tree_perf <- evaluate_model(data$CADD_phred, data$tree_pred)
svm_perf <- evaluate_model(data$CADD_phred, data$svm_pred)

performance <- data.frame(
  Model = c("Linear Regression", "Decision Tree", "SVM"),
  MSE = c(lm_perf$MSE, tree_perf$MSE, svm_perf$MSE),
  R2 = c(lm_perf$R2, tree_perf$R2, svm_perf$R2)
)
print(performance)

# Performance visualization
performance_long <- performance %>%
  pivot_longer(cols = c(MSE, R2), names_to = "Metric", values_to = "Value")

ggplot(performance_long, aes(x = Model, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.1), color = "black", width = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
  labs(title = "Model Performance Comparison",
       x = "Model",
       y = "Value") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")
