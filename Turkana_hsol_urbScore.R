##Turkana continuous measures of urbanicity
setwd("/Users/bm0211/RegEx/Turkana_CVD")
library(tidyverse)
library(rio)
library(here)
library(reshape2)
library(readr)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(stats)
turkana <- import("/Users/bm0211/RegEx/metab_Meta_median_normalized.csv") ##Choose 
dim(turkana)
colnames(turkana)
merge_duplicated_columns <- function(df) {
  # Find duplicated column names
  dup_cols <- names(df)[duplicated(names(df))]
  
  for (col_name in dup_cols) {
    # Identify all columns that are duplicates of the current one
    cols_to_merge <- which(names(df) == col_name)
    
    # Merge the columns
    merged_column <- df[, cols_to_merge[1]]
    if (length(cols_to_merge) > 1) {
      for (i in cols_to_merge[-1]) {
        # If the value is NA in the merged column, take it from the current column
        # If both are not NA and not equal, it keeps the value from the merged column (you can adjust this logic)
        merged_column <- ifelse(is.na(merged_column), df[, i], merged_column)
      }
    }
    
    # Drop the original duplicated columns
    df <- df[, -cols_to_merge]
    
    # Add the merged column back to the dataframe
    df <- cbind(df, setNames(list(merged_column), col_name))
  }
  
  return(df)
}

# Apply the function to your dataframe before the join
turkana <- merge_duplicated_columns(turkana)
## Household style of life index from Gildner 2020: Market integration and soil-transmitted helminth infection among the Shuar of Amazonian Ecuador
turkana$h_sol <- NA
turkana$h_sol <- ifelse(turkana$Rooms > 1, 1, 0) +
  ifelse(turkana$Floor == "Yes", 1, 0) +
  ifelse(turkana$Roof == "Yes", 1, 0) +
  ifelse(turkana$Electricity == "Yes", 1, 0) +
  ifelse(turkana$Flush_toilet == "Yes", 1, 0) +
  ifelse(turkana$Indoor_Water == "Yes", 1, 0)
table(turkana$h_sol)

## Location-based urbanicity score from Novak, 2012: The development and validation of an urbanicity scale in a multi-country study
locs_data <- import("/Users/bm0211/RegEx/Turkana_stuff/Turkana_sampling_loc_coords.csv") ##I updated it with new locations
##
##
locs_data <- locs_data[, 1:5]
head(locs_data)
table(locs_data$Sampling_location)
locs_data <- locs_data %>% distinct(Sampling_location, .keep_all = TRUE)
##turkana <- merge(turkana, locs_data, by = "Sampling_location") # merge by Sampling_location
dim(turkana)
colnames(turkana)
table(turkana$Floor)
table(turkana$Roof)
table(turkana$Electricity)
table(turkana$Flush_toilet)
table(turkana$Mobile_phone)
table(turkana$television)
table(turkana$Flush_toilet)
table(turkana$Meat_frequency)

# Aggregate scores per location
plumbing <- aggregate(turkana$Flush_toilet ~ turkana$Standardized_name, FUN = function(x) length(which(x == "Yes")) / length(which(x == "Yes" | x == "No")))
electricity <- aggregate(turkana$Electricity ~ turkana$Standardized_name, FUN = function(x) length(which(x == "Yes")) / length(which(x == "Yes" | x == "No")))
tv <- aggregate(turkana$television ~ turkana$Standardized_name, FUN = function(x) length(which(x == "Yes")) / length(which(x == "Yes" | x == "No")))
phone <- aggregate(turkana$Mobile_phone ~ turkana$Standardized_name, FUN = function(x) length(which(x == "Yes")) / length(which(x == "Yes" | x == "No")))
not_wage <- aggregate(turkana$Occupation ~ turkana$Standardized_name, FUN = function(x) length(which(x %in% c("Animal Keeping", "Farmer", "Fisherman", "Gathering", "Herding", "Hunting and Gathering"))) / length(which(x != "NA")))
tot <- as.data.frame(aggregate(turkana$Gender ~ turkana$Standardized_name, FUN = function(x) length(x)))
names(tot) <- c("location", "total")
moms <- subset(turkana, Gender == "Female" & Number_of_children > 0)
moms_ed <- as.data.frame(aggregate(moms$h ~ moms$Standardized_name, FUN = function(x) length(which(x != "none" & x != "None")) / length(which(x != "NA"))))
names(moms_ed) <- c("location", "prop_moms_ed")

urban <- as.data.frame(cbind(plumbing[, 1:2], electricity[, 2], tv[, 2], phone[, 2]))
urban <- full_join(urban, not_wage, by = "turkana$Standardized_name")
names(urban) <- c("location", "prop_toilet", "prop_electricity", "prop_tv", "prop_phone", "prop_not_wage")

##urban2 <- merge(urban, turkana_40plus_ed, by = "location")
##urban3 <- merge(urban2, turkana_40less_ed, by = "location")
urban4 <- merge(urban, tot, by = "location")


# Estimate population density per location
density <- import("/Users/bm0211/RegEx/Turkana_stuff/Kenya_density_2020_2pt5_min.txt")
library(geosphere)

# Pull density from closest coordinates
density <- density[density$x < max(locs_data$X_longitude, na.rm = TRUE) &
                    density$x > min(locs_data$X_longitude, na.rm = TRUE), ]
res <- as.data.frame(matrix(ncol = 4, nrow = dim(locs_data[!is.na(locs_data$Y_latitude), ])[1]))
for (i in 1:dim(locs_data[!is.na(locs_data$Y_latitude), ])[1]){
locs_data <- locs_data[!is.na(locs_data$Y_latitude), ]
y1 <- locs_data$Y_latitude[i]
x1 <- locs_data$X_longitude[i]
tmp <- subset(density, x > (x1 - 0.1) & x < (x1 + 0.1) & y > (y1 - 0.1) & y < (y1 + 0.1))
tmp$dist <- 0
for (k in 1:dim(tmp)[1]){
  tmp$dist[k] <- distm(c(tmp$x[k], tmp$y[k]), c(x1, y1), fun = distHaversine)
}
tmp2 <- tmp[(which(tmp$dist == min(tmp$dist))), ]
res[i, 1:4] <- t(as.numeric(tmp2[1, 1:4]))
}
colnames(res) <- c("X_res", "Y_res", "ken_general_2020", "distance")
#
density2 <- as.data.frame(cbind(locs_data[!is.na(locs_data$Y_latitude), ], res))
urban5 <- merge(density2[, c("Sampling_location", "Standardized_name", "ken_general_2020", "Y_latitude", "X_longitude")], urban4, by.y = "location", by.x = "Standardized_name")

urban5$pop_cat <- 1
urban5$pop_cat[which(urban5$ken_general_2020 > 10)] <- 2
urban5$pop_cat[which(urban5$ken_general_2020 > 50)] <- 3
urban5$pop_cat[which(urban5$ken_general_2020 > 100)] <- 4
urban5$pop_cat[which(urban5$ken_general_2020 > 200)] <- 5
urban5$pop_cat[which(urban5$ken_general_2020 > 300)] <- 6
urban5$pop_cat[which(urban5$ken_general_2020 > 400)] <- 7
urban5$pop_cat[which(urban5$ken_general_2020 > 500)] <- 8
urban5$pop_cat[which(urban5$ken_general_2020 > 1000)] <- 9
urban5$pop_cat[which(urban5$ken_general_2020 > 1500)] <- 10

urban5$urb_score <- NA
urban5$urb_score <- urban5$pop_cat + (10 - (10 * urban5$prop_not_wage)) + (5 * urban5$prop_electricity) + (5 * urban5$prop_toilet) + (5 * urban5$prop_tv) + (5 * urban5$prop_phone)
urban5 <- subset(urban5, urb_score != "NaN")
##view(urban5)
# Merge back to full dataset
turkana <- left_join(turkana, urban5[, c("Sampling_location", "urb_score")], by = "Sampling_location")
dim(turkana)
turkana$urb_score <- as.numeric(turkana$urb_score)
turkana$h_sol <- as.numeric(turkana$h_sol)
##
##
##
colnames(turkana)
dim(turkana)
write.csv(turkana, file = "Metabolomics_WITH_Meta.csv", row.names = FALSE)
####Sanity check for the location classifications
table(turkana$Sampling_location)
table(turkana$h_sol[turkana$Sampling_location == "kitale"]) ##Gives you the SCORE of urbanicty
table(turkana$h_sol)
table(turkana$Occupation[turkana$h_sol<2])

##METABOLOMICS Regression AFTER LIFESTYLE CATEGORY
#an empty data frame to store the results
results <- data.frame(Metabolite = character(),
                      Estimate = numeric(),
                      Std.Error = numeric(),
                      t.value = numeric(),
                      Pr..t.. = numeric(),
                      stringsAsFactors = FALSE)

# Loop through each metabolite
for (metabolite in unique(turkana$Metabolites)) {
  # Filter data for the current metabolite
  metabolite_data <- filter(turkana, Metabolites == metabolite)
  
  # Perform linear regression
  fit <- lm(normLogIC ~ h_sol + Age + Sex, data = metabolite_data)
  
  # Get the summary
  summary_fit <- tidy(fit)
  
  # Extract the row corresponding to the 'h_sol' effect
  h_sol_effect <- filter(summary_fit, term == "h_sol")
  
  if (nrow(h_sol_effect) > 0) {
    results <- rbind(results, c(metabolite, h_sol_effect$estimate, h_sol_effect$std.error,
                                h_sol_effect$statistic, h_sol_effect$p.value))
  }
}

# Rename the columns appropriately
names(results) <- c("Metabolite", "Estimate", "Std.Error", "t.value", "P.value")

# Apply Bonferroni correction
results$Adjusted.P.value <- p.adjust(results$P.value, method = "bonferroni")

# Filter for significant results after correction
significant_metabolites <- results[results$Adjusted.P.value < 0.05, ]

# You might want to sort them again based on Adjusted P.value
significant_metabolites <- arrange(significant_metabolites, Adjusted.P.value)

# Print the significant metabolites after Bonferroni correction
write.csv(significant_metabolites, "significant_metabolites_bonferroni.csv", row.names = FALSE)


#####PCA
# Pivot the data so each metabolite is a column
colnames(turkana)
wide_data <- turkana %>%
  dplyr::select(indiv.ID, Metabolites, normLogIC) %>%
  spread(key = Metabolites, value = normLogIC)

# Replace NA values with a very small number
small_number <- 0.00000000000000000000000000000001
wide_data[is.na(wide_data)] <- small_number
colnames(wide_data)
table(wide_data$indiv.ID)
# Extract h_sol, Age, and Sex for merging back after PCA
meta_data <- turkana %>% distinct(indiv.ID, h_sol, Age, Sex, batch)

# Merge metadata back with PCA input data
pca_input <- merge(wide_data, meta_data, by = "indiv.ID")

# Separate features for PCA
features <- pca_input %>% dplyr::select(-indiv.ID, -h_sol, -Age, -Sex, -batch)

# Apply PCA
pca <- prcomp(features, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca$x)

# Proportion of variance explained by each PC
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
percent_explained <- round(var_explained * 100, 2)

# Scree Plot
plot(percent_explained, xlab = "Principal Component", ylab = "Percentage of Variance Explained", 
     type = 'b', pch = 19, main = "Scree Plot", ylim = c(0, max(percent_explained) + 5))
abline(h = 1, col = "red", lty = 2)  # Optional: Line at 1% for reference

# Adding annotations for the first two PCs
text(x = 1:3, y = percent_explained[1:2], labels = paste(percent_explained[1:2], "%"), pos = 4)
# Add 'batch' back for coloring in the plot
pca_df$batch <- pca_input$batch

# Custom color palette with 8 distinct, easy-on-the-eye colors
my_colors <- c("black", "red", "#4DAF4A", "#984EA3", "cyan", "#FFFF33", "#FF7F00", "#F781BF")

# PCA Plot with percentage of variance explained for 'batch'
Batch_PCA<- ggplot(pca_df, aes(x = PC1, y = PC2, color = as.factor(batch))) +
  geom_point() +
  scale_color_manual(values = my_colors) +
  theme_minimal() +
  labs(title = paste("PCA of Metabolites by batch (PC1:", percent_explained[1], "%, PC2:", percent_explained[2], "%)"),
       x = paste("Principal Component 1 (", percent_explained[1], "%)", sep = ""),
       y = paste("Principal Component 2 (", percent_explained[2], "%)", sep = ""))
ggsave("Batch_PCA_plot.png", plot = Batch_PCA, width = 10, height = 8, dpi = 300)
# Assuming 'h_sol' has 8 categories, applying custom colors for 'h_sol'
pca_df$h_sol <- pca_input$h_sol
Hsol_PCA <- ggplot(pca_df, aes(x = PC1, y = PC2, color = as.factor(h_sol))) +
  geom_point() +
  scale_color_manual(values = my_colors) +
  theme_minimal() +
  labs(title = "PCA of Metabolites by h_sol",
       x = "Principal Component 1",
       y = "Principal Component 2")

ggsave("h_sol_PCA_plot.png", plot = Hsol_PCA, width = 10, height = 8, dpi = 300)
getwd()
###########################################################
##################################################################
###############################################################################
##Regression Controlling for PC1 and PC2
library(tidyverse)

# Assuming 'turkana' and 'Market_data' are already loaded in your environment

##### PCA Code (as provided) #####
wide_data <- turkana %>%
  dplyr::select(indiv.ID, Metabolites, normLogIC) %>%
  spread(key = Metabolites, value = normLogIC)

# Replace NA values with a very small number
small_number <- 0.00000000000000000000000000000001
wide_data[is.na(wide_data)] <- small_number

colnames(turkana)
# Extract h_sol, Age, Sex, and other relevant metadata for merging back after PCA
meta_data <- turkana %>% distinct(indiv.ID, h_sol, Age, Sex, Fetch_water, batch)

# Merge metadata back with PCA input data
pca_input <- merge(wide_data, meta_data, by = "indiv.ID")

# Separate features for PCA
features <- pca_input %>% dplyr::select(-indiv.ID, -h_sol, -Age, -Sex, -Fetch_water, -batch)

# Apply PCA
pca <- prcomp(features, center = TRUE, scale. = TRUE)
pca_scores <- as.data.frame(pca$x)

# Add PC1 and PC2 to the meta_data
meta_data$PC1 <- pca_scores$PC1
meta_data$PC2 <- pca_scores$PC2

# Combine PCA scores with original data (assuming 'Market_data' is your main dataset for regression)
Market_data <- as.data.frame(Market_data)
meta_data <- as.data.frame(meta_data)
# Trim spaces just in case
Market_data$indiv.ID <- trimws(Market_data$indiv.ID)
meta_data$indiv.ID <- trimws(meta_data$indiv.ID)

# Attempt merge again
Market_data <- merge(Market_data, meta_data[, c("indiv.ID", "PC1", "PC2")], by = "indiv.ID", all.x = TRUE)


Market_data <- merge(Market_data, meta_data[, c("indiv.ID", "PC1", "PC2")], by = "indiv.ID", all.x = TRUE)

##### Modified Regression Analysis Including PC1 and PC2 #####
# Create a vector of response variable names (assuming these are your metabolites or related outcomes)
response_vars <- c("HDL", "LDL", "Chol", "BMI", "Body_Fat", "Pulse_pressure", "Trig", "cvd_risk")

# Initialize an empty data frame for results
results_df <- data.frame(ResponseVariable = character(), 
                         Coefficient = numeric(), 
                         StdError = numeric(), 
                         tValue = numeric(), 
                         PValue = numeric(),
                         Predictor = character())

# Loop through each response variable
for (response_var in response_vars) {
  # Adjust the formula to include PC1 and PC2
  formula <- as.formula(paste(response_var, "~ Age + Gender + Diet_MI_items + Fetch_water + Alcohol + PC1 + PC2"))
  
  # Fit the linear regression model
  model <- lm(formula, data = Market_data)
  
  # Get summary of the model
  model_summary <- summary(model)
  
  # Extract coefficients, standard errors, t-values, and p-values
  coefficients <- coef(model_summary)
  std_errors <- coefficients[, "Std. Error"]
  t_values <- coefficients[, "t value"]
  p_values <- coefficients[, "Pr(>|t|)"]
  
  # Prepare results for this variable
  results <- data.frame(ResponseVariable = response_var,
                        Coefficient = coefficients[, "Estimate"],
                        StdError = std_errors,
                        tValue = t_values,
                        PValue = p_values,
                        Predictor = rownames(coefficients))
  
  # Append results
  results_df <- rbind(results_df, results)
}

# Display the results
print(results_df)

