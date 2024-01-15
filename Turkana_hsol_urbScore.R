##Turkana continuous measures of urbanicity
setwd("/Users/bm0211/RegEx/Turkana_CVD")
library(tidyverse)
library(rio)
library(here)
library(reshape2)
library(readr)
turkana <- import("Chosen_few.csv") ##Choose 
dim(turkana)
table(turkana$Sampling_location)
turkana <- turkana[!is.na(turkana$Sampling_location), ]
dim(turkana)
colnames(turkana)
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
turkana <- merge(turkana, locs_data, by = "Sampling_location") # merge by Sampling_location
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
turkana <- turkana %>% distinct(Unique.ID, .keep_all = TRUE)
dim(turkana)
turkana$urb_score <- as.numeric(turkana$urb_score)
##
##
##
colnames(turkana)
dim(turkana)
write.csv(turkana, file = "Standardized_Lifestyle_CONT.csv", row.names = FALSE)
####Sanity check for the location classifications
table(turkana$Sampling_location)
table(turkana$h_sol[turkana$Sampling_location == "Ikoree-Katilu"]) ##Gives you the SCORE of urbanicty
table(turkana$Occupation[turkana$h_sol<2])
