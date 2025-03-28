## Load in Packages

library(tidyverse)
library(mice)
library(janitor)
library(readr)
library(ggplot2)
library(visdat)


##Load in Datasheets, these include the metadata found in google drive, aoml's qpcr data, and a key to match the two together
mwk_meta <- read_excel("C:/Users/Robert.Bremer/Downloads/Combined_Mastersheet_GJW.xlsx", sheet = "Master DataSheet") %>%
  clean_names() %>%
  rename(sample_id = sample_id_8)

mwk_qpcr <- read_excel("C:/Users/Robert.Bremer/Documents/mwkBB_qPCR_20250325.xlsx") %>%
  drop_na(sample_name) %>%
  #filter(target_name != "DogBact") %>%
  filter(target_name != "DG37") %>%
  group_by(sample_name,target_name, filename) %>%
  mutate(triplicate = row_number()) %>%
  #pivot_wider(names_from = triplicate, values_from = c(ct, quantity)) %>%
  subset(select= -c(well, reporter, quencher, task, automatic_baseline,baseline_start, highsd, baseline_end, 
                    noamp, outlierrg, expfail, automatic_ct_threshold)) %>%
  #fill(c(quantity_mean,ct_1, ct_2, ct_3, quantity_1, quantity_2, quantity_3), .direction = "downup") %>%
  distinct()

mwk_key <- read_excel("C:/Users/Robert.Bremer/Downloads/MWK_EPA_MST_sample_log.xlsx") %>%
  clean_names() %>%
  rename(sample_name = 2) %>%
  mutate(across(sample_name, str_replace,"mwkBB-","")) %>%
  mutate(sample_name = as.numeric(sample_name))

mwk_data <- mwk_meta %>%
  merge(mwk_key) %>%
  left_join(.,mwk_qpcr, by = "sample_name") %>%
  subset(select = -c(sample_id_29, sample_site_30, sample_date_31, sample_date_38, sample_id_37,
                     comments, block, chemistry, runtime, instrument, reference, offscale,
                     status, sample_date_6, sample_site_1, 
                     tidal_conditions_1_high_2_medium_ebb_3_medium_flood_4_low,
                     salinity_ppt_2)) %>%
  mutate(quantity = as.numeric(quantity))

## Symbioses
# Added .3 to all chlorophyll RFU values
# nMDS based on Euclidean similarity matrix of site-averaged and month-averaged data
# PERMANOVA
# 

#What percentage of qpcr values are missing?
mean(is.na(mwk_data$ct))*100

# MNAR Missing Not at Random
# Dependent variable: quantity

###########################################################
mwk_data2 <- mwk_data %>%
  subset(select = c(sample_name,
          quantity,
           sample_site,
           sample_date,
           target_name,
           #region, 
           #water_type,
           #weather_1_sunny_2_cloudy_3_raining, 
           #current_direction_1_still_2_l_to_r_3_r_to_l,
           #air_temperature, 
           water_temp_c, 
           #do_percent, 
           do_mg_l, 
           salinity_ppt, 
           turb_ntu, 
           #chl_rfu, 
           chl_ug_l 
           #mpn, 
           #tn_m_m, 
           #tp_m_m, 
           #n_n_m_m, 
           #no3_m_m, 
           #nh4_m_m, 
           #srp_m_m,
           #doc_m_m
          )) %>%
  mutate(quantity = as.numeric(quantity)) %>%
  mutate(sample_site = as.factor(sample_site)) %>%
  mutate(target_name = as.factor(target_name))

# Look for correlation, even with site_id or region (or other categorical data)


### Begin the imputation

imp2 <- mice(mwk_data2, maxit = 5,
             predictorMatrix = predM,
             method = "cart", print = F)
head(mwk_data2$quantity)

methods(mice) 

mwk_complete <- mice::complete(imp2,1)
mwk_complete

summary(pool(fitimp))
summary(mwk_data2)
write_xlsx(mwk_data,"C:/Users/Robert.Bremer/Documents/mwkBB_data_20250326.xlsx")
write_xlsx()############################################################
# Set the seed for reproducibility
set.seed(12345)
vis_miss(mwk_data2)
typeof(mwk_data2)

# Perform Multiple Imputation
imp_cart <- mice(mwk_data2, m=200, method='cart', print=T)
complete(imp_cart)
plot(imp_cart)
# Adjust the Prediction Matrix to remove singularity
pred<-imp_cart$pred
pred
pred2 <- pred
pred2[c("sample_name","sample_site","sample_date","target_name","water_temp_c",
       "do_mg_l","salinity_ppt","turb_ntu","chl_ug_l"),"quantity"] <- 0
pred2
# Rerunning the mice algorithm
imp_cart2 <- mice(mwk_data2,meth = "cart", pred = pred2, m = 200)
imp_cart2$meth
plot(imp_cart2)
# Checking convergence
miss2 <- is.na(imp_cart2$data$quantity)
xyplot(imp_cart2, quantity ~ I, na.groups = miss2)


miss <- is.na(imp_cart$data$quantity)
xyplot(imp_cart, quantity ~ I (wgt / (hgt / 100)^2),
       na.groups = miss, cex = c(0.8, 1.2), pch = c(1, 20),
       ylab = "Quantity Imputed", xlab = "Quantity Calculated")
plot(imp_cart, c("quantity"))
pred2<-imp_cart$pred
pred2[c("target_name","water_type"), "quantity"] <- 1
pred2
imp_cart2 <-mice(mwk_data2, meth="cart", pred=pred2, m = 50, print=FALSE)
plot(imp_cart2, c("quantity"))

# Making a 3rd, very small data table
mwk_data3 <- mwk_data2 %>%
  subset(select= c(sample_name,
                   quantity,
                   sample_site,
                   sample_date,
                   target_name))

# Smaller imputation
imp3 <- mice(mwk_data3, maxit = 500)
imp3$meth
plot(imp3)

#Plot some graphs for environmental data


ggplot(mwk_data, aes(x = sample_date, y = water_temp_c))+
  geom_point()

ggplot(mwk_data, aes(x = sample_date, y = salinity_ppt, color = water_type))+
  geom_point()

#Plot some graphs for MST numbers
ggplot(filter(mwk_data, target_name =="DG3"), aes(x = sample_date, y = quantity, color = site_code, shape = region))+
  geom_point()

ggplot(filter(mwk_data, target_name =="HF183"), aes(x = sample_date, y = quantity, color = site_code))+
  geom_point()

ggplot(filter(mwk_data, target_name =="DogBact"), aes(x = sample_date, y = quantity, color = site_code))+
  geom_point()

# Plot some isotope
