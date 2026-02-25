## GBIF bat and rodent data
## maya juman, colleen cronin, alex richardson
## updated february 25, 2026
## LICENSE: MIT

## clean up space
rm(list=ls()) 
graphics.off()

## set up

require(tidyverse)
require(sqldf)
require(here)
require(readxl)
require(foreign)
library(readr)
library(here)

bats = read_delim(here("Chiroptera_2025.7.13_occurrence.txt"))
rods = read_delim(here("Rodentia_2025.7.14_occurrence.txt"))

##select relevant fields

bats2 = bats %>% select("gbifID", 
                        "accessRights", 
                        "bibliographicCitation", 
                        #"identifier",
                        "language",
                        "license",
                        "publisher",
                        "references",
                        "type",
                        "institutionCode",
                        "basisOfRecord",
                        "dynamicProperties",
                        "occurrenceID",
                        "sex",
                        "preparations",
                        "eventDate",
                        "year",
                        "verbatimEventDate",
                        "continent",
                        "countryCode",
                        "stateProvince",
                        "decimalLatitude",
                        "decimalLongitude",
                        "order",
                        "family",
                        "genus",
                        "genericName",
                        "specificEpithet",
                        "taxonRank",
                        "datasetKey",
                        "taxonKey",
                        "species",
                        "verbatimScientificName",
                        "occurrenceRemarks",
                        "materialEntityRemarks"
)

rods2 = rods %>% select("gbifID", 
                        "accessRights", 
                        "bibliographicCitation", 
                        #"identifier",
                        "language",
                        "license",
                        "publisher",
                        "references",
                        "type",
                        "institutionCode",
                        "basisOfRecord",
                        "dynamicProperties",
                        "occurrenceID",
                        "sex",
                        "preparations",
                        "eventDate",
                        "year",
                        "verbatimEventDate",
                        "continent",
                        "countryCode",
                        "stateProvince",
                        "decimalLatitude",
                        "decimalLongitude",
                        "order",
                        "family",
                        "genus",
                        "genericName",
                        "specificEpithet",
                        "taxonRank",
                        "datasetKey",
                        "taxonKey",
                        "species",
                        "verbatimScientificName",
                        "occurrenceRemarks",
                        "materialEntityRemarks"
)

Dataset <- rbind(bats2, rods2) #create overall dataset

Dataset %>% group_by(language) %>% summarize(n=n()) #languages in overall dataset
Dataset %>% group_by(order) %>% summarize(n=n()) #distribution of bats and rodents

Dataset = Dataset %>% filter(order == "Chiroptera" | order == "Rodentia") #only bats and rodents

##search terms to find tissues
terms = c("blood", "sangre", "tiss", "tejid", "serum", "suero", "plasma", 
          "biopsy", "biopsia", "swab", "torunda","heart","kidney","liver","lung",
          "muscle","spleen","colon","ethanol","etanol","EtOH","VTM","RNAlater",
          "shield","lysis","DMSO","EDTA","buffer","froze", "freez", "congelad") %>%
  paste(collapse = "|")

##filter dataset to find tissue samples
Dataset = Dataset %>%
  mutate(tissue_broad = case_when(
    str_detect(tolower(preparations), terms) ~ 1,
    str_detect(tolower(dynamicProperties), terms) ~ 1,
    str_detect(tolower(occurrenceRemarks), terms) ~ 1,
    str_detect(tolower(materialEntityRemarks), terms) ~ 1,
    TRUE ~ 0
  ))

##how many tissue samples
Dataset %>%
  summarize(tissue = sum(as.numeric(tissue_broad), na.rm = T))

##tissue samples per order
table(Dataset$order, Dataset$tissue_broad)

##associating latitude and longitude with countries
extra.countries = data.frame("Country Code" = c("SS", "GG"),
                             "Country" = c("South Sudan", "Guernsey"),
                             "Latitude" = c(6.877, 49.4482),
                             "Longitude" = c(31.307, 2.5895)) %>%
  dplyr::rename("Country Code" = Country.Code)


##import lat/long file from Wang (2022)
#https://github.com/albertyw/avenews/blob/6f2d89a5f68807348400869fb560e0a2bee94809/old/data/average-latitude-longitude-countries.csv
countries = read_csv(here("LatLong.csv")) %>%
  bind_rows(extra.countries)

#cleaening country codes
Dataset$countryCode[Dataset$countryCode == "BQ"] <- "AN"

#matching
Dataset = Dataset %>%
  left_join(countries, by = c("countryCode" = "Country Code"))

##adding unknown coordinates
Dataset = Dataset %>%
  mutate(decimalLatitude = case_when(
    !is.na(decimalLatitude) ~ decimalLatitude,
    is.na(decimalLatitude) ~ Latitude
  ), decimalLongitude = case_when(
    !is.na(decimalLongitude) ~ decimalLongitude,
    is.na(decimalLongitude) ~ Longitude
  ))

##extracting and cleaning dates
Dataset = Dataset %>%
  mutate(year = case_when(
    is.na(year) ~ as.numeric(str_extract(verbatimEventDate, "\\d{4}")),
    TRUE ~ year
  ))

Dataset$year[which(Dataset$year < 1590)] <- NA
Dataset$year[which(Dataset$year > 2025)] <- NA

##generating date frequencies for figure 2
dat_years <- as.data.frame.table(table(Dataset$year))
dat_years$Var1 <- as.numeric(as.character(dat_years$Var1))

dat_years %>%
  write_csv(("dat_years.csv"))

#exporting tissue sample data (filtered dataset)

TissueSamples = Dataset %>% filter(tissue_broad == 1)
TissueSamples %>%
  write_csv(("TissueDataset_10.20.25.csv"))
