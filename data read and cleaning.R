####dataset conversion

require(tidyverse)
require(sqldf)
require(here)
require(readxl)
require(foreign)
library(readr)
# 
bats = read_delim("/Volumes/maya/Cambridge/Colleen paper/Chiroptera_2025.7.13/occurrence.txt")
rods = read_delim("/Volumes/maya/Cambridge/Colleen paper/Rodentia_2025.7.14/occurrence.txt")

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

#load(here("NewPaperEnviro_wRemarks.RData"))

###there are 1142030 bats and 3444815 rodents

Dataset <- rbind(bats2, rods2)

Dataset %>% group_by(language) %>% summarize(n=n())

Dataset %>% group_by(order) %>% summarize(n=n())

Dataset = Dataset %>% filter(order == "Chiroptera" | order == "Rodentia")

Dataset %>% filter(is.na(year) & is.na(eventDate))

#Dataset = read.csv("/Volumes/maya/Cambridge/Colleen paper/Dataset_9.23.25.csv")

#tissue2 = c("blood", "sangre", "tiss", "tejid", "serum", "suero", "plasma", "biopsy", "biopsia", "swab", "torunda") %>%
  #paste(collapse = "|")

tissue_mj = c("blood", "sangre", "tiss", "tejid", "serum", "suero", "plasma", 
              "biopsy", "biopsia", "swab", "torunda","heart","kidney","liver","lung",
              "muscle","spleen","colon") %>%
  paste(collapse = "|")

#buffer2 = c("ethanol", "etanol", "EtOH", "VTM", "RNAlater","shield", "alcohol", "lysis", "DMSO","EDTA") %>%
  #paste(collapse = "|")

buffer_mj = c("ethanol","etanol","EtOH","VTM","RNAlater","shield","alcohol","lysis","DMSO","EDTA","buffer") %>%
  paste(collapse = "|")

frozen2 = c("froze", "freez", "congelad") %>%
  paste(collapse = "|")

###### one big search variable

#removing "alcohol" on oct 16 2025
terms = c("blood", "sangre", "tiss", "tejid", "serum", "suero", "plasma", 
          "biopsy", "biopsia", "swab", "torunda","heart","kidney","liver","lung",
          "muscle","spleen","colon","ethanol","etanol","EtOH","VTM","RNAlater",
          "shield","lysis","DMSO","EDTA","buffer","froze", "freez", "congelad") %>%
  paste(collapse = "|")

#create frozen and tissue variables
Dataset = Dataset %>%
  mutate(frozen = case_when(
    str_detect(tolower(preparations), frozen2) ~ 1,
    str_detect(tolower(dynamicProperties), frozen2) ~ 1,
    str_detect(tolower(occurrenceRemarks), frozen2) ~ 1,
    str_detect(tolower(materialEntityRemarks), frozen2) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(tissue= case_when(
    str_detect(tolower(preparations), tissue_mj) ~ 1,
    str_detect(tolower(dynamicProperties), tissue_mj) ~ 1,
    str_detect(tolower(occurrenceRemarks), tissue_mj) ~ 1,
    str_detect(tolower(materialEntityRemarks), tissue_mj) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(buffer= case_when(
    str_detect(tolower(preparations), buffer_mj) ~ 1,
    str_detect(tolower(dynamicProperties), buffer_mj) ~ 1,
    str_detect(tolower(occurrenceRemarks), buffer_mj) ~ 1,
    str_detect(tolower(materialEntityRemarks), buffer_mj) ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(tissue_broad = case_when(
    str_detect(tolower(preparations), terms) ~ 1,
    str_detect(tolower(dynamicProperties), terms) ~ 1,
    str_detect(tolower(occurrenceRemarks), terms) ~ 1,
    str_detect(tolower(materialEntityRemarks), terms) ~ 1,
    TRUE ~ 0
  ))

Dataset %>%
  summarize(frozen = sum(as.numeric(frozen), na.rm = T), 
            tissue = sum(as.numeric(tissue), na.rm = T),  
            buffer = sum(as.numeric(buffer), na.rm = T))

##create just one variable

Dataset = Dataset %>%
  mutate(tissue_broad = case_when(
    str_detect(tolower(preparations), terms) ~ 1,
    str_detect(tolower(dynamicProperties), terms) ~ 1,
    str_detect(tolower(occurrenceRemarks), terms) ~ 1,
    str_detect(tolower(materialEntityRemarks), terms) ~ 1,
    TRUE ~ 0
  ))

Dataset %>%
  summarize(tissue = sum(as.numeric(tissue_broad), na.rm = T))

table(Dataset$order, Dataset$tissue_broad)

Dataset %>%
  filter(str_detect(occurrenceRemarks, "congelad") | str_detect(dynamicProperties, "congelad") | str_detect(preparations, "congelad")) 
# %>% view()

extra.countries = data.frame("Country Code" = c("SS", "GG"),
                             "Country" = c("South Sudan", "Guernsey"),
                             "Latitude" = c(6.877, 49.4482),
                             "Longitude" = c(31.307, 2.5895)) %>%
  dplyr::rename("Country Code" = Country.Code)

setwd("~/Documents/PhD/colleen MS")
countries = read_csv("LatLong.csv") %>%
  bind_rows(extra.countries)


Dataset$countryCode[Dataset$countryCode == "BQ"] <- "AN"

Dataset = Dataset %>%
  left_join(countries, by = c("countryCode" = "Country Code"))


Dataset %>%
  filter(is.na(decimalLatitude)) %>%
  group_by(countryCode) %>%
  summarize(n = n()) 


Dataset = Dataset %>%
  mutate(decimalLatitude = case_when(
    !is.na(decimalLatitude) ~ decimalLatitude,
    is.na(decimalLatitude) ~ Latitude
  ), decimalLongitude = case_when(
    !is.na(decimalLongitude) ~ decimalLongitude,
    is.na(decimalLongitude) ~ Longitude
  ))

Dataset = Dataset %>%
  mutate(tissUbuff = case_when(
    tissue == 1 ~ 1,
    buffer == 1 ~ 1,
    TRUE ~ 0
  ))

Dataset = Dataset %>%
  mutate(ethanolPrep = case_when(
    str_detect(tolower(preparations), "ethanol|etanol|EtOH|alcohol") ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(ethanolAll= case_when(
    str_detect(tolower(preparations), "ethanol|etanol|EtOH|alcohol") ~ 1,
    str_detect(tolower(dynamicProperties), "ethanol|etanol|EtOH|alcohol") ~ 1,
    str_detect(tolower(occurrenceRemarks), "ethanol|etanol|EtOH|alcohol") ~ 1,
    TRUE ~ 0
  ))

Dataset %>%
  summarize(ethanolPrep = sum(ethanolPrep, na.rm = TRUE),
            ethanolAll = sum(ethanolAll, na.rm = TRUE)) %>%
  mutate(proportion = ethanolPrep/ethanolAll)

Dataset %>%
  group_by(order, species, taxonKey)%>%
  summarize(n = n())%>%
  write_csv("SpeciesFrequency9.23.25.csv")

Dataset %>%
  group_by(order, species,taxonKey)%>%
  summarize(frozen = sum(frozen, na.rm = T),
            tissue = sum(tissue, na.rm = T),
            buffer = sum(buffer, na.rm = T),
            tissUbuff = sum(tissUbuff, na.rm = T))%>%
  write_csv("SpeciesFTBfrequency9.23.25.csv")

Dataset %>%
  group_by(order, species,taxonKey)%>%
  summarize(tissue = sum(tissue_broad, na.rm = T)) %>%
  write_csv("Speciestissuefrequency9.26.25.csv")

#Dataset2 = Dataset %>%
  #select(-eventDate)

##dates

#Dataset = Dataset %>%
#  mutate(year = case_when(
#    is.na(year) ~ as.numeric(str_extract(verbatimEventDate, "\\d{4}")),
#    TRUE ~ year
#  ))

Dataset$year[which(Dataset$year > 2025)] <- NA
summary(Dataset$year)

##get missing species



#exporting data as .csv

setwd("/Volumes/maya/Cambridge/Colleen paper")
Dataset %>%
  write_csv(("Dataset_9.26.25.csv"))

TissueSamples = Dataset %>% filter(tissue_broad == 1)

setwd("/Volumes/maya/Cambridge/Colleen paper")
TissueSamples %>%
  write_csv(("TissueDataset_10.20.25.csv"))

#exporting data as a .dbf
write.dbf(as.data.frame(Dataset2), here("Dataset_2.21.dbf"))

libr_used()
libr_unused()

