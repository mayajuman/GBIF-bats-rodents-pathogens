####tissue analyses

library(tidyr)
library(maptools)
library(car)
library(performance)
library(visreg)
library(emmeans)
library(ggpubr)
library(letsR)
library(raster)
library(dplyr)
library(stringr)
library(plyr)
library(taxize)

###read in

TissueSamples <- read.csv("/Volumes/maya/Cambridge/Colleen paper/TissueDataset_10.20.25.csv")

TissueSamples <- TissueSamples[-which(TissueSamples$order == ""),]

TissueSamples$preparations[108258] <- "Piel, cráneo,  esqueleto y tejido"

for (i in 15:36) {
  TissueSamples[108258,i] <- TissueSamples[108258,i+2]
}

TissueSamples[108258,37:38] <- NA

#Dataset <- read.csv("/Volumes/maya/Cambridge/Colleen paper/Dataset_9.26.25.csv")

####spatial

#for colleen
#percountry <- TissueSamples %>% 
  #group_by(countryCode, order) %>% 
  #dplyr::summarise(n = n()) %>%
  #pivot_wider(names_from = order, values_from = n)

percountry <- TissueSamples %>% 
  group_by(countryCode) %>% 
  dplyr::summarise(n = n())

percountry$countryCode[is.na(percountry$countryCode)] <- "NA"

#load georegion
setwd("~/Documents/PhD/colleen MS")
geo=read.csv("georegion.csv")
geo$countryCode=geo$alpha.2
geo=geo[c("name","region","sub.region","countryCode")]
setdiff(percountry$countryCode,geo$countryCode)

#add missing "country"
geo[250,1] <- "Netherlands Antilles"
geo[250,2] <- "Americas"
geo[250,3] <- "Latin America and the Caribbean"
geo[250,4] <- "AN"

##namibia is not NA
geo$countryCode[which(geo$name == "Namibia")] <- "NA"
geo$sub.region[which(geo$sub.region == "Australia and New Zealand")] <- "Aus, NZ"
geo$sub.region[which(geo$sub.region == "Central Asia")] <- "C. Asia"
geo$sub.region[which(geo$sub.region == "Eastern Asia")] <- "E. Asia"
geo$sub.region[which(geo$sub.region == "Eastern Europe")] <- "E. Europe"
geo$sub.region[which(geo$sub.region == "Latin America and the Caribbean")] <- "L. America, Carib."
geo$sub.region[which(geo$sub.region == "Northern Africa")] <- "N. Africa"
geo$sub.region[which(geo$sub.region == "Northern America")] <- "N. America"
geo$sub.region[which(geo$sub.region == "Northern Europe")] <- "N. Europe"
geo$sub.region[which(geo$sub.region == "South-eastern Asia")] <- "S.E. Asia"
geo$sub.region[which(geo$sub.region == "Southern Asia")] <- "S. Asia"
geo$sub.region[which(geo$sub.region == "Southern Europe")] <- "S. Europe"
geo$sub.region[which(geo$sub.region == "Sub-Saharan Africa")] <- "S.S. Africa"
geo$sub.region[which(geo$sub.region == "Western Asia")] <- "W. Asia"
geo$sub.region[which(geo$sub.region == "Western Europe")] <- "W. Europe"

#remove blanks
percountry <- percountry[-which(percountry$countryCode == "ZZ"),]

#percountry$Total <- rowSums(percountry[,2:3],na.rm=TRUE)
#write.csv(percountry, "Updated country data Oct 16.csv", row.names = FALSE)

setdiff(percountry$countryCode,geo$countryCode)

#merge
gdata=merge(geo,percountry,by="countryCode",all.x=T)
gdata$binsample=ifelse(is.na(gdata$n),0,1)

#remove antarctica
gdata=gdata[!gdata$region=="",]

## GLM for binary sampling
mod1=glm(binsample~region,data=gdata,family=binomial)

## within sampled GLMs
gdata2=gdata[gdata$binsample==1,]
mod2=glm(n~sub.region,data=gdata2,family=poisson)

#range from 1 to 390176
range(gdata2$n)

## Anova
Anova(mod1) #sig for binary effort
Anova(mod2) #vsig for # of samples

## R2
r2_mcfadden(mod1) #0.05
r2_mcfadden(mod2) #0.62

## visreg
visreg(mod1,"region",scale="response",gg=FALSE, xlab="Region", rug=FALSE, partial=FALSE, ylab="Binary sampling of country")
visreg(mod2,"sub.region",scale="response",gg=FALSE, xlab="Subregion",  rug=FALSE, partial=FALSE, ylab="Number of samples collected", cex.axis=0.5)

##spatial biases by institution

setwd("~/Documents/PhD/colleen MS")
mus=read.csv("institutions.csv")
mus <- mus %>% dplyr::select(Institution, Rodentia, Chiroptera, Total, institutionCountry)
mus$Institution <- as.factor(mus$Institution)
mus <- mus %>% dplyr::group_by(institutionCountry) %>% 
  dplyr::summarise(Rodentia = sum(Rodentia), Chiroptera = sum(Chiroptera), Total=sum(Total)) 
mus <- mus[complete.cases(mus),]

inst=mus
inst$countryCode = inst$institutionCountry
inst$n = inst$Total

#geo=read.csv("georegion.csv")
#geo$countryCode=geo$alpha.2
#geo=geo[c("name","region","sub.region","countryCode")]
#setdiff(inst$countryCode,geo$countryCode)


#merge
mdata=merge(geo,inst,by="countryCode",all.x=T)
mdata$binsample=ifelse(is.na(mdata$n),0,1)

#remove antarctica
mdata=mdata[!mdata$region=="",]

## GLM for binary sampling
mmod1=glm(binsample~region,data=mdata,family=binomial)

## within sampled GLMs
mdata2=mdata[mdata$binsample==1,]
mmod2=glm(n~sub.region,data=mdata2,family=poisson)

#range from 1 to 491745; USA
range(mdata2$n)

## Anova
Anova(mmod1) #sig for binary effort
Anova(mmod2) #vsig for # of samples

## R2
r2_mcfadden(mmod1) 
r2_mcfadden(mmod2)


## visreg
visreg(mmod1,"region",scale="response",gg=FALSE, xlab="Region", rug=FALSE, partial=FALSE, ylab="Likelihood of a country having tissue holdings")
visreg(mmod2,"sub.region",scale="response",gg=FALSE, xlab="Subregion", rug=FALSE, partial=FALSE, ylab="Number of samples in holdings", cex.axis=0.8)

##appendix A

mus$`Institution Country`=mus$institutionCountry
mus$`Rodentia Tissue Samples`=mus$Rodentia
mus$`Chiroptera Tissue Samples`=mus$Chiroptera
mus$`Total Tissue Samples`=mus$Total
mus$`Institution Country`=mus$institutionCountry
mus$Total=NULL
mus$Rodentia=NULL
mus$Chiroptera=NULL
mus$institutionCountry=NULL

#write.csv(mus, "Appendix A.csv", row.names = FALSE)

#temporal

##get missing dates from verbatimEventDate

TissueSamples$year <- as.numeric(TissueSamples$year)

TissueSamples = TissueSamples %>%
  mutate(year = case_when(
    is.na(year) ~ as.numeric(str_extract(verbatimEventDate, "\\d{4}")),
    TRUE ~ year
  ))

table(TissueSamples$year, useNA="ifany")

TissueSamples$year[which(TissueSamples$year < 1800)] <- NA
TissueSamples$year[which(TissueSamples$year > 2025)] <- NA
TissueSamples$year[which(TissueSamples$verbatimEventDate == "1864 to 2015")] <- NA

summary(TissueSamples$year)

tis_years <- as.data.frame.table(table(TissueSamples$year))
tis_years$Var1 <- as.numeric(as.character(tis_years$Var1))

##dates

Dataset = Dataset %>%
  mutate(year = case_when(
    is.na(year) ~ as.numeric(str_extract(verbatimEventDate, "\\d{4}")),
    TRUE ~ year
  ))

table(Dataset$year, useNA="ifany")
summary(Dataset$year)

Dataset$year[which(Dataset$year < 1590)] <- NA
Dataset$year[which(Dataset$year > 2025)] <- NA

dat_years <- as.data.frame.table(table(Dataset$year))
dat_years$Var1 <- as.numeric(as.character(dat_years$Var1))

Dataset$tissue_broad <- as.factor(Dataset$tissue_broad)

ggplot(Dataset %>% filter(year > 1806), aes(x = year, color = tissue_broad)) + 
  geom_histogram(fill = "white") +
  xlab("Year") +
  ylab("Count") +
  scale_fill_manual(values = c("lightblue", "red")) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::comma)

###trying gams

##functions from Gibb et al. 2022 Biology Letters (https://doi.org/10.1098/rsbl.2021.0427)

Deriv <- function(mod, n = 200, eps = 1e-7, newdata) {
  if(isTRUE(all.equal(class(mod), "list")))
    mod <- mod$gam
  m.terms <- attr(terms(mod), "term.labels")
  if(missing(newdata)) {
    newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                   function(x) seq(min(x), max(x), length = n))
    names(newD) <- m.terms
  } else {
    newD <- newdata
  }
  X0 <- predict(mod, data.frame(newD), type = "lpmatrix")
  newD <- newD + eps
  X1 <- predict(mod, data.frame(newD), type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xp.r <- NROW(Xp)
  Xp.c <- NCOL(Xp)
  ## dims of bs
  bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
  # number of smooth terms
  t.labs <- attr(mod$terms, "term.labels")
  nt <- length(t.labs)
  ## list to hold the derivatives
  lD <- vector(mode = "list", length = nt)
  names(lD) <- t.labs
  for(i in seq_len(nt)) {
    Xi <- Xp * 0
    want <- grep(t.labs[i], colnames(X1))
    Xi[, want] <- Xp[, want]
    df <- Xi %*% coef(mod)
    df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
    lD[[i]] <- list(deriv = df, se.deriv = df.sd)
    ## Xi <- Xp * 0 ##matrix(0, nrow = Xp.r, ncol = Xp.c)
    ## J <- bs.dims[i]
    ## Xi[,(i-1) * J + 1:J + 1] <- Xp[,(i-1) * J + 1:J +1]
    ## df <- Xi %*% coef(mod)
    ## df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
    ## lD[[i]] <- list(deriv = df, se.deriv = df.sd)
  }
  class(lD) <- "Deriv"
  lD$gamModel <- mod
  lD$eps <- eps
  lD$eval <- newD - eps
  return(lD)
}

confint.Deriv <- function(object, term, alpha = 0.05, ...) {
  l <- length(object) - 3
  term.labs <- names(object[seq_len(l)])
  if(missing(term))
    term <- term.labs
  Term <- match(term, term.labs)
  ##term <- term[match(term, term.labs)]
  if(any(miss <- is.na(Term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  ## if(is.na(term))
  ##     stop("'term' not a valid model term.")
  res <- vector(mode = "list", length = length(term))
  names(res) <- term
  residual.df <- length(object$gamModel$y) - sum(object$gamModel$edf)
  tVal <- qt(1 - (alpha/2), residual.df)
  ## tVal <- qt(1 - (alpha/2), object$gamModel$df.residual)
  for(i in seq_along(term)) {
    upr <- object[[term[i]]]$deriv + tVal * object[[term[i]]]$se.deriv
    lwr <- object[[term[i]]]$deriv - tVal * object[[term[i]]]$se.deriv
    res[[term[i]]] <- list(upper = drop(upr), lower = drop(lwr))
  }
  res$alpha = alpha
  res
}

signifD <- function(x, d, upper, lower, eval = 0) {
  miss <- upper > eval & lower < eval
  incr <- decr <- x
  want <- d > eval
  incr[!want | miss] <- NA
  want <- d < eval
  decr[!want | miss] <- NA
  list(incr = incr, decr = decr)
}

plot.Deriv <- function(x, alpha = 0.05, polygon = TRUE,
                       sizer = FALSE, term, eval = 0, lwd = 3,
                       col = "lightgrey", border = col,
                       ylab, xlab, ...) {
  l <- length(x) - 3
  ## get terms and check specified (if any) are in model
  term.labs <- names(x[seq_len(l)])
  if(missing(term))
    term <- term.labs
  Term <- match(term, term.labs)
  if(any(miss <- is.na(Term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  if(all(is.na(Term)))
    stop("All terms in 'term' not found in model.")
  l <- sum(!miss)
  nplt <- n2mfrow(l)
  ## tVal <- qt(1 - (alpha/2), x$gamModel$df.residual)
  residual.df <- length(x$gamModel$y) - sum(x$gamModel$edf)
  tVal <- qt(1 - (alpha/2), residual.df)
  if(missing(ylab))
    ylab <- expression(italic(hat(f)*"'"*(x)))
  if(missing(xlab)) {
    xlab <- attr(terms(x$gamModel), "term.labels")[Term]
    names(xlab) <- xlab
  }
  layout(matrix(seq_len(l), nrow = nplt[1], ncol = nplt[2]))
  CI <- confint(x, term = term, alpha = alpha)
  for(i in seq_along(term)) {
    ## for(i in seq_len(l)) {
    upr <- CI[[term[i]]]$upper
    lwr <- CI[[term[i]]]$lower
    ylim <- range(upr, lwr)
    plot(x$eval[,term[i]], x[[term[i]]]$deriv, type = "n",
         ylim = ylim, ylab = ylab, xlab = xlab[term[i]], ...)
    if(isTRUE(polygon)) {
      polygon(c(x$eval[,term[i]], rev(x$eval[,term[i]])),
              c(upr, rev(lwr)), col = col, border = border)
    } else {
      lines(x$eval[,term[i]], upr, lty = "dashed")
      lines(x$eval[,term[i]], lwr, lty = "dashed")
    }
    abline(h = 0, ...)
    if(isTRUE(sizer)) {
      lines(x$eval[,term[i]], x[[term[i]]]$deriv, lwd = 1)
      S <- signifD(x[[term[i]]]$deriv, x[[term[i]]]$deriv, upr, lwr,
                   eval = eval)
      lines(x$eval[,term[i]], S$incr, lwd = lwd, col = "blue")
      lines(x$eval[,term[i]], S$decr, lwd = lwd, col = "red")
    } else {
      lines(x$eval[,term[i]], x[[term[i]]]$deriv, lwd = 2)
    }
  }
  layout(1)
  invisible(x)
}


m1 <- mgcv::gam(Freq ~ s(Var1), data=tis_years %>% filter(Var1 < 2025), family="nb", method="REML")
plot(m1)

plot(mgcv::gamm(Freq ~ s(Var1), data=tis_years %>% filter(Var1 < 2025), family="nb", method="REML")$gam)

m2 <- (mgcv::gamm(Freq ~ s(Var1), data=dat_years, method="REML")$gam)
plot(m2)

preds = data.frame(Var1 = seq(min(tis_years$Var1), 2024, length.out=100))
preds = cbind(preds, predict(m1, preds, type = "response", se.fit = TRUE))
preds$upper = (preds$fit + (1.96*preds$se.fit))
preds$lower = (preds$fit - (1.96*preds$se.fit))
preds$fitted = (preds$fit)

m1.d = Deriv(m1, n=length(preds$Var1))
plot(m1.d, sizer=TRUE)

CI = confint(m1.d, alpha = 0.01)
S = signifD(preds$fit, m1.d$Var1$deriv, CI$Var1$upper, CI$Var1$lower, eval = 0)
preds$sig_incr = !is.na(S$incr)
preds$sig_decr = !is.na(S$decr)
preds$sig = !is.na(S$incr) | !is.na(S$decr)

preds$signif_col = NA
preds$signif_col[ preds$sig_incr == TRUE ] = "Increase"
preds$signif_col[ preds$sig_incr == FALSE ] = "Decrease"

preds_inc <- preds[preds$sig_incr == TRUE,]
preds_inc1 <- preds_inc[which(preds_inc$Var1 < 1880),]
preds_inc2 <- preds_inc[which(preds_inc$Var1 > 1880 & preds_inc$Var1 < 1940),]
preds_inc3 <- preds_inc[which(preds_inc$Var1 > 1940),]

preds_dec <- preds[preds$sig_decr == TRUE,]
preds_dec1 <- preds_dec[which(preds_dec$Var1 < 1875),]
preds_dec2 <- preds_dec[which(preds_dec$Var1 > 1900 & preds_dec$Var1 < 1950),]
preds_dec3 <- preds_dec[which(preds_dec$Var1 > 2000),]

fig1 <- ggplot() + 
  geom_bar(mapping = aes(x = Var1, y = Freq), stat = 'identity', color = "grey90", 
           data = dat_years %>% filter(Var1 > 1850)) +
  geom_bar(mapping = aes(x = Var1, y = Freq), stat = 'identity', color = "grey40",
           data = tis_years %>% filter(Var1 > 1850)) +
  xlab("Year") +
  ylab("Count") +
  theme_minimal() +
  geom_ribbon(data=preds, aes(x=Var1, ymin=lower, ymax=upper), alpha=0.2, fill="skyblue4", col=NA, size=0.05) +
  geom_line(data=preds, aes(x=Var1, y=fitted), col="black", size=1.2) +
  geom_line(data=preds_inc1, aes(x=Var1, y=fitted), col="deepskyblue", size=1.2) +
  geom_line(data=preds_inc2, aes(x=Var1, y=fitted), col="deepskyblue", size=1.2) +
  geom_line(data=preds_inc3, aes(x=Var1, y=fitted), col="deepskyblue", size=1.2) +
  geom_line(data=preds_dec1, aes(x=Var1, y=fitted), col="brown", size=1.2) +
  geom_line(data=preds_dec2, aes(x=Var1, y=fitted), col="brown", size=1.2) +
  geom_line(data=preds_dec3, aes(x=Var1, y=fitted), col="brown", size=1.2) + scale_y_sqrt()

ggsave("Fig 1.jpeg", fig1, height=4, width=4)

ggplot() + 
  geom_bar(mapping = aes(x = Var1, y = Freq), stat = 'identity', color = "grey90", 
           data = dat_years %>% filter(Var1 > 1850)) +
  geom_bar(mapping = aes(x = Var1, y = Freq), stat = 'identity', color = "grey40",
           data = tis_years %>% filter(Var1 > 1850)) +
  xlab("Year") +
  ylab("Count") +
  theme_minimal() +
  geom_ribbon(data=preds, aes(x=Var1, ymin=lower, ymax=upper), alpha=0.2, fill="skyblue4", col=NA, size=0.05) +
  geom_line(data=preds, aes(x=Var1, y=fitted), col="black", size=1.2)

##preparations

bloodterms <- c("blood", "sangre", "serum", "suero", "plasma") %>%
  paste(collapse = "|")

tissueterms <- c("tiss", "tejid") %>%
  paste(collapse = "|")

biopsyterms <- c("biopsy", "biopsia") %>%
  paste(collapse = "|")

swabterms <- c("swab", "torunda") %>%
  paste(collapse = "|")

etohterms <- c("ethanol", "etanol", "EtOH") %>%
  paste(collapse = "|")

freezeterms <- c("froze", "freez", "congelad") %>%
  paste(collapse = "|")

bufferterms <- c("DMSO", "EDTA", "VTM", "RNAlater", "shield", "buffer", "lysis") %>%
  paste(collapse = "|")

TissueSamples = TissueSamples %>%
  mutate(bloodterms = case_when(
    str_detect(tolower(preparations), bloodterms) ~ 1,
    str_detect(tolower(dynamicProperties), bloodterms) ~ 1,
    str_detect(tolower(occurrenceRemarks), bloodterms) ~ 1,
    str_detect(tolower(materialEntityRemarks), bloodterms) ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(tissueterms = case_when(
    str_detect(tolower(preparations), tissueterms) ~ 1,
    str_detect(tolower(dynamicProperties), tissueterms) ~ 1,
    str_detect(tolower(occurrenceRemarks), tissueterms) ~ 1,
    str_detect(tolower(materialEntityRemarks), tissueterms) ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(biopsyterms = case_when(
    str_detect(tolower(preparations), biopsyterms) ~ 1,
    str_detect(tolower(dynamicProperties), biopsyterms) ~ 1,
    str_detect(tolower(occurrenceRemarks), biopsyterms) ~ 1,
    str_detect(tolower(materialEntityRemarks), biopsyterms) ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(swabterms = case_when(
    str_detect(tolower(preparations), swabterms) ~ 1,
    str_detect(tolower(dynamicProperties), swabterms) ~ 1,
    str_detect(tolower(occurrenceRemarks), swabterms) ~ 1,
    str_detect(tolower(materialEntityRemarks), swabterms) ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(etohterms = case_when(
    str_detect(tolower(preparations), etohterms) ~ 1,
    str_detect(tolower(dynamicProperties), etohterms) ~ 1,
    str_detect(tolower(occurrenceRemarks), etohterms) ~ 1,
    str_detect(tolower(materialEntityRemarks), etohterms) ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(freezeterms = case_when(
    str_detect(tolower(preparations), freezeterms) ~ 1,
    str_detect(tolower(dynamicProperties), freezeterms) ~ 1,
    str_detect(tolower(occurrenceRemarks), freezeterms) ~ 1,
    str_detect(tolower(materialEntityRemarks), freezeterms) ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(bufferterms = case_when(
    str_detect(tolower(preparations), bufferterms) ~ 1,
    str_detect(tolower(dynamicProperties), bufferterms) ~ 1,
    str_detect(tolower(occurrenceRemarks), bufferterms) ~ 1,
    str_detect(tolower(materialEntityRemarks), bufferterms) ~ 1,
    TRUE ~ 0
  )) %>% 
  mutate(heart = case_when(
    str_detect(tolower(preparations), "heart") ~ 1,
    str_detect(tolower(dynamicProperties), "heart") ~ 1,
    str_detect(tolower(occurrenceRemarks), "heart") ~ 1,
    str_detect(tolower(materialEntityRemarks), "heart") ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(kidney = case_when(
    str_detect(tolower(preparations), "kidney") ~ 1,
    str_detect(tolower(dynamicProperties), "kidney") ~ 1,
    str_detect(tolower(occurrenceRemarks), "kidney") ~ 1,
    str_detect(tolower(materialEntityRemarks), "kidney") ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(liver = case_when(
    str_detect(tolower(preparations), "liver") ~ 1,
    str_detect(tolower(dynamicProperties), "liver") ~ 1,
    str_detect(tolower(occurrenceRemarks), "liver") ~ 1,
    str_detect(tolower(materialEntityRemarks), "liver") ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(lung = case_when(
    str_detect(tolower(preparations), "lung") ~ 1,
    str_detect(tolower(dynamicProperties), "lung") ~ 1,
    str_detect(tolower(occurrenceRemarks), "lung") ~ 1,
    str_detect(tolower(materialEntityRemarks), "lung") ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(muscle = case_when(
    str_detect(tolower(preparations), "muscle") ~ 1,
    str_detect(tolower(dynamicProperties), "muscle") ~ 1,
    str_detect(tolower(occurrenceRemarks), "muscle") ~ 1,
    str_detect(tolower(materialEntityRemarks), "muscle") ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(spleen = case_when(
    str_detect(tolower(preparations), "spleen") ~ 1,
    str_detect(tolower(dynamicProperties), "spleen") ~ 1,
    str_detect(tolower(occurrenceRemarks), "spleen") ~ 1,
    str_detect(tolower(materialEntityRemarks), "spleen") ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(colon = case_when(
    str_detect(tolower(preparations), "colon") ~ 1,
    str_detect(tolower(dynamicProperties), "colon") ~ 1,
    str_detect(tolower(occurrenceRemarks), "colon") ~ 1,
    str_detect(tolower(materialEntityRemarks), "colon") ~ 1,
    TRUE ~ 0
  ))

table(TissueSamples$spleen)
table(TissueSamples$etohterms, TissueSamples$order)

nrow(TissueSamples)
nrow(TissueSamples %>% filter(is.na(species)))

###number of institutions

table(TissueSamples$institutionCode)
length(unique(Dataset$institutionCode))

museums <- TissueSamples %>% 
  group_by(institutionCode, publisher, order) %>% 
  dplyr::summarise(n = n()) %>%
  pivot_wider(names_from = order, values_from = n)

museums$Total <- rowSums(museums[,3:4],na.rm=TRUE)
#write.csv(museums, "Updated museum data Oct 20.csv", row.names = FALSE)

inst2 <- inst2[-which(inst2$institutionCode == "Instituto de Investigación de Recursos Biológicos Alexander von Humboldt (IAvH)"),]
inst2 <- inst2[-which(inst2$institutionCode == "Universidad Icesi (ICESI)"),]

#inst2$institutionCode[which(inst2$institutionCode == "Instituto de Investigación de Recursos Biológicos Alexander von Humboldt (IAvH)")] <- "IAvH"
#inst2$institutionCode[which(inst2$institutionCode == "Universidad Icesi (ICESI)")] <- "ICESI"

length(unique(inst2$institutionCode))

museums <- as.data.frame.matrix(table(TissueSamples$institutionCode, TissueSamples$order))

inst2 <- inst2[,1:4]
inst2[inst2 == 0] <- NA

nrow(inst2[complete.cases(inst2),])
sum(is.na(inst2$Rodentia))

###taxa -- this might be unnecessary

table(TissueSamples$species)

#making a list of unique species keys
keys = TissueSamples %>%
  filter(!duplicated(taxonKey)) %>%
  filter(!is.na(taxonKey)) %>%
  dplyr::select(taxonKey)

#querying GBIF database to match species names with species key
test = taxize::classification(keys$taxonKey, db = "gbif")

#intializing empty dataset
dataset_match = data.frame()

#looping through API response to obtain species names from 
for(n in 1:length(test)){
  temp1 = names(test[n])
  if(length(test[[n]]$name[test[[n]]$rank=="species"])>0){
    temp2 = test[[n]]$name[test[[n]]$rank=="species"]
  } else {
    temp2 = NA
  }
  temp3 = data.frame(taxonKey = temp1,
                     species = temp2)
  dataset_match = bind_rows(dataset_match, temp3)
}
  
dataset_match = dataset_match %>%
  mutate(taxonKey = as.numeric(taxonKey))

TissueSamples = TissueSamples %>% left_join(dataset_match, by = c("taxonKey"= "taxonKey"))

for (i in 1:nrow(TissueSamples)) {
  if(is.na(TissueSamples$species[i])){
    TissueSamples$species[i] <- TissueSamples$species.y[i]
  }
}

###taxonomic biases in samples

## libraries
library(ape)
library(caper)
library(phylofactor)
library(data.table)
library(treeio)
library(ggtree)

sp <- TissueSamples %>% filter(!is.na(species))
#gen <- TissueSamples %>% filter(!is.na(genus))

sp <- sp %>% dplyr::group_by(species) %>% dplyr::summarise(Samples = n())
#gen <- gen %>% group_by(genus) %>% dplyr::summarise(Samples = n())

## load taxonomy and phy
setwd("~/Documents/PhD/1. pteropodidae/phylo")
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')
taxonomy=read.csv("taxonomy_mamPhy_5911species.csv")

taxonomy=taxonomy[taxonomy$ord%in%c("CHIROPTERA","RODENTIA"),]

## fix names
taxonomy$tip=sapply(strsplit(taxonomy$tiplabel,'_'),function(x) paste(x[1],x[2],sep='_'))
sp$tip=gsub(" ","_",sp$species)

#remove any that are just genus
sp <- sp[grepl("_",sp$tip),]

## check names
setdiff(sp$tip,taxonomy$tip) #290 mismatches
diff=setdiff(sp$tip,taxonomy$tip)
#diff=setdiff(gen$genus,taxonomy$gen) #41 mismatches

#manual fixes of genera
gen$genus[which(gen$genus == "Aeorestes")] <- "Lasiurus"
gen$genus[which(gen$genus == "Afronycteris")] <- "Neoromicia"
gen$genus[which(gen$genus == "Alexandromys")] <- "Microtus"
gen$genus[which(gen$genus == "Austronomus")] <- "Tadarida"
gen$genus[which(gen$genus == "Boneia")] <- "Rousettus"
gen$genus[which(gen$genus == "Cassistrellus")] <- "Eptesicus"
gen$genus[which(gen$genus == "Dasypterus")] <- "Lasiurus"
gen$genus[which(gen$genus == "Dipodillus")] <- "Gerbillus"
gen$genus[which(gen$genus == "Doryrhina")] <- "Hipposideros"
gen$genus[which(gen$genus == "Echinoprocta")] <- "Coendou"
gen$genus[which(gen$genus == "Eptesius")] <- "Eptesicus" #typo
gen$genus[which(gen$genus == "Eutamias")] <- "Tamias"
gen$genus[which(gen$genus == "Euxerus")] <- "Xerus"
gen$genus[which(gen$genus == "Gardnerycteris")] <- "Mimon"
gen <- gen[-which(gen$genus == "Hyorhinomys"),] #this is a new genus that isn't included in upham
gen$genus[which(gen$genus == "Lissonycteris")] <- "Myonycteris"
gen$genus[which(gen$genus == "Lyroderma")] <- "Megaderma"
gen$genus[which(gen$genus == "Macronycteris")] <- "Hipposideros"
gen$genus[which(gen$genus == "Micaelamys")] <- "Aethomys"
gen$genus[which(gen$genus == "Micronomus")] <- "Mormopterus"
gen$genus[which(gen$genus == "Mirostrellus")] <- "Pipistrellus"
gen$genus[which(gen$genus == "Myotomys")] <- "Otomys"
gen$genus[which(gen$genus == "Neoeptesicus")] <- "Eptesicus"
gen$genus[which(gen$genus == "Neomicroxus")] <- "Akodon"
gen$genus[which(gen$genus == "Neoplatymops")] <- "Molossops"
gen$genus[which(gen$genus == "Neotamias")] <- "Molossops"
gen$genus[which(gen$genus == "Nesonycteris")] <- "Melonycteris"
gen$genus[which(gen$genus == "Notosciurus")] <- "Sciurus"
gen$genus[which(gen$genus == "Nyctinomus")] <- "Nyctinomops"
gen$genus[which(gen$genus == "Ozimops")] <- "Mormopterus"
gen$genus[which(gen$genus == "Paremballonura")] <- "Emballonura"
gen$genus[which(gen$genus == "Perimyotis")] <- "Pipistrellus"
gen$genus[which(gen$genus == "Pseudoromicia")] <- "Neoromicia"
gen$genus[which(gen$genus == "Salinoctomys")] <- "Tympanoctomys"
gen$genus[which(gen$genus == "Setirostris")] <- "Mormopterus"
gen$genus[which(gen$genus == "Sphiggurus")] <- "Coendou"
gen$genus[which(gen$genus == "Stenonycteris")] <- "Rousettus"
gen$genus[which(gen$genus == "Tanyuromys")] <- "Sigmodontomys"
gen$genus[which(gen$genus == "Thainycteris")] <- "Arielulus"
gen$genus[which(gen$genus == "Vampyriscus")] <- "Vampyressa"
gen$genus[which(gen$genus == "Vansonia")] <- "Pipistrellus"

#setdiff(gen$genus,taxonomy$gen)

## set key for taxonomy -- this may not be necessary if switching to genus above
library(taxize)
library(rentrez)

## set key
set_entrez_key("08bfe1677e98cb10837653b9b0dee88c9c09")
Sys.getenv("ENTREZ_KEY")

## write function
taxtake=function(x){
  y=as.data.frame(x)
  y=tail(y,1)
  return(y)
}

## save records
httr::set_config(httr::config(http_version = 0))
tax=rep(NA,length(diff))

## loop
for(i in 1:length(tax)){
  
  ## classify
  cset=classification(diff[i],db="ncbi",rows=1)
  tax[i]=sapply(cset,function(x) taxtake(as.data.frame(x)))
  #Sys.sleep(1)
  
  # tax[i]=sapply(classification(data$unique_name[i],db="ncbi",rows=1),
  #               function(x) taxtake(as.data.frame(x)))
}

## save tax
tdata=do.call(rbind.data.frame,tax)
names(tdata)="newtip"
tdata$tip=diff

## fix tip
tdata$newtip=gsub(" ","_",tdata$newtip)

## check against taxonomy
diff2 <- setdiff(tdata$newtip,taxonomy$tip)
length(setdiff(tdata$newtip,taxonomy$tip)) #228 mismatches, but NAs on top of this

## merge into data
data=merge(sp,tdata,by="tip",all.x=T)

## setdiff for tip in tree
data$tree=ifelse(data$tip%in%setdiff(data$tip,taxonomy$tip),0,1)

## if 0, use new tips
data$match=ifelse(data$tree==0,data$newtip,data$tip)

## if newtip is 0, use original
data$match=ifelse(is.na(data$match),data$tip,data$match)

## if the tip is in the tree
data$tree2=ifelse(data$match%in%setdiff(data$match,taxonomy$tip),0,1)

## set match as tip
data$tip=data$match
data$newtip=NULL
data$match=NULL
data$tree=data$tree2
data$tree2=NULL

## new setdiff
diff=setdiff(data$tip,taxonomy$tip) #284 mismatches

## export
set=data.frame(old_name=diff)
write.csv(set, "Upham mismatch.csv")

########manual relabeling using MDD
###October 5, 2025
#Current version: v2.3, released September 1, 2025

set$MDD <- NA
set$MDD[which(set$old_name=="Aegialomys_ica")] <- "Aegialomys_xanthaeolus"
set$MDD[which(set$old_name=="Lasiurus_semotus")] <- "Lasiurus_cinereus"
set$MDD[which(set$old_name=="Lasiurus_villosissimus")] <- "Lasiurus_cinereus"
set$MDD[which(set$old_name=="Afronycteris_helios")] <- "Neoromicia_helios"
set$MDD[which(set$old_name=="Afronycteris_nanus")] <- "Neoromicia_nana"
set$MDD[which(set$old_name=="Akodon_caenosus")] <- "Akodon_lutescens"
set$MDD[which(set$old_name=="Alexandromys_middendorffii")] <- "Microtus_middendorffii"
set$MDD[which(set$old_name=="Ammospermophilus_insularis")] <- "Ammospermophilus_leucurus"
set$MDD[which(set$old_name=="Andalgalomys_roigi")] <- "Andalgalomys_olrogi"
set$MDD[which(set$old_name=="Anoura_aequatoris")] <- "Anoura_caudifer"
set$MDD[which(set$old_name=="Anoura_cadenai")] <- "recently described"
set$MDD[which(set$old_name=="Anoura_javieri")] <- "recently described"
set$MDD[which(set$old_name=="Anoura_peruana")] <- "Anoura_geoffroyi"
set$MDD[which(set$old_name=="Artibeus_aequatorialis")] <- "Artibeus_jamaicensis"
set$MDD[which(set$old_name=="Artibeus_anderseni")] <- "Dermanura_anderseni"
set$MDD[which(set$old_name=="Artibeus_aztecus")] <- "Dermanura_aztecus"
set$MDD[which(set$old_name=="Artibeus_cinereus")] <- "Dermanura_cinereus"
set$MDD[which(set$old_name=="Artibeus_glaucus")] <- "Dermanura_glaucus"
set$MDD[which(set$old_name=="Artibeus_ravus")] <- "Dermanura_phaeotis"
set$MDD[which(set$old_name=="Artibeus_rosenbergi")] <- "Dermanura_glaucus"
set$MDD[which(set$old_name=="Artibeus_schwartzi")] <- "Artibeus_jamaicensis"
set$MDD[which(set$old_name=="Artibeus_toltecus")] <- "Dermanura_toltecus"
set$MDD[which(set$old_name=="Artibeus_glaucus_watsoni")] <- "Dermanura_glaucus"
set$MDD[which(set$old_name=="Barbastella_darjelingensis")] <- "Barbastella_leucomelas"
set$MDD[which(set$old_name=="Batomys_hamiguitan")] <- "recently described"
set$MDD[which(set$old_name=="Beamys_major")] <- "Beamys_hindei"
set$MDD[which(set$old_name=="Boneia_bidens")] <- "Rousettus_bidens"
set$MDD[which(set$old_name=="Calomys_achaku")] <- "recently described"
set$MDD[which(set$old_name=="Carollia_brevicaudum")] <- "Carollia_brevicauda"
set$MDD[which(set$old_name=="Cassistrellus_dimissus")] <- "Eptesicus_dimissus"
set$MDD[which(set$old_name=="Cavia_nana")] <- "Cavia_aperea"
set$MDD[which(set$old_name=="Mops_jobimena")] <- "Tadarida_jobimena"
set$MDD[which(set$old_name=="Mops_leucogaster")] <- "recently described"
set$MDD[which(set$old_name=="Chaerephon_pusillus")] <- "recently described"
set$MDD[which(set$old_name=="Chaetodipus_arenarius_siccus")] <- "Chaetodipus_arenarius"
set$MDD[which(set$old_name=="Chilomys_sp._p_NT-2024")] <- "recently described" #Chilomys percequilloi
set$MDD[which(set$old_name=="Chilonatalus_macer")] <- "Chilonatalus_micropus"
set$MDD[which(set$old_name=="Chiroderma_gorgasi")] <- "Chiroderma_trinitatum"
set$MDD[which(set$old_name=="Chiroderma_scopaeum")] <- "Chiroderma_salvini"
set$MDD[which(set$old_name=="Chironax_tumulus")] <- "Chironax_melanocephalus"
set$MDD[which(set$old_name=="Coccymys_shawmayeri")] <- "Coccymys_ruemmleri"
set$MDD[which(set$old_name=="Cricetomys_ansorgei")] <- "Cricetomys_gambianus"
set$MDD[which(set$old_name=="Fukomys_mechowii")] <- "Fukomys_mechowi"
set$MDD[which(set$old_name=="Cynomops_paranus_milleri")] <- "Cynomops_paranus"
set$MDD[which(set$old_name=="Dasymys_alleni")] <- "Dasymys_incomtus"
set$MDD[which(set$old_name=="Dasymys_rwandae")] <- "recently described"
set$MDD[which(set$old_name=="Dasymys_sua")] <- "recently described"
set$MDD[which(set$old_name=="Dendromus_nyasae")] <- "Dendromus_kivu"
set$MDD[which(set$old_name=="Diclidurus_isabella")] <- "Diclidurus_isabellus"
set$MDD[which(set$old_name=="Dipodillus_campestris")] <- "Gerbillus_campestris"
set$MDD[which(set$old_name=="Dipodillus_dasyurus")] <- "Gerbillus_dasyurus"
set$MDD[which(set$old_name=="Dipodillus_harwoodi")] <- "Gerbillus_harwoodi"
set$MDD[which(set$old_name=="Dipodillus_simoni")] <- "Gerbillus_simoni"
set$MDD[which(set$old_name=="Dipodillus_somalicus")] <- "Gerbillus_somalicus"
set$MDD[which(set$old_name=="Dipodomys_ornatus")] <- "Dipodomys_phillipsii"
set$MDD[which(set$old_name=="Doryrhina_camerunensis")] <- "Hipposideros_camerunensis"
set$MDD[which(set$old_name=="Doryrhina_cyclops")] <- "Hipposideros_cyclops"
set$MDD[which(set$old_name=="Doryrhina_muscinus")] <- "Hipposideros_muscinus"
set$MDD[which(set$old_name=="Doryrhina_semoni")] <- "Hipposideros_semoni"
set$MDD[which(set$old_name=="Doryrhina_stenotis")] <- "Hipposideros_stenotis"
set$MDD[which(set$old_name=="Doryrhina_wollastoni")] <- "Hipposideros_wollastoni"
set$MDD[which(set$old_name=="Echimys_guianae")] <- "not in MDD"
set$MDD[which(set$old_name=="Epomophorus_dobsonii")] <- "Epomops_dobsonii"
set$MDD[which(set$old_name=="Epomophorus_minor")] <- "Epomophorus_minimus"
set$MDD[which(set$old_name=="Epomophorus_pusillus")] <- "Micropteropus_pusillus"
set$MDD[which(set$old_name=="Cnephaeus_ognevi")] <- "Eptesicus_bottae"
set$MDD[which(set$old_name=="Erethizon_dorsatus")] <- "Erethizon_dorsatum"
set$MDD[which(set$old_name=="Eumops_nanus")] <- "Eumops_bonariensis"
set$MDD[which(set$old_name=="Gardnerycteris_crenulata")] <- "Mimon_crenulatum"
set$MDD[which(set$old_name=="Gardnerycteris_keenani")] <- "Mimon_crenulatum"
set$MDD[which(set$old_name=="Geomys_bursarius_lutescens")] <- "Geomys_bursarius"
set$MDD[which(set$old_name=="Geoxus_annectens")] <- "Pearsonomys_annectens"
set$MDD[which(set$old_name=="Gerbilliscus_vicinus")] <- "Gerbilliscus_robustus"
set$MDD[which(set$old_name=="Glossophaga_antillarum")] <- "Glossophaga_soricina"
set$MDD[which(set$old_name=="Glossophaga_mutica")] <- "Glossophaga_soricina"
set$MDD[which(set$old_name=="Glossophaga_valens")] <- "Glossophaga_soricina"
set$MDD[which(set$old_name=="Grammomys_poensis")] <- "recently described"
set$MDD[which(set$old_name=="Graomys_chacoensis")] <- "Graomys_centralis"
set$MDD[which(set$old_name=="Heteromys_goldmani")] <- "Heteromys_desmarestianus"
set$MDD[which(set$old_name=="Hipposideros_sp._ARR-2015")] <- "recently described" #Hipposideros cryptovalorona
set$MDD[which(set$old_name=="Hipposideros_gentilis")] <- "Hipposideros_megalotis"
set$MDD[which(set$old_name=="Hipposideros_pendleburyi")] <- "Hipposideros_turpis"
set$MDD[which(set$old_name=="Hipposideros_tephrus")] <- "Hipposideros_caffer"
set$MDD[which(set$old_name=="Holochilus_vulpinus")] <- "Holochilus_brasiliensis"
set$MDD[which(set$old_name=="Hylomyscus_heinrichorum")] <- "recently described"
set$MDD[which(set$old_name=="Hylomyscus_mpungamachagorum")] <- "recently described"
set$MDD[which(set$old_name=="Hylomyscus_pygmaeus")] <- "recently described"
set$MDD[which(set$old_name=="Hylomyscus_stanleyi")] <- "recently described"
set$MDD[which(set$old_name=="Hylomyscus_thornesmithae")] <- "recently described"
set$MDD[which(set$old_name=="Hylopetes_sagitta")] <- "Petinomys_sagitta"
set$MDD[which(set$old_name=="Hyorhinomys_stuempkei")] <- "recently described"
set$MDD[which(set$old_name=="Hypsugo_alaschanicus")] <- "Pipistrellus_savii"
set$MDD[which(set$old_name=="Hypsugo_ariel")] <- "Pipistrellus_ariel"
set$MDD[which(set$old_name=="Hypsugo_cadornae")] <- "Pipistrellus_cadornae"
set$MDD[which(set$old_name=="Hypsugo_imbricatus")] <- "Pipistrellus_imbricatus"
set$MDD[which(set$old_name=="Hypsugo_macrotis")] <- "Pipistrellus_macrotis"
set$MDD[which(set$old_name=="Hypsugo_petersi")] <- "Falsistrellus_petersi"
set$MDD[which(set$old_name=="Hypsugo_pulveratus")] <- "Pipistrellus_pulveratus"
set$MDD[which(set$old_name=="Hypsugo_savii")] <- "Pipistrellus_savii"
set$MDD[which(set$old_name=="Laephotis_capensis")] <- "Neoromicia_capensis"
set$MDD[which(set$old_name=="Laephotis_kirinyaga")] <- "recently described"
set$MDD[which(set$old_name=="Laephotis_malagasyensis")] <- "Neoromicia_malagasyensis"
set$MDD[which(set$old_name=="Laephotis_matroka")] <- "Neoromicia_matroka"
set$MDD[which(set$old_name=="Laephotis_robertsi")] <- "Neoromicia_robertsi"
set$MDD[which(set$old_name=="Laephotis_stanleyi")] <- "recently described"
set$MDD[which(set$old_name=="Lasiurus_frantzii")] <- "Lasiurus_blossevillii"
set$MDD[which(set$old_name=="Leptomys_paulus")] <- "recently described"
set$MDD[which(set$old_name=="Lissonycteris_angolensis")] <- "Myonycteris_angolensis"
set$MDD[which(set$old_name=="Lophostoma_occidentale")] <- "Lophostoma_aequatorialis"
set$MDD[which(set$old_name=="Lophostoma_silvicola")] <- "Lophostoma_silvicolum"
set$MDD[which(set$old_name=="Lophuromys_angolensis")] <- "recently described"
set$MDD[which(set$old_name=="Lophuromys_ansorgei")] <- "Lophuromys_sikapusi"
set$MDD[which(set$old_name=="Lophuromys_aquilus")] <- "Lophuromys_cinereus"
set$MDD[which(set$old_name=="Lophuromys_dudui")] <- "Lophuromys_flavopunctatus"
set$MDD[which(set$old_name=="Lophuromys_laticeps")] <- "not in Upham"
set$MDD[which(set$old_name=="Lophuromys_simensis")] <- "Lophuromys_flavopunctatus"
set$MDD[which(set$old_name=="Lophuromys_verhageni")] <- "Lophuromys_flavopunctatus"
set$MDD[which(set$old_name=="Lophuromys_zena")] <- "not in Upham"
set$MDD[which(set$old_name=="Lyroderma_lyra")] <- "Megaderma_lyra"
set$MDD[which(set$old_name=="Macronycteris_commersonii")] <- "Hipposideros_commersoni"
set$MDD[which(set$old_name=="Macronycteris_gigas")] <- "Hipposideros_gigas"
set$MDD[which(set$old_name=="Macronycteris_vittatus")] <- "Hipposideros_vittatus"
set$MDD[which(set$old_name=="Leiuromys_occasius")] <- "Pattonomys_occasius"
set$MDD[which(set$old_name=="Melanomys_chrysomelas")] <- "Melanomys_caliginosus"
set$MDD[which(set$old_name=="Melanomys_columbianus")] <- "Melanomys_caliginosus"
set$MDD[which(set$old_name=="Micaelamys_namaquensis")] <- "Aethomys_namaquensis"
set$MDD[which(set$old_name=="Micronomus_norfolkensis")] <- "Mormopterus_norfolkensis"
set$MDD[which(set$old_name=="Micronycteris_simmonsae")] <- "recently described"
set$MDD[which(set$old_name=="Micronycteris_tresamici")] <- "recently described"
set$MDD[which(set$old_name=="Miniopterus_africanus")] <- "not in Upham"
set$MDD[which(set$old_name=="Miniopterus_natalensis_arenarius")] <- "Miniopterus_natalensis"
set$MDD[which(set$old_name=="Miniopterus_blepotis")] <- "Miniopterus_schreibersii"
set$MDD[which(set$old_name=="Miniopterus_eschscholtzii")] <- "Miniopterus_schreibersii"
set$MDD[which(set$old_name=="Miniopterus_orianae")] <- "Miniopterus_schreibersii"
set$MDD[which(set$old_name=="Miniopterus_schreibersii_pallidus")] <- "Miniopterus_schreibersii"
set$MDD[which(set$old_name=="Mirostrellus_joffrei")] <- "Pipistrellus_joffrei"
set$MDD[which(set$old_name=="Molossus_currentium_bondae")] <- "Molossus_currentium"
set$MDD[which(set$old_name=="Molossus_melini")] <- "recently described"
set$MDD[which(set$old_name=="Molossus_milleri")] <- "not in Upham"
set$MDD[which(set$old_name=="Molossus_nigricans")] <- "Molossus_rufus"
set$MDD[which(set$old_name=="Molossus_verrilli")] <- "Molossus_molossus"
set$MDD[which(set$old_name=="Mormopterus_kalinowskii")] <- "not in Upham"
set$MDD[which(set$old_name=="Murina_feae")] <- "Murina_cineracea"
set$MDD[which(set$old_name=="Mylomys_cuninghamei")] <- "Mylomys_dybowskii"
set$MDD[which(set$old_name=="Myonycteris_leptodon")] <- "Myonycteris_torquata"
set$MDD[which(set$old_name=="Myotis_annatessae")] <- "recently described"
set$MDD[which(set$old_name=="Myotis_armiensis")] <- "recently described"
set$MDD[which(set$old_name=="Myotis_borneoensis")] <- "Myotis_montivagus"
set$MDD[which(set$old_name=="Myotis_muricola_browni")] <- "Myotis_muricola"
set$MDD[which(set$old_name=="Myotis_caucensis")] <- "Myotis_nigricans"
set$MDD[which(set$old_name=="Myotis_clydejonesi")] <- "recently described"
set$MDD[which(set$old_name=="Myotis_goudotii")] <- "Myotis_goudoti"
set$MDD[which(set$old_name=="Myotis_martiniquensis_nyctor")] <- "Myotis_martiniquensis"
set$MDD[which(set$old_name=="Myotis_peytoni")] <- "Myotis_montivagus"
set$MDD[which(set$old_name=="Myotis_pilosatibialis")] <- "Myotis_keaysi"
set$MDD[which(set$old_name=="Myotis_sibiricus")] <- "Myotis_brandtii"
set$MDD[which(set$old_name=="Neacomys_rosalindae")] <- "recently described"
set$MDD[which(set$old_name=="Neacomys_vargasllosai")] <- "recently described"
set$MDD[which(set$old_name=="Neacomys_xingu")] <- "recently described"
set$MDD[which(set$old_name=="Nectomys_grandis")] <- "Nectomys_magdalenae"
set$MDD[which(set$old_name=="Neomicroxus_bogotensis")] <- "Akodon_bogotensis"
set$MDD[which(set$old_name=="Neoromicia_anchietae")] <- "Hypsugo_anchietae"
set$MDD[which(set$old_name=="Neoromicia_bemainty")] <- "Hypsugo_bemainty"
set$MDD[which(set$old_name=="Neotamias_ruficaudus")] <- "Tamias_ruficaudus"
set$MDD[which(set$old_name=="Neotoma_melanura")] <- "Neotoma_albigula"
set$MDD[which(set$old_name=="Nephelomys_childi")] <- "Nephelomys_albigularis"
set$MDD[which(set$old_name=="Nephelomys_moerex")] <- "Nephelomys_albigularis"
set$MDD[which(set$old_name=="Nephelomys_pectoralis")] <- "Nephelomys_albigularis"
set$MDD[which(set$old_name=="Nephelomys_pirrensis")] <- "Nephelomys_albigularis"
set$MDD[which(set$old_name=="Nesonycteris_fardoulisi")] <- "Melonycteris_fardoulisi"
set$MDD[which(set$old_name=="Niviventer_mekongis")] <- "Niviventer_fulvescens"
set$MDD[which(set$old_name=="Notopteris_macdonaldii")] <- "Notopteris_macdonaldi"
set$MDD[which(set$old_name=="Notopteris_neocaledonicus")] <- "Notopteris_neocaledonica"
set$MDD[which(set$old_name=="Notosciurus_granatensis")] <- "Sciurus_granatensis"
set$MDD[which(set$old_name=="Nycticeinops_crassulus")] <- "Pipistrellus_crassulus"
set$MDD[which(set$old_name=="Nycticeinops_eisentrauti")] <- "Pipistrellus_eisentrauti"
set$MDD[which(set$old_name=="Nycticeinops_grandidieri")] <- "Neoromicia_flavescens"
set$MDD[which(set$old_name=="Nyctophilus_daedalus")] <- "Nyctophilus_bifax"
set$MDD[which(set$old_name=="Nyctophilus_holtorum")] <- "recently described"
set$MDD[which(set$old_name=="Nyctophilus_major")] <- "Nyctophilus_timoriensis"
set$MDD[which(set$old_name=="Octodon_ricardojeda")] <- "Octodon_bridgesi"
set$MDD[which(set$old_name=="Oecomys_franciscorum")] <- "recently described"
set$MDD[which(set$old_name=="Oligoryzomys_guille")] <- "recently described"
set$MDD[which(set$old_name=="Oligoryzomys_mattogrossae")] <- "Oligoryzomys_microtis"
set$MDD[which(set$old_name=="Oligoryzomys_messorius")] <- "Oligoryzomys_fulvescens"
set$MDD[which(set$old_name=="Oligoryzomys_sp._ERS-2016")] <- "Oligoryzomys_longicaudatus" #Oligoryzomys yatesi
set$MDD[which(set$old_name=="Oryzomys_chapmani")] <- "Handleyomys_chapmani"
set$MDD[which(set$old_name=="Otomops_harrisoni")] <- "recently described"
set$MDD[which(set$old_name=="Ozimops_beccarii")] <- "Mormopterus_beccarii"
set$MDD[which(set$old_name=="Ozimops_cobourgianus")] <- "Mormopterus_loriae"
set$MDD[which(set$old_name=="Ozimops_halli")] <- "Mormopterus_halli"
set$MDD[which(set$old_name=="Ozimops_kitcheneri")] <- "Mormopterus_kitcheneri"
set$MDD[which(set$old_name=="Ozimops_loriae")] <- "Mormopterus_loriae"
set$MDD[which(set$old_name=="Ozimops_lumsdenae")] <- "Mormopterus_lumsdenae"
set$MDD[which(set$old_name=="Ozimops_petersi")] <- "Mormopterus_planiceps"
set$MDD[which(set$old_name=="Ozimops_planiceps")] <- "Mormopterus_planiceps"
set$MDD[which(set$old_name=="Ozimops_ridei")] <- "Mormopterus_loriae"
set$MDD[which(set$old_name=="Paratriaenops_furcula")] <- "Paratriaenops_furculus"
set$MDD[which(set$old_name=="Paremballonura_atrata")] <- "Emballonura_atrata"
set$MDD[which(set$old_name=="Paremballonura_tiavato")] <- "Emballonura_tiavato"
set$MDD[which(set$old_name=="Penthetor_lucasii")] <- "Penthetor_lucasi"
set$MDD[which(set$old_name=="Perimyotis_subflavus")] <- "Pipistrellus_subflavus"
set$MDD[which(set$old_name=="Perognathus_mollipilosus")] <- "Perognathus_parvus"
set$MDD[which(set$old_name=="Peromyscus_bakeri")] <- "recently described"
set$MDD[which(set$old_name=="Peromyscus_carolpattonae")] <- "recently described"
set$MDD[which(set$old_name=="Peromyscus_sp._CL-2015b")] <- "recently described" #Peromyscus gardneri
set$MDD[which(set$old_name=="Peromyscus_kilpatricki")] <- "recently described"
set$MDD[which(set$old_name=="Peromyscus_nudipes")] <- "Peromyscus_mexicanus"
set$MDD[which(set$old_name=="Platyrrhinus_incarum")] <- "Platyrrhinus_helleri"
set$MDD[which(set$old_name=="Plecotus_kozlovi")] <- "Plecotus_austriacus"
set$MDD[which(set$old_name=="Plecotus_strelkovi")] <- "recently described"
set$MDD[which(set$old_name=="Plecotus_wardi")] <- "Plecotus_austriacus"
set$MDD[which(set$old_name=="Proechimys_trinitatis")] <- "Proechimys_urichi"
set$MDD[which(set$old_name=="Promops_davisoni")] <- "Promops_centralis"
set$MDD[which(set$old_name=="Pseudoromicia_brunnea")] <- "Neoromicia_brunnea"
set$MDD[which(set$old_name=="Pseudoromicia_kityoi")] <- "recently described"
set$MDD[which(set$old_name=="Pseudoromicia_nyanza")] <- "recently described"
set$MDD[which(set$old_name=="Pseudoromicia_rendalli")] <- "Neoromicia_rendalli"
set$MDD[which(set$old_name=="Pseudoromicia_tenuipinnis")] <- "Neoromicia_tenuipinnis"
set$MDD[which(set$old_name=="Ptenochirus_jagorii")] <- "Ptenochirus_jagori"
set$MDD[which(set$old_name=="Ptenochirus_wetmorei")] <- "Megaerops_wetmorei"
set$MDD[which(set$old_name=="Pteronotus_fulvus")] <- "Pteronotus_davyi"
set$MDD[which(set$old_name=="Pteronotus_parnellii_fuscus")] <- "Pteronotus_parnellii"
set$MDD[which(set$old_name=="Pteronotus_mesoamericanus")] <- "Pteronotus_parnellii"
set$MDD[which(set$old_name=="Pteronotus_parnellii_mexicanus")] <- "Pteronotus_parnellii"
set$MDD[which(set$old_name=="Pteronotus_psilotis")] <- "Pteronotus_personatus"
set$MDD[which(set$old_name=="Pteronotus_pusillus")] <- "Pteronotus_parnellii"
set$MDD[which(set$old_name=="Pteronotus_rubiginosus")] <- "Pteronotus_parnellii"
set$MDD[which(set$old_name=="Pteropus_allenorum")] <- "recently described"
set$MDD[which(set$old_name=="Pteropus_coxi")] <- "recently described"
set$MDD[which(set$old_name=="Pteropus_capistratus_ennisae")] <- "Pteropus_capistratus"
set$MDD[which(set$old_name=="Pteropus_medius")] <- "Pteropus_giganteus"
set$MDD[which(set$old_name=="Pteropus_vetula")] <- "Pteropus_vetulus"
set$MDD[which(set$old_name=="Pygeretmus_shitkovi")] <- "Pygeretmus_zhitkovi"
set$MDD[which(set$old_name=="Reithrodon_physodes")] <- "Reithrodon_auritus"
set$MDD[which(set$old_name=="Rhabdomys_dilectus")] <- "Rhabdomys_pumilio"
set$MDD[which(set$old_name=="Rhinolophus_cornutus")] <- "not in Upham"
set$MDD[which(set$old_name=="Rhinolophus_dobsoni")] <- "Rhinolophus_landeri"
set$MDD[which(set$old_name=="Rhinolophus_macrotis_episcopus")] <- "Rhinolophus_macrotis"
set$MDD[which(set$old_name=="Rhinolophus_lobatus")] <- "Rhinolophus_landeri"
set$MDD[which(set$old_name=="Rhinolophus_microglobosus")] <- "Rhinolophus_stheno"
set$MDD[which(set$old_name=="Rhinolophus_monoceros")] <- "Rhinolophus_pusillus"
set$MDD[which(set$old_name=="Rhinolophus_nippon")] <- "Rhinolophus_ferrumequinum"
set$MDD[which(set$old_name=="Rhinolophus_perditus")] <- "not in Upham"
set$MDD[which(set$old_name=="Rhinolophus_perniger")] <- "Rhinolophus_luctus"
set$MDD[which(set$old_name=="Rhinolophus_refulgens")] <- "Rhinolophus_lepidus"
set$MDD[which(set$old_name=="Rhinolophus_tatar")] <- "Rhinolophus_euryotis"
set$MDD[which(set$old_name=="Rhinopoma_cystops")] <- "Rhinopoma_macinnesi"
set$MDD[which(set$old_name=="Rhogeessa_aenea")] <- "Rhogeessa_aeneus"
set$MDD[which(set$old_name=="Rhogeessa_velilla")] <- "Rhogeessa_io"
set$MDD[which(set$old_name=="Rhynchomys_labo")] <- "recently described"
set$MDD[which(set$old_name=="Rhynchomys_mingan")] <- "recently described"
set$MDD[which(set$old_name=="Scotoecus_albigula")] <- "Scotoecus_hirundo"
set$MDD[which(set$old_name=="Scotoecus_hindei")] <- "Scotoecus_hirundo"
set$MDD[which(set$old_name=="Scotophilus_altilis")] <- "Scotophilus_leucogaster"
set$MDD[which(set$old_name=="Scotophilus_colias")] <- "Scotophilus_nigrita"
set$MDD[which(set$old_name=="Scotophilus_viridis_nigritellus")] <- "Scotophilus_viridis"
set$MDD[which(set$old_name=="Setirostris_eleryi")] <- "Mormopterus_eleryi"
set$MDD[which(set$old_name=="Sigmodon_zanjonensis")] <- "Sigmodon_hispidus"
set$MDD[which(set$old_name=="Soricomys_leonardocoi")] <- "Soricomys_leonardcoi"
set$MDD[which(set$old_name=="Stenonycteris_lanosa")] <- "Rousettus_lanosus"
set$MDD[which(set$old_name=="Sturnira_giannae")] <- "recently described"
set$MDD[which(set$old_name=="Tachyoryctes_ankoliae")] <- "Tachyoryctes_splendens"
set$MDD[which(set$old_name=="Tachyoryctes_audax")] <- "Tachyoryctes_splendens"
set$MDD[which(set$old_name=="Tachyoryctes_daemon")] <- "Tachyoryctes_splendens"
set$MDD[which(set$old_name=="Tachyoryctes_ruandae")] <- "Tachyoryctes_splendens"
set$MDD[which(set$old_name=="Tachyoryctes_ruddi")] <- "Tachyoryctes_splendens"
set$MDD[which(set$old_name=="Tachyoryctes_spalacinus")] <- "Tachyoryctes_splendens"
set$MDD[which(set$old_name=="Tadarida_laticaudata")] <- "Nyctinomops_laticaudatus"
set$MDD[which(set$old_name=="Tamiops_mcclellandii")] <- "Tamiops_macclellandii"
set$MDD[which(set$old_name=="Taterillus_harringtoni")] <- "Taterillus_emini"
set$MDD[which(set$old_name=="Thainycteris_aureocollaris")] <- "Arielulus_aureocollaris"
set$MDD[which(set$old_name=="Thomasomys_dispar")] <- "Thomasomys_cinereiventer"
set$MDD[which(set$old_name=="Thomasomys_fumeus")] <- "Thomasomys_rhoadsi"
set$MDD[which(set$old_name=="Thomasomys_nicefori")] <- "Thomasomys_popayanus"
set$MDD[which(set$old_name=="Thomasomys_princeps")] <- "Thomasomys_aureus"
set$MDD[which(set$old_name=="Tonatia_bakeri")] <- "Tonatia_saurophila"
set$MDD[which(set$old_name=="Tonatia_maresi")] <- "Tonatia_saurophila"
set$MDD[which(set$old_name=="Tylonycteris_fulvida")] <- "Tylonycteris_pachypus"
set$MDD[which(set$old_name=="Uroderma_convexum")] <- "Uroderma_bilobatum"
set$MDD[which(set$old_name=="Uroderma_davisi")] <- "Uroderma_bilobatum"
set$MDD[which(set$old_name=="Vampyriscus_nymphaeus")] <- "Vampyressa_nymphaea"
set$MDD[which(set$old_name=="Vampyrodes_major")] <- "Vampyrodes_caraccioli"
set$MDD[which(set$old_name=="Vansonia_rueppellii")] <- "Pipistrellus_rueppellii"

table(set$MDD)

for (i in 1:nrow(data)) {
  if(data$tree[i] == 0) {
    data$tip[i] <- set$MDD[which(set$old_name == data$tip[i])]
  }
}

setdiff(data$tip,taxonomy$tip)

data <- data[-which(data$tip == "recently described"),] 
data <- data[-which(data$tip == "not in MDD"),] 
data <- data[-which(data$tip == "not in Upham"),] 

setwd("~/Documents/PhD/colleen MS")
#write.csv(data, "harmonized species.csv", row.names = FALSE)

################phylogenetic signal and factorization

setwd("~/Documents/PhD/colleen MS")
data <- read.csv("harmonized species.csv")

## load taxonomy and phy
setwd("~/Documents/PhD/1. pteropodidae/phylo")
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')
taxonomy=read.csv("taxonomy_mamPhy_5911species.csv")

taxonomy=taxonomy[taxonomy$ord%in%c("CHIROPTERA","RODENTIA"),]

## fix names
taxonomy$tip=sapply(strsplit(taxonomy$tiplabel,'_'),function(x) paste(x[1],x[2],sep='_'))

tax_bat=taxonomy[taxonomy$ord%in%c("CHIROPTERA"),]
tax_rod=taxonomy[taxonomy$ord%in%c("RODENTIA"),]

#data$tip <- as.factor(data$ip)
data_sp <- data %>% group_by(tip) %>% dplyr::summarise(Samples = sum(Samples))

## merge
data_sp=merge(data_sp,taxonomy,by="tip") ## sampled

data_rod <- data_sp[which(data_sp$ord == "RODENTIA"),]
data_bat <- data_sp[which(data_sp$ord == "CHIROPTERA"),]

########

## how many in tree but not in data
setdiff(taxonomy$tip,data$tip)

## mark 0/1 for sampled
tdata=taxonomy
tdata$sampled=ifelse(tdata$tip%in%data_sp$tip,1,0)
tdata_bat <- tdata[which(tdata$ord == "CHIROPTERA"),]
tdata_rod <- tdata[which(tdata$ord == "RODENTIA"),]

## subset data
#data$tree=ifelse(data$tip%in%taxonomy$tip,1,0)
table(data$tree)
table(tdata$sampled)
#data=data[data$tree==1,]

## make simple names
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep='_'))

## trim phylo
rtree=tree

tree=keep.tip(rtree,data_sp$tip) ## sampled bats and rodents
tree_bats=keep.tip(rtree,data_bat$tip)
tree_rods=keep.tip(rtree,data_rod$tip)

atree=keep.tip(rtree,taxonomy$tip) ## all bats and rodents
atree_bats=keep.tip(rtree,tax_bat$tip)
atree_rods=keep.tip(rtree,tax_rod$tip)

## match
bdata=data_sp[match(tree$tip.label,data_sp$tip),] ## sampled
bdata_bat=data_bat[match(tree_bats$tip.label,data_bat$tip),] ## sampled
bdata_rod=data_rod[match(tree_rods$tip.label,data_rod$tip),] ## sampled

adata=tdata[match(atree$tip.label,tdata$tip),] ## all
adata_bat=tdata_bat[match(atree_bats$tip.label,tdata_bat$tip),] ## sampled
adata_rod=tdata_rod[match(atree_rods$tip.label,tdata_rod$tip),] ## sampled

## save
bdata$label=bdata$tip
bdata$Species=bdata$tip
adata$label=adata$tip
adata$Species=adata$tip
bdata_bat$label=bdata_bat$tip
bdata_bat$Species=bdata_bat$tip
adata_bat$label=adata_bat$tip
adata_bat$Species=adata_bat$tip
bdata_rod$label=bdata_rod$tip
bdata_rod$Species=bdata_rod$tip
adata_rod$label=adata_rod$tip
adata_rod$Species=adata_rod$tip

## merge
#cdata=comparative.data(phy=tree,data=bdata,names.col=tip,vcv=T,na.omit=F,warn.dropped=T)
#adata=comparative.data(phy=atree,data=adata,names.col=tip,vcv=T,na.omit=F,warn.dropped=T)

cdata_bat=comparative.data(phy=tree_bats,data=bdata_bat,names.col=tip,vcv=T,na.omit=F,warn.dropped=T)
adata_bat=comparative.data(phy=atree_bats,data=adata_bat,names.col=tip,vcv=T,na.omit=F,warn.dropped=T)

cdata_rod=comparative.data(phy=tree_rods,data=bdata_rod,names.col=tip,vcv=T,na.omit=F,warn.dropped=T)
adata_rod=comparative.data(phy=atree_rods,data=adata_rod,names.col=tip,vcv=T,na.omit=F,warn.dropped=T)

## new tip
#cdata$data$tip=cdata$data$label
#adata$data$tip=adata$data$label

cdata_bat$data$tip=cdata_bat$data$label
adata_bat$data$tip=adata_bat$data$label

cdata_rod$data$tip=cdata_rod$data$label
adata_rod$data$tip=adata_rod$data$label

## Holm rejection procedure
HolmProcedure <- function(pf,FWER=0.05){
  
  ## get split variable
  cs=names(coef(pf$models[[1]]))[-1]
  split=ifelse(length(cs)>1,cs[3],cs[1])
  
  ## obtain p values
  if (pf$models[[1]]$family$family%in%c('gaussian',"Gamma","quasipoisson")){
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|t|)'])
  } else {
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|z|)'])
  }
  D <- length(pf$tree$tip.label)
  
  ## this is the line for Holm's sequentially rejective cutoff
  keepers <- pvals<=(FWER/(2*D-3 - 2*(0:(pf$nfactors-1))))
  
  
  if (!all(keepers)){
    nfactors <- min(which(!keepers))-1
  } else {
    nfactors <- pf$nfactors
  }
  return(nfactors)
}

## get species in a clade
cladeget=function(pf,factor){
  spp=pf$tree$tip.label[pf$groups[[factor]][[1]]]
  return(spp)
}

## summarize pf object 
pfsum=function(pf){
  
  ## get formula
  chars=as.character(pf$frmla.phylo)[-1]
  
  ## response
  resp=chars[1]
  
  ## holm
  hp=HolmProcedure(pf)
  
  ## save model
  model=chars[2]
  
  ## set key
  setkey(pf$Data,'Species')
  
  ## make data
  dat=data.frame(pf$Data)
  
  ## make clade columns in data
  for(i in 1:hp){
    
    dat[,paste0(resp,'_pf',i)]=ifelse(dat$Species%in%cladeget(pf,i),'factor','other')
    
  }
  
  ## make data frame to store taxa name, response, mean, and other
  results=data.frame(matrix(ncol=6, nrow = hp))
  colnames(results)=c('factor','taxa','tips','node',"clade",'other')
  
  ## set taxonomy
  taxonomy=dat[c('Species','taxonomy')]
  taxonomy$taxonomy=as.character(taxonomy$taxonomy)
  
  ## loop
  for(i in 1:hp){
    
    ## get taxa
    tx=pf.taxa(pf,taxonomy,factor=i)$group1
    
    ## get tail
    tx=sapply(strsplit(tx,'; '),function(x) tail(x,1))
    
    ## combine
    tx=paste(tx,collapse=', ')
    
    # save
    results[i,'factor']=i
    results[i,'taxa']=tx
    
    ## get node
    tips=cladeget(pf,i)
    node=ggtree::MRCA(pf$tree,tips)
    results[i,'tips']=length(tips)
    results[i,'node']=ifelse(is.null(node) & length(tips)==1,'species',
                             ifelse(is.null(node) & length(tips)!=1,NA,node))
    
    ## get means
    ms=(tapply(dat[,resp],dat[,paste0(resp,'_pf',i)],mean))
    
    ## add in
    results[i,'clade']=ms['factor']
    results[i,'other']=ms['other']
    
  }
  
  ## return
  return(list(set=dat,results=results))
}

## taxonomy
#cdata$data$taxonomy=paste(cdata$data$ord,cdata$data$fam,cdata$data$gen,cdata$data$Species,sep='; ')
#adata$data$taxonomy=paste(adata$data$ord,adata$data$fam,adata$data$gen,adata$data$Species,sep='; ')

cdata_bat$data$taxonomy=paste(cdata_bat$data$ord,cdata_bat$data$fam,cdata_bat$data$gen,cdata_bat$data$Species,sep='; ')
adata_bat$data$taxonomy=paste(adata_bat$data$ord,adata_bat$data$fam,adata_bat$data$gen,adata_bat$data$Species,sep='; ')

cdata_rod$data$taxonomy=paste(cdata_rod$data$ord,cdata_rod$data$fam,cdata_rod$data$gen,cdata_rod$data$Species,sep='; ')
adata_rod$data$taxonomy=paste(adata_rod$data$ord,adata_rod$data$fam,adata_rod$data$gen,adata_rod$data$Species,sep='; ')

## set taxonomy
#taxonomy=data.frame(adata$data$taxonomy)
#names(taxonomy)="taxonomy"
#taxonomy$Species=rownames(adata$data)
#taxonomy=taxonomy[c("Species","taxonomy")]
#taxonomy$taxonomy=as.character(taxonomy$taxonomy)

## D statistic on sampled/not sampled

#set.seed(1)
#dstat=phylo.d(data=adata,binvar=sampled,permut=1000)
#dstat

#bats
set.seed(1)
dstat_bat=phylo.d(data=adata_bat,binvar=sampled,permut=1000)
dstat_bat

#rods
set.seed(1)
dstat_rod=phylo.d(data=adata_rod,binvar=sampled,permut=1000)
dstat_rod

## number of studies and number of bats
hist(log10(cdata_bat$data$Samples))
hist(log10(cdata_rod$data$Samples))

## transform
cdata_bat$data$lSamples=log10(cdata_bat$data$Samples)
cdata_rod$data$lSamples=log10(cdata_rod$data$Samples)

## range
range(cdata_bat$data$Samples)
range(cdata_rod$data$Samples)

## pagel's lambda
pmod_bat=pgls(lSamples~1,data=cdata_bat,lambda="ML") ## lambda = 0.317
pmod_rod=pgls(lSamples~1,data=cdata_rod,lambda="ML") ## lambda = 0.474

## summarize
summary(pmod_bat)
summary(pmod_rod)

#####bats

## GPF for study binary
set.seed(1)
study_pf_bat=gpf(Data=adata_bat$data,tree=adata_bat$phy,
             frmla.phylo=sampled~phylo,
             family=binomial,
             algorithm='phylo',nfactors=25,
             min.group.size = 3)

## summarize
HolmProcedure(study_pf_bat)
study_res_bat=pfsum(study_pf_bat)$results

## lower/greater
study_res_bat$check=ifelse(study_res_bat$clade>study_res_bat$other,"more","less")
#table(study_res$check)

## number of samples
set.seed(1)

cdata_bat$data$lnSamp=log1p(cdata_bat$data$Samples)
cdata_bat$data$lSamp=log10(cdata_bat$data$Samples)

nsamples_pf_bat2=gpf(Data=cdata_bat$data,tree=cdata_bat$phy,
                frmla.phylo=lSamp~phylo,
                family=gaussian,
                algorithm='phylo',nfactors=25,
                min.group.size = 3)

## summarize
HolmProcedure(nsamples_pf_bat2)
nsamples_res_bat2=pfsum(nsamples_pf_bat2)$results

## lower/greater
nsamples_res_bat2$check=ifelse(nsamples_res_bat2$clade>nsamples_res_bat2$other,"more","less")
table(nsamples_res_bat2$check)

dtree_bat=treeio::full_join(as.treedata(adata_bat$phy),adata_bat$data,by="label")
stree_bat=treeio::full_join(as.treedata(cdata_bat$phy),cdata_bat$data,by="label")


#####rods

## GPF for sample binary
set.seed(1)
study_pf_rod=gpf(Data=adata_rod$data,tree=adata_rod$phy,
                 frmla.phylo=sampled~phylo,
                 family=binomial,
                 algorithm='phylo',nfactors=25,
                 min.group.size = 3)

## summarize
HolmProcedure(study_pf_rod)
study_res_rod=pfsum(study_pf_rod)$results

## lower/greater
study_res_rod$check=ifelse(study_res_rod$clade>study_res_rod$other,"more","less")
table(study_res_rod$check)

## number of samples
set.seed(1)

cdata_rod$data$lnSamp=log1p(cdata_rod$data$Samples)
cdata_rod$data$lSamp=log10(cdata_rod$data$Samples)

nsamples_pf_rod2=gpf(Data=cdata_rod$data,tree=cdata_rod$phy,
                    frmla.phylo=lSamp~phylo,
                    family=gaussian,
                    algorithm='phylo',nfactors=25,
                    min.group.size = 3)

## summarize
HolmProcedure(nsamples_pf_rod2)
nsamples_res_rod2=pfsum(nsamples_pf_rod2)$results

## lower/greater
nsamples_res_rod2$check=ifelse(nsamples_res_rod2$clade>nsamples_res_rod2$other,"more","less")
table(nsamples_res_rod2$check)

dtree_rod=treeio::full_join(as.treedata(adata_rod$phy),adata_rod$data,by="label")
stree_rod=treeio::full_join(as.treedata(cdata_rod$phy),cdata_rod$data,by="label")


###############FIG!

##########FIGURE

## set x max
plus=1
pplus=plus+0.75

## fix palette
AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
afun=function(x){
  a=AlberColours[1:x]
  return(a)
}

## make low and high
pcols=afun(2)

## function to loop and add clades
cadd=function(gg,pf,pmax,label="yes"){
  
  ## ifelse
  if(HolmProcedure(pf)==0){
    gg=gg
  }else{
    
    ## make result
    result=pfsum(pf)$results
    
    ## ifelse 
    if(nrow(result)>pmax){
      result=result[1:pmax,]
    }else{
      result=result
    }
    
    ## set tree
    for(i in 1:nrow(result)){
      
      ## highlight clade
      gg=gg+
        geom_hilight(node=result$node[i],
                     alpha=0.25,
                     fill=ifelse(result$clade>
                                   result$other,pcols[2],pcols[1])[i])+
        
        ## add label
        if(label=="yes"){
          geom_cladelabel(node = result$node[i], 
                          label = result$factor[i], 
                          offset = 1,
                          offset.text = 1.5,
                          fontsize=4)
        }    
      
    }
  }
  return(gg)
}

## state pmax
pmax=10

###trees

## make base
base=ggtree(dtree_bat,size=0.05,branch.length='none',layout="circular")
base2=ggtree(stree_bat,size=0.2,branch.length='none',layout="circular")
base3=ggtree(stree_bat,size=0.1)

#batbinary

dtree_bat@extraInfo$sampled <- as.factor(dtree_bat@extraInfo$sampled)
batbin=cadd(base,study_pf_bat,pmax,label="yes") +
  ggtitle("(A) binary sampled bat species") +
  geom_tippoint(aes(color=dtree_bat@extraInfo$sampled), size=0.3) +
  scale_color_manual(values=c("white","black")) +
  theme(legend.position = "none")

#number of batsamples
batsamples=cadd(base2,nsamples_pf_bat,pmax,label="yes") +
  ggtitle("(B) samples per bat species")

base_rod=ggtree(dtree_rod,size=0.05,branch.length='none',layout="circular")
base2_rod=ggtree(stree_rod,size=0.2,branch.length='none',layout="circular")
base3=ggtree(stree_bat,size=0.1)

dtree_rod@extraInfo$sampled <- as.factor(dtree_rod@extraInfo$sampled)
rodbin=cadd(base_rod,study_pf_rod,pmax,label="yes") +
  ggtitle("(C) binary sampled rodent species") +
  geom_tippoint(aes(color=dtree_rod@extraInfo$sampled), size=0.3) +
  scale_color_manual(values=c("white","black")) +
  theme(legend.position = "none")

rodsamples=cadd(base2_rod,nsamples_pf_rod,pmax,label="yes") +
  ggtitle("(D) samples per rodent species")

#number of rod samples
rodsamples=cadd(base2_rod,nsamples_pf_rod,pmax,label="no") +
  geom_cladelabel(node = nsamples_res_rod$node[1], 
                  label = nsamples_res_rod$factor[1], 
                  offset = 1,
                  offset.text = 1,
                  fontsize=4) +
  geom_cladelabel(node = nsamples_res_rod$node[2], 
                  label = nsamples_res_rod$factor[2], 
                  offset = 1,
                  offset.text = 1,
                  fontsize=4) +
  geom_cladelabel(node = nsamples_res_rod$node[3], 
                  label = nsamples_res_rod$factor[3], 
                  offset = 1,
                  offset.text = 4,
                  vjust = 1,
                  fontsize=4) +
  geom_cladelabel(node = nsamples_res_rod$node[4], 
                  label = nsamples_res_rod$factor[4], 
                  offset = 1,
                  offset.text = 1,
                  fontsize=4) +
  geom_cladelabel(node = nsamples_res_rod$node[5], 
                  label = nsamples_res_rod$factor[5], 
                  offset = 1,
                  offset.text = 0.7,
                  vjust = 1,
                  fontsize=4) +
  geom_cladelabel(node = nsamples_res_rod$node[6], 
                  label = nsamples_res_rod$factor[6], 
                  offset = 1,
                  offset.text = 1,
                  fontsize=4) +
  ggtitle("(D) samples per rodent species")


library(patchwork)

#export
setwd("~/Documents/PhD/colleen MS")
png("trees2",width=10,height=12,units="in",res=600)
(batbin|batsamples)/(rodbin|rodsamples)+plot_layout(widths=c(2))
dev.off()


####appendix B

nrow(TissueSamples[which(TissueSamples$specificEpithet == "crenulatum"),])
table((TissueSamples[which(TissueSamples$specificEpithet == "crenulatum"),])$institutionCode)


libr_used()
libr_unused()
