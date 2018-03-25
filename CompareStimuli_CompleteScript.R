
# rm(list = ls())
setwd("/home/caterina/Documents/PhD/PhD_EXPERIMENT_DATA_SETS/CompareStimuli")
options(scipen = 999)



# PACKAGES...............................................................------------------

# Pseudocode for installing necessary packages:
# install.packages( c( "pkg1", "pkg2", as required ) )


# Get this very specfic version of mclust, to preserve results as generated previously:
# install.packages("https://cran.r-project.org/src/contrib/Archive/mclust/mclust_4.3.tar.gz", repos = NULL, type = "source")
# But to be able to use the MclustBootstrap() function, need to upgrade to version 5.2.3

# ! Careful about optmatch however:
# Get this very specfic version of optmatch, on which further analyses depend:
# install.packages("https://cran.r-project.org/src/contrib/Archive/optmatch/optmatch_0.9-6.tar.gz", repos = NULL, type = "source")
# After version 0.9-6 or optmatch,  `pairmatch` and `fullmatch` no longer generate "matched.distances" attributes for their results. To get this information, use `matched.distances`.


# Data viz
library(lattice)
library(ggplot2)
library(GGally)
library(devtools)
# devtools::install_github("ropensci/plotly", force = T) # Most recent version
library(plotly)

# Data manipulation
library(descr)
library(data.table)
library(plyr)
library(dplyr)
library(GenABEL)
library(car)
library(tidyr)
library(stringr)
library(ForeCA)
library(gdata)

# intraclass correlations
library(ICC)

# Case-treatement maching
library(optmatch)

# Normality
library(MVN)
library(mvnormtest)

# Cluster analysis packages
library(GMD)
library(clValid)
library(amap)
library(vegan)
library(vcd)
library(flexclust)
library(mclust)
library(cluster)
library(clue)

# Export result tables
library(stargazer)
library(xtable)
library(texreg)

sessionInfo() # Check that correct versions were loaded.

load("~/Documents/PhD/PhD_EXPERIMENT_DATA_SETS/CompareStimuli/CompareStimuli.RData")



# FUNCTIONS ............................................................-----------

is.outlier = function (x) {
  # See: Davies, P.L. and Gather, U. (1993). "The identification of multiple outliers" (with discussion). J. Amer. Statist. Assoc., 88, 782-801.
  x <- na.omit(x)
  lims <- median(x) + c(-1, 1) * 2.5 * mad(x)
  x < lims[1] | x > lims[2]
}



CIwidth <- function(mean, sd, N){
  se <- sd / sqrt(N);
  upper <- mean + 1.96 * se # Get upper limit for 95% CI, with z=1.96. Or a 99% one with 2.575
  lower <- mean - 1.96 * se # Get lower limit for 95% CI.
  return(upper - lower)
}


cv <- function(ave, std){
  std/ave * 100; 
}


grpdist <- function(X)
{
  require(cluster)
  veg <- as.data.frame(as.factor(X))
  distgr <- daisy(veg,"gower")
  distgr
} 






# DATA MANIPULATION............................................... --------




# First to clean up the individual data(sounds, images, words) that will compose the overall data.


# 1. IAPS -----------------------------------------------------------------

all <- read.table("IAPSAll.txt", header=T, sep="\t") 

# a. Duplicates ----------------------------------------------------------

table(duplicated(all))
table(duplicated(all$IAPS)) # Iaps codes re-occur for 12 cases.

all$ID <- as.numeric(as.factor(all$IAPS)) # This identifies duplicates and gives them the same ID.

doubles <- which(cbind(table(all$ID)) == 2)
duplicates <- all[all$ID %in% doubles, ]

single_record <- !duplicated(duplicates[ , c("desc", "IAPS")])
uniques1 <- duplicates[single_record, c("desc", "IAPS")]

uniques2 <- matrix(ncol=10, nrow=0)
for (i in seq(0, 22, by=2)){
  uniques2 <- rbind(uniques2, colMeans(duplicates[ c(1, 2) + i, 3:12]) )
}

uniques <- data.frame(uniques1, uniques2)

# And now to substitute in real data:
all.unique <- all[! all$ID %in% doubles, ]
all.unique <- rbind(all.unique, uniques)
all.unique <- all.unique[order(all.unique$ID), ]

summary(as.matrix(all.unique[, 3:10]))

IAPS <- data.frame(all.unique)


# b. Dominance missing ---------------------------------------------------

length(which(IAPS$dom1mn == 0)) # 240
IAPS <- subset(IAPS, dom1mn != 0, select = -c(dom2mn, dom2sd, set))


# c. Outliers ------------------------------------------------------------

histkdnc(IAPS$valmn) # non-normal
histkdnc(IAPS$aromn)
histkdnc(IAPS$dom1mn) # yup, non-normal.

table(is.outlier(IAPS$valmn))
table(is.outlier(IAPS$aromn))
table(is.outlier(IAPS$dom1mn)) # 32 here.

IAPS <- IAPS[- which(is.outlier(IAPS$dom1mn)), ]
dim(IAPS)



#  d. 95% CI -------------------------------------------------------------


length(which(CIwidth(IAPS$valmn, IAPS$valsd, 100) > 1)) # 7
length(which(CIwidth(IAPS$aromn, IAPS$arosd, 100) > 1)) # 46
length(which(CIwidth(IAPS$dom1mn, IAPS$dom1sd, 100) > 1)) # 25... but 61 unique, overall.

univar.imprecise <- unique(c(which(CIwidth(IAPS$valmn, IAPS$valsd, 100) > 1),
                             which(CIwidth(IAPS$aromn, IAPS$arosd, 100) > 1),
                             which(CIwidth(IAPS$dom1mn, IAPS$dom1sd, 100) > 1)))

IAPS <- IAPS[ - univar.imprecise, ]

dim(IAPS)



# e. CV --------------------------------------------

length(which(cv(IAPS$valmn, IAPS$valsd) > 30)) # not good, 406.
length(which(cv(IAPS$aromn, IAPS$arosd) > 30)) # 817
length(which(cv(IAPS$dom1mn, IAPS$dom1sd) > 30)) # 730.

length(which(cv(IAPS$dom1mn, IAPS$dom1sd) < 30 & 
               cv(IAPS$aromn, IAPS$arosd) < 30 & 
               cv(IAPS$valmn, IAPS$valsd) < 30)) # Just 1!


# Consistency in chatacter strings:
IAPS$desc <- tolower( as.character( IAPS$desc ) )



with(IAPS, plot( valmn, aromn ) )



# 2. ANEW-----------------------------------------------------------------

ANEW <- read.csv("ANEWall1999.csv", header=T)

m <- as.matrix(ANEW)
m[m=="."] <- NA

ANEW <- data.frame(m)
ANEW[, 2:9] <- apply(ANEW[ , -1], 2, as.numeric)

summary(ANEW)

ANEW <- na.exclude(ANEW)
dim(ANEW)


# a. Word frequency ------------------------------------------------------

cut_points <- seq(0.05, 0.95, by = 0.05)
quantile(ANEW$Word.Frequency, cut_points)
histkdnc(ANEW$Word.Frequency) # Doesn't look like I can extract a group of words that is large enough which also has fairly constant freq... Shall just add it as a covariate in subsequent analyses...



# b. Outliers ------------------------------------------------------------

histkdnc(ANEW$Valence.Mean) # not normal...
histkdnc(ANEW$Arousal.Mean)
histkdnc(ANEW$Dominance.Mean)


table(is.outlier(ANEW$Valence.Mean))
table(is.outlier(ANEW$Arousal.Mean)) # just 3.
table(is.outlier(ANEW$Dominance.Mean)) # just 4.

which(is.outlier(ANEW$Arousal.Mean))
which(is.outlier(ANEW$Dominance.Mean)) # No overlap!

ANEW <- ANEW[- which(is.outlier(ANEW$Arousal.Mean)), ]
ANEW <- ANEW[- which(is.outlier(ANEW$Dominance.Mean)), ]




# c. 95% CI --------------------------------------------------------------


length(which(CIwidth(ANEW$Valence.Mean, ANEW$Valence.SD, 100) > 1)) + #9
  length(which(CIwidth(ANEW$Arousal.Mean, ANEW$Arousal.SD, 100) > 1)) + # wow... 269
  length(which(CIwidth(ANEW$Dominance.Mean, ANEW$Dominance.SD, 100) > 1)) # 89

#just checking
length(unique(ANEW$Word.No.)) == nrow(ANEW) # yup, no duplicates.

univar.imprec <- unique(c(which(CIwidth(ANEW$Valence.Mean, ANEW$Valence.SD, 100) > 1),
                          which(CIwidth(ANEW$Arousal.Mean, ANEW$Arousal.SD, 100) > 1),
                          which(CIwidth(ANEW$Dominance.Mean, ANEW$Dominance.SD, 100) > 1)))

ANEW <- ANEW[- univar.imprec, ]

dim(ANEW)


# d. CV ------------------------------------------------------------------



length(which(cv(ANEW$Dominance.Mean, ANEW$Dominance.SD) < 30 & 
               cv(ANEW$Arousal.Mean, ANEW$Arousal.SD) < 30 & 
               cv(ANEW$Valence.Mean, ANEW$Valence.SD) < 30)) # Just 6!

# Consistency in chatacter strings:
ANEW$Description <- tolower( as.character( ANEW$Description ) )


with(ANEW, plot( valmn, aromn ) )


# 3. IADS -----------------------------------------------------------------

IADS <- read.table("IADS2all.txt", header=T)
dim(IADS)
length(unique(IADS$Number)) #no duplicates
summary(IADS)


# a. Outliers ------------------------------------------------------------

histkdnc(IADS$PlMN) # not normal...
histkdnc(IADS$AroMN)
histkdnc(IADS$DomMN)


table(is.outlier(IADS$PlMN))
table(is.outlier(IADS$AroMN)) 
table(is.outlier(IADS$DomMN)) # No outliers...


# b. 95% CI --------------------------------------------------------------



length(which(CIwidth(IADS$PlMN, IADS$PlSD, 100) > 1)) #2
length(which(CIwidth(IADS$AroMN, IADS$AroSD, 100) > 1)) #1
length(which(CIwidth(IADS$DomMN, IADS$DomSD, 100) > 1)) 

IADS <- IADS[- which(CIwidth(IADS$PlMN, IADS$PlSD, 100) > 1), ]
IADS <- IADS[- which(CIwidth(IADS$AroMN, IADS$AroSD, 100) > 1), ]



# c. Coefficient of variation -------------------------------------------


length(which(cv(IADS$PlMN, IADS$PlSD) <= 30))
length(which(cv(IADS$AroMN, IADS$AroSD) <= 30))
length(which(cv(IADS$DomMN, IADS$DomSD) <= 30)) # Too conservative. Again cannot use this criterion.

length(which(cv(IADS$DomMN, IADS$DomSD) < 30 & 
               cv(IADS$AroMN, IADS$AroSD) < 30 & 
               cv(IADS$PlMN, IADS$PlSD) < 30)) # Just 1 again!



# Consistency in chatacter strings:
IADS$Sound <- tolower( as.character( IADS$Sound ) )



# FINAL CHECKS
summary(IAPS)
summary(ANEW)
summary(IADS)


with(IADS, plot( valmn, aromn ) )



# 4. ANET --------------------------------------------------------------------

ANET <- read.csv("ANETall.csv", sep="\t", header=T)
head( ANET[ , - ncol(ANET) ] )
summary(ANET[, 2:8])
dim(ANET)



# 4a. Outliers ------------------------------------------------------------

histkdnc(ANET$PlMN) # not normal...
histkdnc(ANET$AroMN)
histkdnc(ANET$DomMN)


table(is.outlier(ANET$PlMN))
table(is.outlier(ANET$AroMN)) # 7 cases...
table(is.outlier(ANET$DomMN)) # No outliers...Hmm.

ANET <- ANET[- which(is.outlier(ANET$AroMN)), ]



# 4b. 95% CI --------------------------------------------------------------



length(which(CIwidth(ANET$PlMN, ANET$PlSD, 100) > 1))
length(which(CIwidth(ANET$AroMN, ANET$AroSD, 100) > 1)) # whops... 2!
length(which(CIwidth(ANET$DomMN, ANET$DomSD, 100) > 1)) # whops...1!

ANET <- ANET[-which(CIwidth(ANET$DomMN, ANET$DomSD, 100) > 1),]
ANET <- ANET[-which(CIwidth(ANET$AroMN, ANET$AroSD, 100) > 1),]

# write.csv(ANET, file="OptmatchTexts.csv")


with(ANET, plot( Valence, Arousal ) )



# Binding all modalities --------------------------------------------------------

IAPS <- IAPS[ , -9] # Remove ID
# To keep!
Word.freq <- ANEW[, c(1, 9)] # 9 = Word frequency.
ANEW <- ANEW[, -9]


common_names <-  c("desc", "code", "valmn", "valsd", "aromn", "arosd", "dommn", "domsd")
colnames(IAPS) <- common_names
colnames(IADS) <- common_names
colnames(ANEW) <- common_names

IAPS$type <- "image"
ANEW$type <- "word"
IADS$type <- "sound"

data.frame(dim(IAPS), dim(IADS), dim(ANEW))

AVT <- rbind(IAPS, IADS, ANEW)
head(AVT)
table(AVT$type)

AVT$ID <- paste(AVT$type, AVT$code, AVT$desc, sep = "_")

# write.csv(AVT, file="AVT.csv")

pdf("/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/Graphics/Study2_MatchingStimuli/AVTQuadraticTrends.pdf", height = 7, width = 10)
ggplot( AVT, aes(x = valmn, y = aromn, color = type, pch = type)) +
  geom_point( size = 3.5, alpha = 0.6, stroke = 1.5 ) +
  scale_color_manual( values = c( "image" = "#7fb685",
                                  "word" = "#ef6f6c",
                                  "sound" = "#426a5a")) +
  scale_shape_manual(values = c( 12, 1, 8 ) ) +
  scale_x_continuous( breaks = seq(1, 9, by = 1) ) +
  scale_y_continuous( breaks = seq(1, 9, by = 1) ) +
  labs(x = "Valence", y = "Arousal", color = "Modality", pch = "Modality") +
  ggtitle("Valence-Arousal quadratic trend across modalities") +
  theme_gray( base_size = 17 ) 
dev.off()



#  MATCHING ............................................................... --------




table(AVT$type)

AVT_long <- data.table::melt(AVT, 
                             id.vars = c("ID", "type"), 
                             measure.vars = c("valmn", "aromn", "dommn"))

AVT_long$variable <- mapvalues(AVT_long$variable, 
                               from = c("valmn", "aromn", "dommn"), 
                               to = c("Valence", "Arousal", "Dominance"))

pdf("/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/Graphics/Study2_MatchingStimuli/PADDistribs_BeforeMatch.pdf", 
    width = 9.5, 
    height = 4)
ggplot(AVT_long, aes( x = type, y = value ) ) + 
  geom_point(position = "jitter", alpha = 0.3) + 
  geom_boxplot(alpha = 0.6) + 
  facet_wrap(~ variable) +
  ggtitle("PAD distributions across stimulus types, before matching") +
  xlab("Stimulus type") + ylab("Value")
dev.off()
# interesting. not very similar overall, before matching.





# 1. Sounds with images -----------------------------------------------------



SI <- subset(AVT, type == "sound" | type == "image")
SI$type <- factor(SI$type)

pm1 <- pairmatch(type ~ valmn + aromn + dommn, 
                 controls = 1, 
                 remove.unmatchables = T, 
                 data = SI)
summary(pm1)
stratumStructure(pm1) # perfect match for all sounds, with 0:1 matches for leftover images.

print(pm1, grouped = T, quote = T)


SI$group <- pm1
SI <- SI[with(SI, order(group, type)), ]

# Examples:
SI[which(rownames(SI) %in% c(507, 81)), c("valmn", "aromn", "dommn","ID", "group")]
SI[which(rownames(SI) %in% c(14, 418)), c("valmn", "aromn", "dommn", "ID", "group")]

SI <- na.exclude(SI) # Getting rid of extra images.
table(SI$type) # 164 each.



# 2. Sounds with words ------------------------------------------------------

# Will not work: 
# pm2 <- pairmatch(type ~ valmn + aromn + dommn, remove.unmatchables = T, data = SW)
# Error in pairmatch.matrix(m, controls = controls, data = mfd, remove.unmatchables = remove.unmatchables,: not enough controls in some subclasses
# From ?pairmatch: "In this case matching can still fail, if there is too much competition for certain controls; if you find yourself in that situation you should consider full matching, which necessarily finds a match for everyone with an eligible match somewhere."


SW <- subset(AVT, type == "sound" | type == "word")
SW$type <- factor(SW$type)

fm2 <- fullmatch(type ~ valmn + aromn + dommn, data = SW)
summary(fm2)

sS_fm2 <- stratumStructure(fm2)
no_words <- unlist( str_split( names( sS_fm2 ), ":") )[ seq( 1, length(sS_fm2) * 2, by = 2 ) ]
no_sounds <- unlist( str_split( names( sS_fm2 ), ":") )[ seq( 0, length(sS_fm2) * 2, by = 2 ) ]
sS_fm2 <- data.frame( no_words = as.numeric( no_words ), 
                      no_sounds = as.numeric( no_sounds ), 
                      sS_fm2)

# Reverse sign for odd groups:
sS_fm2[ no_sounds > no_words, "Freq"] <- - sS_fm2[ no_sounds > no_words, "Freq"]
sS_fm2$MatchType <- ifelse(sS_fm2$Freq > 0, 
                           "multiple words : 1 sound", 
                           "1 word : multiple sounds")
sS_fm2$MatchType <- ifelse(sS_fm2$no_words == sS_fm2$no_sounds, "1 word : 1 sound", sS_fm2$MatchType)
sS_fm2$MatchType <- factor(sS_fm2$MatchType , levels = c("multiple words : 1 sound",
                                                         "1 word : 1 sound",
                                                         "1 word : multiple sounds"))

setDT(sS_fm2)
sS_fm2 <- sS_fm2[ order(MatchType, Freq), ]

pdf("/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/Graphics/Study2_MatchingStimuli/FullMatch_StratumStructure.pdf", width = 11, height = 6)
ggplot(sS_fm2, aes(x = Var1, y = Freq, fill = MatchType)) + 
  geom_bar(stat = "identity", position = "identity") +
  scale_x_discrete(limits = sS_fm2$Var1) + # To reorder bars
  ylim( -5, 81) +
  scale_fill_manual(values = c("multiple words : 1 sound" = "#67a9cf",     
                               "1 word : 1 sound" = "#2166ac",
                               "1 word : multiple sounds" = "#b2182b")) +
  geom_text(aes ( x = Var1, y = Freq, label = Var1 ), size = 4,
            hjust = ifelse( sign( sS_fm2$Freq ) > 0, -0.3, 1.3),
            angle = 90 ) +
  labs(x = 'Stimulus matching ratios', y = 'Frequency') +
  ggtitle("Full match results: sound - word groups") +
  theme_grey( base_size = 17 ) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank() ) 
dev.off()


# Add column to SW to code stimulus groupings:
SW$group <- fm2

# Examples:
SW[which(rownames(SW) %in% c(7310, 10010)), c("valmn", "aromn", "dommn","ID", "group")]
SW[which(rownames(SW) %in% c(265, 9401)), c("valmn", "aromn", "dommn","ID", "group")]


# How many groups?
length(table(SW$group)) # 149, as expected. One for each sound matched to one or more words.

# Retrieve distance measures from optmatch structure, for later merge with original data:

# all_distances <- unlist(attr(fm2, 'matched.distances')) # For optmatch versions < 0.9-6
all_distances <- unlist( matched.distances( fm2,
                                            distance = match_on( type ~ valmn + aromn + dommn,
                                                                 data = SW ) ) )
all_distances_frame <- data.frame(UniqueName = as.character(names(all_distances)),
                                  Distance = all_distances,
                                  row.names = NULL,
                                  stringsAsFactors = F)
# So, the distance frame refers just to the 'controls' found for the 'treatments', but the 'treatments' themselves are not listed there, hence we will have to add them back in, and specify a 0 distance for them, since they are the reference point for which a match was found:

head(all_distances_frame)
tail(all_distances_frame)

# Get a unique name to match the format from all_distances_frame:
SW$UniqueName <- paste(SW$group, # Name of matching group
                       rownames(SW), 
                       sep=".")

# Whatever unique names are NOT in the distance frame, were used as reference points, and therefore should get a 0 distance:
SW <- join( SW, all_distances_frame, by = "UniqueName")
SW[ is.na(SW$Distance), "Distance"] <- 0

head(SW)
# Order by ascending distance, so there will be a 0 distance for the treatment, followed by rows of controls, in increasing distance. So, 0 distances belong only to items *for which* pairs were found.

SW <- SW[with(SW, order(fm2, Distance)), ]

head(SW)
head( table(SW$Distance, SW$type) ) # Not sure why the 0s aren't reserved for sounds only.


# Now to grab the best matched word-sound pairs, from the sorted database:
half1 <- which(SW$Distance == 0)
half2 <- half1 + 1

SW <- SW[ c(rbind( half1, half2 ) ), ]
SW <- SW[ with( SW, order(group, type) ), ]

# Some checks:
table(SI$type)
table(SW$type)

length(subset(SW, type == "sound")$code)
length(subset(SI, type == "sound")$code)

table( subset( SI, type == "sound" )$code %in% subset( SW, type == "sound" )$code ) # All ok.

# So: all 164 sounds were matched to images, however only 149 of the sounds found a matching word. 
# So we will further reduce the sample size to 149, to find stimulus trios.



# 3. All 3 modalities together --------------------------------------------


# Now to extract the common sounds from SI and SW, and their respective pairs.

common_sounds <- subset(SI, type == "sound")$code %in% subset(SW, type == "sound")$code 
table(common_sounds)
length(common_sounds)

SI_pairs <- split(SI, as.factor(SI$group))
length(SI_pairs)

relevant_SI_pairs <- SI_pairs[ common_sounds ]
relevant_SI_pairs <- rbindlist(relevant_SI_pairs)

# Check I have the right sounds present in both SI and SW:
SW <- data.table(SW)
bools <- sort( relevant_SI_pairs[ type == "sound", code] ) == sort( SW[ type == "sound", code ] )
table(bools) # Yes, they are.

# Now to give same group codes to SW - SI matches, based on the common sound in them:
relevant_SI_pairs[ , group := as.character(group)]
table( table( relevant_SI_pairs$group ) ) # All ok.




SW[ , UniqueName := NULL]; SW[ , Distance := NULL]; 




all_sounds <- SW[ type == "sound", ID]

trios <- list()
for ( sound in all_sounds ) {
  sound_row <- SW[ ID == sound, ]
  
  group_from_SW <- SW[ ID == sound, ]$group
  word_row <- subset( SW, type == "word" & group == group_from_SW )
  
  group_from_SI <- relevant_SI_pairs[ ID == sound, ]$group
  image_row <- subset( relevant_SI_pairs, type == "image" & group == group_from_SI )
  
  trios[[sound]] <- rbind(sound_row, word_row, image_row)
  
  trios[[sound]]$NewGroup <- which(all_sounds %in% sound)
  
}

table( unlist( lapply(trios, nrow) ) ) # All looks ok.



SIW <- rbindlist(trios)
# SIW[ , desc := tolower(desc)]
SIW[ , ID := NULL]
SIW[ , group := NULL]


# For export:
setnames(SIW, c( "Description", "StimulusDatabaseCode", "Valence", "ValenceSD", "Arousal", "ArousalSD", "Dominance", "DominanceSD", "Type", "Trio"))
setcolorder(SIW, c("Trio", "Description", "StimulusDatabaseCode", "Type", "Valence", "Arousal", "Dominance", "ValenceSD", "ArousalSD", "DominanceSD"))
# capture.output( stargazer( SIW, summary = F ), file = "/home/caterina/Desktop/MyTable.txt")

# write.csv(SIW, "/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/CSVTables/SIW_TrioData.csv")



# Should find ns differences:
summary(aov(Valence ~ Type, data = SIW))
summary(aov(Arousal ~ Type, data = SIW))
summary(aov(Dominance ~ Type, data = SIW))


SIW[ , ID := paste(StimulusDatabaseCode, Description, sep = "_")]
SIW_long <- data.table::melt(SIW, 
                             id.vars = c("ID", "Type"), 
                             measure.vars = c("Valence", "Arousal", "Dominance"))

pdf("/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/Graphics/Study2_MatchingStimuli/PADDistribs_AfterMatch.pdf", 
    width = 9.5, 
    height = 4)
ggplot(SIW_long, aes( x = Type, y = value ) ) + 
  geom_point(position = "jitter", alpha = 0.3) + 
  geom_boxplot(alpha = 0.6) + 
  facet_wrap(~ variable) +
  ggtitle("PAD distributions across stimulus types, after matching") +
  xlab("Stimulus type") + ylab("Value")
dev.off()


cloud(Valence ~ Arousal * Dominance, data = SIW, 
      pch = SIW$Type, col = SIW$Type, type = "p", 
      perspective = T, distance = 0.4, shade = T, lex = 2, 
      xlab = list("Arousal", cex = 1.7), 
      ylab = list("Dominance", cex = 1.7), 
      zlab = list("Valence", cex = 1.7), 
      main = list("Stimulus groups", cex = 2), 
      R.mat = diag(4), screen = list(z = -40, x = -42), aspect = c(1, 1), 
      key = list(text = list( levels( as.factor( SIW$Type ) ), cex = 1.2),
               title = "Stimulus type:", cex.title = 1.3,
               points = list( pch = 1:3, col = 1:nlevels( as.factor( SIW$Type ) ), lwd = 3), 
               columns = 1, 
               lines = list( col = 1:nlevels( as.factor( SIW$Type ) ), lwd = 3) )
)

pdf("/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/Graphics/Study2_MatchingStimuli/ParCoord_BeforeClustering.pdf", width = 12.5, height = 6.5)
# ggparcoord(SIW, columns = c(5, 6, 7), groupColumn = 4) + 
#   scale_colour_manual(values = c("image" = "#ca0020", 
#                                  "word" = "#92c5de", 
#                                  "sound" = "#0571b0"))
ggparcoord(SIW, columns = c(5, 6, 7), alpha = 0.35 ) +
    ggtitle("Parallel coordinates plot of overall PAD relations") +
    xlab("PAD dimension") + ylab("Scaled values") +
    theme_grey( base_size = 17 ) 
dev.off()


# write.csv(SIW, file="/home/caterina/Desktop/NewMatchedstimuli.csv")




# **Compare old matches to current** --------------------------------------



match <- fread("Matchedstimuli.csv", header=T)

length(SIW[ Type == "sound", StimulusDatabaseCode]) == length(match[ type == "sound", code] )

# Words and sounds identical:
table( SIW[ Type == "word", StimulusDatabaseCode] %in%  match[ type == "word", code] )
table( SIW[ Type == "sound", StimulusDatabaseCode] %in%  match[ type == "sound", code] )

# Something going on with one image:
table( SIW[ Type == "image", StimulusDatabaseCode] %in%  match[ type == "image", code] )
table( match[ type == "image", code] %in% SIW[ Type == "image", StimulusDatabaseCode] )

x <- sort(SIW[ Type == "image", paste(Description, StimulusDatabaseCode, sep = "_")]) # unique here: "execution_9414"
y <- sort(match[ type == "image", paste( tolower( desc ), code, sep = "_")]) # unique here: hurtdog_9183

which( table( c( x, y ) ) < 2 )
# Execution 9414 2.060 6.490 3.110 image
# HurtDog   9183 1.690 6.580 2.960 image

# Maybe to prove that it is no biggie that the two images were substituted, we can get a measure of distance between the new image that was matched to its sound + word, vs the old one:
current_duo <- SIW[ Trio == 52 & Type != "image",
                    .SD, .SDcol = c("Valence", "Arousal", "Dominance") ]

correct_image <- SIW[ Trio == 52 & Type == "image",
                      .SD, .SDcol = c("Valence", "Arousal", "Dominance") ]

used_image <- IAPS[ IAPS$code == 9183, c("valmn", "aromn", "dommn") ]

names(used_image) <- c("Valence", "Arousal", "Dominance")


mean( distEuclidean( as.matrix(current_duo),  
                     as.matrix(used_image) ) ) - mean( distEuclidean( as.matrix(current_duo), 
                                                                      as.matrix(correct_image) ) )

# So, pretty much the same.



# MCLUST CLASSIFICATION....................................................-------------------




# Classification itself ---------------------------------------------------


forclus <- SIW[ , .SD, .SDcol = c( 'Valence', 'Arousal', 'Dominance', 'Type', 'ID' ) ]

# For later use:
forclus_3cols <- data.frame( forclus[ , .SD, .SDcol = c( "Valence", "Arousal", "Dominance" ) ] )

fit <- Mclust( forclus_3cols )

pdf("/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/Graphics/Study2_MatchingStimuli/MBC_ClassificationPlot.pdf", width = 7, height = 7)
plot(fit, what = "classification") # plot results
dev.off()

pdf("/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/Graphics/Study2_MatchingStimuli/MBC_UncertaintyPlot.pdf", width = 7, height = 7)
plot(fit, what = "uncertainty") # plot results
dev.off()


summary(fit) # display the best model 
# Mclust EEE (elliposidal, equal volume, shape and orientation) model with 4 components
names(fit)

BICs <- as.matrix( data.frame( fit$BIC ) )

data.frame( sort( BICs, index.return = T ) )
arrayInd( which.max(BICs), .dim = dim(BICs), useNames = T )
arrayInd( 57, .dim = dim(BICs), useNames = T )
arrayInd( which.min(BICs), .dim = dim(BICs), useNames = T )

centroids <- data.frame( t( fit$parameters$mean ), 
            N = as.vector( table( fit$classification ) ),
            `Mixing proportion` = fit$parameters$pro )
xtable(centroids)




# Multivariate normality - actual tests - for the 4 clusters? ------------------------------



C1 <- subset( forclus[ , .SD, .SDcol = 1:3], fit$classification == 1 )
C2 <- subset( forclus[ , .SD, .SDcol = 1:3], fit$classification == 2 )
C3 <- subset( forclus[ , .SD, .SDcol = 1:3], fit$classification == 3 )
C4 <- subset( forclus[ , .SD, .SDcol = 1:3], fit$classification == 4 )

cluster_list <- list( C1, C2, C3, C4 )

mardia_normality <- list()
shapiro_normality <- list()
hz_normality <- list()
royston_normality <- list()
pca_normality <- list()

for ( cluster in 1:length( cluster_list ) ){
  print(cluster)
  mardia_normality[[cluster]] <- mardiaTest( cluster_list[[cluster]], qqplot = F )
  shapiro_normality[[cluster]] <- mshapiro.test( t( cluster_list[[cluster]] ) )   
  hz_normality[[cluster]] <- hzTest( cluster_list[[cluster]], qqplot = F )  
  royston_normality[[cluster]] <- roystonTest( cluster_list[[cluster]], qqplot = F )  
  pca_normality[[cluster]] <- pca( cluster_list[[cluster]] ) 
}

mardia_skew <- unlist( lapply( mardia_normality, function(x) x@p.value.skew ) )
mardia_kurt <- unlist( lapply( mardia_normality, function(x) x@p.value.kurt ) )
shapiro_p <- unlist( lapply( shapiro_normality, function(x) x$p.value ) )
hz_p <- unlist( lapply( hz_normality, function(x) x@p.value ) )
royston_p <- unlist( lapply( royston_normality, function(x) x@p.value ) )
eigen_vals <- matrix( unlist( lapply( pca_normality, function(x) x$eig ) ), nrow = 3, byrow = F)

mv_normality <- round( rbind( mardia_skew, mardia_kurt, shapiro_p, hz_p, royston_p, eigen_vals ), 3)
mv_normality[1:4, ] > 0.05

# Visually inspect normality of clusters, by creating a contour plot:
plot( fit, what = "density" )

contours <- data.frame( forclus[ , .SD, .SDcol = 1:3], 
                        clus = as.factor( paste("Cluster",
                                                fit$classification) ) )

ggplot(contours, aes( x = Valence, y = Arousal ) ) +
  geom_density_2d( aes( group = clus, 
                        color = clus), 
                   binwidth = 0.05 ) +
  geom_point(aes( group = clus, 
                  color = clus), alpha = 0.5) +
  scale_colour_manual(values = c("Cluster 1" = "#2166ac",
                                 "Cluster 2" = "#67a9cf",
                                 "Cluster 3" = "#ef8a62",
                                 "Cluster 4" = "#b2182b")) 

ggplot(contours, aes( x = Valence, y = Dominance ) ) +
  geom_density_2d( aes( group = clus, 
                        color = clus), 
                   binwidth = 0.05 ) +
  geom_point(aes( group = clus, 
                  color = clus), alpha = 0.5) +
  scale_colour_manual(values = c("Cluster 1" = "#2166ac",
                                 "Cluster 2" = "#67a9cf",
                                 "Cluster 3" = "#ef8a62",
                                 "Cluster 4" = "#b2182b")) 

ggplot(contours, aes( x = Arousal, y = Dominance ) ) +
  geom_density_2d( aes( group = clus, 
                        color = clus), 
                   binwidth = 0.05 ) +
  geom_point(aes( group = clus, 
                  color = clus), alpha = 0.5) +
  scale_colour_manual(values = c("Cluster 1" = "#2166ac",
                                 "Cluster 2" = "#67a9cf",
                                 "Cluster 3" = "#ef8a62",
                                 "Cluster 4" = "#b2182b")) 


# What if there is a better way to investigate the influence of non-normality on the results?
# -> Using bootstrap:
# Need a more recent version of mclust for this, however. e.g., v. 5.2.3


# Check for differences using BOOTSTRAPPED params (only safe for sure, for mclust version 5.2.3):
fit_bs <- MclustBootstrap(fit, nboot = 1000, type = "bs")

mean_bs_cis <- summary(fit_bs, what = "ci")$mean
mean_bs_ses <- summary(fit_bs, what = "se")$mean

mean_bs_cis <- data.frame(rbind(mean_bs_cis[ , , 1], 
                                mean_bs_cis[ , , 2],
                                mean_bs_cis[ , , 3],
                                mean_bs_cis[ , , 4]) )

# Reshape the data:
mean_bs_cis$Cluster <- rep( 1:4, each = 2)
mean_bs_cis$Bound <- rep( c("lower_bound_bootstrap_CI", "upper_bound_bootstrap_CI"), 4) 
mean_bs_cis <- melt(mean_bs_cis, measure.vars = names(mean_bs_cis)[1:3] )
mean_bs_cis <- dcast(mean_bs_cis, Cluster + variable ~ Bound, 
      value.var = "value") 


original_params <- data.frame( t( fit$parameters$mean ) )
original_params <- stack(original_params ) 
original_params$Cluster <- rep(1:4, 3)
names(original_params) <- c("estimate", "variable", "Cluster")

mean_bs_cis <- setDT( join(mean_bs_cis, original_params, by = c("Cluster", "variable") ) )
setcolorder(mean_bs_cis, c("Cluster", "variable", "lower_bound_bootstrap_CI", "estimate", "upper_bound_bootstrap_CI") )

# Adding SEs too:
SEs <- data.frame( stack( data.frame( t( mean_bs_ses ) ) ), 
                   Cluster = rep(1:4, 3) )
names(SEs) <- c("error", "variable", "Cluster")
mean_bs_cis <- join(mean_bs_cis, SEs, by = c("Cluster", "variable") )



# Caterpillar plot
mean_bs_cis[ , location := paste( Cluster, variable, sep = "_") ]
mean_bs_cis[ , location := ordered( location, levels = rev( unique( location ) ) ) ]

# pdf("/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/Graphics/Study2_MatchingStimuli/95BootstrappedCIsAndErrors.pdf", width = 11, height = 8)
ggplot(mean_bs_cis, aes(location, estimate, 
                        ymin = lower_bound_bootstrap_CI, 
                        ymax = upper_bound_bootstrap_CI,
                        color = as.factor(Cluster),
                        size = error ) ) + 
  scale_size_continuous( range = c(0.4, 1) ) +
  scale_color_manual(values = c("1" = "#b2182b",     
                                "2" = "#ef8a62",
                                "3" = "#000000",
                                "4" = "#67a9cf") ) + 
  geom_pointrange( show.legend = F ) + coord_flip() +
  labs( y = "Estimate", x = "Cluster and dimension", color = "Cluster") +
  ggtitle( "Bootstrap 95% confidence intervals around the estimated MBC centroids") +
  scale_y_continuous( breaks = seq(1, 9, by = 1), limits = c(1, 9) ) +
  theme_grey( base_size = 17 )
# dev.off()



# PAD corrs - an issue? ---------------------------------------------------


# Correlations:
forclus_wht <- whiten( as.matrix( forclus[ , .SD, .SDcol = c( 'Valence', 'Arousal', 'Dominance') ] ) )$U
ggpairs(forclus_wht)

Mclust(forclus_wht) # Still components, whohoo!
plot( Mclust( forclus_wht ), what = "classification" )
plot( Mclust( forclus_wht ), what = "uncertainty" )




# Select best reps --------------------------------------------------------



forclus$clus <- fit$classification
forclus$uncert <- fit$uncertainty
head(forclus)
forclus <- forclus[ with( forclus, order( clus, uncert ) ), ]

forclus$Cluster <- paste("Cluster", forclus$clus, sep = " ")


# Get first 30 cases from each sorted-by-uncertainty cluster:
best <- forclus %>% 
  group_by(clus) %>% 
  slice(1:30)


xtable( t( table( best$clus, 
                  best$Type ) ) ) # 5 is the minimum number. Therefore, shall get 5 cases/cluster:



best <- best %>% 
  group_by( clus, Type, add = TRUE ) %>% 
  slice(1:5)

table(best$clus, best$Type) # All OK.
xtable( data.frame(best), digits = 3 )



# plot this relative to previous parallel coordinates plot:

best_for_plot <- best %>% 
  group_by( clus, Type, add = TRUE ) %>% 
  slice(1:1)


mylabels <- paste(forclus$Type, forclus$ID, sep = "_")
marking_best_reps <- forclus$ID %in% best_for_plot$ID
marking_best_reps <- ifelse( marking_best_reps == FALSE, NA, marking_best_reps)
mylabels_bestreps <- mylabels[ marking_best_reps ]

isC1C2 <- forclus$clus == 1 | forclus$clus == 2
isC1C2 <- ifelse( isC1C2 == FALSE, NA, isC1C2 )

isC3C4 <- forclus$clus == 3 | forclus$clus == 4
isC3C4 <- ifelse( isC3C4 == FALSE, NA, isC3C4 )

first_axis <- mylabels_bestreps[ isC3C4 ]
middle_axis <- rep("", nrow(forclus) )
last_axis <- mylabels_bestreps[ isC1C2 ]

mylabels <- c(first_axis, middle_axis, last_axis)
myoffset <- rep( c( 1.1, 0, -0.1 ), each = nrow( forclus ) )

# # Careful about running pdf() here. Jitter is random, so might overwrite a plot with nice (but random) label placement:
# pdf("/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/Graphics/Study2_MatchingStimuli/ParCoord_AfterClustering.pdf", width = 12.5, height = 6.5)
ggparcoord(forclus, columns = c(1, 2, 3), groupColumn = 8) +
  geom_line(size = 0.25) + geom_text( aes(label = mylabels, hjust = myoffset ),
                                      nudge_y = jitter( rep( 0, length(myoffset) ),
                                                        factor = 20 ) ) +
  ggtitle("Parallel coordinates plot of PAD relations across clusters") +
  xlab("PAD dimension") + ylab("Scaled values") +
  scale_colour_manual(values = c("Cluster 1" = "#2166ac",
                                 "Cluster 2" = "#67a9cf",
                                 "Cluster 3" = "#ef8a62",
                                 "Cluster 4" = "#b2182b")) +
  theme_grey( base_size = 17 )
# dev.off()



# ADD ARTIFICIAL CLUSTER .............................................................--------


# Add neutral cluster, by first removing stimuli we have already decided to use:
SIW_remainder <- SIW[! SIW$ID %in% best$ID, ]

# So now to extract neutrals from "not.used":
## theoretical neutral = 5 on the 9-point Likert scale
theoretical_model <- as.matrix( data.frame( Valence = 5, Arousal = 5, Dominance = 5 ) )

diffs <- distEuclidean( as.matrix( SIW_remainder[ , .SD, .SDcol = c( "Valence", "Arousal", "Dominance") ] ), 
                        theoretical_model )

SIW_remainder[ , theoretical_diffs := diffs]
SIW_remainder <- SIW_remainder[ order( Type, theoretical_diffs ), ]

artificial_cluster <- SIW_remainder %>% 
  group_by( Type ) %>% 
  slice(1:5)


# Add artificial cluster to 'best' stimuli

best <- data.table(best)
final_75 <- rbind(best[ , .SD, .SDcol = c("Valence", "Arousal", "Dominance", "Type", "ID", "clus") ], 
                  SIW[ SIW$ID %in% artificial_cluster$ID, .SD, .SDcol = c( "Valence", "Arousal", "Dominance", "Type", "ID" ) ],
                  fill = TRUE)
final_75$clus <- ifelse(is.na(final_75$clus), 5, final_75$clus)
final_75 <- final_75[ order(clus, Type), ]


xtable(final_75)


# plot:
final_75$Cluster <- paste("Cluster", final_75$clus)

pdf("/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/Graphics/Study2_MatchingStimuli/2DPlot5Clusters.pdf", width = 7, height = 5)
ggplot(final_75, aes( x = Valence, y = Arousal, pch = substr(Type, 1, 1), color = Cluster ) ) +
  geom_text( aes( label = substr(Type, 1, 1 ) ), size = 6 ) +
  scale_color_manual(values = c("Cluster 1" = "#b2182b",     
                               "Cluster 2" = "#ef8a62",
                               "Cluster 3" = "#000000",
                               "Cluster 4" = "#67a9cf",
                               "Cluster 5" = "#2166ac")) +
  ggtitle("4 data-driven and 1 artificial cluster") +
  theme_grey( base_size = 17 ) 
dev.off()


# Small validation:
ICCest(clus, Valence * Arousal * Dominance, data = final_75, alpha = 0.05, CI.type = "THD")
ICCest(clus, Valence, data = final_75, alpha = 0.05, CI.type = "THD")
# Still fine!!


# lm model comparing all clusters:
final_75$Cluster <- relevel( as.factor( final_75$Cluster ), "Cluster 5")
PAD_model_with_C5 <- lm( I( Valence * Arousal * Dominance ) ~ Cluster + Type, data = final_75 ) # Excellent
summary(PAD_model_with_C5)
texreg(PAD_model_with_C5, single.row = T)




# Difference with closest cluster:
t_v <- t.test( Valence ~ Cluster,  
               data = subset(final_75, 
                             Cluster == "Cluster 5" | Cluster == "Cluster 3") )
t_a <- t.test( Arousal ~ Cluster,  
               data = subset(final_75,
                             Cluster == "Cluster 5" | Cluster == "Cluster 3") )
t_d <- t.test( Dominance ~ Cluster,  
               data = subset(final_75, 
                             Cluster == "Cluster 5" | Cluster == "Cluster 3") )
t_pad <- t.test( Valence * Arousal * Dominance ~ Cluster,  
                 data = subset(final_75, 
                               Cluster == "Cluster 5" | Cluster == "Cluster 3") )

xtable( t_out( t_v,  n.equal = TRUE, 
              welch.df.exact = TRUE, welch.n = NA,
              d.corr = FALSE, print = FALSE) )
xtable( t_out( t_a,  n.equal = TRUE, 
               welch.df.exact = TRUE, welch.n = NA,
               d.corr = FALSE, print = FALSE) )
xtable( t_out( t_d,  n.equal = TRUE, 
               welch.df.exact = TRUE, welch.n = NA,
               d.corr = FALSE, print = FALSE) )
xtable( t_out( t_pad,  n.equal = TRUE, 
               welch.df.exact = TRUE, welch.n = NA,
               d.corr = FALSE, print = FALSE) )



# get group means:

cluster_by_type_means1 <- setDT( aggregate( cbind(Valence, Arousal, Dominance) ~ Type + Cluster, 
                                    FUN = mean, data = final_75) )
cluster_by_type_means2 <- setDT( aggregate( cbind(Valence, Arousal, Dominance) ~ Cluster, 
                                     FUN = mean, data = final_75) )
setnames(cluster_by_type_means2, c( "Cluster", "GrandValence", "GrandArousal", "GrandDominance") )

cluster_by_type_means <- join(cluster_by_type_means1, cluster_by_type_means2, by = "Cluster")
setcolorder(cluster_by_type_means, c( "Cluster", "Type", 
                                      "Valence", "GrandValence", 
                                      "Arousal", "GrandArousal", 
                                      "Dominance",  "GrandDominance" ) )

xtable(cluster_by_type_means)





# 3D plot:
cols <- mapvalues(final_75$Cluster,
          from = unique(final_75$Cluster),
          to = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e") )
# pdf("/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/Graphics/Study2_MatchingStimuli/3DPlot5Clusters.pdf", width = 6, height = 6)
cloud(Dominance ~ Valence + Arousal, 
      # zoom = 1,
      data = final_75, 
      pch = final_75$Type, 
      col = cols, 
      cex = 1,
      lex = 3.5,
      screen = list( z = 140, x = -42 ),
      main = "The multi-modal best representatives from the 5 clusters",
      scales = list(arrows = FALSE),
      key = list( points = list(pch = c( rep(19, 5), NA, 1, 2, 3), lwd = 3.5,
                                col = c("#1b9e77", "#d95f02", "#7570b3", 
                                        "#e7298a", "#66a61e", "white", # Add blank between legend categories
                                        rep( "black", 3 ) ) ),
                  text = list( c( unique( final_75$Cluster ), "", unique( as.character( final_75$Type ) ) ) ),
                  space = 'right',
                  column = 1 )
)
# dev.off()



# # **Compare with old classification** 
# 
# old_selection <- read.csv("/home/caterina/Documents/PhD/PhD_Experiment_Archives/RandomizeStimsOriginal.csv")
# 
# old_selection_1_4 <- subset(old_selection, Type != "film" & Type != "VE" & Cluster != 5)
# 
# old_selection_1_4 <- sort( paste( old_selection_1_4$Type,
#                         old_selection_1_4$DatabaseCode,
#                         tolower(old_selection_1_4$Descr), 
#                         sep = "_" ) )
# 
# current_selection <- sort( paste(best$Type,
#                            best$ID,
#                            sep = "_") )
# 
# # Remove already-known difference:
# old_selection_1_4 <- old_selection_1_4[! old_selection_1_4 == "image_9183_hurtdog"]
# current_selection <- current_selection[! current_selection == "image_9414_execution"]
# 
# table( old_selection_1_4 == current_selection ) # All OK.
# 
# 
# # Artificial cluster
# 
# old_selection_5 <- subset(old_selection, Type != "film" & Type != "VE" & Cluster == 5)
# 
# old_selection_5 <- sort( paste( old_selection_5$Type,
#                                 old_selection_5$DatabaseCode,
#                                 tolower(old_selection_5$Descr), 
#                                 sep = "_" ) )
# artificial_selection <- sort( paste(artificial_cluster$Type,
#                                     artificial_cluster$ID,
#                                     sep = "_") )
# 
# table( old_selection_5 == artificial_selection ) # All OK.


# MBC validation ----------------------------------------------------------


half1 <- forclus_3cols[sample( rownames( forclus_3cols ), 
                               as.integer( nrow(forclus_3cols) / 2 ),
                               replace = FALSE), ]

half2 <- forclus_3cols[ setdiff( rownames( forclus_3cols ), 
                                 rownames( half1 ) ), ]

# Must make half1 and half 2 exactly equal in dim
to_exclude <- sample( rownames( half2 ), 1 )
half2 <- forclus_3cols[ setdiff( rownames(half2), to_exclude ), ]

MBC_h1 <- Mclust(half1, G = 4)
MBC_h2 <- Mclust(half2, G = 4)


# Predicting cluster structure in one half based on the other half of the data:

# h1 -> h2

pred_h2_from_h1 <- cl_predict(MBC_h1, 
                              newdata = half2, 
                              type = "class_ids")

plot( half2, col = pred_h2_from_h1 ) # Looks very reasonable

pred_h2_from_h1_cramer <- assocstats( table( pred_h2_from_h1,
                                             MBC_h2$classification)  )$cramer


# h2 -> h1

pred_h1_from_h2 <- cl_predict(MBC_h2, 
                              newdata = half1, 
                              type = "class_ids")

plot( half1, col = pred_h1_from_h2 ) # Looks very reasonable

pred_h1_from_h2_cramer <- assocstats( table( pred_h1_from_h2, 
                                             MBC_h1$classification ) )$cramer

mean( c( pred_h2_from_h1_cramer, pred_h1_from_h2_cramer) )



# Overlap

# For h2 -> h1
`rand_h2 -> h1` <- clusteval::cluster_similarity(MBC_h1$classification, 
                                                 cl_predict(MBC_h2, newdata = half1), 
                                                 similarity = "rand")

clusteval::comembership_table(MBC_h1$classification, 
                              cl_predict(MBC_h2, newdata = half1))

`adjusted_rand_h2 -> h1` <- mclust::adjustedRandIndex(MBC_h1$classification, 
                                                      cl_predict(MBC_h2, newdata = half1))


# For h1 -> h2
`rand_h1 -> h2` <- clusteval::cluster_similarity(MBC_h2$classification, 
                                                 cl_predict(MBC_h1, newdata = half2), 
                                                 similarity = "rand")

clusteval::comembership_table(MBC_h2$classification, 
                              cl_predict(MBC_h1, newdata = half2))

`adjusted_rand_h1 -> h2` <-mclust::adjustedRandIndex(MBC_h2$classification, 
                                                     cl_predict(MBC_h1, newdata = half2))


mean( c(`rand_h2 -> h1`, `rand_h1 -> h2`) )
mean( c(`adjusted_rand_h2 -> h1`, `adjusted_rand_h1 -> h2`) )


# Variation of Information

`vi_h2 -> h1` <- igraph::compare(MBC_h1$classification, 
                                 cl_predict(MBC_h2, newdata = half1), 
                                 method = c("vi")) 

`vi_h1 -> h2` <- igraph::compare(MBC_h2$classification, 
                                 cl_predict(MBC_h1, newdata = half2), 
                                 method = c("vi")) 


# max = log(nclust)*2 
`vi_h2 -> h1` / ( log(4) * 2 ) # OR:
`normalised_vi_h2 -> h1` <- `vi_h2 -> h1` / log( nrow( half1 ) ) 

`vi_h1 -> h2` / ( log(4) * 2 ) # OR
`normalised_vi_h1 -> h2` <- `vi_h1 -> h2`  / log( nrow( forclus ) ) 

mean( c(`normalised_vi_h2 -> h1`, `normalised_vi_h1 -> h2`) )


# Removing 10% and reassessing number of clusters

ideal_G <- vector()

for (rep in 1:1000) {

  print(rep)

  cases_to_exclude <- sample( rownames( forclus_3cols ),
                              as.integer(10 * nrow(forclus_3cols) / 100 ),
                              replace = F )
  test_data <- forclus_3cols[ ! rownames(forclus_3cols) %in% cases_to_exclude, ]

  ideal_G <- append(ideal_G,
                    Mclust(test_data)$G )

}

prop.table( table( ideal_G) )



# intraclass correlations

icc_val <- unlist( ICCest(clus, Valence, data = forclus, alpha = 0.05, CI.type = "THD") )
# Smith gives upper limits over 1 which is impossible.
icc_aro <- unlist( ICCest(clus, Arousal, data = forclus, alpha = 0.05, CI.type = "THD") )
icc_dom <- unlist( ICCest(clus, Dominance, data = forclus, alpha = 0.05, CI.type = "THD") )
icc_pad <- unlist( ICCest(clus, Valence * Arousal * Dominance, data = forclus, alpha = 0.05, CI.type = "THD") )

icc_SIW_4clusters <- rbind( icc_val, icc_aro, icc_dom, icc_pad )
icc_SIW_4clusters <- round( icc_SIW_4clusters, 3 )

xtable(icc_SIW_4clusters)






# ADD ANET TEXTS TO FINAL SELECTION OF STIMULI.............................-----------------------------

ANET$ID <- paste( ANET$ANET, ANET$Sentence, sep = "_")
ANET$clus <- NA
ANET$Type <- "text"
ANET <- data.table(ANET)

setnames(ANET, c( "ANET", "Valence", "ValenceSD", "Arousal", "ArousalSD", "Dominance", 
                  "DominanceSD", "Set", "Sentence", "ID", "clus", "Type" ) )

final_75_plus_ANET <- rbind( final_75, ANET[ , .SD, .SDcol = names(final_75) ] )
final_75_plus_ANET[ , MatchingStatus := ifelse( Type == "text", "ExplanatoryText", "ActualStimulus" ) ]

# Full match without the max.controls constraint:
final_fm <- fullmatch( MatchingStatus ~ Valence + Arousal + Dominance, 
                       # max.controls = 1,
                       data = final_75_plus_ANET )
summary(final_fm)
stratumStructure(final_fm)

final_75_plus_ANET[, sentence_groups := final_fm ]
table( final_75_plus_ANET$sentence_groups, 
       final_75_plus_ANET$Type )



# Add matching distances back into the data:
text_distances <- unlist( matched.distances( final_fm,
                                             distance = match_on( MatchingStatus ~ Valence + Arousal + Dominance, 
                                                                  # max.controls = 1,
                                                                  data = final_75_plus_ANET ) ) )
text_distances_frame <- data.frame(UniqueName = as.character(names(text_distances)),
                                   Distance = text_distances,
                                   row.names = NULL,
                                   stringsAsFactors = F)

# Get a unique name to match the format from all_distances_frame:
final_75_plus_ANET$UniqueName <- paste(final_75_plus_ANET$sentence_groups, # Name of matching group
                                       rownames(final_75_plus_ANET), 
                                       sep = ".")

# Whatever unique names are NOT in the distance frame, were used as reference points, and therefore should get a 0 distance:
final_75_plus_ANET <- join( final_75_plus_ANET, text_distances_frame, by = "UniqueName")
final_75_plus_ANET[ is.na(final_75_plus_ANET$Distance), "Distance"] <- 0

final_75_plus_ANET <- final_75_plus_ANET[ order(sentence_groups, MatchingStatus, Distance), ]

# Breaking down the types of matches into several types, handled differently:
# First, to deal with 'proper' matches, i.e., 1 actual stimulus taken as reference point:
final_75_plus_ANET$MatchType <- ifelse(final_75_plus_ANET$MatchingStatus == "ActualStimulus" & 
                                         final_75_plus_ANET$Distance == 0,
                                       "ProperMatchAnchor", 
                                       ifelse(final_75_plus_ANET$MatchingStatus == "ExplanatoryText" & 
                                                final_75_plus_ANET$Distance == 0,
                                              "InverseMatchAnchor", 
                                              "MatchesForAnchor"))

sentence_groups_list <- split(final_75_plus_ANET, final_75_plus_ANET$sentence_groups)
proper_matches_index <- unlist( lapply( sentence_groups_list, 
                                        function(x) any( which(x$MatchType == "ProperMatchAnchor") ) ) )
final_75_plus_ANET_proper_matches <- rbindlist( sentence_groups_list[ proper_matches_index ] )

half1 <- which(final_75_plus_ANET_proper_matches$Distance == 0)
half2 <- half1 + 1

final_75_plus_ANET_proper_matches <- final_75_plus_ANET_proper_matches[ c(rbind( half1, half2 ) ), ]


# And now to tackle second type of matches, i.e., 'improper' matches, where ANET text served as reference point:
inverse_matches_index <- unlist( lapply( sentence_groups_list, 
                                        function(x) ! any( which(x$MatchType == "ProperMatchAnchor") ) ) )

final_75_plus_ANET_inverse_matches <- rbindlist( sentence_groups_list[ inverse_matches_index ] )

# Improper match subtypes:
# Subtype 1: inversed pairs, where just 1 reference-point text was matched to 1 actual stimulus:
inversed_pairs <- names( which( table( final_75_plus_ANET_inverse_matches$sentence_groups ) == 2 ) )

final_75_plus_ANET_inverse_pairs <- final_75_plus_ANET_inverse_matches[ final_75_plus_ANET_inverse_matches$sentence_groups %in% inversed_pairs, ]


# Subtype 2: improper matches where 1 reference-point text was matched to MULTIPLE actual stimuli:
inversed_full_matches <- names( which( table( final_75_plus_ANET_inverse_matches$sentence_groups ) > 2 ) )

final_75_plus_ANET_inverse_full_matches <- final_75_plus_ANET_inverse_matches[ final_75_plus_ANET_inverse_matches$sentence_groups %in% inversed_full_matches, ]

repetition_no <- table( factor( final_75_plus_ANET_inverse_full_matches$sentence_groups) ) - 2 
# -2 because: 
# a. we do not want to include the Inverse Match Anchors themselves, just the other stimuli, and 
# b. we need one less when adding back to original data.

# Have to repeat ANET sentences, since they are shared by multiple actual stimuli:
inverse_anchors_to_duplicate <- which(final_75_plus_ANET_inverse_full_matches$MatchType == "InverseMatchAnchor")

final_75_plus_ANET_inverse_full_matches_expanded <- 
  final_75_plus_ANET_inverse_full_matches[ rep( inverse_anchors_to_duplicate, repetition_no), ]

full_inverse_matches <- rbind(final_75_plus_ANET_inverse_full_matches_expanded, final_75_plus_ANET_inverse_full_matches)


# Prepare for rbind and sorting of all 3 mini-datasets, by:

## One subtype: filling in NAs in the clus column, by using multiple lag values:
full_inverse_matches <- full_inverse_matches[ order(sentence_groups, MatchingStatus), ]

lagged_clusters <- vector()
for (gr in unique(full_inverse_matches$sentence_groups) ) {
  gr_subset <- subset(full_inverse_matches, sentence_groups == gr) 
  lagged_clusters <- append(lagged_clusters, lag(gr_subset$clus, ( nrow( gr_subset ) / 2 ) ) ) 
}

lagged_clusters[ is.na(lagged_clusters) ] <- 0
full_inverse_matches$clus[ is.na(full_inverse_matches$clus) ] <- 0
full_inverse_matches$clus <- lagged_clusters + full_inverse_matches$clus

## Another subtype: 
final_75_plus_ANET_inverse_pairs <- final_75_plus_ANET_inverse_pairs[ order(sentence_groups, MatchingStatus), ]
final_75_plus_ANET_inverse_pairs[ , clus := fill(final_75_plus_ANET_inverse_pairs, 
                                                 clus, 
                                                 .direction = "down")$clus ]

## And another subtype: 
final_75_plus_ANET_proper_matches <- final_75_plus_ANET_proper_matches[ order(sentence_groups, MatchingStatus), ]
final_75_plus_ANET_proper_matches[ , clus := fill(final_75_plus_ANET_proper_matches, 
                                                 clus, 
                                                 .direction = "down")$clus ]

# Bind all 3 together, and then sorting resulting dataset:
ANET_matches_to_75_stims <- rbind(final_75_plus_ANET_proper_matches, 
                                  final_75_plus_ANET_inverse_pairs, 
                                  full_inverse_matches)
ANET_matches_to_75_stims <- ANET_matches_to_75_stims[ order(clus, sentence_groups, MatchingStatus ), ]



# Now to interleave texts with actual stimuli, when those texts are shared with multiple stims in the same matching group:
ANET_matches_to_75_stims <- interleave(ANET_matches_to_75_stims[ANET_matches_to_75_stims$MatchingStatus=="ActualStimulus", ], 
                                       ANET_matches_to_75_stims[ANET_matches_to_75_stims$MatchingStatus=="ExplanatoryText", ])

ANET_matches_to_75_stims$sorting_aid <- ifelse( ANET_matches_to_75_stims$Type == "text", 
                                                   NA, 
                                                   as.character(ANET_matches_to_75_stims$Type))

ANET_matches_to_75_stims$sorting_aid <- fill(ANET_matches_to_75_stims, 
                                             sorting_aid, 
                                             .direction = "down")$sorting_aid

ANET_matches_to_75_stims <- ANET_matches_to_75_stims[ order(clus, sorting_aid), ]

table(ANET_matches_to_75_stims$clus, ANET_matches_to_75_stims$Type) # Looks ok!


# Some cleaning:
ANET_matches_to_75_stims$sentence_groups <- as.character( ANET_matches_to_75_stims$sentence_groups )
ANET_matches_to_75_stims$sentence_groups <- mapvalues(ANET_matches_to_75_stims$sentence_groups, 
                                                      unique( ANET_matches_to_75_stims$sentence_groups ), 
                                                      1:57 )
ANET_matches_to_75_stims[ , sorting_aid := NULL]
ANET_matches_to_75_stims[ , UniqueName := NULL]
ANET_matches_to_75_stims[ , Distance := NULL]
ANET_matches_to_75_stims[ , MatchType := NULL]
setcolorder( ANET_matches_to_75_stims,  c("sentence_groups", "clus", "ID", "Type", "MatchingStatus", "Valence", "Arousal", "Dominance"  ) )


# Export for APPENDIX:
ANET_matches_to_75_stims_for_export <- ANET_matches_to_75_stims
IDs <- unlist( str_split( ANET_matches_to_75_stims_for_export$ID, "_" ) )
Code <- IDs[ seq( 1, 300, by = 2 ) ]
Descr <- IDs[ seq( 0, 300, by = 2 ) ]
ANET_matches_to_75_stims_for_export$ID <- NULL
ANET_matches_to_75_stims_for_export[ , Code := Code]
ANET_matches_to_75_stims_for_export[ , Descr := Descr]
setcolorder( ANET_matches_to_75_stims_for_export, c( "sentence_groups", "Type", "MatchingStatus", "Code", "Descr", "Valence", "Arousal", "Dominance", "clus" ) )
ANET_matches_to_75_stims_for_export[ , Descr := str_replace(Descr, ";", ",")]
# write.csv2(ANET_matches_to_75_stims_for_export, file = "/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/CSVTables/Chapter3/ANETMatchesTo75Stims_Updated.csv", quote = F)
# write.table(ANET_matches_to_75_stims_for_export, file = "/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/CSVTables/Chapter3/ANETMatchesTo75Stims_Updated.txt", quote = F, sep = "\t", dec = ".")

# Verification
icc_val <- unlist( ICCest(clus, Valence, data = ANET_matches_to_75_stims, alpha = 0.05, CI.type = "THD") )
icc_aro <- unlist( ICCest(clus, Arousal, data = ANET_matches_to_75_stims, alpha = 0.05, CI.type = "THD") )
icc_dom <- unlist( ICCest(clus, Dominance, data = ANET_matches_to_75_stims, alpha = 0.05, CI.type = "THD") )
icc_pad <- unlist( ICCest(clus, Valence * Arousal * Dominance, data = ANET_matches_to_75_stims, alpha = 0.05, CI.type = "THD") )

icc_ANET_matches_5clusters <- rbind( icc_val, icc_aro, icc_dom, icc_pad )
icc_ANET_matches_5clusters <- round( icc_ANET_matches_5clusters, 3 )

xtable(icc_ANET_matches_5clusters)


ANET_matches_to_75_stims$clus <- paste( "Cluster", ANET_matches_to_75_stims$clus )
ANET_matches_to_75_stims_long <- data.table::melt(ANET_matches_to_75_stims, 
                             id.vars = c("ID", "Type", "clus", "MatchingStatus"), 
                             measure.vars = c("Valence", "Arousal", "Dominance"))


pdf("/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/Graphics/Study2_MatchingStimuli/BoxplotsAfterANETMatching.pdf", width = 12, height = 13)
ggplot(ANET_matches_to_75_stims_long, aes( x = Type, y = value ) ) + 
  geom_point(alpha = 0.3) + # geom_jitter(alpha = 0.3) + show.legend = F +
  geom_boxplot( aes(fill = MatchingStatus), alpha = 0.6, lwd = 0.4 ) + 
  facet_grid( clus ~ variable ) +
  scale_fill_manual(values = c("ActualStimulus" = "#ca0020",     
                               "ExplanatoryText" = "#0571b0")) +
  ggtitle("PAD distributions across stimulus types and clusters, after ANET matching") +
  xlab("Stimulus type") + ylab("PAD value") +
  scale_y_continuous( name = "PAD value", limits = c( 1, 9 ), breaks = seq( 0, 9, by = 2 ) ) +
  theme_grey( base_size = 17 )  
dev.off()


`3D` <- plot_ly(ANET_matches_to_75_stims, 
                type="scatter3d",
                mode = 'markers',
                x = ~Valence, y = ~Arousal, z = ~Dominance,
                symbol = ~Type,
                symbols = c("circle", "square", "diamond", "x"),
                color = ~clus,
                colors = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e'),
                text = ~ID,
                hoverinfo = 'text', textposition = "bottom right",
                marker = list(size = 5.5, opacity = 0.7) ) %>%
  layout(scene = list(xaxis = list(title = 'Valence' ),
                      yaxis = list(title = 'Arousal', range = c(2, 8) ),
                      zaxis = list(title = 'Dominance' ) ),
         legend = list(orientation = 'h',
                       # bgcolor = "#E2E2E2",
                       bordercolor = "#E2E2E2",
                       borderwidth = 2,
                       font = list(size = 12) ),
         title = 'PAD 3D space: clusters and stimulus types')

Sys.setenv("plotly_username" = "LexKonstantine")
Sys.setenv("plotly_api_key" = "PRFUP3jLgZWbfzdZdAZ6")
# plotly_IMAGE(`3D`, width = 1000, height = 500, format = "png", scale = 1, out_file = "plotly-test-image.png")
plotly_POST(`3D`, filename = "test.png")







# COMPARISON WITH MY IAPS PAPER...............................................-------------------------------

# Getting and editing old data associated with IAPS publication:
IAPS_paper_results <- fread("https://raw.githubusercontent.com/CaterinaC/IAPSClustering2016/master/IAPSClassification_Constantinescuetal2016.csv")
IAPS_paper_results[ , isIAPSOnlyCluster := TRUE ]
IAPS_paper_results[ , Description := tolower( Description ) ]
IAPS_paper_results[ , ImageCode := as.numeric( ImageCode ) ]
IAPS_paper_results[ , Under75thQuantile := NULL ]
setnames(IAPS_paper_results, "Cluster", "IAPSCluster")

# Get and annotate first 5 cases from each sorted-by-uncertainty cluster,
# and count these as the IAPS best reps:
IAPS_paper_results_annotated <- IAPS_paper_results %>%
  group_by(IAPSCluster) %>%
  slice(1:5)
setDT( IAPS_paper_results_annotated )[ , isIAPSOnlyClusterRep := TRUE]

# Import back to the full IAPS data, the info about which cases count as best reps:
IAPS_paper_results <- join(IAPS_paper_results, 
                           IAPS_paper_results_annotated[ , .SD, .SDcols = c( "ImageCode", 
                                                                             "Description", 
                                                                             "isIAPSOnlyClusterRep" )], 
                           by = c("ImageCode", "Description"))


# Getting and editing CURRENT data on image reps, after matching the modalities:
matched_image_reps <- final_75[ Type == "image", ] # Only images relevant for IAPS comparison
matched_image_reps[ , isMatchedClusterRep := TRUE ]
setnames(matched_image_reps, "clus", "ClusterOfMatchedReps")

ImageCode <- unlist( str_split(matched_image_reps$ID, "_") )[ seq( 1, 
                                                                   nrow(matched_image_reps) * 2, 
                                                                   by = 2 ) ]
matched_image_reps[ , ImageCode := as.numeric( ImageCode ) ]

Description <- unlist( str_split(matched_image_reps$ID, "_") )[ seq( 0,  
                                                                     nrow(matched_image_reps) * 2, 
                                                                     by = 2 ) ]
matched_image_reps[ , Description := Description]



IAPS_paper_and_matching_results <- join(IAPS_paper_results, 
                                        matched_image_reps[ , .SD, .SDcols = c("ImageCode",
                                                                               "Description", 
                                                                               "isMatchedClusterRep",
                                                                               "ClusterOfMatchedReps") ],
                                        by = c( "ImageCode", "Description") )


# Finally, getting info on how ALL the images were clustered, right after being matched:
matched_images <- forclus[ Type == "image", ]

ImageCode <- unlist( str_split(matched_images$ID, "_") )[ seq( 1, 
                                                               nrow(matched_images) * 2, 
                                                               by = 2 ) ]
Description <- unlist( str_split(matched_images$ID, "_") )[ seq( 0, 
                                                                 nrow(matched_images) * 2,
                                                                 by = 2 ) ]

matched_images[ , ImageCode := as.numeric( ImageCode ) ]
matched_images[ , Description := Description]

setnames(matched_images, "Cluster", "MatchedCluster")
matched_images[ , MatchedCluster := str_replace(MatchedCluster, "Cluster ", "") ]

IAPS_paper_and_matching_results <- join(IAPS_paper_and_matching_results, 
                                        matched_images[ , .SD, .SDcols = c("ImageCode",
                                                                           "Description", 
                                                                           "MatchedCluster") ],
                                        by = c( "ImageCode", "Description") )

names(IAPS_paper_and_matching_results)



# All 25 images that were used in my MATCHED data are accounted for:
with(IAPS_paper_and_matching_results, table(isMatchedClusterRep, isIAPSOnlyCluster) )


# Assessing overlap / co-membership for the 25 cases:
`25only` <- IAPS_paper_and_matching_results[ isMatchedClusterRep == TRUE, ]

with(`25only`, table(IAPSCluster, MatchedCluster) )

# clusteval::comembership_table(`25only`$IAPSCluster, 
#                               `25only`$MatchedCluster)

# clusteval::cluster_similarity(`25only`$IAPSCluster, 
#                               `25only`$MatchedCluster, 
#                               similarity = "rand") # This crashes R for whatever reason.

mclust::adjustedRandIndex(`25only`$IAPSCluster, 
                          `25only`$MatchedCluster) # Not ideal.

# How about overlap between all the IAPS-only clusters, and all the images from the matched clusters?
with(IAPS_paper_and_matching_results, table(MatchedCluster, IAPSCluster) )
mclust::adjustedRandIndex(IAPS_paper_and_matching_results$IAPSCluster, 
                          IAPS_paper_and_matching_results$MatchedCluster) 




# Where are the CURRENT best reps from the MATCHED dataset, relative to the IAPS clusters?

pdf("/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/Graphics/Study2_MatchingStimuli/ClusterDrift.pdf", height = 8, width = 12)
# Put all the info together in a layered plot :
ggplot() +
  # images from matched data:
  geom_point( data = IAPS_paper_and_matching_results[ ! is.na( MatchedCluster ), ], 
              aes( x = Valence, y = Arousal, color = as.factor( MatchedCluster ) ), 
              size = 4, alpha = 0.2, pch = 12, stroke = 0.9 ) +
  # IAPS-only:
  geom_point( data = IAPS_paper_and_matching_results[ isIAPSOnlyCluster == TRUE, ], 
              aes( x = Valence, y = Arousal, fill = as.factor( IAPSCluster ) ), 
              size = 4, alpha = 0.075, pch = 24 ) +
  ## Best reps layers:
  # Layer for matched data classification reps:
  geom_point(data = IAPS_paper_and_matching_results[ isMatchedClusterRep == TRUE, ], 
             aes( x = Valence, y = Arousal, color = as.factor( ClusterOfMatchedReps ) ), 
             size = 3, alpha = 1, pch = 12, stroke = 1.7 ) +
  scale_color_manual(values = c("1" = "#f92c2c",     
                                "2" = "#389bcd",
                                "3" = "#f08b00",
                                "4" = "#003663",
                                "5" = "#3c9e60")) + 
  # Another layer for IAPS-only classification reps:
  geom_point(data = IAPS_paper_and_matching_results[ isIAPSOnlyClusterRep == TRUE, ], 
             aes( x = Valence, y = Arousal, fill = as.factor(IAPSCluster) ), 
             size = 3.5,  alpha = 1, pch = 24 ) +
  scale_fill_manual(values = c( "1" = "#FFC300",     
                                "2" = "#FF5733",
                                "3" = "#C70039",
                                "4" = "#900C3F",
                                "5" = "#581845" ) ) + 
  ggtitle("IAPS cluster drift pre- and post-matching") +
  labs(fill = "Original IAPS\nclusters", col = "Image clusters\npost-matching") +
  theme_grey( base_size = 18 ) 
dev.off()





# COMPARISON WITH OTHER PAPERS...............................................--------------------------------



# Larsen, J. T., Norris, C. J., & Cacioppo, J. T. (2003). Effects of positive and negative affect on electromyographic activity over zygomaticus major and corrugator supercilii. Psychophysiology, 40(5), 776-785.

# "Sixty-six affective pictures, sounds, and words were chosen from standardized stimulus sets. [...] Stimuli were selected to fall into 11 valence categories (each containing six stimuli) that spanned the bipolar valence dimension  and  varied  widely  in  content. Pictures  were representative of most of the affective categories described by Bradley, Codispoti, Cuthbert, and Lang (2001) and included mutilated bodies, accidents, household objects, people, foods, sports, and so forth. "

# "In our sample, the most negatively rated stimuli were a picture of a victim of severe burns (IAPS #2800), the sound of a fight (IADS #290), and the word war (ANEW #482). The most positively rated stimuli were a picture of a trio of puppies (IAPS #1710), music from Beethoven (IADS #810), and the word sweetheart (ANEW #424)."

# "averages were constructed for stimuli rated as very negative (P = 0, N = 3–4), moderately negative (P = 0, N = 1–2), moderately positive (P = 1–2, N = 0), and very positive (P = 3–4, N = 0). For comparison, a fifth set of averages was constructed for stimuli rated as neutral (N = 0, P = 0)."


rm( list = ls() )
load("~/Documents/PhD/PhD_EXPERIMENT_DATA_SETS/CompareStimuli/CompareStimuli.RData")


IAPS_codes <- c("1052", "1121", "1274", "1660", "1710", "1750", "1930", "2050", "2057", "2208", "2340", "2661", "2800", "3100", "3130", "3160", "3250", "3261", "3400", "4000", "4004", "4230", "4302", "4520", "4534", "4571", "4598", "4641", "4656", "4680", "4770", "4800", "5300", "5731", "5910", "5920", "5971", "5972", "6150", "6260", "6561", "6610", "6930", "7002", "7006", "7010", "7050", "7217", "7430", "7460", "7600", "8060", "8080", "8185", "8210", "8260", "8311", "8370", "8460", "8501", "9102", "9500", "9520", "9560", "9570", "9810")

IADS_codes <- c("105", "106", "109", "110", "111", "112", "113", "116", "120", "130", "132", "133", "152", "171", "201", "202", "205", "206", "215", "216", "220", "221", "226", "230", "251", "252", "261", "262", "280", "290", "292", "310", "319", "320", "322", "325", "351", "353", "358", "360", "361", "380", "400", "403", "420", "422", "425", "500", "501", "502", "602", "625", "700", "701", "702", "706", "708", "709", "711", "722", "725", "730", "802", "810", "815", "826")

ANEW_codes <- c("005", "016", "069", "077", "112", "152", "167", "198", "206", "212", "241", "261", "301", "335", "337", "385", "391", "393", "394", "403", "423", "424", "433", "435", "437", "438", "447", "456", "472", "478", "482", "486", "498", "503", "530", "541", "549", "570", "571", "573", "638", "644", "664", "675", "677", "682", "683", "731", "734", "743", "757", "758", "759", "772", "777", "829", "845", "854", "878", "884", "890", "897", "908", "958", "964", "979")

all.unique$ID <- NULL
original_IAPS <- all.unique
# Replace Dominance 0 values with NA?
original_IAPS[ original_IAPS$dom1mn == 0, "dom1mn"] <- NA

# Also get ANEW and IADS
original_ANEW <- read.csv("ANEWall1999.csv", header=T)
original_IADS <- read.table("IADS2all.txt", header=T)


setDT(original_IAPS); setDT(original_ANEW); setDT(original_IADS)
setnames(original_ANEW, c("Descr", "Code", "Valence", "Valence.SD", "Arousal", "Arousal.SD",
                          "Dominance", "Dominance.SD", "Word.Freq") )
setnames(original_IADS, c("Descr", "Code", "Valence", "Valence.SD", "Arousal", "Arousal.SD",
                          "Dominance", "Dominance.SD") )
setnames(original_IAPS, c("Descr", "Code", "Valence", "Valence.SD", "Arousal", "Arousal.SD",
                          "Dominance", "Dominance.SD", "DomExtra", "DomExtra.SD", "Set") )


# Mark separately, for each modality, which stimuli were used in this paper:
# (Separately because the codes may not be unique between modalities)

# IAPS
original_IAPS[ , `Larsen et al., 2003` := ifelse( as.character(Code) %in% IAPS_codes, 
                                   "Larsen et al., 2003", 
                                   NA) ]
original_IAPS[ , Descr := tolower(Descr) ]
original_IAPS[ , ID := paste(Code, Descr, sep = "_") ]
original_IAPS[ , Type := "image" ]


# IADS
original_IADS[ , `Larsen et al., 2003` := ifelse( as.character(Code) %in% IADS_codes, 
                          "Larsen et al., 2003", 
                          NA) ]
original_IADS[ , Descr := tolower(Descr) ]
original_IADS[ , ID := paste(Code, Descr, sep = "_") ]
original_IADS[ , Type := "sound" ]


# ANEW
original_ANEW[ , `Larsen et al., 2003` := ifelse( as.character(Code) %in% ANEW_codes, 
                          "Larsen et al., 2003", 
                          NA) ]
original_ANEW[ , Descr := tolower(Descr) ]
original_ANEW[ , ID := paste(Code, Descr, sep = "_") ]
original_ANEW[ , Type := "word" ]


comparing_selection_methods <- rbind(original_IAPS, original_IADS, original_ANEW, fill = TRUE)
table(comparing_selection_methods$`Larsen et al., 2003`)
# I should have 75, and Larsen 66 * 3 = 198



# Adding info on which stimuli were sampled as part of our trios:
SIW[ , `Our method` := "Our matched trios"]
comparing_selection_methods <- join(comparing_selection_methods, 
                                    SIW[ , .SD, .SDcol = c("ID", "Our method")],
                                    by = "ID")

# Adding info on which stimuli were sampled as part of our best reps:
comparing_selection_methods <- join(comparing_selection_methods,
                                    final_75[ , .SD, .SDcol = c("ID", "clus")],
                                    by = "ID")

table(comparing_selection_methods$`Larsen et al., 2003`,
      comparing_selection_methods$Type)

table(comparing_selection_methods$`Our method`,
      comparing_selection_methods$Type)




# How to merge the info in the Author and Source columns ?
# Particular care for cases that were sampled by Larsen AND me at the same time, therefore, will be different to show independent boxplots without some hackery..
# So, what we have to do is duplicate the entries that both Larsen and WE have used, and thus allow use in different groups:

merging_aid <- data.table::melt(comparing_selection_methods, 
                 id.vars = c("ID"), 
                 measure.vars = c("Our method", "Larsen et al., 2003"))

merging_aid <- merging_aid[ , .SD, .SDcols = c("ID", "value")]
merging_aid <- na.exclude(merging_aid[ order(ID), ])
setnames(merging_aid, "value", "Source")


# Merge back into the previous dataset:
comparing_selection_methods <- join(comparing_selection_methods, 
                                    merging_aid,
                                    by = "ID")
# So the source column codes who sampled a stimulus, and deals with situations where a stimulus was sampled by multiple authors by duplicating that stimulus 
# Get a measure of how many stimuli were sampled by both me and Larsen:
table( duplicated( comparing_selection_methods$ID) ) # Considerable amount

comparing_selection_methods$Status <- ifelse( is.na(comparing_selection_methods$Source), 
                                              "Not sampled", 
                                              "Sampled" ) 
table( comparing_selection_methods$Status )

comparing_selection_methods$Source_complete <- ifelse( is.na( comparing_selection_methods$Source ),
                                                       "Not sampled", 
                                                       comparing_selection_methods$Source )
table( comparing_selection_methods$Source_complete, 
       comparing_selection_methods$clus, 
       useNA = "always" )



comparing_selection_methods[ , Source_complete2 := ifelse( Source_complete == "Our matched trios" & 
                                                             ! is.na( clus ),
                                                           as.character(clus), 
                                                           as.character(Source_complete) )]

comparing_selection_methods[ , Source_complete := relevel( factor( Source_complete ), "Not sampled" ) ]
comparing_selection_methods[ , Source_complete2 := relevel( factor( Source_complete2 ), "Not sampled" ) ]
comparing_selection_methods[ , alphas := ifelse( Source_complete2 == "Not sampled",
                                                 "transparent",
                                                 "opaque") ]

comparing_selection_methods$Source_complete2 <- mapvalues(comparing_selection_methods$Source_complete2,
                                                          from = c("Not sampled", 
                                                                   "Our matched trios",
                                                                   "Larsen et al., 2003", 
                                                                   "1", "2", "3", "4", "5"),
                                                          to =  c("Not sampled", 
                                                                  "Our matched trios", 
                                                                  "Larsen et al., 2003", 
                                                                  "Cluster 1", "Cluster 2", 
                                                                  "Cluster 3", "Cluster 4",
                                                                  "Cluster 5"))


# Some plots, based on all the parameters created above:
# Valence - Arousal
ggplot(comparing_selection_methods, 
       aes( x = Valence, y = Arousal, 
            pch = Type, 
            color = Source_complete2,
            alpha = alphas, 
            size = Source_complete2) ) + 
  geom_point() + 
  scale_color_manual(values = c("black",  
                                c("green", "orange", "hotpink", "yellow", "cyan"), 
                                "red", "blue" ) ) +
  scale_alpha_discrete(range = c(1, 0.1) ) +
  scale_size_discrete( range = c(4.5, 2) ) 


# Turning to 3D plotting in lattice now:
cols <- mapvalues(comparing_selection_methods$Source_complete2, 
          from = levels( comparing_selection_methods$Source_complete2 ),
          to = c("black",  c("blue", "purple", "hotpink", "yellow", "cyan"), 
                 "red", "green" ) ) 

comparing_selection_methods$colors <- with(comparing_selection_methods, 
                                           ifelse(Source_complete2 == "Not sampled", 
                                                  adjustcolor(cols, alpha.f = 0.05), 
                                                  adjustcolor(cols, alpha.f = 1) ) )

comparing_selection_methods$sizes <- with(comparing_selection_methods,
                                          ifelse(Source_complete2 == "Not sampled", 
                                                 2, 
                                                 0.5 ) )

# Going back in to increase sizes of actual cluster reps a bit more:
comparing_selection_methods$sizes <- with(comparing_selection_methods,
                                          ifelse(Source_complete2 == "Cluster 1" |
                                                 Source_complete2 == "Cluster 2" |
                                                 Source_complete2 == "Cluster 3" |
                                                 Source_complete2 == "Cluster 4" |
                                                 Source_complete2 == "Cluster 5",
                                                 1, 
                                                 comparing_selection_methods$sizes ) )

# General 3D plot, with everything in it:
cloud(Dominance ~ Valence + Arousal, 
      zoom = 1.1,
      data = comparing_selection_methods, 
      pch = 19, 
      col = comparing_selection_methods$colors, 
      cex = comparing_selection_methods$sizes,
      screen = list( z = 170, x = -60 ), # list( z = 140, x = -42 )
      key = list(points = list(pch = 19, 
                               col =  c("black",  "blue", "orange", "hotpink", 
                                        "yellow", "cyan", "red", "green" ) ), 
                 text = list( levels( comparing_selection_methods$Source_complete2 ) ), 
                 space = 'right', 
                 columns = 1 ),
      main = "Comparison with other authors")



##### Subplots:

stimulus_of_interest_name <- "image"
# stimulus_of_interest_name <- "word"
# stimulus_of_interest_name <- "sound"

stimulus_of_interest <- comparing_selection_methods[ comparing_selection_methods$Type == stimulus_of_interest_name, ]

# Editing graphical params from before:
## Going to stop highlighting our matched trios, since we're already going to make the point about Larsen's cases not being matches as well as ours in the ggplot2 boxplots below:
stimulus_of_interest$Source_complete2 <- ordered( stimulus_of_interest$Source_complete2,
                                                  levels = c("Not sampled", "Our matched trios", 
                                                             "Cluster 1", "Cluster 2", "Cluster 3", 
                                                             "Cluster 4", "Cluster 5" , 
                                                             "Larsen et al., 2003") )

# Set color mapping:
cols <- mapvalues(stimulus_of_interest$Source_complete2, 
                  from = levels( stimulus_of_interest$Source_complete2 ),
                  to = c("black", "black", "#e6ab02", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "black") ) 

# Making both non-sampled and trio data faded out:
stimulus_of_interest$colors <- with(stimulus_of_interest, 
                                    ifelse(Source_complete2 == "Not sampled" |
                                             Source_complete2 == "Our matched trios", 
                                           adjustcolor(cols, alpha.f = 0.025), 
                                           adjustcolor(cols, alpha.f = 1) ) )

# Set sizes:
stimulus_of_interest$sizes <- with(stimulus_of_interest,
                                   ifelse(Source_complete2 == "Not sampled"  |
                                            Source_complete2 == "Our matched trios" , 
                                          2, 
                                          0.5 ) )

stimulus_of_interest$sizes <- with(stimulus_of_interest,
                                   ifelse(Source_complete2 == "Cluster 1" |
                                            Source_complete2 == "Cluster 2" |
                                            Source_complete2 == "Cluster 3" |
                                            Source_complete2 == "Cluster 4" |
                                            Source_complete2 == "Cluster 5",
                                          1.5, 
                                          stimulus_of_interest$sizes ) )

# Set symbol mapping:
symbs <- mapvalues(stimulus_of_interest$Source_complete2,
                   from = levels( stimulus_of_interest$Source_complete2 ),
                   to = c(19, 19, 7, 9, 10, 12, 8, 20) ) 
symbs <- as.numeric( as.character( symbs ) )

pdf(paste("/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/Graphics/Study2_MatchingStimuli/3D_Larsen2003_vs_Clusters_",
          stimulus_of_interest_name, ".pdf", 
          sep = ""), 
    width = 10, height = 10)
cloud(Dominance ~ Valence + Arousal, 
      zoom = 0.85,  scales = list(arrows = FALSE, cex = 1.5 ),
      xlab = list( label = "Valence", cex = 1.5, rot = -5),
      ylab = list( label = "Arousal", cex = 1.5, rot = 80),
      zlab = list( label = "Dominance", cex = 1.5, rot = 95),
      data = stimulus_of_interest, 
      pch = symbs, 
      lex = 3,
      col = stimulus_of_interest$colors, 
      cex = stimulus_of_interest$sizes * 1.5,
      screen = list( z = 170, x = -60 ), # list( z = 140, x = -42 )
      key = list(points = list(pch = c(19, NA, 7, 9, 10, 12, 8, NA, 20),
                               alpha = c(0.2, 0, 1, 1, 1, 1, 1, 0, 1),
                               lwd = 2, cex = 2,
                               col =  c("black", "white", "#e6ab02", "#d95f02", "#7570b3",
                                        "#e7298a", "#66a61e", "white", "black" ) ),
                 text = list( c("Not sampled", "",
                                "Cluster 1", "Cluster 2", "Cluster 3",
                                "Cluster 4", "Cluster 5",
                                "", "Larsen et al., 2003"),
                              cex = 1.5 ),
                 space = 'right',
                 columns = 1 ),
      main = list( paste( tools::toTitleCase( stimulus_of_interest_name ), 
                   " sampling locations",
                   sep = ""),
                   cex = 1.5) )
dev.off()


# All together:
# library(latticeExtra)
# P1$main[[1]] <- "Stimulus sampling locations"
# P1$legend$right$args$key$space <- "bottom"
# P1$legend$bottom <- P1$legend$right
# P1$legend$right <- NULL
# P1$legend$bottom$args$key$columns <- 9
# P1$bottom$args$key$columns <- 9
# pdf("/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/Graphics/Study2_MatchingStimuli/3D_Larsen2003_vs_Clusters_IAPSANEWIADS_2.pdf", width = 18, height = 7)
# c(Images = P1, Words = P2, Sounds = P3)
# dev.off()





# Long format:
comparing_selection_methods_long <- data.table::melt(comparing_selection_methods, 
                                                     id.vars = c("ID", "Type", "Source"), 
                                                     measure.vars = c("Valence", "Arousal", "Dominance"))


## Boxplots showing that our matching procedure is superior to 'manual' methods:
# pdf("/home/caterina/Documents/PhD/PhDThesis/ThesisTEX/Graphics/Study2_MatchingStimuli/Larsen2003ANEWIAPSIADS.pdf", width = 15.5, height = 4.5)
ggplot( na.exclude( comparing_selection_methods_long ), 
        aes( y = value, fill = Type, x = Source ) ) + 
  geom_boxplot( position = position_dodge(width = 0.85) ) +
  facet_grid( ~ variable) +
  ggtitle("Comparing sampling strategies") +
  ylab("Values") +
  theme_grey( base_size = 18 ) +
  scale_fill_manual(values = c("image" = "#67a9cf",     
                               "sound" = "#ca0020",
                               "word" = "#f4a582"))
# dev.off()




plot_ly(comparing_selection_methods, 
        type="scatter3d",
        mode = 'markers',
        x = ~Valence, y = ~Arousal, z = ~Dominance,
        symbol = ~Type,
        symbols = c("circle", "square", "diamond", "x"),
        color = ~Source_complete2,
        colors = c('#1b9e77', '#d95f02', '#FFC300', '#C70039', '#581845', '#DAF7A6', '#349b88', '#f5006c'),
        text = ~ID,
        hoverinfo = 'text', textposition = "bottom right",
        marker = list(size = 5.5, opacity = 0.3) ) %>%
  layout(scene = list(xaxis = list(title = 'Valence' ),
                      yaxis = list(title = 'Arousal', range = c(2, 8) ),
                      zaxis = list(title = 'Dominance' ) ),
         legend = list(orientation = 'h',
                       # bgcolor = "#E2E2E2",
                       bordercolor = "#E2E2E2",
                       borderwidth = 2,
                       font = list(size = 12) ),
         title = 'PAD 3D space: clusters and stimulus types')






# HIERARCHICAL CLUSTERING..................................................---------------------


# > 1) Ward-type plot -----------------------------------------------------

wss <- vector()
for (i in 2:15) { 
  wss[i] <- sum( kmeans( forclus_3cols, centers = i )$withinss, iter.max = 100 )
}

wss_data <- data.frame( ClusterNumber = 1:15, WSS = wss )

ggplot(wss_data, aes(x = ClusterNumber, y = WSS ) ) + geom_point( size = 2 ) +
  geom_line(alpha = 0.3) + theme_grey( base_size = 17 ) +
  xlab("Number of Clusters") + ylab("Within-groups sum of squares (WSS)") +
  ggtitle( "Number of clusters vs. reduction in WSS" )

# ggtitle( expression( paste( "Number of ", italic("k"), "-means clusters vs. reduction in WSS" ) ) )



# > 2) The "Elbow" Method for Clustering Evaluation -----------------------


Step1 <- dist( forclus_3cols, method = "euclidean", diag = FALSE, upper = FALSE )
Step2 <- hclust( Step1, method = "ward.D", members= NULL )
Step3 <- css.hclust( Step1, hclust.obj = Step2 )
Step4 <- elbow.batch( Step3, 0.05, 0.80, precision = 3 )
print(Step4) # It says a good k = 5!

# And:

ggplot(Step3, aes(x = Step3$k, y = Step3$ev ) ) + geom_point( size = 2 ) +
  geom_line(alpha = 0.3) + theme_grey( base_size = 17 ) +
  xlab("Number of Clusters") + ylab("Explained variance") +
  ggtitle( "Variance explained by clusters" )



# > 3) Cluster validation indices:: clValid -------------------------------


hierk1 <- clValid( as.matrix( forclus_3cols ), 2:8, 
                  clMethods = c("hierarchical"), 
                  validation = c("internal","stability"),
                  maxitems = 10000, 
                  metric = "euclidean")
optimalScores(hierk1) # mostly 2 or 8... for hierarchical. 



# > 4) Cophenetic correlations --------------------------------------------
# And
# > 5) Gower (1983) distance ----------------------------------------------

# Same loop for both:


hc_solutions <- c( list(), list() )

for ( method in c( "euclidean", "correlation" ) ) {
  for ( link in c( "complete", "mcquitty", "ward", "average" ) ) {
    
    print( paste( method, "+", link) )
    
    hc_solution <- hcluster( forclus_3cols,
                             method = method, 
                             link = link )
    hc_solutions[[ method ]][[ link ]] <- hc_solution
    
    print( paste( "Cophenetic correlation =",
                  cor( Dist( forclus_3cols, 
                             method = method ), 
                       cophenetic( hc_solutions[[ method ]][[ link ]] ) ),
                  "Gower distance =",
                  sum( ( Dist( forclus_3cols, method = method ) - 
                           cophenetic( hc_solutions[[ method ]][[ link ]] ) ) ^ 2 )
                )
         )
    cat("\n")
  }
}

layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
layout.show(4)
lapply(hc_solutions[[1]], plot)
lapply(hc_solutions[[2]], plot)
layout(matrix(c(1, 1, 1, 1), 2, 2, byrow = TRUE))
layout.show(1)

# "correlation" + "ward" #  *** Best.
# "correlation" + "average" #  *** Second best.

# Gower: smallest value! Average linkage + correlation = best match





# > 6) Average silhouette widths ------------------------------------------

k_value <- vector()
asw <- vector()

for ( k in 2 : ( nrow( forclus_3cols ) - 1 ) ) {
  sil <- silhouette( cutree( hc_solutions[[ 2 ]][[ 4 ]], k = k), 
                     Dist(forclus_3cols, method = "correlation" ) )
  k_value[k] <- k
  asw[k] <- summary( sil )$avg.width
}

k.best <- which.max( asw )

# windows(title="Silhouettes - Average - k = 2 to n-1")
plot(k_value, asw, 
     type = "h", xlim = c( 0, 20 ), 
     main = "Silhouette-optimal number of clusters", 
     xlab = "k (number of groups)", ylab= "Average silhouette width" ) 

axis(1, k.best, 
     paste("optimum", k.best, sep = "\n"), 
     col = "red", font = 2, col.axis = "red")

points( k.best, max( asw, na.rm = T ), pch = 16, col = "red", cex = 1.5)

cat("Silhouette-optimal number of clusters k =", k.best, "\n", 
    "with an average silhouette width of", max( asw, na.rm = T ), "\n")



# > 7) Shepard-like diagrams ----------------------------------------------


plot.new()
par(mfrow = c( 2, 2 ) )

# With correlation distances:

for ( link in names( hc_solutions[[2]] ) ) {
  
  plot( Dist( forclus_3cols, method = "correlation" ), 
        cophenetic( hc_solutions[[2]][[link]] ), 
        xlab = "Stimulus distance", ylab = "Cophenetic distance", 
        asp = 1, xlim = c( 0, 3 ), ylim = c( 0, 3 ), 
        main = c( link, paste("Cophenetic correlation =", 
                              round(cor( Dist( forclus_3cols, 
                                               method = "correlation" ), 
                                         cophenetic( hc_solutions[[2]][[link]] ) ), 
                                    3 ) ) ) )
  abline(0,1)
  lines( lowess( Dist( forclus_3cols, 
                       method = "correlation" ), 
                 cophenetic( hc_solutions[[2]][[link]] ) ), 
         col = "red" )
  
}




# > 8) Graphs of fusion level values --------------------------------------


plot.new()
par(mfrow = c( 2, 2 ) )
# hc_solutions[[2]], i.e. all the linkage methods with correlation distances

for ( link in names( hc_solutions[[2]] ) ) {
  
  # Plot the fusion level values of the single linkage clustering
  plot(hc_solutions[[2]][[link]]$height, 
       nrow( forclus_3cols ) : 2, type = "S", 
       main = paste("Fusion levels - Linkage:", link), 
       ylab = "k (number of clusters)", xlab = "h (node height)",
       col = "grey")
  
  text( hc_solutions[[2]][[link]]$height, 
       nrow( forclus_3cols ) : 2, 
       nrow( forclus_3cols ) : 2, 
       col = "red", cex = 0.8)
  
}


par(mfrow = c( 1, 1 ) )

# Choose a common number of groups:  where at least a small jump is present in all four graphs of fusion levels
# Maybe k = 5 here!

# So:

# gr <- cutree( hc_solutions[[2]][['ward']], 5 )
gr <- cutree( hc_solutions[[2]][['average']], 5 )

cloud(Valence ~ Arousal * Dominance, data = forclus_3cols, 
      pch = gr, col = gr, type = "p", 
      perspective = T, distance = 0.4, shade = T, lex = 2, 
      xlab = list( "Arousal", cex = 1.7 ), 
      ylab = list( "Dominance", cex = 1.7 ), 
      zlab = list( "Valence", cex = 1.7 ), 
      # main = list( "Clusters", cex = 2 ), 
      R.mat = diag(4), screen = list( z = -40, x = -42 ), 
      aspect = c( 1, 1 ),
      key = list( text = list( levels( as.factor( gr ) ), cex = 1.2 ),
                  title = "Cluster number:", cex.title = 1.3,
                  points = list( pch = 1 : 5, col = 1 : nlevels( as.factor( gr ) ), lwd = 3 ), 
                  columns = 2, 
                  lines = list( col = 1 : nlevels( as.factor( gr ) ), lwd = 3 ) )
)



# > 9) Mantel Optimality --------------------------------------------------

max_k <- 30
kt <- data.frame( k = 1 : max_k, r = NA )

for ( i in 2 : max_k ) {
  
  gr <- cutree( hc_solutions[[2]][['average']], i )
  # gr <- cutree( hc_solutions[[2]][['ward']], i ) # Use interchangeably
  
  distgr <- grpdist( gr )
  mt <- cor( Dist( forclus_3cols, method = "correlation" ),
             distgr, method = "pearson")
  kt[ i, 2 ] <- mt
}

k.best <- which.max( kt$r )

# windows(title="Optimal number of clusters - Mantel")
plot(kt$k, kt$r, type = "h", xlim = c(0,30),
     main = "Mantel-optimal number of clusters - Average", 
     xlab = "k (number of groups)", 
     ylab="Pearson's correlation")

axis(1, k.best, 
     paste("optimum", k.best, sep = "\n"), 
     col = "red", font = 2, col.axis = "red")

points(k.best, max( kt$r, na.rm = T ), 
       pch = 16, col = "red", cex = 1.5)

cat("Mantel-optimal number of clusters k =", k.best, "\n", 
    "with a matrix linear correlation of", max( kt$r, na.rm = T ), "\n")

# Hmmm, suggests k = 2 for Ward, and k = 3 for Average linkage.




# K-MEANS................................................................... -----------------------------


# > 1) SSI ----------------------------------------------------------------

cas1 <- cascadeKM( forclus_3cols, 2, 8, iter = 3000, criterion = "ssi" )
plot(cas1, sortg = T) # suggests 8

cas2 <- cascadeKM( forclus_3cols, 2, 8, iter = 3000, criterion = "calinski" )
plot( cas2, sortg = T ) # suggests 2


# > 2) ClValid ------------------------------------------------------------


hierk2 <- clValid(forclus_3cols, 2 : 8,
                  clMethods = c("kmeans"), 
                  validation = c("internal", "stability"), 
                  maxitems = 10000, metric = "euclidean")

optimalScores(hierk2) # 2, 6 or 8... for kmeans.




# > 3) Build 2 random halves ----------------------------------------------


half1 <- forclus_3cols[sample( rownames( forclus_3cols ), 
                              as.integer( nrow(forclus_3cols) / 2 ),
                              replace = FALSE), ]

half2 <- forclus_3cols[ setdiff( rownames( forclus_3cols ), 
                                 rownames( half1 ) ), ]

# Must make half1 and half 2 exactly equal in dim
to_exclude <- sample( rownames( half2 ), 1 )
half2 <- forclus_3cols[ setdiff( rownames(half2), to_exclude ), ]


# Compare ideal cluster no. for both halves
s1_ssi <- cascadeKM(half1, 2, 8, iter = 1000, criterion = "ssi")
plot(s1_ssi, sortg = T) # says 6 or 8 , usually

s2_ssi <- cascadeKM(half2, 2, 8, iter = 1000,  criterion = "ssi")
plot(s2_ssi, sortg = T) # 6 or 8 again


s1_c <- cascadeKM(half1, 2, 8, iter = 1000, criterion = "calinski")
plot(s1_c, sortg = T) # 3

s2_c <- cascadeKM(half2, 2, 8, iter = 1000,  criterion = "calinski")
plot(s2_c, sortg = T) # 2


# Compare clusterings with a chi-square.

for ( distance in c( "correlation", "euclidean" ) ) {
  for ( k in 2 : 8 ) {

    half1_KM <- Kmeans(as.matrix( half1 ), k, iter.max = 2000, nstart = 1000, method = distance)
    half2_KM <- Kmeans(as.matrix( half2 ), k, iter.max = 2000, nstart = 1000, method = distance)
    
    cramer <- assocstats( table( half1_KM$cluster, half2_KM$cluster ) )$cramer 
    
    print( paste( "For distance metric =", distance,
                  ", and k =", k, 
                  ", Cramer's V =", cramer ) )
    cat("\n")
    
  }
}

# Poor results, which only tend to improve with k - not a great sign.




