####################################################################
# Machine Learning - Project
#
# Javier de Jesús Flores Herrera
# Nejada Karriqi
####################################################################

####################################################################
# SECTION 1: LOADING NEEDED LIBRARIES & SETUP SOME VARIABLES
####################################################################

library(chemometrics)
library(corrplot)
library(DMwR)
library(factoextra)
library(FactoMineR)
library(MASS)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(purrr) 
library(tidyr) 

# Color palette
colors <- sample(c('#000000', '#FF0000', '#808000', '#00FF00', '#008000', 
                   '#00FFFF', '#0000FF', '#FF00FF', '#800080', '#ffa500'))

####################################################################
# SECTION 1: LOADING THE DATA
####################################################################

# Read red wine csv file
red = read.delim("winequality-red.csv", sep=";", header=TRUE);

# Read white wine csv file
white = read.delim("winequality-white.csv", sep=";", header=TRUE);


####################################################################
# SECTION 1: PREPROCESSING
####################################################################



#Check number of observations, number of variables
dim(red) 
dim(white)

#Check types of the variables for each wine
str(red)
str(white)

## Turn quality variable into factor
red$quality = as.factor(red$quality)
white$quality = as.factor(white$quality)


# There are not missing values
sum(is.na(red))
sum(is.na(white))

# Check how quality scores are distributed

table(red$quality)
table(white$quality)


####################################################################
# SECTION 7: GAUSSIANITY AND TRANSFORMATIONS
####################################################################

##########----------------------RED-----------------------##########
# Plot distribution of variables
red %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value,fill=key)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram(bins=sqrt(nrow(red))) +
  theme(legend.position="none") 
# Plot distribution of variables using log10
red %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value,fill=key)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram(bins=sqrt(nrow(red))) +
  theme(legend.position="none") +
  scale_x_continuous(trans='log10')

##########---------------------WHITE----------------------##########
white %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value,fill=key)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram(bins=sqrt(nrow(white))) +
  theme(legend.position="none") 
# Plot distribution of variables using log10
white %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value,fill=key)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram(bins=sqrt(nrow(white))) +
  theme(legend.position="none") +
  scale_x_continuous(trans='log10')


#*******TODO: use boxcox to make variables gaussian********

####################################################################
# SECTION 1:  Outlier detection
####################################################################

#index for quality variable
indexQuality <- 12

##########----------------------RED-----------------------##########
# create a new dataframe without the quality variable
numRed <- red[,-indexQuality]

out.scores <- data.frame(lof=lofactor(numRed, k=5))
out.scores$index <- rownames(out.scores)
rownames(out.scores) <- rownames(numRed)
plot(lof ~ index, out.scores, main="Multivariate outlier detection with Local Outlier Factor",
     ylab="Local Outlier Factor", xlab="Index")
abline(h=1.4, col="red")
#identify which cases are outliers 
red_outliers.lof <- subset(out.scores, lof > 1.4, row.names=rownames(out.scores))
text(lof~index, red_outliers.lof, labels=rownames(red_outliers.lof), pos=1)

#Remove outliers, if the outlier
table(red$quality)
redOutliers <- red[as.numeric(rownames(red_outliers.lof)),]
red <- red[-as.numeric(rownames(redOutliers[which(redOutliers$quality != 8 & redOutliers$quality != 7),])),]
table(red$quality)

plotRedOut <- ggplot(redOutliers, aes(x=quality, fill = quality)) +
  geom_bar(stat="count") +
  geom_text(position = "stack", stat='count',aes(label=..count..), vjust = -0.5)+
  labs(y="Num of Observations", x="Red wine Quality") +
  theme(legend.position="none")

##########---------------------WHITE----------------------##########

# create a new dataframe without the quality variable
numWhite <- white[,-indexQuality]

out.scores <- data.frame(lof=lofactor(numWhite, k=5))
out.scores$index <- rownames(out.scores)
rownames(out.scores) <- rownames(numWhite)
plot(lof ~ index, out.scores, main="Multivariate outlier detection with Local Outlier Factor",
     ylab="Local Outlier Factor", xlab="Index")
abline(h=1.4, col="red")
#identify which cases are outliers 
white_outliers.lof <- subset(out.scores, lof > 1.4, row.names=rownames(out.scores))
text(lof~index, white_outliers.lof, labels=rownames(white_outliers.lof), pos=1)


#Remove outliers
table(white$quality)
whiteOutliers <- white[as.numeric(rownames(white_outliers.lof)),]
white <- white[-as.numeric(rownames(whiteOutliers[which(whiteOutliers$quality != 8 & whiteOutliers$quality != 9),])),]
table(white$quality)

plotWhiteOut <- ggplot(whiteOutliers, aes(x=quality, fill = quality)) +
  geom_bar(stat="count") +
  geom_text(position = "stack", stat='count',aes(label=..count..), vjust = -0.5)+
  labs(y="Num of Observations", x="White wine Quality") +
  theme(legend.position="none")


## Comparison of # of outliers in each dataset.
grid.arrange(plotRedOut, plotWhiteOut, top = "Outliers",nrow = 1, ncol=2)

####################################################################
# SECTION 1:  perform PCA
####################################################################


###########-------------Correlation plot-----------################

red1 <- red
red1$quality <- as.numeric(red$quality)
white1 <- white
white1$quality <- as.numeric(white$quality)
par(mfrow=c(1,2))
corrplot(cor(red1,use="complete.obs"),tl.cex = 0.75, title = "Red wines")
corrplot(cor(white1,use = "complete.obs"), tl.cex = 0.75, title = "White wines") 
par(mfrow=c(1,1))
##########---------------------RED----------------------##########
pca.red <- PCA(red1, quali.sup = c(12), scale.unit=T, graph = F)

fviz_pca_var(pca.red,col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), ggtheme = theme_minimal())

fviz_pca_ind(pca.red,
             label = "none", # hide individual labels
            # habillage = red1$quality, # color by region
            col.ind.sup = "blue",
             palette = colors,
             invisible = "ind.sup",
             show.legend = FALSE
);

# using raiser rule 4, elbow suggest 5 components. We take 5 since it's almost 80% of the variance
plot(pca.red$eig[,1],type="b",col='blue',ylab="Eigenvalue",xlab="Component Number", main= "Cumulative Variance Explained of Eigen Values")
abline(h=1,lty=2,col="red")
text(pca.red$eig[,1], labels = round(pca.red$eig[,3], digits = 2), col = 4, pos = 1)

##########---------------------WHITE----------------------##########
pca.white <- PCA(white,quali.sup = c(12), scale.unit=T, graph = F)

fviz_pca_var(pca.white,col.var="contrib",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), ggtheme = theme_minimal())

fviz_pca_ind(pca.white,
            label = "none", # hide individual labels
            #habillage = white$quality, # color by region
            palette = colors,
            invisible = "ind.sup",
            show.legend = FALSE
)
# using raiser rule 4, elbow suggest 5 components. We take 5 since it's almost 80% of the variance
plot(pca.white$eig[,1],type="b",col='blue',ylab="Eigenvalue",xlab="Component Number", main= "Cumulative Variance Explained of Eigen Values")
abline(h=1,lty=2,col="red")
text(pca.white$eig[,1], labels = round(pca.white$eig[,3], digits = 2), col = 4, pos = 1)
####################################################################
# SECTION 1:  Clustering
####################################################################


##########----------------------RED-----------------------##########

## first, perform a hierarchical clustering
# significant dimensions
signi <- 5
# if we perform on clustering on pca, get factorial coordinates of individuals
df.pca <- pca.red$ind$coord[, 1:signi]
cluster<- HCPC(df.pca,nb.clust=-1,method="ward",consol=TRUE)
##########---------------------WHITE----------------------##########


