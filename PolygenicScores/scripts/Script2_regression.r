################################################################################################################
# Step 4: Testing if polygenic scores are associated with phenotype 
################################################################################################################
# Change your working directory (you may need to change the filepath in next line)
setwd("~/genetic-epidemiology-practicals/PolygenicScores/data")

# removes anything already in memory
rm(list=ls())

# We have already created a file for you that contains all the schizophrenia polygenic scores together with 
# the phenotypic data (generated using SCZ_score_files_merge.r). Open this file:
data <- read.table("pgc_scz_scores_with_outcomes.txt", header=T)
# Explore the dataset
names(data)
head(data)
dim(data)

# Test if polygenic scores for schizophrenia are associated with BMI (x = z_score, y = BMI)
bmiassoc <- summary(glm(bmi ~ stdSCORE_S6, data=data))

#	Calculate the 95% confidence interval
lci <- bmiassoc$coefficients[2,1] - (1.96*bmiassoc$coefficients[2,2])
uci <- bmiassoc$coefficients[2,1] + (1.96*bmiassoc$coefficients[2,2])

# results:
bmiassoc
lci
uci


################################################################################################################
# (Optional) Step 5: Interpreting results across all polygenic score p-value thresholds (RStudio)
################################################################################################################

# create an empty dataframe with 10 columns and 12 rows
all.results <-data.frame(matrix(ncol=10,nrow=12))
# give each column a name
colnames(all.results) <- c("PRS","beta","se","t","p","lci","uci","r2","N","n_SNPs")
# create a list of the PRS p value thresholds (this is only used to name the rows of the PRS column) and 
# set j as 1 (to be used later)
thresh <-list("0.5", "0.4", "0.3", "0.2", "0.1", "0.05", "0.01", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7")
j=1

# loop to perform 12 linear regressions and save results to all.results data.frame
for (i in 1:12) {
results<-glm(data$bmi ~ data[, paste0("stdSCORE_S",i)], na.action="na.exclude")
data$pred <- predict(results)
r2 <- cor(data$bmi,data$pred,use="complete.obs")^2

all.results[i,1] <- thresh[j]                                   #PRS
all.results[i,2] <- summary(results)$coefficients[2,1]          #beta
all.results[i,3] <- summary(results)$coefficients[2,2]          #se
all.results[i,4] <- summary(results)$coefficients[2,3]          #t
all.results[i,5] <- summary(results)$coefficients[2,4]          #p
all.results[i,6] <- all.results[i,2]-(1.96*all.results[i,3])    #lci
all.results[i,7] <- all.results[i,2]+(1.96*all.results[i,3])    #uci
all.results[i,8] <- r2                                          #r2
all.results[i,9] <- nobs(results)                               #N
all.results[i,10] <- max(data[, paste0("CNT_S",i)])/2           #n_snps
j=j+1
}

# View the results
all.results

# To visualise the r2 of across all association analyses use:
jpeg(file="~/genetic-epidemiology-practicals/PolygenicScores/results/bmi_sczPRS_r2_plot.jpeg")
barplot(all.results$r2, ylim=c(0,0.02), ylab="BMI ~ standardised PRS (r2)", col="lightblue",
        xlab="Schizophrenia PRS p value thresholds", names.arg = all.results$PRS,
        main="r2 values for associations between SCZ PRS and BMI \nacross multiple PRS thresholds")
dev.off()

# Note that you can also use ggplot to visualise association results across all PRS thresholds:
# Install and load ggplot
# install.packages("ggplot2")
library(ggplot2)

ggplot(all.results, aes(x = PRS, y = beta)) +
  geom_point(size = 4) +
  ylim(-1,1) +
  geom_hline(yintercept = 0, col="red") +
  theme(panel.background = element_rect(fill = "grey97"), axis.text=element_text(size=12), axis.title=element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits=c("0.5", "0.4", "0.3", "0.2", "0.1", "0.05", "0.01", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7")) +
  geom_errorbar(aes(ymax = uci, ymin = lci), width=0.3) + 
  xlab("Schizophrenia PRS p value thresholds") + 
  ylab("Unit change in BMI per s.d. increase in schziophrenia PRS") +
  ggtitle("Associations between SCZ PRS and BMI across multiple PRS thresholds")
ggsave("~/genetic-epidemiology-practicals/PolygenicScores/results/bmi_scz_associations.jpeg")
