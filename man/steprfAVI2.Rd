% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steprfAVI2.R
\name{steprfAVI2}
\alias{steprfAVI2}
\title{Select predictive variables for random forest by AVI and accuracy
in a stepwise algorithm}
\usage{
steprfAVI2(
  trainx,
  trainy,
  cv.fold = 10,
  ntree = 500,
  rpt = 20,
  predacc = "VEcv",
  importance = TRUE,
  maxk = c(4),
  nsim = 100,
  min.n.var = 2,
  corr.threshold = 0.5,
  rseed = 1234,
  ...
)
}
\arguments{
\item{trainx}{a dataframe or matrix contains columns of predictor variables.}

\item{trainy}{a vector of response, must have length equal to the number of
rows in trainx.}

\item{cv.fold}{integer; number of folds in the cross-validation. if > 1,
then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
validation that is recommended.}

\item{ntree}{number of trees to grow. This should not be set to too small a
number, to ensure that every input row gets predicted at least a few times.
By default, 500 is used.}

\item{rpt}{iteration of cross validation.}

\item{predacc}{"VEcv" for vecv for numerical data, or "ccr" (i.e., correct
classification rate) or "kappa" for categorical data.}

\item{importance}{imprtance of predictive variables.}

\item{maxk}{maxk split value. By default, 4 is used.}

\item{nsim}{iteration number. By default, 100 is used.}

\item{min.n.var}{minimum number of predictive variables remained in the final
predictive model the default is 2. If 1 is used, then warnings: 'invalid mtry:
 reset to within valid range' will be issued, which should be ignored.}

\item{corr.threshold}{correlation threshold and the defaults value is 0.5.}

\item{rseed}{random seed. By default, 1234 is used.}

\item{...}{other arguments passed on to randomForest.}
}
\value{
A list with the following components: 1) variable.removed:
variable removed based on AVI, 2) predictive.accuracy: averaged predictive
 accuracy of the model after excluding the variable.removed, 3)
 delta.accuracy: contribution to accuracy by each variable.removed, and 4)
  predictive.accuracy2: predictive accuracy matrix of the model after
   excluding the variable.removed for each iteration.
}
\description{
This function is similar to 'steprfAVI'; the only difference
is that 'set.seed()' is added before each code line that involves
 randomness and such additions alter the results considerably.
}
\examples{
\donttest{
library(spm)
data(petrel)
steprf1 <- steprfAVI2(trainx = petrel[, c(1,2, 6:9)], trainy = petrel[, 5],
 rpt = 2, predacc = "VEcv", importance = TRUE, nsim = 3, min.n.var = 2)
steprf1

#plot steprf1 results
library(reshape2)
pa1 <- as.data.frame(steprf1$predictive.accuracy2)
names(pa1) <- steprf1$variable.removed
pa2 <- melt(pa1, id = NULL)
names(pa2) <- c("Variable","VEcv")
library(lattice)
par (font.axis=2, font.lab=2)
with(pa2, boxplot(VEcv~Variable, ylab="VEcv (\%)", xlab="Predictive variable removed"))

barplot(steprf1$delta.accuracy, col = (1:length(steprf1$variable.removed)),
names.arg = steprf1$variable.removed, main = "Predictive accuracy vs variable removed",
font.main = 4, cex.names=1, font=2, ylab="Increase in VEcv (\%)")
}

}
\references{
Li, J. (2022). Spatial Predictive Modeling with R. Boca Raton,
Chapman and Hall/CRC.

Li, J. 2013. Predicting the spatial distribution of seabed
gravel content using random forest, spatial interpolation methods and
their hybrid methods. Pages 394-400  The International Congress on
Modelling and Simulation (MODSIM) 2013, Adelaide.

Li, J., Alvarez, B., Siwabessy, J., Tran, M., Huang, Z., Przeslawski, R.,
Radke, L., Howard, F. and Nichol, S. (2017). "Application of random forest,
 generalised linear model and their hybrid methods with geostatistical
 techniques to count data: Predicting sponge species richness."
 Environmental Modelling & Software 97: 112-129.

Li, J., Siwabessy, J., Huang, Z., Nichol, S. (2019). "Developing an optimal
 spatial predictive model for seabed sand content using machine learning,
 geostatistics and their hybrid methods." Geosciences 9 (4):180.

Li, J., Siwabessy, J., Tran, M., Huang, Z. and Heap, A., 2014. Predicting
 Seabed Hardness Using Random Forest in R. Data Mining Applications with R. Y.
  Zhao and Y. Cen. Amsterdam, Elsevier: 299-329.

Liaw, A. and M. Wiener (2002). Classification and Regression by
randomForest. R News 2(3), 18-22.

Smith, S.J., Ellis, N., Pitcher, C.R., 2011. Conditional variable
importance in R package extendedForest. R vignette
[http://gradientforest.r-forge.r-project.org/Conditional-importance.pdf].

Chang, W. 2021. Cookbook for R. http://www.cookbook-r.com/.
}
\author{
Jin Li
}
