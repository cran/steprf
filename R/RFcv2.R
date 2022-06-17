#' @title Cross validation, n-fold for random forest (RF)
#'
#' @description This function is a cross validation function
#' for random forest. It is for functions 'steprf', 'steprfAVI', ect.
#'
#' @param trainx a dataframe or matrix contains columns of predictor variables.
#' @param trainy a vector of response, must have length equal to the number of
#' rows in trainx.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param mtry a function of number of remaining predictor variables to use as
#' the mtry parameter in the randomForest call.
#' @param ntree number of trees to grow. This should not be set to too small a
#' number, to ensure that every input row gets predicted at least a few times.
#' By default, 500 is used.
#' @param predacc "VEcv" for vecv for numerical data, or "ccr" (i.e., correct
#'  classification rate) or "kappa" for categorical data.
#' @param ... other arguments passed on to randomForest.
#'
#' @return A list with the following component:
#' vecv for numerical data: ; or
#'  ccr (correct classification rate) for categorical data: .
#'
#' @note This function is largely based on rf.cv (see Li et al. 2013) and
#' rfcv in randomForest.
#'
#' @references Li, J., J. Siwabessy, M. Tran, Z. Huang, and A. Heap. 2013.
#' Predicting Seabed Hardness Using Random Forest in R. Pages 299-329 in Y.
#' Zhao and Y. Cen, editors. Data Mining Applications with R. Elsevier.
#'
#' Li, J. 2013. Predicting the spatial distribution of seabed gravel content
#' using random forest, spatial interpolation methods and their hybrid methods.
#' Pages 394-400  The International Congress on Modelling and Simulation
#' (MODSIM) 2013, Adelaide.
#'
#' Liaw, A. and M. Wiener (2002). Classification and Regression by
#' randomForest. R News 2(3), 18-22.
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(spm)
#' data(hard)
#' data(petrel)
#'
#' rfcv1 <- RFcv2(petrel[, c(1,2, 6:9)], petrel[, 5], predacc = "VEcv")
#' rfcv1
#'
#' rfcv2 <- RFcv2(hard[, -c(1, 17)], hard[, 17], predacc = "ccr")
#' rfcv2
#'
#' rfcv3 <- RFcv2(hard[, -c(1, 17)], hard[, 17], predacc = "kappa")
#' rfcv3
#'
#' n <- 10 # number of iterations, 60 to 100 is recommended.
#' VEcv <- NULL
#' for (i in 1:n) {
#' rfcv1 <- RFcv2(petrel[, c(1,2,6:9)], petrel[, 5], predacc = "VEcv")
#' VEcv [i] <- rfcv1
#' }
#' plot(VEcv ~ c(1:n), xlab = "Iteration for RF", ylab = "VEcv (%)")
#' points(cumsum(VEcv) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(VEcv), col = 'blue', lwd = 2)
#'
#' n <- 10 # number of iterations, 60 to 100 is recommended.
#' measures <- NULL
#' for (i in 1:n) {
#' rfcv1 <- RFcv2(hard[, c(4:6)], hard[, 17], predacc = "ccr")
#' measures <- rbind(measures, rfcv1)
#' }
#' plot(measures ~ c(1:n), xlab = "Iteration for RF", ylab = "Correct
#' classification Rate  (%)")
#' points(cumsum(measures) / c(1:n) ~ c(1:n), col = 2)
#' abline(h = mean(measures), col = 'blue', lwd = 2)
#' }
#'
#' @export
RFcv2 <- function (trainx, trainy, cv.fold = 10, mtry = if (!is.null(trainy) &&
  !is.factor(trainy)) max(floor(ncol(trainx) / 3), 1) else
  floor(sqrt(ncol(trainx))), ntree = 500, predacc = "VEcv", ...) {
  classRF <- is.factor(trainy)
  n <- nrow(trainx)
  if (classRF) {
    f <- trainy
  } else {if (dim(table(trainy)) <= 4) {
    f <- trainy
      } else {
    f <- cut(trainy, c(-Inf, stats::quantile(trainy, 1:4/5), Inf))
    }
  }

  nlvl <- table(f)
  idx <- numeric(n)
  for (i in 1:length(nlvl)) {
    idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold, length = nlvl[i]))
  }
  cv.pred <- NULL
  for (i in 1:cv.fold) {
    all.rf <- randomForest::randomForest(trainx[idx != i, , drop = FALSE],
      trainy[idx != i], trainx[idx == i, , drop = FALSE],
      trainy[idx == i], mtry = mtry, ntree=ntree)
    cv.pred[idx == i] <- all.rf$test$predicted
  }
  predictive.accuracy <- NULL
  if (predacc == "VEcv") {predictive.accuracy = spm::vecv(trainy, cv.pred)} else (
  if (predacc == "ccr") {predictive.accuracy = spm2::ccr(trainy, cv.pred)} else
    if (predacc == "kappa") {predictive.accuracy = psy::ckappa(as.data.frame(cbind(cv.pred, trainy)))$kappa} else (
  stop ("This measure is not supported in this version!")))
  predictive.accuracy
}
