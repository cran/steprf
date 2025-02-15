#' @title Select predictive variables for random forest by  various variable importance methods and predictive accuracy
#' in a stepwise algorithm
#'
#' @description This function is to select predictive variables for random forest
#' by various variable importance methods (i.e., AVI, Knowledge informed AVI
#' (KIAVI), KIAVI2) and predictive accuracy. It is implemented via the functions 'steprfAVI'
#' and 'steprfAVIPredictors'.
#'
#' @param trainx a dataframe or matrix contains columns of predictor variables.
#' @param trainy a vector of response, must have length equal to the number of
#' rows in trainx.
#' @param method a variable selection method for 'RF'; can be: "AVI", "KIAVI"
#'  and "KIAVI2". If "AVI" is used, it would produce tha same results as
#'   'steprfAVI'. By default, "KIAVI" is used.
#' @param cv.fold integer; number of folds in the cross-validation. if > 1,
#' then apply n-fold cross validation; the default is 10, i.e., 10-fold cross
#' validation that is recommended.
#' @param ntree number of trees to grow. This should not be set to too small a
#' number, to ensure that every input row gets predicted at least a few times.
#' By default, 500 is used.
#' @param rpt iteration of cross validation.
#' @param predacc "VEcv" for vecv for numerical data, or "ccr" (i.e., correct
#'  classification rate) or "kappa" for categorical data.
#' @param importance imprtance of predictive variables.
#' @param maxk maxk split value. By default, 4 is used.
#' @param nsim iteration number. By default, 100 is used.
#' @param delta.predacc minimum changes between the accuracies of two consecutive
#' predictive models.
#' @param min.n.var minimum number of predictive variables remained in the final
#' predictive model the default is 2. If 1 is used, then warnings: 'invalid mtry:
#'  reset to within valid range' will be issued, which should be ignored.
#' @param corr.threshold correlation threshold and the defaults value is 0.5.
#' @param ... other arguments passed on to randomForest.
#'
#' @return A list with the following components: 1) steprfPredictorsFinal: the
#'  variables selected for the last RF model, whether it is of the highest
#'  predictive accuracy need to be confirmed using 'max.predictive.accuracy'
#'  that is listed next; 2) max.predictive.accuracy: the predictive accuracy
#'  of the most accurate RF model for each run of 'steprfAVI', which can be used
#'  to confirm the model with the highest accuracy, 3) numberruns: number of runs
#'  of 'steprfAVI'; 4) laststepAVI: the outpouts of last run of 'steprfAVI'; 5)
#'   steprfAVIOutputsAll: the outpouts of all 'steprfAVI' produced during the
#'   variable selection process; 6) steprfPredictorsAll: the outpouts of
#'   'steprfAVIPredictors' for all 'steprfAVI' produced during the variable
#'    selection process; 7) KIAVIPredictorsAll: predictors used for all
#'    'steprfAVI' produced during the variable selection process; for a method
#'    "AVI", if the variables are different from those in the traning dataset,
#'    it suggests that these variables should be tested if the predictive
#'    accuracy can be further improved.
#'
#' @note In 'steprf', 'steprfAVI' is used instead of 'steprfAVI1' and
#' 'steprfAVI2'. This is because: 1) 'avi' is expected to change with the
#' removal of each predictor, but in 'steprfAVI1' the averaged variable
#' importance is calculated only once and is from the full model only, so
#' its use is expected to produce a less optimal model, hence not used; and
#' 2) the 'steprf' would lead to the same set of predictors as that for
#'  'steprfAVI2' if 'steprfAVI2' is used, so it is not used either.
#'
#' @references Li, J. (2022). Spatial Predictive Modeling with R. Boca Raton,
#' Chapman and Hall/CRC.
#'
#' Li, J. (2019). "A critical review of spatial predictive modeling process
#' in environmental sciences with reproducible examples in R." Applied
#' Sciences 9: 2048.
#'
#' Li, J. 2013. Predicting the spatial distribution of seabed
#' gravel content using random forest, spatial interpolation methods and
#' their hybrid methods. Pages 394-400  The International Congress on Modelling
#'  and Simulation (MODSIM) 2013, Adelaide.
#'
#' Li, J., Alvarez, B., Siwabessy, J., Tran, M., Huang, Z., Przeslawski, R.,
#' Radke, L., Howard, F. and Nichol, S. (2017). "Application of random forest,
#'  generalised linear model and their hybrid methods with geostatistical
#'  techniques to count data: Predicting sponge species richness."
#'  Environmental Modelling & Software 97: 112-129.
#'
#' Li, J., Siwabessy, J., Huang, Z., and Nichol, S. (2019). "Developing an
#'  optimal spatial predictive model for seabed sand content using machine
#'  learning, geostatistics and their hybrid methods." Geosciences 9 (4):180.
#'
#' Li, J., Siwabessy, J., Tran, M., Huang, Z. and Heap, A., 2014. Predicting
#'  Seabed Hardness Using Random Forest in R. Data Mining Applications with R. Y.
#'   Zhao and Y. Cen. Amsterdam, Elsevier: 299-329.
#'
#' Li, J., Tran, M. and Siwabessy, J., 2016. Selecting optimal random forest
#' predictive models: a case study on predicting the spatial distribution of
#'  seabed hardness. PLOS ONE 11(2): e0149089.
#'
#' Liaw, A. and M. Wiener (2002). Classification and Regression by
#' randomForest. R News 2(3), 18-22.
#'
#' Smith, S.J., Ellis, N., Pitcher, C.R., 2011. Conditional variable
#' importance in R package extendedForest. R vignette
#' [http://gradientforest.r-forge.r-project.org/Conditional-importance.pdf].
#'
#' @author Jin Li
#' @examples
#' \donttest{
#' library(spm)
#' data(petrel)
#' set.seed(1234)
#' steprf1 <- steprf(trainx = petrel[, c(1,2, 6:9)], trainy =
#' petrel[, 5], method = "KIAVI", rpt = 2, predacc = "VEcv", importance = TRUE,
#'  nsim = 3, delta.predacc = 0.01)
#' names(steprf1)
#' steprf1$steprfPredictorsFinal$variables.most.accurate
#' steprf1$max.predictive.accuracy
#' }
#'
#' @export
steprf <- function (trainx, trainy, method = "KIAVI", cv.fold = 10, ntree = 500, rpt = 20, predacc = "VEcv", importance = TRUE, maxk = c(4), nsim = 100, delta.predacc = 0.001, min.n.var = 2, corr.threshold = 0.5, ...) {
  n <- nrow(trainx); p <- ncol(trainx)
  predictive.accuracy1 <- NULL
  predictive.accuracy2 <- NULL

  #predictive.accuracy2 <- vector((p - min.n.var), mode = "list")
  steprfAVI1 <- list()
  steprfAVIPredictors1 <- list()
  kiavi.predictors <- list()
  max.predictive.accuracy <- NULL

  if (method == "AVI") {
    steprf1 <- steprfAVI(trainx = trainx, trainy = trainy, ntree = ntree, cv.fold = cv.fold, predacc = predacc,  rpt = rpt, importance = TRUE, nsim = nsim, min.n.var = min.n.var, corr.threshold = corr.threshold)

    i <- 1
    steprfAVI1[[i]] <- steprf1

    predictors.selected <- steprfAVIPredictors(steprf1, trainx = trainx)
    steprfAVIPredictors1[[i]] <- predictors.selected
    predictive.accuracy1 <- predictors.selected$max.predictive.accuracy
    max.predictive.accuracy[i] <- predictive.accuracy1
    kiavi.predictors[[i]] <- union(predictors.selected$variables.most.accurate, predictors.selected$PABV)

  } else (
  if (method == "KIAVI") {
  steprf1 <- steprfAVI(trainx = trainx, trainy = trainy, ntree = ntree, cv.fold = cv.fold, predacc = predacc,  rpt = rpt, importance = TRUE, nsim = nsim, min.n.var = min.n.var, corr.threshold = corr.threshold)

  i <- 1
  steprfAVI1[[i]] <- steprf1

  predictors.selected <- steprfAVIPredictors(steprf1, trainx = trainx)
  steprfAVIPredictors1[[i]] <- predictors.selected
  predictive.accuracy1 <- predictors.selected$max.predictive.accuracy
  max.predictive.accuracy[i] <- predictive.accuracy1
  kiavi.predictors[[i]] <- union(predictors.selected$variables.most.accurate, predictors.selected$PABV)

  trainx1 <- trainx[, kiavi.predictors[[i]], drop = FALSE]
  p <- ncol(trainx1)

  for (i in 2:(p - min.n.var)) {
    steprf1 <- NULL
    steprf1 <- steprfAVI(trainx = trainx1, trainy = trainy, ntree = ntree, cv.fold = cv.fold, predacc = predacc,  rpt = rpt, importance = TRUE, nsim = nsim, min.n.var = min.n.var, corr.threshold = corr.threshold)

    steprfAVI1[[i]] <- steprf1

    predictors.selected <- steprfAVIPredictors(steprf1, trainx = trainx1)
    steprfAVIPredictors1[[i]] <- predictors.selected
    predictive.accuracy2 <- predictors.selected$max.predictive.accuracy
    max.predictive.accuracy[i] <- predictive.accuracy2
    kiavi.predictors[[i]] <- union(predictors.selected$variables.most.accurate, predictors.selected$PABV)

    if ((predictive.accuracy2 - predictive.accuracy1) <= delta.predacc)
    {break
      } else
      {predictive.accuracy1 <- predictive.accuracy2
      trainx1 <- trainx[, kiavi.predictors[[i]], drop = FALSE]
      }
  }
  } else (
    if (method == "KIAVI2") {
      steprf1 <- steprfAVI(trainx = trainx, trainy = trainy, ntree = ntree, cv.fold = cv.fold, predacc = predacc,  rpt = rpt, importance = TRUE, nsim = nsim, min.n.var = min.n.var, corr.threshold = corr.threshold)

      i <- 1
      steprfAVI1[[i]] <- steprf1

      predictors.selected <- steprfAVIPredictors(steprf1, trainx = trainx)
      steprfAVIPredictors1[[i]] <- predictors.selected
      predictive.accuracy1 <- predictors.selected$max.predictive.accuracy
      max.predictive.accuracy[i] <- predictive.accuracy1
      kiavi.predictors[[i]] <- predictors.selected$PABV

      trainx1 <- trainx[, kiavi.predictors[[i]], drop = FALSE]
      p <- ncol(trainx1)

      for (i in 2:(p - min.n.var)) {
        steprf1 <- NULL
        steprf1 <- steprfAVI(trainx = trainx1, trainy = trainy, ntree = ntree, cv.fold = cv.fold, predacc = predacc, rpt = rpt, importance = TRUE, nsim = nsim, min.n.var = min.n.var, corr.threshold = corr.threshold)

        steprfAVI1[[i]] <- steprf1

        predictors.selected <- steprfAVIPredictors(steprf1, trainx = trainx1)
        steprfAVIPredictors1[[i]] <- predictors.selected
        predictive.accuracy2 <- predictors.selected$max.predictive.accuracy
        max.predictive.accuracy[i] <- predictive.accuracy2
        kiavi.predictors[[i]] <- predictors.selected$PABV

        if ((predictive.accuracy2 - predictive.accuracy1) <= delta.predacc)
        {break
        } else
        {predictive.accuracy1 <- predictive.accuracy2
          trainx1 <- trainx[, kiavi.predictors[[i]], drop = FALSE]
          }
      }
    } else (
      stop ("This variable selection method is not supported!")
    )))
  list(steprfPredictorsFinal = predictors.selected, max.predictive.accuracy = max.predictive.accuracy, numberruns = i, laststepAVI = steprf1, steprfAVIOutputsAll = steprfAVI1, steprfPredictorsAll = steprfAVIPredictors1, KIAVIPredictorsAll = kiavi.predictors)
}
