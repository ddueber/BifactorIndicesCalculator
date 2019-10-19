test_that("ECV_SS Works", {
  Lambda <- matrix(c(  0, .82, .10,   0,
                       0, .77, .35,   0,
                       0, .79, .32,   0,
                       0, .66, .39,   0,
                       0, .51,   0, .71,
                       0, .56,   0, .43,
                       0, .68,   0, .13,
                       0, .60,   0, .50,
                     .47, .83,   0,   0,
                     .27, .60,   0,   0,
                     .28, .78,   0,   0,
                     .75, .55,   0,   0),
                   ncol = 4, byrow = TRUE)
  ECV_SS_Expected <- c(.3229228, .7138154, .1426991, .4064201)
  ## First with no names
  expect_equal(ECV_SS(Lambda), ECV_SS_Expected, tolerance = .000001)
  colnames(Lambda) = c("A", "B", "C", "D")
  names(ECV_SS_Expected) = c("A", "B", "C", "D")
  ## Now with names
  expect_equal(ECV_SS(Lambda), ECV_SS_Expected, tolerance = .000001)

  ## Now let's do one with no general factor
  Lambda2 <- matrix(c(.3,  0,  0,
                      .4,  0,  0,
                      .5,  0,  0,
                      0, .4,  0,
                      0, .6,  0,
                      0, .7,  0,
                      0,  0, .4,
                      0,  0, .5,
                      0,  0, .3),
                    ncol = 3, byrow = TRUE)
  ECV_SS_Expected2 <- c(1, 1, 1)
  expect_equal(ECV_SS(Lambda2), ECV_SS_Expected2, tolerance = .000001)

})

test_that("ECV_SG Works", {
  Lambda <- matrix(c(  0, .82, .10,   0,
                       0, .77, .35,   0,
                       0, .79, .32,   0,
                       0, .66, .39,   0,
                       0, .51,   0, .71,
                       0, .56,   0, .43,
                       0, .68,   0, .13,
                       0, .60,   0, .50,
                     .47, .83,   0,   0,
                     .27, .60,   0,   0,
                     .28, .78,   0,   0,
                     .75, .55,   0,   0),
                   ncol = 4, byrow = TRUE)
  ECV_SG_Expected <- c(.11744676, .71381542, .04862725, .12011057)
  ## First with no names
  expect_equal(ECV_SG(Lambda), ECV_SG_Expected, tolerance = .000001)
  colnames(Lambda) = c("A", "B", "C", "D")
  names(ECV_SG_Expected) = c("A", "B", "C", "D")
  ## Now with names
  expect_equal(ECV_SG(Lambda), ECV_SG_Expected, tolerance = .000001)

  ## Now let's do one with no general factor
  Lambda2 <- matrix(c(.3,  0,  0,
                      .4,  0,  0,
                      .5,  0,  0,
                       0, .4,  0,
                       0, .6,  0,
                       0, .7,  0,
                       0,  0, .4,
                       0,  0, .5,
                       0,  0, .3),
                      ncol = 3, byrow = TRUE)
  expect_equal(ECV_SG(Lambda2), NULL)
})

test_that("ECV_GS Works", {
  Lambda <- matrix(c(  0, .82, .10,   0,
                       0, .77, .35,   0,
                       0, .79, .32,   0,
                       0, .66, .39,   0,
                       0, .51,   0, .71,
                       0, .56,   0, .43,
                       0, .68,   0, .13,
                       0, .60,   0, .50,
                       .47, .83,   0,   0,
                       .27, .60,   0,   0,
                       .28, .78,   0,   0,
                       .75, .55,   0,   0),
                   ncol = 4, byrow = TRUE)
  ECV_GS_Expected <- c(.6770772, .71381542, .8573009, .5935799)
  ## First with no names
  expect_equal(ECV_GS(Lambda), ECV_GS_Expected, tolerance = .000001)
  colnames(Lambda) = c("A", "B", "C", "D")
  names(ECV_GS_Expected) = c("A", "B", "C", "D")
  ## Now with names
  expect_equal(ECV_GS(Lambda), ECV_GS_Expected, tolerance = .000001)

  ## Now let's do one with no general factor
  Lambda2 <- matrix(c(.3,  0,  0,
                      .4,  0,  0,
                      .5,  0,  0,
                      0, .4,  0,
                      0, .6,  0,
                      0, .7,  0,
                      0,  0, .4,
                      0,  0, .5,
                      0,  0, .3),
                    ncol = 3, byrow = TRUE)
  expect_equal(ECV_GS(Lambda2), NULL)

})

test_that("IECV Works", {
  Lambda <- matrix(c(  0, .82, .10,   0,
                       0, .77, .35,   0,
                       0, .79, .32,   0,
                       0, .66, .39,   0,
                       0, .51,   0, .71,
                       0, .56,   0, .43,
                       0, .68,   0, .13,
                       0, .60,   0, .50,
                       .47, .83,   0,   0,
                       .27, .60,   0,   0,
                       .28, .78,   0,   0,
                       .75, .55,   0,   0),
                   ncol = 4, byrow = TRUE)
  IECV_Expected <- c(0.9853458, 0.8287671, 0.8590502, 0.7411945, 0.3403559,
                       0.6290873, 0.9647402, 0.5901639, 0.7571994, 0.8316008,
                       0.8858474, 0.3497110)
  ## First with no names
  expect_equal(IECV(Lambda), IECV_Expected, tolerance = .000001)
  rownames(Lambda) = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L")
  names(IECV_Expected) = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L")
  ## Now with names
  expect_equal(IECV(Lambda), IECV_Expected, tolerance = .000001)

  ## Now let's do one with no general factor
  Lambda2 <- matrix(c(.3,  0,  0,
                      .4,  0,  0,
                      .5,  0,  0,
                      0, .4,  0,
                      0, .6,  0,
                      0, .7,  0,
                      0,  0, .4,
                      0,  0, .5,
                      0,  0, .3),
                    ncol = 3, byrow = TRUE)
  expect_equal(IECV(Lambda2), NULL)

})

test_that("Omega_S Works", {
  Lambda <- matrix(c(  0, .82, .10,   0,
                       0, .77, .35,   0,
                       0, .79, .32,   0,
                       0, .66, .39,   0,
                       0, .51,   0, .71,
                       0, .56,   0, .43,
                       0, .68,   0, .13,
                       0, .60,   0, .50,
                       .47, .83,   0,   0,
                       .27, .60,   0,   0,
                       .28, .78,   0,   0,
                       .75, .55,   0,   0),
                   ncol = 4, byrow = TRUE)
  Theta <- rep(1, nrow(Lambda)) - rowSums(Lambda^2)
  Omega_S_Expected <- c(.9067561, .9482359, .8915387, .8400528)
  ## First with no names
  expect_equal(Omega_S(Lambda, Theta), Omega_S_Expected, tolerance = .000001)
  colnames(Lambda) = c("A", "B", "C", "D")
  names(Omega_S_Expected) = c("A", "B", "C", "D")
  ## Now with names
  expect_equal(Omega_S(Lambda, Theta), Omega_S_Expected, tolerance = .000001)

  ## Now let's do one with no general factor
  Lambda2 <- matrix(c(.3,  0,  0,
                      .4,  0,  0,
                      .5,  0,  0,
                      0, .4,  0,
                      0, .6,  0,
                      0, .7,  0,
                      0,  0, .4,
                      0,  0, .5,
                      0,  0, .3),
                    ncol = 3, byrow = TRUE)
  Theta2 <- rep(1, nrow(Lambda2)) - rowSums(Lambda2^2)
  Omega_S_Expected2 <- c(.3654822, .5922131, .3654822)
  expect_equal(Omega_S(Lambda2, Theta2), Omega_S_Expected2,  tolerance = .000001)

})

test_that("Omega_H Works", {
  Lambda <- matrix(c(  0, .82, .10,   0,
                       0, .77, .35,   0,
                       0, .79, .32,   0,
                       0, .66, .39,   0,
                       0, .51,   0, .71,
                       0, .56,   0, .43,
                       0, .68,   0, .13,
                       0, .60,   0, .50,
                       .47, .83,   0,   0,
                       .27, .60,   0,   0,
                       .28, .78,   0,   0,
                       .75, .55,   0,   0),
                   ncol = 4, byrow = TRUE)
  Theta <- rep(1, nrow(Lambda)) - rowSums(Lambda^2)
  Omega_H_Expected <- c(.2642460, .8507481, .1133118, .3040647)
  ## First with no names
  expect_equal(Omega_H(Lambda, Theta), Omega_H_Expected, tolerance = .000001)

  ## Now with names
  colnames(Lambda) = c("A", "B", "C", "D")
  names(Omega_H_Expected) = c("A", "B", "C", "D")
  expect_equal(Omega_H(Lambda, Theta), Omega_H_Expected, tolerance = .000001)

  ## Now let's do one with no general factor
  Lambda2 <- matrix(c(.3,  0,  0,
                      .4,  0,  0,
                      .5,  0,  0,
                      0, .4,  0,
                      0, .6,  0,
                      0, .7,  0,
                      0,  0, .4,
                      0,  0, .5,
                      0,  0, .3),
                    ncol = 3, byrow = TRUE)
  Theta2 <- rep(1, nrow(Lambda2)) - rowSums(Lambda2^2)
  Omega_H_Expected2 <- c(.3654822, .5922131, .3654822)
  expect_equal(Omega_H(Lambda2, Theta2), Omega_H_Expected2,  tolerance = .000001)


})

test_that("PUC Works", {
  Lambda <- matrix(c(  0, .82, .10,   0,
                       0, .77, .35,   0,
                       0, .79, .32,   0,
                       0, .66, .39,   0,
                       0, .51,   0, .71,
                       0, .56,   0, .43,
                       0, .68,   0, .13,
                       0, .60,   0, .50,
                     .47, .83,   0,   0,
                     .27, .60,   0,   0,
                     .28, .78,   0,   0,
                     .75, .55,   0,   0),
                   ncol = 4, byrow = TRUE)
  expect_equal(PUC(Lambda), .7272727272727, tolerance = .000001)

  ## Now let's do one with no general factor
  Lambda2 <- matrix(c(.3,  0,  0,
                      .4,  0,  0,
                      .5,  0,  0,
                      0, .4,  0,
                      0, .6,  0,
                      0, .7,  0,
                      0,  0, .4,
                      0,  0, .5,
                      0,  0, .3),
                    ncol = 3, byrow = TRUE)
  expect_equal(PUC(Lambda2), NULL)

})

test_that("ARPB Works", {
  Lambda <- matrix(c(  0, .82, .10,   0,
                       0, .77, .35,   0,
                       0, .79, .32,   0,
                       0, .66, .39,   0,
                       0, .51,   0, .71,
                       0, .56,   0, .43,
                       0, .68,   0, .13,
                       0, .60,   0, .50,
                     .47, .83,   0,   0,
                     .27, .60,   0,   0,
                     .28, .78,   0,   0,
                     .75, .55,   0,   0),
                   ncol = 4, byrow = TRUE)
  UniLambda <- c(.78, .84, .82, .77, .69, .62, .69, .66, .82, .56, .74, .65)
  ## First with no names
  ARPB_Expected <- list(ARPB = c(.102578),
                        AbsRelBias = c(0.04878049, 0.09090909, 0.03797468, 0.16666667,
                                       0.35294118, 0.10714286, 0.01470588, 0.10000000,
                                       0.01204819, 0.06666667, 0.05128205, 0.18181818))
  expect_equal(ARPB(Lambda, UniLambda), ARPB_Expected, tolerance = .000001)

  # Now, with names
  rownames(Lambda) = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L")
  names(ARPB_Expected$AbsRelBias) = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L")
  expect_equal(ARPB(Lambda, UniLambda), ARPB_Expected, tolerance = .000001)

  ## Now let's do one with no general factor
  Lambda2 <- matrix(c(.3,  0,  0,
                      .4,  0,  0,
                      .5,  0,  0,
                      0, .4,  0,
                      0, .6,  0,
                      0, .7,  0,
                      0,  0, .4,
                      0,  0, .5,
                      0,  0, .3),
                    ncol = 3, byrow = TRUE)
  UniLambda2 <- c(.78, .84, .82, .77, .69, .62, .69, .66, .82)
  expect_equal(ARPB(Lambda2, UniLambda2), NULL)


})

test_that("getGen Works", {
  Lambda <- matrix(c(  0, .82, .10,   0,
                       0, .77, .35,   0,
                       0, .79, .32,   0,
                       0, .66, .39,   0,
                       0, .51,   0, .71,
                       0, .56,   0, .43,
                       0, .68,   0, .13,
                       0, .60,   0, .50,
                     .47, .83,   0,   0,
                     .27, .60,   0,   0,
                     .28, .78,   0,   0,
                     .75, .55,   0,   0),
                   ncol = 4, byrow = TRUE)
  expect_equal(getGen(Lambda), 2)


  ## Now let's do one with no general factor
  Lambda2 <- matrix(c(.3,  0,  0,
                      .4,  0,  0,
                      .5,  0,  0,
                      0, .4,  0,
                      0, .6,  0,
                      0, .7,  0,
                      0,  0, .4,
                      0,  0, .5,
                      0,  0, .3),
                    ncol = 3, byrow = TRUE)
  expect_equal(getGen(Lambda2), NULL)


})

test_that("isBifactor Works", {
  Lambda <- matrix(c(  0, .82, .10,   0,
                       0, .77, .35,   0,
                       0, .79, .32,   0,
                       0, .66, .39,   0,
                       0, .51,   0, .71,
                       0, .56,   0, .43,
                       0, .68,   0, .13,
                       0, .60,   0, .50,
                     .47, .83,   0,   0,
                     .27, .60,   0,   0,
                     .28, .78,   0,   0,
                     .75, .55,   0,   0),
                   ncol = 4, byrow = TRUE)
  expect_equal(isBifactor(Lambda), TRUE)

  ## Now let's do one with no general factor
  Lambda2 <- matrix(c(.3,  0,  0,
                      .4,  0,  0,
                      .5,  0,  0,
                      0, .4,  0,
                      0, .6,  0,
                      0, .7,  0,
                      0,  0, .4,
                      0,  0, .5,
                      0,  0, .3),
                    ncol = 3, byrow = TRUE)
  expect_equal(isBifactor(Lambda2), FALSE)

})

test_that("bifactorIndices Works", {
  ##bifactor from matrices
  Lambda <- matrix(c(  0, .82, .10,   0,
                       0, .77, .35,   0,
                       0, .79, .32,   0,
                       0, .66, .39,   0,
                       0, .51,   0, .71,
                       0, .56,   0, .43,
                       0, .68,   0, .13,
                       0, .60,   0, .50,
                     .47, .83,   0,   0,
                     .27, .60,   0,   0,
                     .28, .78,   0,   0,
                     .75, .55,   0,   0),
                   ncol = 4, byrow = TRUE)
  UniLambda <- c(.78, .84, .82, .77, .69, .62, .69, .66, .82, .56, .74, .65)
  expect_equal(bifactorIndices(Lambda, UniLambda = UniLambda), readRDS("bindices_from_matrix.rds"), tolerance = .000001)
  expect_error(bifactorIndices(Lambda, UniLambda = UniLambda, standardized = FALSE), "Not enough information is provided to compute indicator residual variances. Either provide indicator residual variances or use a standardized solution.")

  ## bifactor from lavaan

  ## bifactor from mirt

})

test_that("bifactorIndicesMplus Works", {
  cont_output <- MplusAutomation::readModels("continuous.out")
  cat_output <- MplusAutomation::readModels("categorical.out")

  expect_error(bifactorIndicesMplus(cont_output), "You must request standardized output from Mplus when standardized = TRUE")
  expect_equal(bifactorIndicesMplus(cont_output, standardized = FALSE), readRDS("cont_unst.rds"), tolerance = .000001)
  expect_equal(bifactorIndicesMplus(cat_output), readRDS("cat_stdyx.rds"), tolerance = .000001)
  expect_error(bifactorIndicesMplus(cat_output, standardized = FALSE), "Bifactor indices based on unstandardized coefficients with categorical variables is not available")

})