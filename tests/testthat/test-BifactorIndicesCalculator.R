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
  expect_equal(ECV_SS(Lambda), ECV_SS_Expected, tolerance = .0001)

  ## Now with names
  colnames(Lambda) = c("A", "B", "C", "D")
  names(ECV_SS_Expected) = c("A", "B", "C", "D")
  expect_equal(ECV_SS(Lambda), ECV_SS_Expected, tolerance = .0001)

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
  expect_equal(ECV_SS(Lambda2), ECV_SS_Expected2, tolerance = .0001)

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
  expect_equal(ECV_SG(Lambda), ECV_SG_Expected, tolerance = .0001)
  colnames(Lambda) = c("A", "B", "C", "D")
  names(ECV_SG_Expected) = c("A", "B", "C", "D")
  ## Now with names
  expect_equal(ECV_SG(Lambda), ECV_SG_Expected, tolerance = .0001)

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
  expect_equal(ECV_GS(Lambda), ECV_GS_Expected, tolerance = .0001)
  colnames(Lambda) = c("A", "B", "C", "D")
  names(ECV_GS_Expected) = c("A", "B", "C", "D")
  ## Now with names
  expect_equal(ECV_GS(Lambda), ECV_GS_Expected, tolerance = .0001)

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
  expect_equal(IECV(Lambda), IECV_Expected, tolerance = .0001)
  rownames(Lambda) = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L")
  names(IECV_Expected) = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L")
  ## Now with names
  expect_equal(IECV(Lambda), IECV_Expected, tolerance = .0001)

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
  expect_equal(Omega_S(Lambda, Theta), Omega_S_Expected, tolerance = .0001)
  colnames(Lambda) = c("A", "B", "C", "D")
  names(Omega_S_Expected) = c("A", "B", "C", "D")
  ## Now with names
  expect_equal(Omega_S(Lambda, Theta), Omega_S_Expected, tolerance = .0001)

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
  expect_equal(Omega_S(Lambda2, Theta2), Omega_S_Expected2,  tolerance = .0001)

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
  expect_equal(Omega_H(Lambda, Theta), Omega_H_Expected, tolerance = .0001)

  ## Now with names
  colnames(Lambda) = c("A", "B", "C", "D")
  names(Omega_H_Expected) = c("A", "B", "C", "D")
  expect_equal(Omega_H(Lambda, Theta), Omega_H_Expected, tolerance = .0001)

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
  expect_equal(Omega_H(Lambda2, Theta2), Omega_H_Expected2,  tolerance = .0001)


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
  expect_equal(PUC(Lambda), .7272727272727, tolerance = .0001)

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
  expect_equal(ARPB(Lambda, UniLambda), ARPB_Expected, tolerance = .0001)

  # Now, with names
  rownames(Lambda) = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L")
  names(ARPB_Expected$AbsRelBias) = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L")
  expect_equal(ARPB(Lambda, UniLambda), ARPB_Expected, tolerance = .0001)

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

test_that("FD Works", {
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
  Phi <- diag(nrow = 4)
  FD_Expected <- c(0.8764811, 0.9500041, 0.6280990, 0.8438307)
  ## First with no names
  expect_equal(FD(Lambda, Phi), FD_Expected, tolerance = .0001)
  colnames(Lambda) = c("A", "B", "C", "D")
  names(FD_Expected) = c("A", "B", "C", "D")
  ## Now with names
  expect_equal(FD(Lambda, Phi), FD_Expected, tolerance = .0001)

  ## Now let's do one with correlated traits
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
  Phi2 <- matrix(c(1, .3, .4, .3, 1, .5, .4, .5, 1), nrow = 3)
  FD_Expected_2 <- c(0.6489614, 0.8052508, 0.6808390)
  expect_equal(FD(Lambda2, Phi2), FD_Expected_2, tolerance = .0001)

})

test_that("H Works", {
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
  H_Expected <- c(0.6340948, 0.9282446, 0.3070802, 0.6144805)
  ## First with no names
  expect_equal(H(Lambda), H_Expected, tolerance = .0001)
  colnames(Lambda) = c("A", "B", "C", "D")
  names(H_Expected) = c("A", "B", "C", "D")
  ## Now with names
  expect_equal(H(Lambda), H_Expected, tolerance = .0001)

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
  H_Expected_2 <- c(0.3837472, 0.6315076, 0.3837472)
  expect_equal(H(Lambda2), H_Expected_2, tolerance = .0001)

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
  expect_equal(bifactorIndices(Lambda, UniLambda = UniLambda), readRDS("bindices_from_matrix.rds"), tolerance = .0001)
  expect_error(bifactorIndices(Lambda, UniLambda = UniLambda, standardized = FALSE), "Not enough information is provided to compute indicator residual variances. Either provide indicator residual variances or use a standardized solution.")

  # correlated traits example
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
  Phi2 <- matrix(c(1, .5, .4, .5, 1, .6, .4, .6, 1), nrow = 3)
  expect_equal(bifactorIndices(Lambda2, Phi = Phi2), readRDS("Lambda2_indices.rds"),  tolerance = .0001)

  ## bifactor from lavaan
  SRS_UnidimensionalModel <-
    "SRS =~ SRS_1  + SRS_2  + SRS_3  + SRS_4  + SRS_5  +
            SRS_6  + SRS_7  + SRS_8  + SRS_9  + SRS_10 +
            SRS_11 + SRS_12 + SRS_13 + SRS_14 + SRS_15 +
            SRS_16 + SRS_17 + SRS_18 + SRS_19 + SRS_20"

SRS_BifactorModel <-
   "SRS =~ SRS_1  + SRS_2  + SRS_3  + SRS_4  + SRS_5  +
           SRS_6  + SRS_7  + SRS_8  + SRS_9  + SRS_10 +
           SRS_11 + SRS_12 + SRS_13 + SRS_14 + SRS_15 +
           SRS_16 + SRS_17 + SRS_18 + SRS_19 + SRS_20
  Function     =~ SRS_5  + SRS_9  + SRS_12 + SRS_15 + SRS_18
  Pain         =~ SRS_1  + SRS_2  + SRS_8  + SRS_11 + SRS_17
  SelfImage    =~ SRS_4  + SRS_6  + SRS_10 + SRS_14 + SRS_19
  MentalHealth =~ SRS_3  + SRS_7  + SRS_13 + SRS_16 + SRS_20"

  SRS_Unidimensional <- lavaan::cfa(SRS_UnidimensionalModel,
                                  SRS_data,
                                  ordered = paste0("SRS_", 1:20),
                                  orthogonal = TRUE)

  SRS_bifactor <- lavaan::cfa(SRS_BifactorModel,
                            SRS_data,
                            ordered = paste0("SRS_", 1:20),
                            orthogonal = TRUE)

  expect_equal(bifactorIndices(SRS_bifactor, UniLambda = SRS_Unidimensional), readRDS("lav_indices.rds"), tolerance = .0001)

  ## Two tier from lavaan
  MTMM_model <- "
  Trait1  =~ T1M1_1+T1M1_2+T1M1_3+T1M2_1+T1M2_2+T1M2_3+T1M3_1+T1M3_2+T1M1_3
  Trait2  =~ T2M1_1+T2M1_2+T2M1_3+T2M2_1+T2M2_2+T2M2_3+T2M3_1+T2M3_2+T2M1_3
  Trait3  =~ T3M1_1+T3M1_2+T3M1_3+T3M2_1+T3M2_2+T3M2_3+T3M3_1+T3M3_2+T3M1_3
  Method1 =~ T1M1_1+T1M1_2+T1M1_3+T2M1_1+T2M1_2+T2M1_3+T3M1_1+T3M1_2+T3M1_3
  Method2 =~ T1M2_1+T1M2_2+T1M2_3+T2M2_1+T2M2_2+T2M2_3+T3M2_1+T3M2_2+T3M2_3
  Method3 =~ T1M3_1+T1M3_2+T1M3_3+T2M3_1+T2M3_2+T2M3_3+T3M3_1+T3M3_2+T3M3_3

  Trait1 ~~ 0*Method1
  Trait1 ~~ 0*Method2
  Trait1 ~~ 0*Method3
  Trait2 ~~ 0*Method1
  Trait2 ~~ 0*Method2
  Trait2 ~~ 0*Method3
  Trait3 ~~ 0*Method1
  Trait3 ~~ 0*Method2
  Trait3 ~~ 0*Method3

  Method1 ~~ 0*Method2
  Method1 ~~ 0*Method3
  Method2 ~~ 0*Method3"

  MTMM_fit <- lavaan::cfa(MTMM_model, MTMM_data)
  expect_equal(bifactorIndices(MTMM_fit), readRDS("mtmm.rds"), tolerance = .0001)


  ## bifactor from mirt -- these lines commented out because they take too long for the R CMD check on CRAN
  #bi_data <- read.csv("bifactorData.csv")
  #colnames(bi_data) <- c(paste0("x", 1:24))
  #specific <- c(1, 1, 2, 1, 2, 2, 2, 1, 3, 3, 1, 1, 2, 3, 2, 2, 1, 2, 3, 3, 3, 1, 3, 3)
  #bi_fit_mirt <- mirt::bfactor(bi_data, specific)
  #expect_equal(bifactorIndices(bi_fit_mirt), readRDS("mirt_indices.rds"), tolerance = .0001)

})

test_that("bifactorIndices_expl Works", {
  library(psych)
  SRS_BEFA <- fa(SRS_data, nfactors = 5, rotate = "bifactor")

  ItemsBySF = list(MR4 = paste0("SRS_", c(5, 9, 12, 15, 18)),
                   MR2 = paste0("SRS_", c(1, 2, 8, 11, 17)),
                   MR3 = paste0("SRS_", c(4, 6, 10, 14, 19)),
                   MR5 = paste0("SRS_", c(3, 7, 13, 16, 20)))

  expect_equal(bifactorIndices_expl(SRS_BEFA), readRDS("exploratory_bindices_SRS.rds"), tolerance = .0001)
  expect_equal(bifactorIndices_expl(SRS_BEFA, ItemsBySF), readRDS("exploratory_bindices_SRS_fixed.rds"), tolerance = .0001)

})

test_that("bifactorIndicesMplus Works", {
  cont_output <- MplusAutomation::readModels("continuous.out")
  cat_output <- MplusAutomation::readModels("categorical.out")
  cont_output_facvar <- MplusAutomation::readModels("bifactor_continuous_wrongfacvar.out")

  expect_error(bifactorIndicesMplus(cont_output), "You must request STDYX output for computing bifactor indices based on standardized coefficients.")
  expect_equal(bifactorIndicesMplus(cont_output, standardized = FALSE), readRDS("cont_unst.rds"), tolerance = .0001)
  expect_equal(bifactorIndicesMplus(cat_output), readRDS("cat_stdyx.rds"), tolerance = .0001)
  expect_message(bifactorIndicesMplus(cat_output, standardized = FALSE), "Only bifactor indices based on standardized coefficients make sense for categorical indicators.")
  expect_error(bifactorIndicesMplus(cont_output_facvar, standardized = FALSE), NA)
})

test_that("bifactorIndicesMplus_expl Works", {
  efa_output <- MplusAutomation::readModels("bifactor_efa.out")

  expect_error(bifactorIndicesMplus_expl(efa_output), "MplusAutomation does not support EFA output yet, but should soon!")

})

test_that("bifactorIndicesMplus_ESEM Works", {
  std_output <- MplusAutomation::readModels("bifactor_esem.out")
  nostd_output <- MplusAutomation::readModels("bifactor_esem_nostd.out")

  expect_equal(bifactorIndicesMplus_ESEM(std_output), readRDS("ESEM.rds"), tolerance = .0001)
  expect_error(bifactorIndicesMplus_ESEM(nostd_output), "You must request standardized output from Mplus when standardized = TRUE")

})
