### library(metafor); library(testthat); Sys.setenv(NOT_CRAN="true"); Sys.setenv(RUN_VIS_TESTS="true")

context("Checking misc: location-scale models")

source("settings.r")

dat <- dat.bangertdrowns2004

test_that("location-scale model works correctly for an intercept-only model", {

   res1 <- rma(yi, vi, data=dat)
   res2 <- rma.mv(yi, vi, random = ~ 1 | id, data=dat, sparse=.sparse)
   res3 <- rma(yi, vi, data=dat, scale = ~ 1)
   res4 <- rma(yi, vi, data=dat, scale = res3$Z)

   expect_equivalent(res1$tau2, res2$sigma2, tolerance=.tol[["var"]])
   expect_equivalent(res1$tau2, exp(res3$alpha[1]), tolerance=.tol[["var"]])
   expect_equivalent(res1$tau2, exp(res4$alpha[1]), tolerance=.tol[["var"]])

})

test_that("location-scale model works correctly for two subgroups with different tau^2 values", {

   res1 <- rma.mv(yi, vi, data=dat, random = ~ factor(meta) | id, struct="DIAG", subset=!is.na(meta), cvvc="transf", sparse=.sparse)
   expect_warning(res2 <- rma(yi, vi, data=dat, scale = ~ meta))
   expect_warning(res3 <- rma(yi, vi, data=dat, scale = res2$Z.f))

   expect_equivalent(res1$tau2, c(exp(res2$alpha[1]), exp(res2$alpha[1] + res2$alpha[2])), tolerance=.tol[["var"]])
   expect_equivalent(res1$tau2, c(exp(res3$alpha[1]), exp(res3$alpha[1] + res3$alpha[2])), tolerance=.tol[["var"]])

   expect_warning(res4 <- rma(yi, vi, data=dat, scale = ~  0 + factor(meta)))

   expect_equivalent(unname(sqrt(diag(res1$vvc))), res4$se.alpha, tolerance=.tol[["se"]])

   expect_warning(res5 <- rma(yi, vi, data=dat, scale = ~  0 + factor(meta), link="identity"))
   expect_equivalent(res1$tau2, res5$alpha, tolerance=.tol[["var"]])

   skip_on_cran()

   conf1 <- confint(res1)
   conf5 <- confint(res5, control=list(vc.min=0, vc.max=.5))
   expect_equivalent(conf1[[1]]$random[1,], conf5[[1]]$random, tolerance=.tol[["var"]])
   expect_equivalent(conf1[[2]]$random[1,], conf5[[2]]$random, tolerance=.tol[["var"]])

})

test_that("profile() and confint() work correctly for location-scale models", {

   skip_on_cran()

   png(filename="images/test_misc_rma_ls_profile_1_test.png", res=200, width=1800, height=1600, type="cairo")

   par(mfrow=c(2,2))

   res1  <- rma(yi, vi, data=dat)
   prof1 <- profile(res1, progbar=FALSE, cline=TRUE, xlim=c(.01,.15))
   conf1 <- confint(res1, type="PL")
   abline(v=conf1$random[1,2:3], lty="dotted")

   res2  <- rma.mv(yi, vi, random = ~ 1 | id, data=dat, sparse=.sparse)
   prof2 <- profile(res2, progbar=FALSE, cline=TRUE, xlim=c(.01,.15))
   conf2 <- confint(res2)
   abline(v=conf2$random[1,2:3], lty="dotted")

   res3  <- rma(yi, vi, data=dat, scale = ~ 1)
   prof3 <- profile(res3, progbar=FALSE, cline=TRUE, xlim=log(c(.01,.15)))
   conf3 <- confint(res3)
   abline(v=conf3$random[1,2:3], lty="dotted")

   expect_equivalent(prof1$ll[c(1,20)], prof3$ll[c(1,20)], tolerance=.tol[["fit"]])
   expect_equivalent(conf1$random[1,], exp(conf3$random), tolerance=.tol[["var"]])

   res4 <- rma(yi, vi, data=dat, scale = ~ 1, link="identity")
   prof4 <- profile(res4, progbar=FALSE, cline=TRUE, xlim=c(.01,.15))
   conf4 <- confint(res4, control=list(vc.max=.2))
   abline(v=conf4$random[1,2:3], lty="dotted")

   dev.off()

   expect_true(.vistest("images/test_misc_rma_ls_profile_1_test.png", "images/test_misc_rma_ls_profile_1.png"))

   expect_equivalent(prof1$ll, prof2$ll, tolerance=.tol[["fit"]])
   expect_equivalent(conf1$random[1,], conf2$random[1,], tolerance=.tol[["var"]])

   expect_equivalent(prof1$ll, prof4$ll, tolerance=.tol[["fit"]])
   expect_equivalent(conf1$random[1,], conf4$random, tolerance=.tol[["var"]])

})

test_that("location-scale model works correctly for a continuous predictor", {

   skip_on_cran()

   res1 <- rma(yi, vi, data=dat, scale = ~ grade)
   expect_equivalent(res1$beta, 0.2220791, tolerance=.tol[["coef"]])
   expect_equivalent(res1$alpha, c(-3.10513013522415, 0.041361925354706), tolerance=.tol[["coef"]])

   res2 <- rma(yi, vi, data=dat, scale = ~ grade, link="identity")
   expect_equivalent(res2$alpha, c(0.042926535, 0.002729234), tolerance=.tol[["coef"]])
   #expect_equivalent(res1$tau2, res2$tau2, tolerance=.tol[["var"]]) # not true

   res3 <- rma.mv(yi, vi, data=dat, random = ~ sqrt(grade) | id, rho=0, struct="GEN", cvvc=TRUE, sparse=.sparse)
   expect_equivalent(c(res2$alpha), diag(res3$G), tolerance=.tol[["coef"]])
   expect_equivalent(diag(res2$M),  diag(res3$M), tolerance=.tol[["var"]])
   expect_equivalent(unname(sqrt(diag(res3$vvc))), res2$se.alpha, tolerance=.tol[["se"]])

   conf11 <- confint(res1, alpha=1)
   expect_equivalent(conf11$random, c(-3.10513, -5.25032, -1.21713), tolerance=.tol[["var"]])
   conf12 <- confint(res1, alpha=2, xlim=c(-1,1))
   expect_equivalent(conf12$random, c( 0.04136, -0.65819,  0.69562), tolerance=.tol[["var"]])

   conf21 <- confint(res2, alpha=1, control=list(vc.min=-0.4, vc.max=0.3))
   conf22 <- confint(res2, alpha=2, control=list(vc.min=-0.1, vc.max=0.05))
   conf2  <- list(conf21, conf22)
   class(conf2) <- "list.confint.rma"
   expect_equivalent(conf2[[1]]$random, c(0.04293, -0.00137, 0.23145), tolerance=.tol[["var"]])
   expect_equivalent(conf2[[2]]$random, c(0.00273, -0.04972, 0.04411), tolerance=.tol[["var"]])

   conf3 <- confint(res3)
   expect_equivalent(conf3[[1]]$random[1,], c(0.04291, 0.00000, 0.11333), tolerance=.tol[["var"]])
   expect_equivalent(conf3[[2]]$random[1,], c(0.00273, 0.00000, 0.04062), tolerance=.tol[["var"]])

   # conf2 and conf3 are not the same because in res3 the two components must
   # be >= 0 while this restriction does not apply to res2 (and when profiling
   # or getting the CIs, fixing a particular component can lead to the other
   # component becoming negative)

   png(filename="images/test_misc_rma_ls_profile_2_test.png", res=200, width=1800, height=2200, type="cairo")

   par(mfrow=c(3,2))

   profile(res1, alpha=1, progbar=FALSE, cline=TRUE)
   abline(v=conf11$random[2:3], lty="dotted")
   profile(res1, alpha=2, progbar=FALSE, cline=TRUE)
   abline(v=conf12$random[2:3], lty="dotted")

   profile(res2, alpha=1, progbar=FALSE, cline=TRUE, xlim=c(0,0.3))
   abline(v=conf2[[1]]$random[2:3], lty="dotted")
   profile(res2, alpha=2, progbar=FALSE, cline=TRUE, xlim=c(-0.1,0.05))
   abline(v=conf2[[2]]$random[2:3], lty="dotted")

   profile(res3, tau2=1, progbar=FALSE, cline=TRUE, xlim=c(0,.3))
   abline(v=conf3[[1]]$random[1,2:3], lty="dotted")
   profile(res3, tau2=2, progbar=FALSE, cline=TRUE, xlim=c(0,.05))
   abline(v=conf3[[2]]$random[1,2:3], lty="dotted")

   dev.off()

   expect_true(.vistest("images/test_misc_rma_ls_profile_2_test.png", "images/test_misc_rma_ls_profile_2.png"))

})

test_that("location-scale model works correctly for multiple predictors", {

   skip_on_cran()

   expect_warning(res1 <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni)))
   expect_equivalent(res1$beta, 0.1110317, tolerance=.tol[["coef"]])
   expect_equivalent(res1$alpha, c(-1.08826059, -0.03429344, 2.09197456, -0.28439165), tolerance=.tol[["coef"]])

   expect_warning(res2 <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni), control=list(scaleZ=FALSE)))
   expect_equivalent(res2$beta, 0.1110317, tolerance=.tol[["coef"]])
   expect_equivalent(res2$alpha, c(-1.08826210, -0.03429332, 2.09197501, -0.28439156), tolerance=.tol[["coef"]])

   out <- capture.output(print(res1))

   expect_warning(res2  <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni), control=list(optimizer="Nelder-Mead")))
   expect_warning(res3  <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni), control=list(optimizer="BFGS")))
   expect_warning(res4  <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni), control=list(optimizer="bobyqa")))
   expect_warning(res5  <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni), control=list(optimizer="nloptr")))
   expect_warning(res6  <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni), control=list(optimizer="hjk")))
   expect_warning(res7  <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni), control=list(optimizer="nmk")))
   expect_warning(res8  <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni), control=list(optimizer="mads")))
   expect_warning(res9  <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni), control=list(optimizer="ucminf")))
   expect_warning(res10 <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni), control=list(optimizer="lbfgsb3c")))
   expect_warning(res11 <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni), control=list(optimizer="subplex")))
   expect_warning(res12 <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni), control=list(optimizer="BBoptim")))
   expect_warning(res13 <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni), control=list(optimizer="Rcgmin")))
   expect_warning(res14 <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni), control=list(optimizer="Rvmmin")))

   expect_equivalent(res1$alpha,  c(-1.08826059, -0.03429344, 2.09197456, -0.28439165), tolerance=.tol[["coef"]])
   expect_equivalent(res2$alpha,  c(-1.08879415, -0.03426271, 2.09166227, -0.28432946), tolerance=.tol[["coef"]])
   expect_equivalent(res3$alpha,  c(-1.08791095, -0.03439789, 2.09179476, -0.28438389), tolerance=.tol[["coef"]])
   expect_equivalent(res4$alpha,  c(-1.08826099, -0.03429340, 2.09197460, -0.28439162), tolerance=.tol[["coef"]])
   expect_equivalent(res5$alpha,  c(-1.09036615, -0.03393392, 2.09205708, -0.28429889), tolerance=.tol[["coef"]])
   expect_equivalent(res6$alpha,  c(-1.08825599, -0.03429422, 2.09197166, -0.28439180), tolerance=.tol[["coef"]])
   expect_equivalent(res7$alpha,  c(-1.08867491, -0.03415188, 2.09213170, -0.28436838), tolerance=.tol[["coef"]])
   expect_equivalent(res8$alpha,  c(-1.08825988, -0.03429568, 2.09198084, -0.28439174), tolerance=.tol[["coef"]])
   expect_equivalent(res9$alpha,  c(-1.08826216, -0.03429383, 2.09197932, -0.28439198), tolerance=.tol[["coef"]])
   expect_equivalent(res10$alpha, c(-1.08825730, -0.03429256, 2.09197369, -0.28439170), tolerance=.tol[["coef"]])
   expect_equivalent(res11$alpha, c(-1.08826074, -0.03429341, 2.09197437, -0.28439162), tolerance=.tol[["coef"]])
   expect_equivalent(res12$alpha, c(-1.08823316, -0.03429494, 2.09194049, -0.28439102), tolerance=.tol[["coef"]])
   expect_equivalent(res13$alpha, c(-1.08826085, -0.03429338, 2.09197445, -0.28439162), tolerance=.tol[["coef"]])
   expect_equivalent(res14$alpha, c(-1.08826091, -0.03429340, 2.09197450, -0.28439161), tolerance=.tol[["coef"]])

})

test_that("permutation tests work correctly for a location-scale model", {

   skip_on_cran()

   expect_warning(res <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni)))

   set.seed(1234)
   sav <- permutest(res, iter=100, progbar=FALSE)

   out <- capture.output(print(sav))

   expect_equivalent(sav$pval, 0.01, tolerance=.tol[["pval"]])
   expect_equivalent(sav$pval.alpha, c(0.81, 0.95, 0.02, 0.04), tolerance=.tol[["coef"]])

   png(filename="images/test_misc_rma_ls_permutest_light_test.png", res=200, width=1800, height=1800, type="cairo")
   plot(sav, QS=TRUE, alpha=1:4)
   dev.off()

   expect_true(.vistest("images/test_misc_rma_ls_permutest_light_test.png", "images/test_misc_rma_ls_permutest_light.png"))

   png(filename="images/test_misc_rma_ls_permutest_dark_test.png", res=200, width=1800, height=1800, type="cairo")
   setmfopt(theme="dark")
   plot(sav, QS=TRUE, alpha=1:4)
   setmfopt(theme="default")
   dev.off()

   expect_true(.vistest("images/test_misc_rma_ls_permutest_dark_test.png", "images/test_misc_rma_ls_permutest_dark.png"))

})

test_that("predict() works correctly for location-scale models", {

   skip_on_cran()

   expect_warning(res <- rma(yi, vi, data=dat, mods = ~ meta, scale = ~ meta))
   res0 <- rma(yi, vi, data=dat, subset=meta==0)
   res1 <- rma(yi, vi, data=dat, subset=meta==1)

   pred  <- predict(res, addx=TRUE, addz=TRUE)
   pred0 <- predict(res0)
   pred1 <- predict(res1)

   expect_equivalent(pred$pred[1:2],  c(pred1$pred,  pred0$pred), tolerance=.tol[["pred"]])
   expect_equivalent(pred$se[1:2] ,   c(pred1$se,    pred0$se),   tolerance=.tol[["pred"]])
   expect_equivalent(pred$ci.lb[1:2], c(pred1$ci.lb, pred0$ci.lb), tolerance=.tol[["pred"]])
   expect_equivalent(pred$ci.ub[1:2], c(pred1$ci.ub, pred0$ci.ub), tolerance=.tol[["pred"]])
   expect_equivalent(pred$pi.lb[1:2], c(pred1$pi.lb, pred0$pi.lb), tolerance=.tol[["pred"]])
   expect_equivalent(pred$pi.ub[1:2], c(pred1$pi.ub, pred0$pi.ub), tolerance=.tol[["pred"]])

   pred <- predict(res, newmods=0:1)
   expect_equivalent(pred$pred, c(pred0$pred, pred1$pred), tolerance=.tol[["pred"]])

   pred2 <- predict(res, newmods=cbind(1,0:1))
   expect_equivalent(pred, pred2)

   pred <- predict(res, newmods=0:1, newscale=0:1)

   expect_equivalent(pred$pred,  c(pred0$pred,  pred1$pred), tolerance=.tol[["pred"]])
   expect_equivalent(pred$se ,   c(pred0$se,    pred1$se),   tolerance=.tol[["pred"]])
   expect_equivalent(pred$ci.lb, c(pred0$ci.lb, pred1$ci.lb), tolerance=.tol[["pred"]])
   expect_equivalent(pred$ci.ub, c(pred0$ci.ub, pred1$ci.ub), tolerance=.tol[["pred"]])
   expect_equivalent(pred$pi.lb, c(pred0$pi.lb, pred1$pi.lb), tolerance=.tol[["pred"]])
   expect_equivalent(pred$pi.ub, c(pred0$pi.ub, pred1$pi.ub), tolerance=.tol[["pred"]])

   pred2 <- predict(res, newmods=cbind(1,0:1), newscale=0:1)
   expect_equivalent(pred, pred2)
   pred2 <- predict(res, newmods=0:1, newscale=cbind(1,0:1))
   expect_equivalent(pred, pred2)
   pred2 <- predict(res, newmods=cbind(1,0:1), newscale=cbind(1,0:1))
   expect_equivalent(pred, pred2)

   pred <- predict(res, newscale=0:1, transf=exp)
   expect_equivalent(pred$pred, c(res0$tau2, res1$tau2), tolerance=.tol[["var"]])

   expect_warning(res <- rma(yi, vi, data=dat, mods = ~ meta, scale = ~ meta, link="identity"))

   pred <- predict(res, newscale=0:1)
   expect_equivalent(pred$pred, c(res0$tau2, res1$tau2), tolerance=.tol[["var"]])

})

test_that("anova() works correctly for location-scale models", {

   skip_on_cran()

   expect_warning(res1 <- rma(yi, vi, data=dat, mods = ~ factor(grade) + meta + sqrt(ni), scale = ~ factor(grade) + meta + sqrt(ni)))
   expect_warning(res0 <- rma(yi, vi, data=dat, mods = ~ factor(grade) + meta + sqrt(ni), scale = ~ 1))

   sav <- anova(res1, res0)
   expect_equivalent(sav$LRT,  3.146726, tolerance=.tol[["test"]])
   expect_equivalent(sav$pval, 0.6773767, tolerance=.tol[["pval"]])

   sav <- anova(res1, btt=2:4)
   expect_equivalent(sav$QM,  5.286715, tolerance=.tol[["test"]])
   expect_equivalent(sav$QMp, 0.1519668, tolerance=.tol[["pval"]])

   sav <- anova(res1, att=2:4)
   expect_equivalent(sav$QS,  2.030225, tolerance=.tol[["test"]])
   expect_equivalent(sav$QSp, 0.5661571, tolerance=.tol[["pval"]])

   expect_error(anova(res1, btt=2:4, att=2:4))

   sav <- anova(res1, X=c(0,1,-1,0,0,0))
   expect_equivalent(sav$QM,  4.463309, tolerance=.tol[["test"]])
   expect_equivalent(sav$QMp, 0.03463035, tolerance=.tol[["pval"]])
   tmp <- predict(res1, newmods=c(1,-1,0,0,0), intercept=FALSE)
   expect_equivalent(sav$Xb[1,1], tmp$pred, tolerance=.tol[["test"]])
   tmp <- predict(res1, newmods=cbind(0,1,-1,0,0,0))
   expect_equivalent(sav$Xb[1,1], tmp$pred, tolerance=.tol[["test"]])

   sav <- anova(res1, Z=c(0,1,-1,0,0,0))
   expect_equivalent(sav$QS,  0.3679934, tolerance=.tol[["test"]])
   expect_equivalent(sav$QSp, 0.5441001, tolerance=.tol[["pval"]])
   tmp <- predict(res1, newscale=c(1,-1,0,0,0), intercept=FALSE)
   expect_equivalent(sav$Za[1,1], tmp$pred, tolerance=.tol[["test"]])
   tmp <- predict(res1, newscale=cbind(0,1,-1,0,0,0))
   expect_equivalent(sav$Za[1,1], tmp$pred, tolerance=.tol[["test"]])

   expect_error(anova(res1, X=c(0,1,-1,0,0,0), Z=c(0,1,-1,0,0,0)))

})

test_that("vif() works correctly for location-scale models", {

   skip_on_cran()

   expect_warning(res <- rma(yi, vi, data=dat, scale = ~ grade + meta + sqrt(ni)))
   sav <- round(vif(res)$vifs, 4)
   expect_equivalent(sav, c(grade = 1.3087, meta = 1.06, `sqrt(ni)` = 1.2847))

})

rm(list=ls())
