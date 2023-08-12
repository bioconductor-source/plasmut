test_that("correct bayes factor computation", {

  param.list <- list(ctc=list(a=1, b=9999),
                      ctdna=list(a=1, b=9),
                      chip=list(a=1, b=9),
                      montecarlo.samples=50e3,
                      prior.weight=0.1)

   dat <- data.frame(y=c(4, 1),
                 n=c(1000, 1000),
                 analyte=c("plasma", "buffy coat"),
                 mutation="mutA",
                 sample_id="id1")

  output <- importance_sampler(dat, param.list)


  expect_equal(output$bayesfactor$bayesfactor, 0.0410561, tolerance = 0.05)


})