test_that("calculate_nca_tcell works", {
    result <- calculate_nca_tcell(time = c(-34, 3, 7, 14, 29, 56, 85, 168), conc = c(50.0, 18998.6, 29261.4, 62113.0, 10920.7, 729.9, 184.9, 561.2), dose = 1, metrics = c("Lambdaz","AucD28"))
    expect_true("AucD28" %in% names(result))
    expect_equal(result$AucD28, 874996.883)
    expect_equal(result$Lambdaz, 0.015903226011427)
})
