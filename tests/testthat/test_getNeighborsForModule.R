testthat::context("Testing getting the neighbors for a module")

testthat::test_that("Test simple graph",{
    
    ig <- graph_from_data_frame(data.frame(from=c(1, 2, 3), to=c(2, 4, 5)), directed = FALSE)
    
    testthat::expect_equal(getNeighborsForModule(ig, c(1,2)), c(4))
})