testthat::context("Testing the identifaction of functional modules")

testthat::test_that("Testing on example data set",{
    data(qmdiab.data)
    data(qmdiab.annos)
    data(qmdiab.phenos)
    
    data <- qmdiab.data[, 1:75]
    annotations <- qmdiab.annos[1:75]
    
    net.graph <- generateNetwork(data = data, annotations = annotations)
    
    mods <- identifyModules(
           graph = net.graph, data = data, annotations =
             annotations, phenotype = qmdiab.phenos$T2D, alpha = 0.05)
    
    
    testthat::expect_equal(unname(unlist(mods$modules)),
                           c(1, 7.276806e-09, 5.861762e-01, 5.457604e-07), 
                           tolerance = .002)
    testthat::expect_equal(mods$nodes[moduleID==1, name], 
                           c("P::2-hydroxybutyrate (AHB)","P::3-hydroxyisobutyrate"))
    testthat::expect_equal(mods$nodes[moduleID==1, score.after.adding], 
                           c(2.985106e-08, 7.276806e-09), 
                           tolerance = .002)
})