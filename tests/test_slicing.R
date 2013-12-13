# Testing the slicing functions

#library(testthat)
#source("../slicing/length_slicing_funcs.R")
#test_file("C:/Work/flr/R4Med/R/tests/test_slicing.R")
#test_dir("C:/Work/flr/R4Med/R/tests")

context("FAO proportional slicing")

# get_prop(L1, L2, Latage1, Latage2)
# Returns the proportion of L1-L2 that is inside VBL1-VBL2
test_that("get_prop calculations correct proportions for different length combinations", {
    # Six situations
    # Inside
    expect_that(get_prop(1,2,0,3), is_identical_to(1))
    # Envelopes
    expect_that(get_prop(0,3,1,2), is_identical_to(1/3))
    # Outside left
    expect_that(get_prop(0,1,2,3), is_identical_to(0))
    # Outside right
    expect_that(get_prop(2,3,1,2), is_identical_to(0))
    # Overlap left
    expect_that(get_prop(0,2,1,4), is_identical_to(0.5))
    # Overlap right
    expect_that(get_prop(3,6,2,5), is_identical_to(2/3))
    # Acceptable but unlikely
    expect_that(get_prop(0,0,1,1), is_identical_to(0))
    # Ls the wrong way round
    expect_that(get_prop(6,3,2,5), throws_error())
    expect_that(get_prop(3,6,5,2), throws_error())
    expect_that(get_prop(-1,0,2,5), throws_error())
    expect_that(get_prop(0,-1,2,5), throws_error())
    expect_that(get_prop(1,2,-1,2), throws_error())
    expect_that(get_prop(1,2,5,-1), throws_error())
})
