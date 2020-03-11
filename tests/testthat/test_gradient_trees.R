context("gradient tree functions")

library(gbex)

X = data.frame(X1 = 1:10)
rr = 1:10
rr2 = 1:10

test_that("gradient_tree returns a complete gradient_tree object", {
  depth = 2
  min_leaf_size = 1
  grad_tree = gradient_tree(X,rr,rr2,depth,min_leaf_size)
  expect_s3_class(grad_tree,"gradient_tree")
  expect_s3_class(grad_tree$tree,"rpart")
  expect_s3_class(grad_tree$update_table,"data.frame")
})

test_that("gradient_tree fits a stump when depth = 0", {
  depth = 0
  min_leaf_size = 1
  grad_tree = gradient_tree(X,rr,rr2,depth,min_leaf_size)
  expect_equal(nrow(grad_tree$update_table),1)
  expect_equal(length(unique(grad_tree$tree$where)),1)
})


