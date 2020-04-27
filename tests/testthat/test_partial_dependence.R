context("partial dependence")

library(gbex)

data_positive = data.frame(x1 = c(-1,-1,1,1),
                  x2 = c(-1,1,-1,1),
                  y = c(-3,-1,1,3))
data_negative = data.frame(x1 = c(-1,-1,1,1),
                           x2 = c(-1,1,-1,1),
                           y = c(3,1,-1,-3))

test_that("PD for stump tree with lower to left child", {
  ctrl = rpart::rpart.control(maxdepth = 1, minsplit=1, cp=0, maxcompete = 0,maxsurrogate = 0, minbucket = 1)
  tree = rpart::rpart(y~.,data=data_positive, method='anova',control=ctrl)
  update_table = data.frame(leaf = c(2,3),
                            update = c(-2,2))
  tree = list(tree = tree, update_table = update_table)
  class(tree) = "gradient_tree"

  expect_equal(PD_tree(tree,"x1",c(-1,1)),c(-2,2))
  expect_equal(PD_tree(tree,"x2",c(-1,1)),c(0,0))
})

test_that("PD for stump tree with lower to right child", {
  ctrl = rpart::rpart.control(maxdepth = 1, minsplit=1, cp=0, maxcompete = 0,maxsurrogate = 0, minbucket = 1)
  tree = rpart::rpart(y~.,data=data_negative, method='anova',control=ctrl)
  update_table = data.frame(leaf = c(2,3),
                            update = c(-2,2))
  tree = list(tree = tree, update_table = update_table)
  class(tree) = "gradient_tree"

  expect_equal(PD_tree(tree,"x1",c(-1,1)),c(2,-2))
  expect_equal(PD_tree(tree,"x2",c(-1,1)),c(0,0))
})

test_that("PD for double split tree with lower to left child", {
  ctrl = rpart::rpart.control(maxdepth = 2, minsplit=1, cp=0, maxcompete = 0,maxsurrogate = 0, minbucket = 1)
  tree = rpart::rpart(y~.,data=data_positive, method='anova',control=ctrl)
  update_table = data.frame(leaf = c(3,4,6,7),
                            update = c(-3,-1,1,3))
  tree = list(tree = tree, update_table = update_table)
  class(tree) = "gradient_tree"

  expect_equal(PD_tree(tree,"x1",c(-1,1)),c(-2,2))
  expect_equal(PD_tree(tree,"x2",c(-1,1)),c(-1,1))
})

test_that("PD for double split tree with lower to right child", {
  ctrl = rpart::rpart.control(maxdepth = 2, minsplit=1, cp=0, maxcompete = 0,maxsurrogate = 0, minbucket = 1)
  tree = rpart::rpart(y~.,data=data_negative, method='anova',control=ctrl)
  update_table = data.frame(leaf = c(3,4,6,7),
                            update = c(-3,-1,1,3))
  tree = list(tree = tree, update_table = update_table)
  class(tree) = "gradient_tree"

  expect_equal(PD_tree(tree,"x1",c(-1,1)),c(2,-2))
  expect_equal(PD_tree(tree,"x2",c(-1,1)),c(1,-1))
})

