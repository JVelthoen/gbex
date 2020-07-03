#' Fit gradient Tree random splits
#'
#' Estimate a gradient tree then calculate for the leafnodes the optimal gradient based on a single newton raphson step
#'
#' @param X the covariates a data.frame
#' @param rr the derivative of likelihood or powerdivergence
#' @param rr2 the second derivative of likelihood or powerdivergence
#' @param depth A value indicating the maximum depth
#' @param min_leaf_size A value indicating the minimum leafsize
#' @return A list with tree object and a named vector with for each leafnode the update value
#' @export
gradient_tree_random_split <- function(X,rr,rr2,depth,min_leaf_size,nsplit = 5){
  splits_df = data.frame(index=1,depth = 0,split_point = NA, split_var = NA,child_left = NA, child_right = NA)
  index_list = list(1:nrow(X))
  row_index = 1
  while(row_index <= nrow(splits_df)){
    if(splits_df$depth[row_index] < depth & length(index_list[[row_index]]) > 2*min_leaf_size){
      # Make a single split on node row_index
      split = find_optimal_split(X[index_list[[row_index]],],
                         rr[index_list[[row_index]]],
                         min_leaf_size,
                         nsplit)
      splits_df$split_var[row_index] = split$var
      splits_df$split_point[row_index] = split$split
      splits_df$child_left[row_index] = nrow(splits_df) + 1
      splits_df$child_right[row_index] = nrow(splits_df) + 2
      # Add the resulting nodes to splits data frame
      new_lines = data.frame(index = c(splits_df$child_left[row_index],splits_df$child_right[row_index]),
                             depth = splits_df$depth[row_index] + 1,
                             split_point = NA,
                             split_var = NA,
                             child_left = NA,
                             child_right = NA)
      splits_df = rbind(splits_df,new_lines)
      # Save indices of data in each of the nodes
      index_list[[splits_df$child_left[row_index]]] = which(X[index_list[[row_index]],split$var] < split$split)
      index_list[[splits_df$child_right[row_index]]] = which(X[index_list[[row_index]],split$var] >= split$split)
      }
    row_index = row_index + 1
  }


  predict_table = data.frame(index = splits_df$index[is.na(splits_df$split_var)])
  predict_table$prediction = sapply(index_list[predict_table$index],function(indices){
    prediction = sum(rr[indices])/sum(rr2[indices])
    return(ifelse(abs(prediction)>1,sign(prediction),prediction))
  })


  gradient_tree = list(splits=splits_df, predict_table = predict_table)
  class(gradient_tree) = "Rgradient_tree"
  gradient_tree
  return(gradient_tree)
}


#' predict with gradient Tree random splits
#'
#' Estimate a gradient tree then calculate for the leafnodes the optimal gradient based on a single newton raphson step
#'
#' @param X the covariates a data.frame
#' @param rr the derivative of likelihood or powerdivergence
#' @param rr2 the second derivative of likelihood or powerdivergence
#' @param depth A value indicating the maximum depth
#' @param min_leaf_size A value indicating the minimum leafsize
#' @return A list with tree object and a named vector with for each leafnode the update value
#' @export
predict.Rgradient_tree = function(object,newdata){
    predict_table = object$predict_table
    splits_df = object$splits

    row_index = 1
    index_list = list(1:nrow(newdata))
    while(row_index < nrow(splits_df)){
      if(!is.na(splits_df$split_var[row_index])){
        index_list[[splits_df$child_left[row_index]]] =
          which(newdata[index_list[[row_index]],splits_df$split_var[row_index]] <
                  splits_df$split_point[row_index])
        index_list[[splits_df$child_right[row_index]]] =
          which(newdata[index_list[[row_index]],splits_df$split_var[row_index]] >=
                  splits_df$split_point[row_index])
      }
      row_index = row_index + 1
    }

    leaf_loc_list = lapply(predict_table$index, function(index){
      if(length(index_list[[index]]) > 0){
        data.frame(index = index_list[[index]],leaf_index = index)
      } else{
        NULL
      }
    })
    leaf_loc = do.call("rbind",leaf_loc_list)
    leaf_loc = leaf_loc[order(leaf_loc$index),]
    prediction = predict_table$prediction[match(leaf_loc$leaf_index,predict_table$index)]
    return(prediction)
  }

#' helper function for gradien tree with random splits
#'
#' Estimate a gradient tree then calculate for the leafnodes the optimal gradient based on a single newton raphson step
#'
#' @param X the covariates a data.frame
#' @param rr the derivative of likelihood or powerdivergence
#' @param rr2 the second derivative of likelihood or powerdivergence
#' @param depth A value indicating the maximum depth
#' @param min_leaf_size A value indicating the minimum leafsize
#' @return A list with tree object and a named vector with for each leafnode the update value
#' @export
find_optimal_split = function(X,rr,min_leaf_size,nsplit){

  split_point_per_var = lapply(colnames(X),function(var){
    possible_splits = sort(X[,var])[min_leaf_size:(nrow(X)-min_leaf_size)]
    return(data.frame(var = var,
                      split = sample(sort(possible_splits),min(nsplit,length(possible_splits))),
                      stringsAsFactors=F))
  })
  split_points = do.call("rbind",split_point_per_var)
  mse_split = apply(split_points,1,function(split){
    var= split[1]
    point = split[2]
    bool_left = X[var]<point
    left_leaf_se = (mean(rr[bool_left]) - rr[bool_left])^2
    right_leaf_se = (mean(rr[!bool_left]) - rr[!bool_left])^2
    return(mean(c(left_leaf_se,right_leaf_se)))
  })
  chosen_split = split_points[which(mse_split == min(mse_split)),]
  return(chosen_split[1,])
}
