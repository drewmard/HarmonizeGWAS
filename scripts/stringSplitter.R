stringSplitter = function(data_input,
                          column_to_use="Compare_nr-vs-r",
                          divider="_",
                          cols_to_use=c("A1","A2"),
                          idx_to_use=c(2,3),
                          custom_idx_to_use=NA,
                          # remove_custom_idx_to_use=FALSE,
                          removeCol=FALSE,
                          idCol=FALSE) { 
  
  data_input_subset = data_input[,column_to_use]
  y = strsplit(data_input_subset,divider)
  
  if (!anyNA(custom_idx_to_use)) {
    # custom_idx_to_use_param = custom_idx_to_use
    # custom_idx_to_use = as.numeric(df1[idx,custom_idx_to_use])
    data_input[,cols_to_use[1]] = unlist(lapply(1:length(custom_idx_to_use),function(i) y[[i]][custom_idx_to_use[i]]))
    # if (remove_custom_idx_to_use) {
    #   data_input = data_input[,colnames(data_input)!=custom_idx_to_use_param]
    # }
  } else {
    for (i in 1:length(cols_to_use)) {
      data_input[,cols_to_use[i]] = unlist(lapply(y,function(x)x[[idx_to_use[i]]]))
    }
  }
  
  if (idCol) {
    data_input$id = data_input_subset
  }
  if (removeCol) {
    data_input = data_input[,colnames(data_input)!=column_to_use]
  }
  return(data_input)
}
