ladd = function(L1, L2){
  n = length(L1)
  outL = L1
  for(i in 1:n){
    outL[[i]] = L1[[i]] + L2[[i]]
  }
  return(outL)
}

lsubtract = function(L1, L2){
  n = length(L1)
  outL = L1
  for(i in 1:n){
    outL[[i]] = L1[[i]] - L2[[i]]
  }
  return(outL)
}

labsdiff = function(L1, L2){
  n = length(L1)
  outL = L1
  for(i in 1:n){
    outL[[i]] = abs(L1[[i]] - L2[[i]])
  }
  return(outL)
}
