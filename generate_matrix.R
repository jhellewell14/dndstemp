get_mut <- function(c1, c2){
  cc1 <- strsplit(x = as.character(c1), split = "")[[1]]
  cc2 <- strsplit(x = as.character(c2), split = "")[[1]]
  
  n_mut <- sum(cc1 != cc2)
  if(n_mut > 1){
    return("0")
  }
  if(n_mut == 0){
    return("diag")
  }
  
  mut <- c(cc1[cc1 != cc2], cc2[cc1 != cc2])
  tv <- ifelse(all(mut %in% c("C", "T")) | all(mut %in% c("A", "G")), "transition", "transversion")
  return(tv)
}

get_mut2 <- function(c1, c2){
  cc1 <- strsplit(x = as.character(c1), split = "")[[1]]
  cc2 <- strsplit(x = as.character(c2), split = "")[[1]]
  
  n_mut <- sum(cc1 != cc2)
  if(n_mut > 1){
    return("0")
  }
  if(n_mut == 0){
    return("diag")
  }
  
  tv <- 0
  for(i in 1:3){
    tv <- ifelse(cc1[i] != cc2[i], paste0(cc1[i], cc2[i], collapse = ""), tv)  
  }
  
  return(sort(tv))
}



prot <- unlist(geneticCode)[-c(11, 12, 15)]

mdt <- as.data.table(expand.grid(tripletNames_noSTO, tripletNames_noSTO))
mdt <- mdt[, .(c1 = Var1, c2 = Var2)]
mdt$p1 <- prot[mdt$c1]
mdt$p2 <- prot[mdt$c2]

mdt[, trs := get_mut(c1, c2), c("c1", "c2")]
mdt$j <- match(mdt$c1, tripletNames_noSTO)
mdt$i <- match(mdt$c2, tripletNames_noSTO)

mdt[, mut := get_mut2(c1, c2), c("c1", "c2")]

mdt[, el := ifelse(mut == "AG", "alpha", NA)]
mdt[, el := ifelse(mut == "GA", "alpha", el)]

mdt[, el := ifelse(mut == "AC", "beta", el)]
mdt[, el := ifelse(mut == "CA", "beta", el)]

mdt[, el := ifelse(mut == "AT", "gamma", el)]
mdt[, el := ifelse(mut == "TA", "gamma", el)]

mdt[, el := ifelse(mut == "GC", "delta", el)]
mdt[, el := ifelse(mut == "CG", "delta", el)]

mdt[, el := ifelse(mut == "GT", "epsilon", el)]
mdt[, el := ifelse(mut == "TG", "epsilon", el)]

mdt[, el := ifelse(mut == "CT", "eta", el)]
mdt[, el := ifelse(mut == "TC", "eta", el)]

mdt[, el := ifelse(mut == "diag", 0, el)]

mdt[, ns := ifelse((c1 != c2) & (p1 == p2), "S", 0)]
mdt[, ns := ifelse((c1 != c2) & (p1 != p2), "NS", ns)]

mdt[, el := ifelse(ns == "NS", paste0(el, " * omega"), el)]

# mdt[, el := ifelse((trs == "transversion" & ns == "S"), "1", "0")]
# mdt[, el := ifelse(trs == "transversion" & ns == "NS", "omega", el)]
# mdt[, el := ifelse(trs == "transition" & ns == "S", "kappa", el)]
# mdt[, el := ifelse(trs == "transition" & ns == "NS", "omega * kappa", el)]
#
res <- mdt[trs != "0" & el != 0]
#
sink("output.txt")
# # M[1,2] = kappa;
for(i in 1:nrow(res)){
  cat(paste0("M[",res$i[i],",",res$j[i],"] = ",res$el[i], ";"), append = TRUE, fill = TRUE)
}
sink()

# mout <- matrix(0, 61, 61)
# for(i in 1:nrow(res)){
#   mout[res$i[i], res$j[i]] = res$ns[i]
# }
# plot(mout)
#



