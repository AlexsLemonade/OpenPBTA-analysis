multi_defined <- function() {
  return(1)
}

multi_defined <- function() {
  return(2)
}

single_lhs_ld_rhs_defined_ret1 <- function() {
  return(1)
}


assign("single_assign_defined_ret2", function() { return(2) })

`<-`(single_ld_lhs_rhs_defined_ret3, function() { return(3) })
