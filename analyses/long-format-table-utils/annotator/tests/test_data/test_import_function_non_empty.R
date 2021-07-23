multi_defined <- function() {
  return(1)
}

multi_defined <- function() {
  return(2)
}
# Test comment
multi_defined_multi_ways1 <- function() {
  return(1)
}

# Test comment
assign("multi_defined_multi_ways1", function() { return(1) })
`<-`(multi_defined_multi_ways2, function() { return(1) })
# Test comment
assign("multi_defined_multi_ways2", function() { return(1) })

# Test comment
# single_lhs_ld_rhs_defined_ret1 <- function() {
#   return(1)
# }
single_lhs_ld_rhs_defined_ret1 <- function() {
  return(1)
}

assign("single_assign_defined_ret2", function() { return(2) })
`<-`(single_ld_lhs_rhs_defined_ret3, function() { return(3) })

single_lhs_ld_rhs_defined_ret2 <- function() {
  nested_single_lhs_ld_rhs_defined_ret1 <- function() {
    return(1)
  }
  return(nested_single_lhs_ld_rhs_defined_ret1() + 1)
}
{
  nested_single_lhs_ld_rhs_defined_ret2 <- function() {
    return(2)
  }
}
nested_single_lhs_ld_rhs_defined_ret3 <- {
  function() {
    return(3)
  }
}
# test comment
# test comment
# test comment
single_lhs_ld_rhs_defined_ret3 <- function() {
  assign("nested_assign_defined_ret3", function() { return(3) })
  nested_assign_defined_ret3()
}
# test comment
# test comment

single_lhs_ld_rhs_defined_ret4 <- function() {
  `<-`("nested_single_ld_lhs_rhs_defined_ret4", function() { return(4) })
  return(nested_single_ld_lhs_rhs_defined_ret4())
}
# test comment
# test comment

# test comment
# test comment <-
single_lhs_ld_rhs_defined_ret5 <- function(
  nested_single_ld_lhs_rhs_defined_ret5 = 1) {

  nested_single_ld_lhs_rhs_defined_ret5 <- function() { return(5) }
  return(nested_single_ld_lhs_rhs_defined_ret5())
}

single_lhs_ld_rhs_defined_ret6 <- function() {
  return(single_lhs_ld_rhs_defined_ret5() + 1)
}
