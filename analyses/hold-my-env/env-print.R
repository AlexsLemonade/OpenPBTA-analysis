# Testing environment variables in CI
# Casey Greene (please do not merge this though)
#
# This script doesn't do anything useful

# Declare command line options
option_list <- list(
  optparse::make_option(
    c("-p", "--print"),
    type = "character",
    default = "42",
    help = "Something to print",
  )
)

# Read the arguments passed
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

print(opt$print)
