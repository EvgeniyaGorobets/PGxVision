# System setup for development
if (sys.nframe() == 0) {
  # Initial setup:
  # Run this every time you're about to start working on the package
  install.packages(c("devtools", "roxygen2", "testthat", "knitr", "usethis"))
  devtools::has_devel()

  # To update documenation:
  devtools::document()

  # To check the package:
  devtools::test()
  devtools::check()

  # TO load the package:
  devtools::load_all()

  # To update the README:
  devtools::build_readme()
}

# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   Load Package:              'Ctrl + Shift + L'

# [END]
