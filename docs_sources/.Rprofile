# .Rprofile

## To use tinytext with quarto
## quarto install tinytex --update-path

# This script is executed when R starts in this project directory.
# It's used here to activate the renv environment.

# Define the relative path to the renv activation script.
# Based on your information, 'activate.R' is in '../renv/activate.R'
# relative to your project root (where this .Rprofile and _quarto.yml are).
renv_activate_path <- "../renv/activate.R"

# Check if the renv activation script exists at the specified path
if (file.exists(renv_activate_path)) {
  # If it exists, source it to activate the renv environment
  source(renv_activate_path)
  # Provide a confirmation message in the console
  message("Successfully sourced renv activation script: ", normalizePath(renv_activate_path))
  
  # Optional: You can add other project-specific R startup commands here.
  # For example, setting global knitr options:
  # renv::load()
  # if (requireNamespace("knitr", quietly = TRUE)) {
  #   knitr::opts_chunk$set(
  #     echo = TRUE,
  #     fig.width = 7,
  #     fig.height = 5,
  #     warning = FALSE,
  #     message = FALSE
  #   )
  #   message("Default knitr chunk options set.")
  # }
  
} else {
  # If the script doesn't exist, print a warning message.
  # This helps in diagnosing path issues.
  warning("renv activation script NOT found at: '", renv_activate_path,
          "'. Current working directory is: '", getwd(),
          "'. Please ensure the path in .Rprofile is correct relative to the project root.")
}

# It's good practice for .Rprofile to end with a newline character.


