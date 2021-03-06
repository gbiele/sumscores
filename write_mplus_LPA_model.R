#' Estimate parameters for profiles for a specific solution (requires purchasing and installing MPlus to use)
#' @details Creates an mplus model (.inp) and associated data file (.dat)
#' @param data_filename name of data file to prepare; defaults to d.dat
#' @param script_filename name of script to prepare; defaults to i.inp
#' @param output_filename name of the output; defaults to o.out
#' @param savedata_filename name of the output for the save data (with the original data conditional probabilities); defaults to o-mod.out
#' @param the_title title of the model; defaults to test
#' @param starts number of initial stage starts and number of final stage optimizations; defaults to c(20, 4); can be set to be more conservative to c(500, 50)
#' @param m_iterations number of iterations for the EM algorithm; defaults to 500
#' @param st_iterations the number of initial stage iterations; defaults to 10; can be set more to be more conservative to 50
#' @param convergence_criterion convergence criterion for the Quasi-Newton algorithm for continuous outcomes; defaults to 1E-6 (.000001); can be set more conservatively to 1E-7 (.0000001)
#' @param remove_tmp_files whether to remove data, script, and output files; defaults to TRUE
#' @param print_input_file whether to print the input file to the console
#' @param return_save_data whether to return the save data (with the original data and the posterior probabilities for the classes and the class assignment) as a data.frame along with the MPlus output; defaults to TRUE
#' @param optseed random seed for analysis
#' @param include_LMR whether to include the Lo-Mendell-Rubin likelihood-ratio test; defaults to TRUE
#' @param include_BLRT whether to include the bootstrapped LRT; defaults to FALSE because of the time this takes to run
#' @param n_processors = 1
#' @inheritParams estimate_profiles
#' @import dplyr
#' @import tidyr
#' @importFrom tibble tibble
#' @examples
#' \dontrun{
#' m <- estimate_profiles_mplus(iris,
#'                             Sepal.Length, Sepal.Width, Petal.Length, Petal.Width,
#'                             n_profiles = 2,
#'                             model = 1)
#' }
#' @return either a tibble or a ggplot2 plot of the BIC values for the explored models
#' @export

write_mplus_LPA_model <- function(df,
                                    ...,
                                    n_profiles,
                                    the_title = "test",
                                    data_filename = "d.dat",
                                    script_filename = "i.inp",
                                    output_filename = "i.out",
                                    savedata_filename = "d-mod.dat",
                                    model = 1,
                                    starts = c(20, 4),
                                    m_iterations = 500,
                                    st_iterations = 10,
                                    convergence_criterion = 1E-6,
                                    remove_tmp_files = TRUE,
                                    print_input_file = FALSE,
                                    return_save_data = TRUE,
                                    optseed = NULL,
                                    n_processors = 1,
                                    include_LMR = TRUE,
                                    include_BLRT = FALSE) {
  message("Note that this and other functions that use MPlus are at the experimental stage! Please provide feedback at https://github.com/jrosen48/tidyLPA")
  
  d <- select_ancillary_functions_mplus(df, ...)
  x <- utils::capture.output(suppressWarnings(MplusAutomation::prepareMplusData(d, data_filename, inpfile = FALSE)))
  
  unquoted_variable_name <- paste0(names(d), collapse = " ")
  
  var_list <- list()
  for (i in 1:length(names(d))) {
    var_list[[i]] <- names(d)[i]
  }
  
  TITLE <- paste0("TITLE: ", the_title)
  
  DATA <- paste0("DATA: File is ", data_filename, ";")
  
  VARIABLE_line0 <- "VARIABLE:"
  VARIABLE_line1 <- paste0("Names are ", unquoted_variable_name, ";")
  VARIABLE_line2 <- paste0("Classes = c(", n_profiles, ");")
  
  ANALYSIS_line0 <- "ANALYSIS:"
  ANALYSIS_line1 <- "Type is mixture;"
  ANALYSIS_line2 <- paste0("start = ", starts[1], " ", starts[2], ";")
  ANALYSIS_line3 <- paste0("miterations = ", m_iterations, ";")
  ANALYSIS_line4 <- paste0("stiterations = ", st_iterations, ";")
  ANALYSIS_line5 <- paste0("convergence = ", convergence_criterion, ";")
  
  if (is.null(optseed)) {
    ANALYSIS_line6 <- NULL
  } else {
    ANALYSIS_line6 <- paste0("optseed = ", optseed, ";")
  }
  
  ANALYSIS_line7 <- paste0("processors = ", n_processors, ";")
  
  MODEL_overall_line000 <- paste0("! model specified is: ", model)
  MODEL_overall_line00 <- paste0("MODEL:")
  MODEL_overall_line0 <- paste0("%overall%")
  MODEL_overall_line1 <- paste0("[", unquoted_variable_name, "];")
  MODEL_overall_line2 <- paste0(unquoted_variable_name, ";")
  
  if (include_LMR == TRUE) {
    OUTPUT_line0 <- "OUTPUT: TECH1 TECH11;"
    if (include_BLRT == TRUE) {
      OUTPUT_line0 <- "OUTPUT: TECH1 TECH11 TECH14;"
    }
  } else {
    OUTPUT_line0 <- "OUTPUT: TECH1;"
    if (include_BLRT == TRUE) {
      OUTPUT_line0 <- "OUTPUT: TECH1  TECH14;"
    }
  }
  
  SAVEDATA_line0 <- paste0("SAVEDATA: File is ", savedata_filename, ";")
  SAVEDATA_line1 <- "SAVE = CPROBABILITIES;"
  
  if (model == 1) {
    overall_collector <- list()
    for (j in 1:length(var_list)) {
      for (k in j:length(var_list)) {
        if (var_list[[j]] != var_list[[k]]) {
          the_index <- length(overall_collector)
          overall_collector[[the_index + 1]] <- paste0(var_list[[j]], " WITH ", var_list[[k]], "@0;")
        }
      }
    }
    the_index <- 0
    class_collector <- list()
    for (i in 1:n_profiles) {
      if (the_index != 0) {
        the_index <- the_index + 1
      }
      class_collector[[the_index + 1]] <- paste0("%c#", i, "%")
      class_collector[[the_index + 2]] <- paste0("[", unquoted_variable_name, "];")
      class_collector[[the_index + 3]] <- paste0(unquoted_variable_name, "(", 1, "-", length(var_list), ");")
      for (j in 1:length(var_list)) {
        for (k in j:length(var_list)) {
          if (var_list[[j]] != var_list[[k]]) {
            the_index <- length(class_collector)
            class_collector[[the_index + 1]] <- paste0(var_list[[j]], " WITH ", var_list[[k]], "@0;")
          }
        }
      }
    }
  } else if (model == 2) {
    overall_collector <- list()
    for (j in 1:length(var_list)) {
      for (k in j:length(var_list)) {
        if (var_list[[j]] != var_list[[k]]) {
          the_index <- length(overall_collector)
          overall_collector[[the_index + 1]] <- paste0(var_list[[j]], " WITH ", var_list[[k]], ";")
        }
      }
    }
    the_index <- 0
    class_collector <- list()
    for (i in 1:n_profiles) {
      if (the_index != 0) {
        the_index <- the_index + 1
      }
      class_collector[[the_index + 1]] <- paste0("%c#", i, "%")
      class_collector[[the_index + 2]] <- paste0("[", unquoted_variable_name, "];")
      class_collector[[the_index + 3]] <- paste0(unquoted_variable_name, "(", 1, "-", length(var_list), ");")
      
      temp_index <- 0
      for (j in 1:length(var_list)) {
        for (k in j:length(var_list)) {
          if (var_list[[j]] != var_list[[k]]) {
            the_index <- length(class_collector)
            class_collector[[the_index + 1]] <- paste0(var_list[[j]], " WITH ", var_list[[k]], "(", length(var_list) + temp_index + 1, ");")
            temp_index <- (temp_index + 1)
          }
        }
      }
    }
  } else if (model == 3) {
    overall_collector <- list()
    for (j in 1:length(var_list)) {
      for (k in j:length(var_list)) {
        if (var_list[[j]] != var_list[[k]]) {
          the_index <- length(overall_collector)
          overall_collector[[the_index + 1]] <- paste0(var_list[[j]], " WITH ", var_list[[k]], "@0;")
        }
      }
    }
    the_index <- 0
    class_collector <- list()
    for (i in 1:n_profiles) {
      if (the_index != 0) {
        the_index <- the_index + 1
      }
      class_collector[[the_index + 1]] <- paste0("%c#", i, "%")
      class_collector[[the_index + 2]] <- paste0("[", unquoted_variable_name, "];")
      class_collector[[the_index + 3]] <- paste0(unquoted_variable_name, ";")
      for (j in 1:length(var_list)) {
        for (k in j:length(var_list)) {
          if (var_list[[j]] != var_list[[k]]) {
            the_index <- length(class_collector)
            class_collector[[the_index + 1]] <- paste0(var_list[[j]], " WITH ", var_list[[k]], "@0;")
          }
        }
      }
    }
  } else if (model == 4) {
    overall_collector <- list()
    for (j in 1:length(var_list)) {
      for (k in j:length(var_list)) {
        if (var_list[[j]] != var_list[[k]]) {
          the_index <- length(overall_collector)
          overall_collector[[the_index + 1]] <- paste0(var_list[[j]], " WITH ", var_list[[k]], ";")
        }
      }
    }
    the_index <- 0
    class_collector <- list()
    for (i in 1:n_profiles) {
      if (the_index != 0) {
        the_index <- the_index + 1
      }
      class_collector[[the_index + 1]] <- paste0("%c#", i, "%")
      class_collector[[the_index + 2]] <- paste0("[", unquoted_variable_name, "];")
      class_collector[[the_index + 3]] <- paste0(unquoted_variable_name, ";")
      
      temp_index <- 0
      for (j in 1:length(var_list)) {
        for (k in j:length(var_list)) {
          if (var_list[[j]] != var_list[[k]]) {
            the_index <- length(class_collector)
            class_collector[[the_index + 1]] <- paste0(var_list[[j]], " WITH ", var_list[[k]], "(", temp_index + 1, ");")
            temp_index <- (temp_index + 1)
          }
        }
      }
    }
  } else if (model == 5) {
    overall_collector <- list()
    for (j in 1:length(var_list)) {
      for (k in j:length(var_list)) {
        if (var_list[[j]] != var_list[[k]]) {
          the_index <- length(overall_collector)
          overall_collector[[the_index + 1]] <- paste0(var_list[[j]], " WITH ", var_list[[k]], ";")
        }
      }
    }
    the_index <- 0
    class_collector <- list()
    for (i in 1:n_profiles) {
      if (the_index != 0) {
        the_index <- the_index + 1
      }
      class_collector[[the_index + 1]] <- paste0("%c#", i, "%")
      class_collector[[the_index + 2]] <- paste0("[", unquoted_variable_name, "];")
      class_collector[[the_index + 3]] <- paste0(unquoted_variable_name, "(", 1, "-", length(var_list), ");")
      
      temp_index <- 0
      for (j in 1:length(var_list)) {
        for (k in j:length(var_list)) {
          if (var_list[[j]] != var_list[[k]]) {
            the_index <- length(class_collector)
            class_collector[[the_index + 1]] <- paste0(var_list[[j]], " WITH ", var_list[[k]], ";")
            temp_index <- (temp_index + 1)
          }
        }
      }
    }
  } else if (model == 6) {
    overall_collector <- list()
    for (j in 1:length(var_list)) {
      for (k in j:length(var_list)) {
        if (var_list[[j]] != var_list[[k]]) {
          the_index <- length(overall_collector)
          overall_collector[[the_index + 1]] <- paste0(var_list[[j]], " WITH ", var_list[[k]], ";")
        }
      }
    }
    the_index <- 0
    class_collector <- list()
    for (i in 1:n_profiles) {
      if (the_index != 0) {
        the_index <- the_index + 1
      }
      class_collector[[the_index + 1]] <- paste0("%c#", i, "%")
      class_collector[[the_index + 2]] <- paste0("[", unquoted_variable_name, "];")
      class_collector[[the_index + 3]] <- paste0(unquoted_variable_name, ";")
      
      temp_index <- 0
      for (j in 1:length(var_list)) {
        for (k in j:length(var_list)) {
          if (var_list[[j]] != var_list[[k]]) {
            the_index <- length(class_collector)
            class_collector[[the_index + 1]] <- paste0(var_list[[j]], " WITH ", var_list[[k]], ";")
            temp_index <- (temp_index + 1)
          }
        }
      }
    }
  }
  
  readr::write_lines(
    c(
      TITLE,
      DATA,
      VARIABLE_line0, VARIABLE_line1, VARIABLE_line2,
      MODEL_overall_line00, MODEL_overall_line0, MODEL_overall_line1, MODEL_overall_line2,
      overall_collector,
      class_collector,
      ANALYSIS_line0, ANALYSIS_line1, ANALYSIS_line2, ANALYSIS_line3, ANALYSIS_line4, ANALYSIS_line5, ANALYSIS_line6, ANALYSIS_line7,
      OUTPUT_line0,
      SAVEDATA_line0,
      SAVEDATA_line1
    ),
    script_filename
  )
  

}

#' Extract summary statistics from an Mplus model
#' @details Extract log likelihood, BIC, and entropy statistics from an Mplus model
#' @param x an mplus model
#' @return a tibble with summary statistics
#' @export

extract_mplus_summary <- function(x) {
  x$summaries[c("LL", "BIC", "Entropy")]
}

check_list <- function(x, check) {
  x[1] == check
}

check_warnings <- function(x, check) {
  if (any(purrr::map_lgl(x$warnings, check_list, check = check))) {
    return(stringr::str_c("Warning: ", "The best loglikelihood was not replicated"))
  } else {
    return("No warning")
  }
}

check_errors <- function(x, check) {
  if (any(purrr::map_lgl(x$errors, check_list, check = check))) {
    return(stringr::str_c("Error: ", "Convergence issue"))
  } else {
    return("No error")
  }
}

select_ancillary_functions_mplus <- function(df, ...) {
  if (!is.data.frame(df)) stop("df must be a data.frame (or tibble)")
  df <- tibble::as_tibble(df)
  df_ss <- dplyr::select(df, ...)
  cases_to_keep <- stats::complete.cases(df_ss) # to use later for comparing function to index which cases to keep
  d <- df_ss[cases_to_keep, ] # removes incomplete cases
  names(d) <- stringr::str_replace(names(d), "\\.", "_")
  return(d)
}