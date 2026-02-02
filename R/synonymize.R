#' Translate plant names to accepted names using multiple synonym sources
#'
#' This function standardizes scientific names by translating synonyms to
#' accepted names using a combination of:
#' \itemize{
#'   \item User-supplied lookup tables (LUTs)
#'   \item Built-in synonym tables (NatureServe, SEINet, USDA, WCVP)
#'   \item Optional fuzzy and exact matching via the World Checklist of Vascular Plants (rWCVP)
#' }
#'
#' The function proceeds in priority order:
#' \enumerate{
#'   \item Direct matches to the checklist
#'   \item User-supplied LUTs (in order)
#'   \item Built-in LUTs (in order given by \code{synonym_sources})
#'   \item Optional WCVP exact matches
#'   \item Optional WCVP fuzzy matches
#'   \item Optional binomial-only pass (if \code{ssp_mods = TRUE})
#' }
#'
#' A new column \code{acceptedName} is added to the data, along with
#' \code{translationSource} indicating how the name was translated.
#'
#' @param input_df A data frame containing taxon names to be standardized.
#' @param name_col Name of the column in \code{input_df} containing scientific names.
#' @param checklist Optional data frame of accepted names. If \code{NULL} or missing,
#'   the built-in NatureServe Network tracheophyte checklist is used. This may contain names which are not correct in your region.
#' @param checklist_name_col Column in \code{checklist} containing accepted names.
#' @param synonym_LUTs Optional list of user-supplied data frames with columns
#'   \code{inputName} and \code{outputName}.
#' @param synonym_sources Character vector giving the order of built-in synonym sources
#'   to apply. Any of \code{"NatureServe"}, \code{"SEINet"}, \code{"USDA"}, \code{"WCVP"}.
#' @param synonym_sources_rerun Logical. If TRUE, reload synonym sources (ex: with a new target checklist).
#' @param fuzzy Logical. If TRUE, attempt exact and fuzzy matching using \pkg{rWCVP}.
#' @param wcvp_rerun Logical. If TRUE, ignored cached \pkg{rWCVP} results and rerun. Note that you can pass saved LUTs from previous runs if desired in synonym_LUTs.
#' @param ssp_mods Logical. If TRUE, perform a second pass using only binomial names
#'   (genus + species) for unresolved infraspecific names.
#'
#' @details
#' The function never overwrites an already-resolved \code{acceptedName}.
#' Earlier sources in the priority chain take precedence over later ones.
#'
#' Built-in synonym tables are shipped with the package and loaded from
#' \code{inst/extdata}.
#'
#' @return The input data frame with two additional columns:
#' \itemize{
#'   \item \code{acceptedName} – the resolved accepted name (or NA if unresolved)
#'   \item \code{translationSource} – source used for the translation
#' }
#'
#' @seealso \code{\link[rWCVP]{wcvp_match_names}}
#'
#' @examples
#' \dontrun{
#' df <- data.frame(scientificName = c("Pinus murrayana", "Quercus alba"))
#' synonymize(df)
#' }
#'
#' @export


synonymize <- function(input_df,
                       name_col = "scientificName",
                       checklist = NA,
                       checklist_name_col = "SNAME",
                       synonym_LUTs = list(),
                       synonym_sources = c("NatureServe", "SEINet", "USDA", "WCVP"),
                       synonym_sources_rerun = FALSE,
                       fuzzy = FALSE,
                       wcvp_rerun = FALSE,
                       ssp_mods = FALSE) {


  # -----------------------------------------
  # Load built in synonym sources
  # -----------------------------------------
  get_builtin_LUTs <- function(sources, full_source = FALSE, synonym_sources_rerun = FALSE) {

    # Make sure cache exists
    if (!exists(".synon_cache", envir = .GlobalEnv)) {
      print("Creating .synon_cache")
      assign(".bulkcat_cache", new.env(parent = emptyenv()), envir = .GlobalEnv)
    }

    # Only load LUTs if not already in cache
    if (!exists("builtin_LUTs", envir = .bulkcat_cache, inherits = FALSE) | synonym_sources_rerun) {
      LUTs <- list()
      for (source_name in sources) {
        lut_path <- switch(source_name,
                           "NatureServe" = system.file("extdata", "NatureServe.csv", package = "synon"),
                           "SEINet"      = system.file("extdata", "SEINet.csv", package = "synon"),
                           "USDA"        = system.file("extdata", "USDA.csv", package = "synon"),
                           "WCVP"        = system.file("extdata", "WCVP.csv", package = "synon"),
                           NULL
        )
        if (!is.null(lut_path) && nzchar(lut_path)) {
          cat("Loading built in lookup table from: ", source_name, "\n")
          LUTs[[source_name]] <- utils::read.csv(lut_path, stringsAsFactors = FALSE)

          # save unfiltered rWCVP version for the rWCVP package
          if (source_name == "WCVP"){
            LUTs[["rWCVP"]] <- LUTs[[source_name]]
          }
          # Filter out trivial rows and keep only relevant names
          LUTs[[source_name]] <- LUTs[[source_name]] %>%
            dplyr::filter(
              inputName != outputName,
              inputName %in% checklist[[checklist_name_col]] |
                outputName %in% checklist[[checklist_name_col]]
            )
          #  Swap input/output if inputName is in checklist
          LUTs[[source_name]] <- LUTs[[source_name]] %>%
            dplyr::mutate(
              swap_flag = inputName %in% checklist[[checklist_name_col]],
              input_orig = inputName,          # store original inputName
              inputName  = dplyr::if_else(swap_flag, outputName, inputName),
              outputName = dplyr::if_else(swap_flag, input_orig, outputName)
            ) %>%
            dplyr::select(-swap_flag, -input_orig) %>%
            dplyr::select(inputName, outputName)

          #write.csv(LUTs[[source_name]], paste0("filtered_built_ins_", source_name, ".csv"), row.names=FALSE)
        }
      }
      assign("builtin_LUTs", LUTs, envir = .bulkcat_cache)
      cat("Saved builtin_LUTs to cache!\n")
    } else {
      cat("Using cached builtin_LUTs\n")
    }

    return(.bulkcat_cache$builtin_LUTs)
  }
  builtin_LUTs <- get_builtin_LUTs(synonym_sources, synonym_sources_rerun)
  # ------------------------------------------
  # Create unique_df from input and initialize acceptedName and translationSource
  # ------------------------------------------

  unique_df <- input_df %>%
    dplyr::distinct(.data[[name_col]]) %>%
    dplyr::mutate(
      acceptedName = NA_character_,
      translationSource = NA_character_
    ) %>%
    # rename the column back to the clean name
    dplyr::rename(!!name_col := .data[[name_col]])


  if (missing(checklist) || is.null(checklist)) {
    message("No checklist supplied. Using NatureServe Network Tracheophyta checklist...")
    checklist <- utils::read.csv(system.file("extdata", "NatureServe.csv", package = "BulkCAT"), stringsAsFactors = FALSE)
    checklist_name_col <- "outputName"
  }

  # ------------------------------------------
  # 2️⃣ Direct matches to checklist
  # ------------------------------------------
    print("Finding direct matches to checklist...")
    direct_idx <- unique_df[[name_col]] %in% checklist[[checklist_name_col]]
    unique_df$acceptedName[direct_idx] <- unique_df[[name_col]][direct_idx]
    unique_df$translationSource[direct_idx] <- "Direct"


  # ------------------------------------------
  # 3️⃣ Function to process one LUT
  # ------------------------------------------
    process_LUT <- function(LUT, df, checklist, source_name = "LUT", name_col = "scientificName") {

        # Filter out trivial rows and keep only relevant names
        LUT <- LUT %>%
          dplyr::filter(
            inputName != outputName,
            inputName %in% checklist[[checklist_name_col]] |
              outputName %in% checklist[[checklist_name_col]]
          )

        # 2️⃣ Swap input/output if inputName is in checklist
        LUT <- LUT %>%
          dplyr::mutate(
            swap_flag = inputName %in% checklist[[checklist_name_col]],
            input_orig = inputName,          # store original inputName
            inputName  = dplyr::if_else(swap_flag, outputName, inputName),
            outputName = dplyr::if_else(swap_flag, input_orig, outputName)
          ) %>%
          dplyr::select(-swap_flag, -input_orig)

        # 3️⃣ Keep unique inputName, prioritizing outputName in checklist
        LUT <- LUT %>%
          dplyr::mutate(in_checklist = outputName %in% checklist[[checklist_name_col]]) %>%
          dplyr::group_by(inputName) %>%
          dplyr::slice_max(in_checklist, n = 1, with_ties = FALSE) %>%
          dplyr::ungroup() %>%
          dplyr::select(-in_checklist)

        # 4️⃣ Join to main df

        df <- df %>%
          dplyr::left_join(LUT, by = setNames("inputName", name_col)) %>%
          dplyr::mutate(
            acceptedName = ifelse(is.na(acceptedName) & !is.na(outputName),
                                  outputName, acceptedName),
            translationSource = ifelse(is.na(translationSource) & !is.na(outputName),
                                       source_name, translationSource)
          ) %>%
          dplyr::select(-outputName)

      return(df)
    }

  # ------------------------------------------
  # 4️⃣ Run synonyms function
  # ------------------------------------------
  run_synonyms <- function(unique_df, name_col) {

    # 4a apply user-provided LUTs in order
    print("Applying user supplied LUTs...")
    if (length(synonym_LUTs) > 0) {
      for (i in seq_along(synonym_LUTs)) {
        LUT <- synonym_LUTs[[i]]
        source_name <- ifelse(i <= length(synonym_sources),
                              synonym_sources[i],
                              paste0("LUT", i))
        unique_df <- process_LUT(LUT, unique_df, checklist, source_name, name_col)
      }
    }

    # 4b apply built-in LUTs (already loaded)
    print("Applying built-in LUTs...")

    if (length(builtin_LUTs) > 0) {
      for (source_name in names(builtin_LUTs)) {
        if (source_name != "rWCVP") {
          print(paste0("Applying built-in LUTs for ", source_name))
          unique_df <- process_LUT(
            builtin_LUTs[[source_name]],
            unique_df,
            checklist,
            source_name,
            name_col
          )
        }
      }
    }

    return(unique_df)
  }

  # Run synonyms on full names
  print("############# Running synonyms...")
  unique_df <- run_synonyms(unique_df, name_col = name_col)

  # ------------------------------------------
  # 5️⃣ Repeat the process for binomial names only (ssp_mods)
  # ------------------------------------------
  if (ssp_mods) {
    unique_df <- unique_df %>%
      dplyr::mutate(binomialName = stringr::word(.data[[name_col]], 1, 2))

    cat("Finding direct matches for binomial names...\n")

    # Direct updates for binomial matches
    direct_idx <- is.na(unique_df$acceptedName) &
      unique_df$binomialName %in% checklist[[checklist_name_col]]

    unique_df$acceptedName[direct_idx] <- unique_df$binomialName[direct_idx]
    unique_df$translationSource[direct_idx] <- "Direct binomial"

    print("############# Running binomial synonyms...")
    unique_df <- run_synonyms(unique_df, name_col = "binomialName")
  }

###################################
  # --------------------------------------------
    # fuzzy matching using rWCVP
    # ------------------------------------------
    wcvp_run <- function(unique_df, name_col, rerun = FALSE, is_binomial = FALSE){

      if (!requireNamespace("rWCVP", quietly = TRUE)) {
        warning("rWCVP not installed; skipping WCVP translation.")
        return(NA)
      }
    dir.create("rWCVP", showWarnings = FALSE)



    query_df <- unique_df %>%
        dplyr::filter(is.na(acceptedName)) %>%
        dplyr::distinct(.data[[name_col]]) %>%
        dplyr::rename(query = .data[[name_col]])

    if (nrow(query_df) == 0) {
      return(NA)
    }

    # check cache for results from previous run, unless rerun desired.
    cache_tag <- if (is_binomial) "binomial" else "full"
    fuzzy_key <- paste0("rWCVP/wcvp_fuzzy_", cache_tag, "_", name_col)
    cache_file <- paste0(fuzzy_key, ".rds")

    if (file.exists(cache_file) && !rerun) {
      wcvp_fuzzy_LUT <- readRDS(cache_file)
      cat("Found cached rWCVP results at: ", cache_file, "loading stored results.", "\n")
    } else {
      cat("Cached rWCVP results not found or rerun selected, running rWCVP for: ", fuzzy_key)

      # match using rWCVP
      wcvp_matches <- rWCVP::wcvp_match_names(
        query_df,
        name_col = "query",
        fuzzy = TRUE
      )

      # ---------------------------------------
      # Build FUZZY LUT
      # ---------------------------------------
      # Assume builtin_LUTs[["WCVP"]] exists and is preloaded
      wcvp_fuzzy_LUT <- wcvp_matches %>%
        dplyr::left_join(
          builtin_LUTs[["rWCVP"]],
          by = c("wcvp_accepted_id" = "accepted_plant_name_id")) %>%
          dplyr::mutate(
          inputName = query,
          outputName = outputName
        ) %>%
        dplyr::select(inputName, outputName) %>%
        dplyr::filter(inputName != outputName)

      # Ensure one output per inputName, preferring checklist matches
      wcvp_fuzzy_LUT <- wcvp_fuzzy_LUT %>%
        dplyr::filter(!is.na(outputName)) %>%
        dplyr::group_by(inputName) %>%
        dplyr::arrange(dplyr::desc(outputName %in% checklist[[checklist_name_col]])) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()

      # ---------------------------------------
      # Save to cache
      # ---------------------------------------

      saveRDS(wcvp_fuzzy_LUT, cache_file)
      # ---------------------------------------
      # Also write CSVs (for review and use as custom user supplied LUTs to avoid retranslating the same names.
      # ---------------------------------------

      ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

      if (is_binomial) {
        file_name <- paste0("rWCVP/rWCVP_fuzzy_binomial_", ts, ".csv")
        utils::write.csv(wcvp_fuzzy_LUT, file_name, row.names = FALSE)
      } else {
        file_name <- paste0("rWCVP/rWCVP_fuzzy_", ts, ".csv")
        utils::write.csv(wcvp_fuzzy_LUT, file_name, row.names = FALSE)
      }
    }
    return(wcvp_fuzzy_LUT)
  }
  #######################################

  if (fuzzy)
  {
    print("Running rWCVP...")
    rWCVP_lut <- wcvp_run(unique_df, name_col = name_col, rerun = wcvp_rerun, is_binomial = FALSE)
    if (!(length(rWCVP_lut) == 1 && is.na(rWCVP_lut))) {
        unique_df <- process_LUT(rWCVP_lut, unique_df, checklist, source_name = "rWCVP", name_col = name_col)
      }

    if (ssp_mods) {
      print("Running WCVP binomial")
      rWCVP_lut <- wcvp_run(unique_df, name_col = "binomialName", rerun = wcvp_rerun, is_binomial = TRUE)
      if (!(length(rWCVP_lut) == 1 && is.na(rWCVP_lut))) {
        unique_df <- process_LUT(rWCVP_lut, unique_df, checklist, source_name = "rWCVP_binomial", name_col = "binomialName")
      }
    }
  }

  # apply translated names to the input_df
  output_df <- input_df %>% dplyr::left_join(unique_df, by = name_col)

  return(output_df)
}




