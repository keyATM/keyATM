#' Run Seeded LDA
#'
#' @param data_folder_path Full path to the data folder.
#' @param seed_csv_path Full path to seed.
#' @param num_regular_topic a number of regular topic
#' @param iter_num iteration number.
#' @param show_words_num Number of words to show in the output object.
#' @param full_output Return topic assignments and indicators for all words.
#' @param seed seed for random number generators (not working well??).
#' @return Seeded LDA results.
#' @examples
#' \dontrun{seededlda(data_folder_path="data/is/stored/here/", seed_path="path/to/seed.txt")}
#'
#' @export
seededlda <- function(data_folder_path,
                      seed_csv_path,
                      num_regular_topic = 0,
                      iter_num = 10,
                      show_words_num = 10,
                      full_output = F,
                      rand_seed = 0) {
  # Check variables
  if (!iter_num >= 0) {
    stop("Number of iteration is negative. Check iter_num argument.")
  }
  if (!show_words_num >= 0) {
    stop("Number of words shown in the output is negative. Check show_words_num argument.")
  }
  if (rand_seed < 0) {
    stop("Seed is negative. Check seed argument.")
  }
  if (rand_seed == 0) {
    message("Seed is not set. Randomly set seed.")
  }

  # Modify
  nc <- nchar(data_folder_path)
  lastcharacter <- substr(data_folder_path, nc, nc)
  if (lastcharacter != "/") {
    data_folder_path <- paste0(data_folder_path, "/")
  }

  # Create CSV file
  seed_vec <- read_seed_file(seed_csv_path)
  create_seed_file(seed_vec, num_regular_topic, save_file_name = "seedwords")
  seed_full_path <- normalizePath(".", "seedwords.txt")

  # Run
  res <- RunSeededLDA(
    datafolder = data_folder_path,
    seed_path = seed_full_path,
    show_words_num = show_words_num,
    iter_num = iter_num,
    full_output = full_output,
    seed = rand_seed
  )
  class(res) <- c("seededlda", class(res))
  return(res)
}

#' Create Seed file
#'
#' @param seeds_words a vector of seed words.
#' @param num_regular_topic a number of regular topic
#'
create_seed_file <- function(seed_words,
                             num_regular_topic = 0,
                             save_file_name = "seedwords") {
  savefilepath <- paste0(save_file_name, ".txt")
  seeds_words <- c(seed_words, rep(" ", num_regular_topic))

  fileConn <- file(savefilepath)
  writeLines(seeds_words, fileConn)
  close(fileConn)
}

#' Read Seed file
#'
#' @param file_path a full path to the seed file (csv)
#' @return a vector of seeds
#'
read_seed_file <- function(file_path) {
  seed_csv <- as.matrix(read.csv(file_path,
                                 stringsAsFactors = F)) # skip first line
  seed_vec <- c()
  for (i in 1:nrow(seed_csv)) {
    seed_vec <- c(seed_vec, paste(seed_csv[i, ], collapse = " "))
  }
  return(seed_vec)
}
