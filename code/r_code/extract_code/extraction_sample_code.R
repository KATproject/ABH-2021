# Illustration code to extract search index 
# from Google Trends API
#################################################

options(stringsAsFactors = F)

## Required packages

# Download package-managing package 'pacman'
if (!require("pacman")) install.packages("pacman")

# We need to use the github version of 'gtrendsR' as CRAN packages are outdated
if (!require("devtools")) install.packages("devtools")
devtools::install_github("PMassicotte/gtrendsR", force = T)

# Load all the necessary packages
pacman::p_load(data.table, magrittr, zoo, stringr, gtrendsR)

## Directory
directory <- "/home/ltmb/Dropbox/Sebastien-Leo-Sol/time writing/code/r_code/extract_code/"

#############################################################################################
# Obtain search volume at default Google Trends API resolution  for the same queries as above
# For a search period lasting a few years, resolution will be at monthly level
#############################################################################################


## Create table containing searched keywords and corresponding URL topic codes
search_examples.dt <- data.table(search_term = c("2015", "2016", "april", "september", "valentine's day", "eid al-fitr", 
                                                 "sunscreen", "donald trump" , "cat"),
                                 url_keyword = c("%2Fg%2F11b77d3kyf", "%2Fg%2F11b77d6b6d", "%2Fm%2F0lkm", "%2Fm%2F06vkl", "%2Fm%2F018y5m", "eid al-fitr",
                                                 "%2Fm%2F01r0j5", "%2Fm%2F0cqt90", "%2Fm%2F01yrx"))


## Extract the data from Google Trends API and store in 'placebo_comparison.txt'
for(i in 1:9){
  # Decode the URL to inject into the extracting function
  search_term <- URLdecode(search_examples.dt[i, ]$url_keyword)
  # Extract the search volume from 2008 to 2018
  records <- gtrends(keyword = search_term[1] , time = "2008-01-01 2018-12-31") 
  # Store the relevant variables
  interest <- records$interest_over_time[, c("date", "hits", "keyword", "geo")] %>% as.data.table()
  # Rename the 'keyword' for convenience (search term instead of URL code)
  interest[, keyword := search_examples.dt[i, 1]]
  # Save the resulting database
  fwrite(interest, paste0(directory, "placebo_comparison.txt"), append = T)
  # Clean the intermediary objects
  rm(records, interest)
}

## Read the data we just downloaded
placebo.dt <- fread(paste0(directory, "placebo_comparison.txt"))

## Make sure the date is read as such
placebo.dt[, date :=  substr(date, 1, 10) %>% as.Date()]

## Convert 'hits' variable into numeric
placebo.dt[hits == "<1", hits := "0"]
placebo.dt[, hits := as.numeric(hits)]

###################################################################################
# Obtain daily-level search volume for the same queries as above but at daily level
###################################################################################

## Create table containing searched keywords and corresponding URL topic codes
search_examples.dt <- data.table(search_term = c("2016", "2015", "september",  "april",  "valentine's day", "diwali", 
                                                 "sunscreen", "donald trump" , "cat"),
                                 url_keyword = c("%2Fg%2F11b77d6b6d", "%2Fg%2F11b77d3kyf", "%2Fm%2F06vkl", "%2Fm%2F0lkm",  "%2Fm%2F018y5m", "%2Fm%2F019fpj",
                                                 "%2Fm%2F01r0j5", "%2Fm%2F0cqt90", "%2Fm%2F01yrx"))

## Extracting function that downloads searches at daily resolution (i.e. 270 days windows)
DailygTrendsR <- function(searched_term, begin_date = "2008-01-01", end_date = "2018-12-31", bridge_size = 5){
  # Function that downloads seach interest for target at the year level (e.g. "2016")
  # Resulting interest data is at a daily resolution
  # The function rescales search interest over various time windows
  ## Inputs:
  # searched_year: The year being searched, should be a number, e.g. 2016
  # size_tails:    The number of days outside of the searched_year to consider. For instance, size_tails = 300
  #                would mean that the output search window would include 300 before the start of searched_year
  #                and 300 days after the end of searched_year. Default is 240
  # bridge_size:   The number of overlapping days between each searching time period of 270 days. Default is 5 days.
  # geo:           The localization to use for the search. Default is "US".
  ## Output:
  # Daily search interest over the period extending from [searched_year - size_tails] to [searched_year + size_tails]
  # for searched_year. Maximum is rescaled to be 100, minimum to be zero, rounding to the next unit.
  require(stringr)
  require(zoo)
  ## Print the country and year being downloaded:
  ## URL decode the search term
  searched_term %<>% URLdecode()
  ## Print the country and year being downloaded:
  print(paste("Downloading worldwide data for the searched term", searched_term))
  
  size.period <- (as.Date(end_date) - as.Date(begin_date) + 1) %>% as.numeric()  
  # Create a full vector of dates to be used after the search results have been stitched together
  time.vector <- seq(as.IDate(begin_date), as.IDate(end_date), by = "1 day") %>% as.character()
  # We are going to create the search boundaries recursively, starting with the beginning date
  search.boundaries <- as.Date(begin_date)
  
  # We add the intermediary boundaries
  for (i in 1:15){
    search.boundaries <- c(search.boundaries, 
                           as.Date(begin_date) + i*269 - (i - 1)*(bridge_size - 1), 
                           as.Date(begin_date) + i*269 - i*(bridge_size - 1))
  } 
  
  # We now add the end date
  search.boundaries <- c(search.boundaries, as.Date(end_date))
  
  # Convert time boundary vector to dataframe
  search.boundaries <- search.boundaries %>% matrix(ncol = 2, byrow = T) %>% as.data.table()
  # Rename
  setnames(search.boundaries, c("beginning", "end"))
  # Convert to date format
  search.boundaries[, `:=` (beginning = as.Date(beginning) %>% as.character(), 
                            end = as.Date(end) %>% as.character())]
  
  # Create the matrix that will hold the search results
  results <- matrix(data = NA, nrow = 270, ncol = nrow(search.boundaries))
  
  for(i in 1: nrow(search.boundaries)){
    b.period <- search.boundaries[i, beginning]
    e.period <- search.boundaries[i, end]
    date.window <- paste(b.period, e.period, collapse = " ")
    records <- gtrends(keyword = searched_term[1], 
                       time = date.window[1]) 
    # Fill the entries of the matrix that holds all the split search results
    results[, i] <- c(records$interest_over_time$hits, 
                      rep(NA, 270 - length(records$interest_over_time$hits)))
    
  }
  # Convert columns of 'results' to numeric. Whenever "hits == <1", replace by "hits == 0"
  results %<>% as.data.frame()
  results <- data.frame(lapply(results, function(x) {gsub("<1", "0.5", x)})) %>% sapply(., as.numeric) %>% as.matrix()
  # Rescale the sub-period search volume to match that of the first sub-period
  return(results)
}

## Function that rescales and stitches together the search volume for the different search periods using 
## the overlapping days (i.e. 'bridge_size' in DailyTendsR()) 
RescaleSearches <- function(matrix_to_rescale , bridge_size = 5){
  Ncol <- ncol(matrix_to_rescale)
  for(a in Ncol:2){
    rescale.factor <- (matrix_to_rescale[1:bridge_size, a]/matrix_to_rescale[((270 - bridge_size + 1):270), (a-1)]) %>% mean(na.rm = T)
    matrix_to_rescale[ , a:Ncol] <- matrix_to_rescale[ , a:Ncol]/rescale.factor
  }
  # Stitch the rescaled sub-period search volume
  stitched.vec <- matrix_to_rescale[1:270, 1]
  for(i in 2:ncol(matrix_to_rescale)){stitched.vec <- c(stitched.vec, matrix_to_rescale[((bridge_size + 1):270), i])}
  
  # Remove the missing observations. The final vector should have the same size as the date vector for the full period
  stitched.vec <- stitched.vec[!is.na(stitched.vec)]
  
  # Rescale so that the index extends from 0 to 100
  # stitched.vec <- round(stitched.vec * 100 / max(stitched.vec), 0)
  return(stitched.vec)
}

## Function that englobes both DailyTrendsR() and RescaleSearches() defined above
PlaceboDailyFunction <- function(Nrow, begin_date = "2008-01-01", end_date = "2018-12-31", bridge_size = 5){
  search_term <- search_examples.dt$url_keyword[Nrow]
  time.vector <- seq(as.Date(begin_date), as.Date(end_date), 1)
  keyword <- search_examples.dt$search_term[Nrow]
  hits.matrix <- DailygTrendsR(searched_term = search_term)
  hits.matrix.rescaled <- RescaleSearches(hits.matrix, bridge_size = bridge_size)
  results.hits <- data.table(date = time.vector, hits = hits.matrix.rescaled, keyword = keyword)
  fwrite(results.hits, paste0(directory, "placebo_comparison_new.txt"), append = T)
  
}

## Run the function for all the 9 search terms contained in 'search_examples.dt'
for(a in 1:nrow(search_examples.dt)){
  PlaceboDailyFunction(a)
}

## Read the data we just downloaded
placebo_daily.dt <- fread(paste0(directory, "placebo_comparison_new.txt"))

## Convert the date
placebo_daily.dt[, date :=  substr(date, 1, 10) %>% as.Date()]

## Convert the search volume to numeric
placebo_daily.dt[hits == "<1", hits := "0"]
placebo_daily.dt[, hits := as.numeric(hits)]
