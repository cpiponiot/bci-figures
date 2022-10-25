word_file <- readtext::readtext("D:/projets-autres/bci-50ha/bibliography/Bibliography for AGB chapter.docx")

citations <- strsplit(word_file$text, "ADDIN EN.REFLIST ")[[1]][2]

bib_input <- unlist(data.table::tstrsplit(citations, split = "\n"))

authors <- strsplit(bib_input, ",|[0-9]")

author1 <- gsub(" ", "", sapply(authors, data.table::first))

authors <- sapply(authors, function(x) {
  x <- gsub("and |. $", "", x[1:(which(x=="")[1] - 1)])
  x_first <- paste(c(x[1], x[2]), collapse = ",")
  if (length(x) > 2) {
    x <- c(x_first, 
           lapply(strsplit(x[3:length(x)], " "), function(i) {
             i <- paste(
               c(paste(c(i[length(i)], ","), collapse = ""), 
                 i[-c(1, length(i))]), collapse = " ")
           }))
  } else {x <- x_first}
  
  x <- paste(x, collapse = " and ")
})

year <- sapply(strsplit(bib_input, "\\."), function(x) {
  x <- as.numeric(gsub("a|b| ", "", x))
  x <- x[!is.na(x) & x > 1950 & x < 2025][1]
})

key <- paste0(author1, year)
key[duplicated(key)] <- paste0(key[duplicated(key)], "a")
## remove accents from key
key <- gsub("ü", "u", key)

full_title <- sapply(strsplit(bib_input, "\\. "), function(x) {
  x_year <- as.numeric(gsub("a|b| ", "", x))
  n_title <- which(!is.na(x_year) & x_year > 1950 & x_year < 2025)[1] + 1
  paste0(x[n_title], ".")
})

journal <- sapply(strsplit(bib_input, "\\."), function(x) {
  x_year <- as.numeric(gsub("a|b| ", "", x))
  n_journal <- which(!is.na(x_year) & x_year > 1950 & x_year < 2025)[1] + 2
  x <- data.table::tstrsplit(x[n_journal], "\\, ")[[1]]
  x <- gsub("^ ", "", x)
})

journal[journal == "Philos Trans R Soc Lond B Biol Sci"] <- "Philosophical Transactions Of The Royal Society Of London Series B-Biological Sciences"

volume <- sapply(strsplit(bib_input, "\\."), function(x) {
  x_year <- as.numeric(gsub("a|b| ", "", x))
  n_journal <- which(!is.na(x_year) & x_year > 1950 & x_year < 2025)[1] + 2
  if (grepl("\\, ", x[n_journal])) {
    x <- data.table::tstrsplit(x[n_journal], "\\, ")[[2]] 
    x <- data.table::tstrsplit(x, "\\: ")[[1]] 
  } else 
    x <- ""
})
volume[volume=="ed"|volume=="TX"] <- ""

page <- sapply(strsplit(bib_input, "\\."), function(x) {
  x_year <- as.numeric(gsub("a|b| ", "", x))
  n_journal <- which(!is.na(x_year) & x_year > 1950 & x_year < 2025)[1] + 2
  if (grepl("\\, ", x[n_journal])) {
    x <- data.table::tstrsplit(x[n_journal], "\\, ")[[2]] 
    if (grepl("\\: ", x)) {
      x <- data.table::tstrsplit(x, "\\: ")[[2]]
    } else {x <- ""}
  } else {x <- ""}
})

page[key == "Muller-Landau2014"] <- "381-415"

doi <- sapply(strsplit(bib_input, "\\. "), function(x) {
  if (grepl("doi\\.org", x[length(x)])) {
    return(gsub("https://doi.org/", "", x[length(x)]))
  } else return("")
})

type <- rep("article", length(key))
type[key == "Condit1998"] <- "book"

## book chapters
type[grepl("^In ", journal) & !grepl("DRYAD", journal)] <- "incollection"

booktitle <- rep("", length(key))
booktitle[type == "incollection"] <- gsub("^In ", "", journal[type == "incollection"])
journal[type %in% c("book", "incollection")] <- ""

editor <- rep("", length(key))
editor[grep("The First 100 Years", booktitle)] <- "H. C. Muller-Landau and S. J. Wright"
editor[booktitle == "Forests and Global Change"] <- "H. C. Muller-Landau and S. J. Wright"

editor <- sapply(strsplit(bib_input, "\\, |: "), function(x) {
  x <- gsub("^ed\\. ", "", x[grep("^ed.", x)])
  if (length(x) == 0)  x <- ""
  return(x)
})

publisher <- rep("", length(key))
publisher[key == "Condit1998"] <- "Springer-Verlag, Berlin, and R. G. Landes Company"
publisher[grep("Barro Colorado Island", booktitle)] <- "Smithsonian Institution Scholarly Press"
publisher[grep("Global", booktitle)] <- "Cambridge University Press"


address <- rep("", length(key))
address[key == "Condit1998"] <- "Georgetown, TX, USA"
address[key == "Muller-Landau2014"] <- "Cambridge, England"


bib_file <- sapply(1:length(key), function(i) {
  paste0(
    c("@", type[i], "{", key[i], ",\n", 
      "author = {", authors[i], "},\n",
      "year = {", year[i], "},\n",  
      "title = {{", full_title[i], "}},\n", 
      "journal = {", journal[i], "},\n",
      "volume = {", volume[i], "},\n", 
      "pages = {", page[i], "},\n", 
      "address = {", address[i], "},\n", 
      "publisher = {", publisher[i], "},\n", 
      "editor = {", editor[i], "},\n", 
      "booktitle = {", booktitle[i], "},\n", 
      "doi = {", doi[i], "}\n",
      "}"
    ), collapse = "")
})


bib_file <- paste(bib_file, collapse = "\n\n")

### add some refs

martinezcano <- "\n@article{MartinezCano2019, 
author = {{Mart{\'{i}}nez Cano}, Isabel and Muller-Landau, Helene C. and Wright, S. Joseph and Bohlman, Stephanie A. and Pacala, Stephen W.},
doi = {10.5194/bg-16-847-2019},
journal = {Biogeosciences},
number = {4},
pages = {847--862},
title = {{Tropical tree height and crown allometries for the Barro Colorado Nature Monument, Panama: a comparison of alternative hierarchical models incorporating interspecific variation in relation to life history traits}},
volume = {16},
year = {2019}
}
"

chave2005 <- "\n@article{Chave2005,
author = {Chave, J. and Andalo, C. and Brown, S. and Cairns, M. a. and Chambers, J. Q. and Eamus, D. and F{\"{o}}lster, H. and Fromard, F. and Higuchi, N. and Kira, T. and Lescure, J-P P. and Nelson, B. W. and Ogawa, H. and Puig, H. and Ri{\'{e}}ra, B. and Yamakura, T.},
doi = {10.1007/s00442-005-0100-x},
journal = {Oecologia},
number = {1},
pages = {87--99},
title = {{Tree allometry and improved estimation of carbon stocks and balance in tropical forests}},
volume = {145},
year = {2005}
}
"


Thomas2012 <-"\n@article{Thomas2012,
author = {Thomas, Sean C. and Martin, Adam R.},
doi = {10.3390/f3020332},
journal = {Forests},
number = {2},
pages = {332--352},
title = {{Carbon content of tree tissues: A synthesis}},
volume = {3},
year = {2012}
}
"
bib_file <- paste(bib_file, chave2005, martinezcano, Thomas2012, collapse = "\n\n")


## some corrections

## accents  and special characters
bib_file <- gsub("'\\{", "\\\\'\\{", bib_file)
bib_file <- gsub("&Thinsp;", "-", bib_file)

bib_file <- gsub("doi: ","", bib_file)

## remove empty info 
for (info in c("key", "author", "year", "title", "journal", "volume", 
               "pages", "address", "publisher", "editor", "booktitle", "doi")){
  bib_file <- gsub(paste0(",\\\n", info, " = \\{\\}"),"", bib_file)
}
  
fileConn<-file("common/bibliography.bib")
writeLines(bib_file, fileConn)
close(fileConn)

