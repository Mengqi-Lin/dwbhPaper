
install.packages("pdftools")
library(pdftools)
text <- pdftools::pdf_text(pdf = "GRE填空机经1000题难度分级版（第1版）.pdf")
text <- str_replace_all(text, "[\\t\\n\\r]+", " ")
text1 <- sapply(text, function(x){
              str_extract_all(x, "\\w+\\w+\\w+\\w+\\w+")[[1]]
              })
text2 <- unlist(text1, F, F)
text3 <- sort(table(text2), decreasing = T)
write.table(text3, file = "GRE1000.txt", sep = "")
text4 <- which(text3>1 & text3<6)
write.table(text4, file = "GRE1000频率2-5.txt", sep = "")


what_to_Eat <- function(x = NULL, n = 1) {
  if (is.null(x)) {
    return(sample(c("太哼", "松涛", "砖", "樊阿姨", "清真","学五"), 1, replace = T))
  } else {
    return(sample(x, size = n, replace = T))
  }
}

what_to_Eat(c("太哼", "松涛", "砖"))
