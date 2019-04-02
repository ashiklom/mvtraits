mprint <- function(...) {
  paste(capture.output(print(...)), collapse = "\n")
}
