terminal.twigs <- function (topology, times) {
  n <- length(topology);
  return(seq(1,n)[times == max(times)]);
}

