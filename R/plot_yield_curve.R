#
plot_yield_curve <- function(lh.params, z, MNPL_in){
  p1 <- pop_vs_yield(z.vec = z, lh.params = lh.params, ggp = TRUE)
  p2 <-   p1 + geom_vline(xintercept = MNPL_in, colour="grey", lty=2, lwd=1.1)
  return(p2)
}
