library(mvtraits)

set.seed(12345)

r <- 0.7
cor_group_diff <- list(
  matrix(c(1, r, r, 1), 2, 2),
  diag(2),
  matrix(c(1, -r, -r, 1), 2, 2)
)
cor_group_same <- rep(list(matrix(c(1, 0.6, 0.6, 1), 2, 2)), 3)
sdd <- diag(c(0.2, 0.2))
sigma_group_same <- lapply(cor_group_same, function(x) sdd %*% x %*% sdd)
sigma_group_diff <- lapply(cor_group_diff, function(x) sdd %*% x %*% sdd)
rand_same <- random_data_hier(sigma_group = sigma_group_same)
rand_diff <- random_data_hier(sigma_group = sigma_group_diff)

abplot <- function() {
  par(mar = c(5, 5, 1, 1), cex = 0.7, cex.axis = 1.2, cex.lab = 1.5)
  plot(0, 0, type = "n", xlab = "Trait A", ylab = "Trait B",
    xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))
}

covplot <- function(rand) {
  abplot()
  points(rand$dat_all[, 1], rand$dat_all[, 2], col = rand$groups + 1, pch = 19)
  legend("topleft", sprintf("PFT %d", 1:3), pch = 19, col = 2:4, cex = 1.5)
}

mypng <- function(filename) {
  png(filename = filename, width = 6, height = 6, units = "in", res = 300)
}
mypng("figures/conceptual_groupDiff.png")
covplot(rand_diff)
dev.off()
mypng("figures/conceptual_groupSame.png")
covplot(rand_same)
dev.off()

## With an ellipse around them
ellipse_diff <- Map(
  ellipse_axes,
  mean = rand_diff$mu_group,
  cov = rand_diff$sigma_group
) %>%
  purrr::map2(1:3, ~dplyr::mutate(.x, group = .y)) %>%
  purrr::map(tidyr::unnest)

ellipse_global <- ellipse_axes(
  mean = rand_diff$mu_global,
  cov = rand_diff$sigma_global,
  prob = 0.9
)

mypng("figures/conceptual_modelHier.png")
covplot(rand_diff)
purrr::walk(1:3,
  ~lines(
    ellipse_diff[[.]]$ellipse_x,
    ellipse_diff[[.]]$ellipse_y,
    col = ellipse_diff[[.]]$group + 1,
    lwd = 4
  )
)
purrr::map(rand_diff$mu_group,
  ~points(
    .[1], .[2],
    col = "black",
    pch = "+",
    cex = 8
  )
)
lines(ellipse_global[[1, "ellipse_x"]], ellipse_global[[1, "ellipse_y"]],
  col = "black", lty = "dashed", lwd = 2)
dev.off()

mypng("figures/conceptual_modelMulti.png")
abplot()
points(rand_diff$dat_all[, 1], rand_diff$dat_all[, 2], pch = 19)
lines(ellipse_global[[1, "ellipse_x"]], ellipse_global[[1, "ellipse_y"]],
  col = "black", lty = "solid", lwd = 2)
dev.off()
