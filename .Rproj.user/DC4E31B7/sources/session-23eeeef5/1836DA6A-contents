################################################################################
#### Cluster (delete-one-cluster) jackknife covariance for compositional regression
####
#### beta_(g) = estimate with cluster g removed,  g = 1, ..., G
#### covb     = (G - 1) / G * sum_g (beta_(g) - c)(beta_(g) - c)'
#### with c = mean of the beta_(g) (center = "mean", the classic jackknife, CV3)
#### or   c = beta.hat from the full sample (center = "full", slightly larger).
#### With id = NULL every observation is its own cluster (delete-one jackknife).
#### Same types and fitters as clust.boot.
#### Tsagris Michail 11/2026
################################################################################
clust.jack <- function(y, x, type = "kl", a = 1, id = NULL, center = "mean", ncores = 1) {

  if ( type == "lm" ) {
    y <- Compositional::alr(y)
    pkg <- NULL                                   # base R only
    fitfun <- function(y, x) {
      X <- cbind(1, x)
      solve( crossprod(X), crossprod(X, y) )
    }
  } else if ( type == "kl" )  {
    pkg <- "Compositional"
    fitfun <- function(y, x)  Compositional::kl.compreg(y, x)$be
  } else if ( type == "diri" )  {
    pkg <- "Compositional"
    fitfun <- function(y, x)  Compositional::diri.reg.irls(y, x)$be
  } else if ( type == "zadr" ) {
    pkg <- "CompositionalZADR"
    fitfun <- function(y, x)  CompositionalZADR::zadr.irls(y, x)$be
  } else if ( type == "areg" ) {
    pkg <- "CompositionalSR"
    fitfun <- function(y, x)  CompositionalSR::areg(y, x, a = a)$be
  } else if ( type == "cln" ) {
    pkg <- "Compositionalcln"
    fitfun <- function(y, x)  Compositionalcln::cln.reg(y, x)$beta
  } else stop("Unknown type.")

  x <- as.matrix(x)
  n <- dim(y)[1]
  fit <- fitfun(y, x)
  k <- length(fit)
  if ( is.null(id) )  id <- 1:n                  # no clusters: delete-one jackknife
  rows <- split(seq_along(id), id)
  G <- length(rows)
  if ( G < 2 )  stop("At least two clusters are needed.")

  onedel <- function(g, yy, xx, rows, k, fitfun) {
    ind <- rows[[g]]
    tryCatch(as.vector(fitfun(yy[-ind, , drop = FALSE], xx[-ind, , drop = FALSE])),
             error = function(e) rep(NA_real_, k))
  }

  if ( ncores <= 1 ) {                           # serial
    bb <- sapply(1:G, onedel, yy = y, xx = x, rows = rows, k = k, fitfun = fitfun)
  } else {                                       # parallel
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    if ( !is.null(pkg) )                         # load the right package on every worker
      parallel::clusterCall(cl, library, pkg, character.only = TRUE)
    bb <- parallel::parSapply(cl, 1:G, onedel, yy = y, xx = x, rows = rows,
                              k = k, fitfun = fitfun)
  }

  bb <- t( matrix(bb, nrow = k) )                # G x k, row g = vec(be) without cluster g
  rn <- rownames(fit)  ;  cn <- colnames(fit)
  if ( is.null(rn) )  rn <- c("constant", paste0("x", 1:(nrow(fit) - 1)))
  rn[rn == ""] <- "constant"                     # unnamed intercept of cbind(1, x)
  if ( is.null(cn) )  cn <- paste0("Y", 1 + 1:ncol(fit))
  colnames(bb) <- paste(rep(cn, each = length(rn)), rn, sep = ":")

  ok <- complete.cases(bb)
  if ( sum(ok) < G )  warning( paste(G - sum(ok), "of the", G, "leave-one-cluster-out fits failed and were dropped.") )
  bb <- bb[ok, , drop = FALSE]
  Gs <- dim(bb)[1]

  if ( center == "mean" ) {
    cc <- colMeans(bb)
  } else if ( center == "full" ) {
    cc <- as.vector(fit)
  } else stop("center must be 'mean' or 'full'.")
  e <- t( t(bb) - cc )
  covb <- (Gs - 1) / Gs * crossprod(e)           # k x k covariance of vec(be)

  list(be = fit, covb = covb, G = Gs)
}
