clust.boot <- function(y, x, type = "kl", a = 1, id = NULL, B = 999, ncores = 1) {

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
  }

  x <- as.matrix(x)
  n <- dim(y)[1]
  fit <- fitfun(y, x)
  k <- length(fit)
  if ( is.null(id) )  id <- 1:n                  # no clusters: ordinary bootstrap
  rows <- split(seq_along(id), id)
  G <- length(rows)

  onerep <- function(b, yy, xx, rows, G, k, fitfun) {
    ind <- unlist(rows[sample.int(G, G, replace = TRUE)], use.names = FALSE)
    tryCatch(as.vector(fitfun(yy[ind, , drop = FALSE], xx[ind, , drop = FALSE])),
             error = function(e) rep(NA_real_, k))
  }

  if ( ncores <= 1 ) {                           # serial
    bb <- sapply(1:B, onerep, yy = y, xx = x, rows = rows, G = G, k = k, fitfun = fitfun)
  } else {                                       # parallel
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    if ( !is.null(pkg) )                         # load the right package on every worker
      parallel::clusterCall(cl, library, pkg, character.only = TRUE)
    bb <- parallel::parSapply(cl, 1:B, onerep, yy = y, xx = x, rows = rows,
                              G = G, k = k, fitfun = fitfun)
  }

  bb <- t(bb)
  rn <- rownames(fit)  ;  cn <- colnames(fit)
  if ( is.null(rn) )  rn <- c("constant", paste0("x", 1:(nrow(fit) - 1)))
  if ( is.null(cn) )  cn <- paste0("Y", 1 + 1:ncol(fit))
  colnames(bb) <- paste(rep(cn, each = length(rn)), rn, sep = ":")

  covb <- cov(bb, use = "complete.obs")          # k x k covariance of vec(be)
  list(be = fit, covb = covb)
}
