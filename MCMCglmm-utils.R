vif.MCMCglmm <- function (fit, intercept.columns = c(1)) {
    nF <- fit$Fixed$nfl
    v <- cov(as.matrix(fit$X[,1:nF]))
    nam <- colnames(fit$Sol[,1:nF])

    v <- v[-intercept.columns, -intercept.columns, drop = FALSE]
    nam <- nam[-intercept.columns]
    
    d <- diag(v)^0.5
    v <- diag(solve(v/(d %o% d)))
    names(v) <- nam
    v
}

kappa.MCMCglmm <- function (fit, add.intercept = TRUE, scale = TRUE, intercept.columns = c(1)) {
    nF <- fit$Fixed$nfl
    X <- fit$X[,1:nF]
    if (with.intercept & scale) {
        kappa(cBind(rep(1), scale(X[, -intercept.columns])))
    } else if (with.intercept & !scale) {
        kappa(X)
    } else if (!add.intercept & scale) {
        kappa(scale(X[,-intercept.columns]))
    } else {
        kappa(X[,-intercept.columns])
    }
}
