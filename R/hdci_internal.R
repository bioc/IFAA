# Internalized routines adapted from HDCI 1.0-2.
# HDCI is licensed under GPL-2.
# These functions are retained internally because HDCI was archived from CRAN.

# ---- bootLOPR ----
bootLOPR <- function (x, y, lambda2, B = 500, type.boot = "residual", thres = 0.5, 
    alpha = 0.05, OLS = TRUE, cv.method = "cv", nfolds = 10, 
    foldid, cv.OLS = TRUE, tau = 0, parallel = FALSE, standardize = TRUE, 
    intercept = TRUE, parallel.boot = FALSE, ncores.boot = 1, 
    ...) 
{
    x <- as.matrix(x)
    y <- as.numeric(y)
    n <- dim(x)[1]
    p <- dim(x)[2]
    if ((type.boot != "residual") & (type.boot != "paired")) {
        stop("type.boot should take value of 'residual' or 'paired'.")
    }
    if (missing(lambda2)) {
        lambda2 <- 1/n
    }
    selectset <- rep(0, p)
    Beta <- rep(0, p)
    Beta.LPR <- rep(0, p)
    globalfit <- glmnet(x, y, standardize = standardize, intercept = intercept, 
        ...)
    lambda <- globalfit$lambda
    cvfit <- escv.glmnet(x, y, lambda = lambda, nfolds = nfolds, 
        tau = tau, cv.OLS = cv.OLS, parallel = parallel, standardize = standardize, 
        intercept = intercept, ...)
    if (cv.method == "cv") {
        lambda.opt <- cvfit$lambda.cv
    }
    if (cv.method == "escv") {
        lambda.opt <- cvfit$lambda.escv
    }
    if (cv.method == "cv1se") {
        lambda.opt <- cvfit$lambda.cv1se
    }
    fitlasso <- predict(globalfit, type = "coefficients", s = lambda.opt)
    fit.value <- predict(globalfit, newx = x, s = lambda.opt)
    betalasso <- fitlasso[-1]
    Beta <- betalasso
    selectvar <- betalasso != 0
    if (OLS & sum(selectvar) > 0) {
        mls.obj <- mls(x[, selectvar, drop = FALSE], y, tau = tau, 
            standardize = standardize, intercept = intercept)
        Beta[selectvar] <- mls.obj$beta
        fit.value <- mypredict(mls.obj, newx = x[, selectvar, 
            drop = FALSE])
    }
    PR.obj <- PartRidge(x = x, y = y, lambda2 = lambda2, varset = selectvar, 
        standardize = standardize, intercept = intercept)
    Beta.LPR <- PR.obj$beta
    if (type.boot == "residual") {
        fit <- fit.value
        residual <- y - fit
        residual_center <- residual - mean(residual)
        Beta.boot <- matrix(0, nrow = B, ncol = p)
        Beta.boot.LPR <- matrix(0, nrow = B, ncol = p)
        out <- list()
        if (!parallel.boot) {
            for (i in 1:B) {
                out[[i]] <- list()
                resam <- sample(1:n, n, replace = TRUE)
                ystar <- fit + residual_center[resam]
                if (!OLS) {
                  boot.obj <- Lasso(x, ystar, lambda = lambda.opt, 
                    standardize = standardize, intercept = intercept, 
                    ...)
                }
                if (OLS) {
                  boot.obj <- LassoOLS(x, ystar, lambda = lambda.opt, 
                    standardize = standardize, intercept = intercept, 
                    ...)
                }
                out[[i]][[1]] <- boot.obj$beta
                boot.selectvar <- boot.obj$beta != 0
                boot.PR.obj <- PartRidge(x = x, y = ystar, lambda2 = lambda2, 
                  varset = boot.selectvar, standardize = standardize, 
                  intercept = intercept)
                out[[i]][[2]] <- boot.PR.obj$beta
            }
        }
        else {
            resams <- matrix(sample(1:n, n * B, replace = TRUE), 
                nrow = n)
            out <- foreach(i = 1:B) %dopar% {
                resam <- resams[, i]
                ystar <- fit + residual_center[resam]
                if (!OLS) {
                  boot.obj <- Lasso(x, ystar, lambda = lambda.opt, 
                    standardize = standardize, intercept = intercept, 
                    ...)
                }
                if (OLS) {
                  boot.obj <- LassoOLS(x, ystar, lambda = lambda.opt, 
                    standardize = standardize, intercept = intercept, 
                    ...)
                }
                beta1 <- boot.obj$beta
                boot.selectvar <- boot.obj$beta != 0
                boot.PR.obj <- PartRidge(x = x, y = ystar, lambda2 = lambda2, 
                  varset = boot.selectvar, standardize = standardize, 
                  intercept = intercept)
                beta2 <- boot.PR.obj$beta
                list(beta1 = beta1, beta2 = beta2)
            }
        }
        for (i in 1:B) {
            Beta.boot[i, ] <- out[[i]][[1]]
            Beta.boot.LPR[i, ] <- out[[i]][[2]]
            out[[i]] <- 0
        }
        out <- NULL
        interval <- ci(Beta, Beta.boot, alpha = alpha, type = "basic")
        interval.LPR <- ci(Beta.LPR, Beta.boot.LPR, alpha = alpha, 
            type = "basic2", Beta2 = Beta)
        var.thres <- apply(Beta.boot != 0, 2, mean) > thres
        interval.LOPR <- interval.LPR
        interval.LOPR[, var.thres] <- interval[, var.thres]
    }
    if (type.boot == "paired") {
        Beta.boot <- matrix(0, nrow = B, ncol = p)
        Beta.boot.LPR <- matrix(0, nrow = B, ncol = p)
        out <- list()
        if (!parallel.boot) {
            for (i in 1:B) {
                out[[i]] <- list()
                resam <- sample(1:n, n, replace = TRUE)
                rx <- x[resam, ]
                ry <- y[resam]
                if (!OLS) {
                  boot.obj <- Lasso(rx, ry, lambda = lambda.opt, 
                    standardize = standardize, intercept = intercept, 
                    ...)
                }
                if (OLS) {
                  boot.obj <- LassoOLS(rx, ry, lambda = lambda.opt, 
                    standardize = standardize, intercept = intercept, 
                    ...)
                }
                out[[i]][[1]] <- boot.obj$beta
                boot.selectvar <- boot.obj$beta != 0
                boot.PR.obj <- PartRidge(x = rx, y = ry, lambda2 = lambda2, 
                  varset = boot.selectvar, standardize = standardize, 
                  intercept = intercept)
                out[[i]][[2]] <- boot.PR.obj$beta
            }
        }
        else {
            resams <- matrix(sample(1:n, n * B, replace = TRUE), 
                nrow = n)
            out <- foreach(i = 1:B) %dopar% {
                resam <- resams[, i]
                rx <- x[resam, ]
                ry <- y[resam]
                if (!OLS) {
                  boot.obj <- Lasso(rx, ry, lambda = lambda.opt, 
                    standardize = standardize, intercept = intercept, 
                    ...)
                }
                if (OLS) {
                  boot.obj <- LassoOLS(rx, ry, lambda = lambda.opt, 
                    standardize = standardize, intercept = intercept, 
                    ...)
                }
                beta1 <- boot.obj$beta
                boot.selectvar <- boot.obj$beta != 0
                boot.PR.obj <- PartRidge(x = rx, y = ry, lambda2 = lambda2, 
                  varset = boot.selectvar, standardize = standardize, 
                  intercept = intercept)
                beta2 <- boot.PR.obj$beta
                list(beta1 = beta1, beta2 = beta2)
            }
        }
        for (i in 1:B) {
            Beta.boot[i, ] <- out[[i]][[1]]
            Beta.boot.LPR[i, ] <- out[[i]][[2]]
            out[[i]] <- 0
        }
        out <- NULL
        interval <- ci(Beta, Beta.boot, alpha = alpha, type = "quantile")
        interval.LPR <- ci(Beta.LPR, Beta.boot.LPR, alpha = alpha, 
            type = "quantile")
        var.thres <- apply(Beta.boot != 0, 2, mean) > thres
        interval.LOPR <- interval.LPR
        interval.LOPR[, var.thres] <- interval[, var.thres]
    }
    object <- list(lambda.opt = lambda.opt, Beta = Beta, Beta.LPR = Beta.LPR, 
        interval = interval, interval.LPR = interval.LPR, interval.LOPR = interval.LOPR)
    object
}

# ---- escv.glmnet ----
escv.glmnet <- function (x, y, lambda = NULL, nfolds = 10, foldid, cv.OLS = FALSE, 
    tau = 0, parallel = FALSE, standardize = TRUE, intercept = TRUE, 
    ...) 
{
    if (!is.null(lambda) && length(lambda) < 2) {
        stop("Need more than one value of lambda for escv.glmnet")
    }
    n <- nrow(x)
    p <- ncol(x)
    y <- drop(y)
    glmnet.call <- match.call(expand.dots = TRUE)
    which <- match(c("nfolds", "foldid"), names(glmnet.call), 
        F)
    if (any(which)) {
        glmnet.call <- glmnet.call[-which]
    }
    glmnet.call[[1]] <- as.name("glmnet")
    glmnet.object <- glmnet(x, y, lambda = lambda, standardize = standardize, 
        intercept = intercept, ...)
    lambda <- glmnet.object$lambda
    glmnet.object$call <- glmnet.call
    if (missing(foldid)) {
        foldid <- sample(rep(seq(nfolds), length = n))
    }
    else {
        nfolds <- max(foldid)
    }
    if (nfolds < 3) {
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    }
    out <- list()
    if (!parallel) {
        for (k in 1:nfolds) {
            test <- foldid == k
            train <- foldid != k
            obj <- glmnet(x[train, , drop = FALSE], y[train], 
                lambda = lambda, standardize = standardize, intercept = intercept, 
                ...)
            fitmat <- predict(obj, newx = x)
            predtest <- predict(obj, newx = x[test, , drop = FALSE])
            residmat <- apply((y[test] - predtest)^2, 2, mean)
            out[[k]] <- list(residmat = residmat, fitmat = fitmat)
        }
    }
    else {
        out <- foreach(k = 1:nfolds, .packages = c("glmnet")) %dopar% 
            {
                test <- foldid == k
                train <- foldid != k
                obj <- glmnet(x[train, , drop = FALSE], y[train], 
                  lambda = lambda, standardize = standardize, 
                  intercept = intercept, ...)
                fitmat <- predict(obj, newx = x)
                predtest <- predict(obj, newx = x[test, , drop = FALSE])
                residmat <- apply((y[test] - predtest)^2, 2, 
                  mean)
                list(residmat = residmat, fitmat = fitmat)
            }
    }
    residmat <- matrix(0, length(lambda), nfolds)
    residmates <- matrix(0, length(lambda), nfolds)
    fitmat <- array(0, dim = c(n, length(lambda), nfolds))
    for (k in 1:nfolds) {
        residmat[, k] <- out[[k]]$residmat
        fitmat[, , k] <- out[[k]]$fitmat
        out[[k]]$residmat <- NULL
        out[[k]]$fitmat <- NULL
    }
    meanfit <- apply(fitmat, c(1, 2), mean)
    meanfit2 <- apply(meanfit^2, 2, sum)
    for (k in 1:nfolds) {
        residmates[, k] <- apply((fitmat[, , k] - meanfit)^2, 
            2, sum)/meanfit2
    }
    residmates[is.na(residmates)] <- Inf
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/nfolds)
    es <- apply(residmates, 1, mean)
    es.error <- sqrt(apply(residmates, 1, var)/nfolds)
    indcv <- which.min(cv)
    lambda.cv <- lambda[indcv]
    indcv0 <- indcv
    cv1se <- cv
    cv1se[cv <= (cv[indcv] + cv.error[indcv])] <- cv[indcv] + 
        cv.error[indcv]
    indcv1se <- which.min(cv1se)
    lambda.cv1se <- lambda[indcv1se]
    indescv <- which.min(es[1:indcv])
    lambda.escv <- lambda[indescv]
    if (cv.OLS) {
        out <- list()
        if (!parallel) {
            for (k in 1:nfolds) {
                test <- foldid == k
                train <- foldid != k
                obj <- glmnet(x[train, , drop = FALSE], y[train], 
                  lambda = lambda[1:indcv], standardize = standardize, 
                  intercept = intercept, ...)
                fitmat <- predict(obj, newx = x)
                predtest <- predict(obj, newx = x[test, , drop = FALSE])
                selectset0 <- rep(0, p)
                for (i in 1:indcv) {
                  selectset <- abs(obj$beta[, i]) > 0
                  if (sum(selectset) > 0) {
                    if (sum(abs(selectset - selectset0)) > 0) {
                      mls.obj <- mls(x[train, selectset, drop = FALSE], 
                        y[train], tau = tau, standardize = standardize, 
                        intercept = intercept)
                      fitmat[, i] <- mypredict(mls.obj, newx = x[, 
                        selectset, drop = FALSE])
                      predtest[, i] <- mypredict(mls.obj, newx = x[test, 
                        selectset, drop = FALSE])
                    }
                    else {
                      fitmat[, i] <- fitmat[, i - 1]
                      predtest[, i] <- predtest[, i - 1]
                    }
                    selectset0 <- selectset
                  }
                }
                residmat <- apply((y[test] - predtest)^2, 2, 
                  mean)
                out[[k]] <- list(residmat = residmat, fitmat = fitmat)
            }
        }
        else {
            out <- foreach(k = 1:nfolds) %dopar% {
                test <- foldid == k
                train <- foldid != k
                obj <- glmnet(x[train, , drop = FALSE], y[train], 
                  lambda = lambda[1:indcv], standardize = standardize, 
                  intercept = intercept, ...)
                fitmat <- predict(obj, newx = x)
                predtest <- predict(obj, newx = x[test, , drop = FALSE])
                selectset0 <- rep(0, p)
                for (i in 1:indcv) {
                  selectset <- abs(obj$beta[, i]) > 0
                  if (sum(selectset) > 0) {
                    if (sum(abs(selectset - selectset0)) > 0) {
                      mls.obj <- mls(x[train, selectset, drop = FALSE], 
                        y[train], tau = tau, standardize = standardize, 
                        intercept = intercept)
                      fitmat[, i] <- mypredict(mls.obj, newx = x[, 
                        selectset, drop = FALSE])
                      predtest[, i] <- mypredict(mls.obj, newx = x[test, 
                        selectset, drop = FALSE])
                    }
                    else {
                      fitmat[, i] <- fitmat[, i - 1]
                      predtest[, i] <- predtest[, i - 1]
                    }
                    selectset0 <- selectset
                  }
                }
                residmat <- apply((y[test] - predtest)^2, 2, 
                  mean)
                list(residmat = residmat, fitmat = fitmat)
            }
        }
        residmat <- matrix(0, indcv, nfolds)
        residmates <- matrix(0, indcv, nfolds)
        fitmat <- array(0, dim = c(n, indcv, nfolds))
        for (k in 1:nfolds) {
            residmat[, k] <- out[[k]]$residmat
            fitmat[, , k] <- out[[k]]$fitmat
            out[[k]]$residmat <- NULL
            out[[k]]$fitmat <- NULL
        }
        meanfit <- apply(fitmat, c(1, 2), mean)
        meanfit2 <- apply(meanfit^2, 2, sum)
        for (k in 1:nfolds) {
            residmates[, k] <- apply((fitmat[, , k] - meanfit)^2, 
                2, sum)/meanfit2
        }
        residmates[is.na(residmates)] <- Inf
        cv <- apply(residmat, 1, mean)
        cv.error <- sqrt(apply(residmat, 1, var)/nfolds)
        es <- apply(residmates, 1, mean)
        es.error <- sqrt(apply(residmates, 1, var)/nfolds)
        indcv <- which.min(cv)
        lambda.cv <- lambda[indcv]
        cv1se <- cv
        cv1se[cv <= (cv[indcv] + cv.error[indcv])] <- cv[indcv] + 
            cv.error[indcv]
        indcv1se <- which.min(cv1se)
        lambda.cv1se <- lambda[indcv1se]
        indescv <- which.min(es[1:indcv0])
        lambda.escv <- lambda[indescv]
    }
    object <- list(lambda = lambda, glmnet.fit = glmnet.object, 
        cv = cv, cv.error = cv.error, es = es, es.error = es.error, 
        lambda.cv = lambda.cv, lambda.cv1se = lambda.cv1se, lambda.escv = lambda.escv)
    object
}

# ---- mls ----
mls <- function (x, y, tau = 0, standardize = TRUE, intercept = TRUE) 
{
    x <- as.matrix(x)
    n <- nrow(x)
    m <- ncol(x)
    one <- rep(1, n)
    if (intercept) {
        meanx <- drop(one %*% x)/n
        x <- scale(x, meanx, FALSE)
        mu <- mean(y)
        y <- drop(y - mu)
    }
    else {
        meanx <- rep(0, m)
        mu <- 0
        y <- drop(y)
    }
    if (standardize) {
        normx <- sqrt(drop(one %*% (x^2)))
        x <- scale(x, FALSE, normx)
    }
    else {
        normx <- rep(1, m)
    }
    if (abs(tau) <= 1e-08) {
        obj <- lm(y ~ . - 1, data = data.frame(y = y, x = x))
        beta <- coef(obj)
        beta[is.na(beta)] <- 0
    }
    else {
        s <- svd(x)
        Ux <- s$u
        Vx <- s$v
        Dx <- s$d/sqrt(n)
        indi <- Dx > tau
        D <- rep(0, length(Dx))
        D[indi] <- 1/Dx[indi]
        beta <- 1/sqrt(n) * Vx %*% diag(D, length(Dx), length(Dx)) %*% 
            t(Ux) %*% y
    }
    beta <- drop(scale(t(beta), FALSE, normx))
    object <- list()
    object$beta <- beta
    object$beta0 <- mu - drop(meanx %*% beta)
    object$meanx <- meanx
    object$mu <- mu
    object$normx <- normx
    object$tau <- tau
    object
}

# ---- mypredict ----
mypredict <- function (object, newx) 
{
    p <- length(object$meanx)
    drop(scale(newx, object$meanx, FALSE) %*% matrix(object$beta, 
        nrow = p)) + object$mu
}

# ---- PartRidge ----
PartRidge <- function (x, y, lambda2 = 0, varset, standardize = TRUE, intercept = TRUE) 
{
    x <- as.matrix(x)
    n <- nrow(x)
    m <- ncol(x)
    one <- rep(1, n)
    if (intercept) {
        meanx <- drop(one %*% x)/n
        x <- scale(x, meanx, FALSE)
        mu <- mean(y)
        y <- drop(y - mu)
    }
    else {
        meanx <- rep(0, m)
        mu <- 0
        y <- drop(y)
    }
    if (standardize) {
        normx <- sqrt(drop(one %*% (x^2)))
        x <- scale(x, FALSE, normx)
    }
    else {
        normx <- rep(1, m)
    }
    penalty.factor <- rep(1, m)
    if (sum(varset) > 0) {
        penalty.factor[varset] <- 0
    }
    xy <- t(x) %*% y/n
    xx <- t(x) %*% x/n
    diag(xx) <- diag(xx) + lambda2 * penalty.factor
    beta <- solve(xx, xy, tol = 1e-64)
    beta <- drop(beta/normx)
    object <- list()
    object$beta <- beta
    object$beta0 <- mu - drop(meanx %*% beta)
    object$meanx <- meanx
    object$mu <- mu
    object$normx <- normx
    object$lambda2 <- lambda2
    object
}

# ---- Lasso ----
Lasso <- function (x, y, lambda = NULL, fix.lambda = TRUE, cv.method = "cv", 
    nfolds = 10, foldid, cv.OLS = FALSE, tau = 0, parallel = FALSE, 
    standardize = TRUE, intercept = TRUE, ...) 
{
    x <- as.matrix(x)
    n <- dim(x)[1]
    p <- dim(x)[2]
    if (is.null(lambda) & fix.lambda) {
        stop("Should give a value of lambda for fix.lambda=TRUE")
    }
    if (length(lambda) > 1 & fix.lambda) {
        stop("The length of lambda should be 1 if fix.lambda=TRUE")
    }
    if (fix.lambda) {
        globalfit <- glmnet(x, y, standardize = standardize, 
            intercept = intercept, ...)
        fitlasso <- predict(globalfit, type = "coefficients", 
            s = lambda)
        beta0 <- fitlasso[1]
        beta <- fitlasso[-1]
    }
    else {
        globalfit <- glmnet(x, y, lambda = lambda, standardize = standardize, 
            intercept = intercept, ...)
        cvobj <- escv.glmnet(x = x, y = y, lambda = lambda, nfolds = nfolds, 
            foldid = foldid, tau = tau, cv.OLS = cv.OLS, parallel = parallel, 
            standardize = standardize, intercept = intercept, 
            ...)
        if (cv.method == "cv") {
            lambda.opt <- cvobj$lambda.cv
        }
        if (cv.method == "cv1se") {
            lambda.opt <- cvobj$lambda.cv1se
        }
        if (cv.method == "escv") {
            lambda.opt <- cvobj$lambda.escv
        }
        fitlasso <- predict(globalfit, type = "coefficients", 
            s = lambda.opt)
        beta0 <- fitlasso[1]
        beta <- fitlasso[-1]
        lambda <- lambda.opt
    }
    if (intercept) {
        meanx <- apply(x, 2, mean)
        mu <- mean(y)
    }
    else {
        meanx <- rep(0, p)
        mu <- 0
    }
    object <- list()
    object$beta0 <- beta0
    object$beta <- beta
    object$lambda <- lambda
    object$meanx <- meanx
    object$mu <- mu
    object
}

# ---- LassoOLS ----
LassoOLS <- function (x, y, OLS = TRUE, lambda = NULL, fix.lambda = TRUE, 
    cv.method = "cv", nfolds = 10, foldid, cv.OLS = TRUE, tau = 0, 
    parallel = FALSE, standardize = TRUE, intercept = TRUE, ...) 
{
    x <- as.matrix(x)
    n <- dim(x)[1]
    p <- dim(x)[2]
    if (is.null(lambda) & fix.lambda) {
        stop("Should given a value of lambda for fix.lambda=TRUE")
    }
    if (length(lambda) > 1 & fix.lambda) {
        stop("The length of lambda should be 1 if fix.lambda=TRUE")
    }
    if (fix.lambda) {
        globalfit <- glmnet(x, y, standardize = standardize, 
            intercept = intercept, ...)
        fitlasso <- predict(globalfit, type = "coefficients", 
            s = lambda)
        betalasso <- fitlasso[-1]
        selectvar <- betalasso != 0
        beta0 <- fitlasso[1]
        beta <- fitlasso[-1]
        if (OLS & sum(selectvar) > 0) {
            ls.obj <- mls(x[, selectvar, drop = FALSE], y, tau, 
                standardize, intercept)
            beta0 <- ls.obj$beta0
            beta[selectvar] <- ls.obj$beta
        }
    }
    else {
        globalfit <- glmnet(x, y, lambda = lambda, standardize = standardize, 
            intercept = intercept, ...)
        cvobj <- escv.glmnet(x = x, y = y, lambda = lambda, nfolds = nfolds, 
            foldid = foldid, tau = tau, cv.OLS = cv.OLS, parallel = parallel, 
            standardize = standardize, intercept = intercept, 
            ...)
        if (cv.method == "cv") {
            lambda.opt <- cvobj$lambda.cv
        }
        if (cv.method == "cv1se") {
            lambda.opt <- cvobj$lambda.cv1se
        }
        if (cv.method == "escv") {
            lambda.opt <- cvobj$lambda.escv
        }
        fitlasso <- predict(globalfit, type = "coefficients", 
            s = lambda.opt)
        betalasso <- fitlasso[-1]
        selectvar <- betalasso != 0
        beta0 <- fitlasso[1]
        beta <- fitlasso[-1]
        if (OLS & sum(selectvar) > 0) {
            ls.obj <- mls(x[, selectvar, drop = FALSE], y, tau, 
                standardize, intercept)
            beta0 <- ls.obj$beta0
            beta[selectvar] <- ls.obj$beta
        }
    }
    if (intercept) {
        meanx <- apply(x, 2, mean)
        mu <- mean(y)
    }
    else {
        meanx <- rep(0, p)
        mu <- 0
    }
    object <- list()
    object$beta0 <- beta0
    object$beta <- beta
    object$lambda <- lambda
    object$meanx <- meanx
    object$mu <- mu
    object$tau <- tau
    object
}

# ---- ci ----
ci <- function (Beta, Beta_bootstrap, alpha = 0.05, type = c("basic", 
    "quantile", "bca", "basic2"), a, Beta2) 
{
    p <- dim(Beta_bootstrap)[2]
    B <- dim(Beta_bootstrap)[1]
    interval <- matrix(0, 2, p)
    if (type == "basic") {
        bound.percentile <- apply(Beta_bootstrap, 2, function(u) {
            quantile(u, prob = c(1 - alpha/2, alpha/2))
        })
        interval[1, ] <- 2 * Beta - bound.percentile[1, ]
        interval[2, ] <- 2 * Beta - bound.percentile[2, ]
    }
    if (type == "basic2") {
        bound.percentile <- apply(Beta_bootstrap, 2, function(u) {
            quantile(u, prob = c(1 - alpha/2, alpha/2))
        })
        interval[1, ] <- Beta + Beta2 - bound.percentile[1, ]
        interval[2, ] <- Beta + Beta2 - bound.percentile[2, ]
    }
    if (type == "quantile") {
        bound.percentile <- apply(Beta_bootstrap, 2, function(u) {
            quantile(u, prob = c(alpha/2, 1 - alpha/2))
        })
        interval[1, ] <- bound.percentile[1, ]
        interval[2, ] <- bound.percentile[2, ]
    }
    if (type == "bca") {
        if (missing(a)) {
            a <- 0
        }
        fstar <- Beta_bootstrap < (matrix(rep(1, B), B, 1) %*% 
            matrix(Beta, 1, B))
        pstar <- apply(fstar, 2, mean)
        xi <- qnorm(pstar)
        for (j in 1:p) {
            if (xi[j] == -Inf || xi[j] == Inf) {
                alpha1 <- alpha/2
                alpha2 <- 1 - alpha/2
            }
            else {
                alpha1 <- pnorm(xi[j] + (xi[j] + qnorm(alpha/2))/(1 - 
                  a * (xi[j] + qnorm(alpha/2))))
                alpha2 <- pnorm(xi[j] + (xi[j] + qnorm(1 - alpha/2))/(1 - 
                  a * (xi[j] + qnorm(1 - alpha/2))))
            }
            interval[, j] <- quantile(Beta_bootstrap[, j, drop = FALSE], 
                prob = c(alpha1, alpha2))
        }
    }
    return(interval)
}

