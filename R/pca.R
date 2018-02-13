# 
# # MJM NOTE: likely not necessary, just figure out how to properly import the
# # function
# 
# 
# # this is just a copied and pasted version of the "pca" function from the pcaMethods package.
# # I have no idea at all why it works when I do this, but it does, so, I guess I will use it.
# pca_sj <- function (object, method, nPcs = 2, scale = c("none", "pareto", 
#                                                         "vector", "uv"), center = TRUE, completeObs = TRUE, subset = NULL, 
#                     cv = c("none", "q2"), ...) 
# {
#   if (inherits(object, "data.frame")) {
#     num <- vapply(object, is.numeric, logical(1))
#     if (sum(num) < 2) 
#       stop("no numeric data in supplied data.frame")
#     Matrix <- as.matrix(object[, num])
#   }
#   else if (inherits(object, "ExpressionSet")) {
#     Matrix <- t(exprs(object))
#   }
#   else Matrix <- as.matrix(object, rownames.force = TRUE)
#   if (!is.null(subset)) 
#     Matrix <- Matrix[, subset]
#   cv <- match.arg(cv)
#   scale <- match.arg(scale)
#   if (nPcs > ncol(Matrix)) {
#     warning("more components than matrix columns requested")
#     nPcs <- min(dim(Matrix))
#   }
#   if (nPcs > nrow(Matrix)) {
#     warning("more components than matrix rows requested")
#     nPcs <- min(dim(Matrix))
#   }
#   if (!checkData(Matrix, verbose = interactive())) 
#     stop("Invalid data format.", "Run checkData(data, verbose=TRUE) for details")
#   missing <- is.na(Matrix)
#   if (missing(method)) {
#     if (any(missing)) 
#       method <- "nipals"
#     else method <- "svd"
#   }
#   if (any(missing) & method == "svd") {
#     warning("data has missing values using nipals instead of user requested svd")
#     method <- "nipals"
#   }
#   method <- match.arg(method, choices = listPcaMethods())
#   prepres <- prep(Matrix, scale = scale, center = center, simple = FALSE, 
#                   ...)
#   switch(method, svd = {
#            res <- svdPca(prepres$data, nPcs = nPcs)
#                   }, nipals = {
#                     res <- nipalsPca(prepres$data, nPcs = nPcs)
#                   }, rnipals = {
#                     res <- RnipalsPca(prepres$data, nPcs = nPcs)
#                   }, bpca = {
#                     res <- bpca(prepres$data, nPcs = nPcs)
#                   }, ppca = {
#                     res <- ppca(prepres$data, nPcs = nPcs)
#                   }, svdImpute = {
#                     res <- svdImpute(prepres$data, nPcs = nPcs)
#                   }, robustPca = {
#                     res <- robustPca(prepres$data, nPcs = nPcs)
#                   }, nlpca = {
#                     res <- nlpca(prepres$data, nPcs = nPcs)
#                   })
#   nPcs <- ncol(res@scores)
#   if (is.null(scores(res)) | is.null(loadings(res)) | is.null(R2cum(res)) | 
#       is.null(method(res))) 
#     stop(paste("bad result from pca method", method))
#   colnames(res@scores) <- paste("PC", 1:nPcs, sep = "")
#   rownames(res@scores) <- rownames(Matrix)
#   if (all(dim(loadings(res)) == c(ncol(Matrix), nPcs))) {
#     colnames(res@loadings) <- paste("PC", 1:nPcs, sep = "")
#     rownames(res@loadings) <- colnames(Matrix)
#   }
#   if (!is.null(subset)) 
#     res@subset <- subset
#   res@missing <- missing
#   res@nPcs <- nPcs
#   res@nObs <- nrow(Matrix)
#   res@nVar <- ncol(Matrix)
#   res@sDev <- apply(scores(res), 2, sd)
#   res@center <- prepres$center
#   res@centered <- center
#   res@scale <- prepres$scale
#   res@scaled <- scale
#   res@R2 <- res@R2cum[1]
#   if (length(res@R2cum) > 1) 
#     res@R2 <- c(res@R2, diff(res@R2cum))
#   if (completeObs) {
#     cObs <- Matrix
#     if (method %in% listPcaMethods("nonlinear")) 
#       cObs[missing] <- fitted(res, Matrix, pre = TRUE, 
#                               post = TRUE)[missing]
#     else cObs[missing] <- fitted(res, post = TRUE)[missing]
#     res@completeObs <- cObs
#   }
#   if (cv == "q2") 
#     res@cvstat <- Q2(res, Matrix, nruncv = 1, ...)
#   return(res)
# }
# 
