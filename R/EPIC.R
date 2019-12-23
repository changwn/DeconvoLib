# code from EPIC r package. 
#signature_epic <- EPIC::TRef$refProfiles
#gene_epic <- EPIC::TRef$sigGenes

EPIC <- function (bulk, reference = NULL, mRNA_cell = NULL, mRNA_cell_sub = NULL, 
          sigGenes = NULL, scaleExprs = TRUE, withOtherCells = TRUE, 
          constrainedSum = TRUE, rangeBasedOptim = FALSE) 
{
  if (!is.matrix(bulk) && !is.data.frame(bulk)) 
    stop("'bulk' needs to be given as a matrix or data.frame")
  with_w <- TRUE
  if (is.null(reference)) {
    reference <- EPIC::TRef
  }
  else if (is.character(reference)) {
    if (reference %in% prebuiltRefNames) {
      reference <- get(reference, pos = "package:EPIC")
    }
    else stop("The reference, '", reference, "' is not part of the allowed ", 
              "references:", paste(prebuiltRefNames, collapse = ", "))
  }
  else if (is.list(reference)) {
    refListNames <- names(reference)
    if ((!all(c("refProfiles", "sigGenes") %in% 
              refListNames)) || (("refProfiles" %in% refListNames) && 
                                 !is.null(sigGenes))) 
      stop("Reference, when given as a list needs to contain at least the ", 
           "fields 'refProfiles' and 'sigGenes' (sigGenes could also be ", 
           "given as input to EPIC instead)")
    if (!is.matrix(reference$refProfiles) && !is.data.frame(reference$refProfiles)) 
      stop("'reference$refProfiles' needs to be given as a matrix or data.frame")
    if (!("refProfiles.var" %in% refListNames)) {
      warning("'refProfiles.var' not defined; using identical weights ", 
              "for all genes")
      with_w <- FALSE
    }
    else if (!is.matrix(reference$refProfiles.var) && !is.data.frame(reference$refProfiles.var)) {
      stop("'reference$refProfiles.var' needs to be given as a matrix or ", 
           "data.frame when present.")
    }
    else if (!identical(dim(reference$refProfiles.var), dim(reference$refProfiles)) || 
             !identical(dimnames(reference$refProfiles.var), dimnames(reference$refProfiles))) 
      stop("The dimensions and dimnames of 'reference$refProfiles' and ", 
           "'reference$refProfiles.var' need to be the same")
  }
  else {
    stop("Unknown format for 'reference'")
  }
  bulk <- merge_duplicates(bulk, in_type = "bulk samples")
  refProfiles <- merge_duplicates(reference$refProfiles, in_type = "reference profiles")
  if (with_w) {
    refProfiles.var <- merge_duplicates(reference$refProfiles.var, 
                                        warn = F)
  }
  else {
    refProfiles.var <- 0
  }
  nSamples <- NCOL(bulk)
  samplesNames <- colnames(bulk)
  if (is.null(samplesNames)) {
    samplesNames <- 1:nSamples
    colnames(bulk) <- samplesNames
  }
  nRefCells <- NCOL(refProfiles)
  refCellsNames <- colnames(refProfiles)
  bulk_NA <- apply(is.na(bulk), MARGIN = 1, FUN = all)
  if (any(bulk_NA)) {
    warning(sum(bulk_NA), " genes are NA in all bulk samples, removing these.")
    bulk <- bulk[!bulk_NA, ]
  }
  bulkGenes <- rownames(bulk)
  refGenes <- rownames(refProfiles)
  commonGenes <- intersect(bulkGenes, refGenes)
  if (is.null(sigGenes)) 
    sigGenes <- unique(reference$sigGenes)
  sigGenes <- sigGenes[sigGenes %in% commonGenes]
  nSigGenes <- length(sigGenes)
  if (nSigGenes < nRefCells) 
    stop("There are only ", nSigGenes, " signature genes", 
         " matching common genes between bulk and reference profiles,", 
         " but there should be more signature genes than reference cells")
  if (scaleExprs) {
    if (length(commonGenes) < 2000) 
      warning("there are few genes in common between the bulk samples and ", 
              "reference cells:", length(commonGenes), 
              ", so the data scaling ", "might be an issue")
    bulk <- scaleCounts(bulk, sigGenes, commonGenes)$counts
    temp <- scaleCounts(refProfiles, sigGenes, commonGenes)
    refProfiles <- temp$counts
    if (with_w) 
      refProfiles.var <- scaleCounts(refProfiles.var, sigGenes, 
                                     normFact = temp$normFact)$counts
  }
  else {
    bulk <- bulk[sigGenes, ]
    refProfiles <- refProfiles[sigGenes, ]
    if (with_w) 
      refProfiles.var <- refProfiles.var[sigGenes, ]
  }
  if (is.null(mRNA_cell)) 
    mRNA_cell <- EPIC::mRNA_cell_default
  if (!is.null(mRNA_cell_sub)) {
    if (is.null(names(mRNA_cell_sub)) || !is.numeric(mRNA_cell_sub)) 
      stop("When mRNA_cell_sub is given, it needs to be a named numeric vector")
    mRNA_cell[names(mRNA_cell_sub)] <- mRNA_cell_sub
  }
  minFun <- function(x, A, b, w) {
    return(sum((w * (A %*% x - b)^2), na.rm = TRUE))
  }
  minFun.range <- function(x, A, b, A.var) {
    val.max <- (A + A.var) %*% x - b
    val.min <- (A - A.var) %*% x - b
    cErr <- rep(0, length(b))
    outOfRange <- (sign(val.max) * sign(val.min) == 1)
    cErr[outOfRange] <- pmin(abs(val.max[outOfRange]), abs(val.min[outOfRange]))
    return(sum(cErr, na.rm = TRUE))
  }
  if (with_w && !rangeBasedOptim) {
    w <- rowSums(refProfiles/(refProfiles.var + 1e-12), na.rm = TRUE)
    med_w <- stats::median(w[w > 0], na.rm = TRUE)
    w[w > 100 * med_w] <- 100 * med_w
  }
  else w <- 1
  if (withOtherCells) {
    cMin <- 0
  }
  else {
    cMin <- 0.99
  }
  cMax <- 1
  ui <- diag(nRefCells)
  ci <- rep(0, nRefCells)
  if (constrainedSum) {
    ui <- rbind(ui, rep(1, nRefCells), rep(-1, nRefCells))
    ci <- c(ci, cMin, -cMax)
  }
  cInitProp <- (min(1, cMax) - 1e-05)/nRefCells
  tempPropPred <- lapply(1:nSamples, FUN = function(cSample) {
    b <- bulk[, cSample]
    if (!rangeBasedOptim) {
      fit <- stats::constrOptim(theta = rep(cInitProp, 
                                            nRefCells), f = minFun, grad = NULL, ui = ui, 
                                ci = ci, A = refProfiles, b = b, w = w)
    }
    else {
      fit <- stats::constrOptim(theta = rep(cInitProp, 
                                            nRefCells), f = minFun.range, grad = NULL, ui = ui, 
                                ci = ci, A = refProfiles, b = b, A.var = refProfiles.var)
    }
    fit$x <- fit$par
    if (!withOtherCells) 
      fit$x <- fit$x/sum(fit$x, na.rm = TRUE)
    b_estimated <- refProfiles %*% fit$x
    if (nSigGenes > 2) {
      suppressWarnings(corSp.test <- stats::cor.test(b, 
                                                     b_estimated, method = "spearman"))
      corPear.test <- stats::cor.test(b, b_estimated, method = "pearson")
    }
    else {
      corSp.test <- corPear.test <- list()
      corSp.test$estimate <- corSp.test$p.value <- corPear.test$estimate <- corPear.test$p.value <- NA
    }
    regLine <- stats::lm(b_estimated ~ b)
    regLine_through0 <- stats::lm(b_estimated ~ b + 0)
    if (!rangeBasedOptim) {
      rmse_pred <- sqrt(minFun(x = fit$x, A = refProfiles, 
                               b = b, w = w)/nSigGenes)
      rmse_0 <- sqrt(minFun(x = rep(0, nRefCells), A = refProfiles, 
                            b = b, w = w)/nSigGenes)
    }
    else {
      rmse_pred <- sqrt(minFun.range(x = fit$x, A = refProfiles, 
                                     b = b, A.var = refProfiles.var)/nSigGenes)
      rmse_0 <- sqrt(minFun.range(x = rep(0, nRefCells), 
                                  A = refProfiles, b = b, A.var = refProfiles.var)/nSigGenes)
    }
    gof <- data.frame(fit$convergence, ifelse(is.null(fit$message), 
                                              "", fit$message), rmse_pred, rmse_0, corSp.test$estimate, 
                      corSp.test$p.value, corPear.test$estimate, corPear.test$p.value, 
                      regLine$coefficients[2], regLine$coefficients[1], 
                      regLine_through0$coefficients[1], sum(fit$x), stringsAsFactors = FALSE)
    return(list(mRNAProportions = fit$x, fit.gof = gof))
  })
  mRNAProportions <- do.call(rbind, lapply(tempPropPred, function(x) x$mRNAProportions))
  dimnames(mRNAProportions) <- list(samplesNames, refCellsNames)
  fit.gof <- do.call(rbind, lapply(tempPropPred, function(x) x$fit.gof))
  dimnames(fit.gof) <- list(samplesNames, c("convergeCode", 
                                            "convergeMessage", "RMSE_weighted", "Root_mean_squared_geneExpr_weighted", 
                                            "spearmanR", "spearmanP", "pearsonR", 
                                            "pearsonP", "regline_a_x", "regline_b", 
                                            "regline_a_x_through0", "sum_mRNAProportions"))
  if (any(fit.gof$convergeCode != 0)) 
    warning("The optimization didn't fully converge for some samples:\n", 
            paste(samplesNames[fit.gof$convergeCode != 0], collapse = "; "), 
            "\n - check fit.gof for the convergeCode and convergeMessage")
  if (withOtherCells) 
    mRNAProportions <- cbind(mRNAProportions, otherCells = 1 - 
                               rowSums(mRNAProportions))
  tInds <- match(colnames(mRNAProportions), names(mRNA_cell))
  if (anyNA(tInds)) {
    defaultInd <- match("default", names(mRNA_cell))
    if (is.na(defaultInd)) {
      tStr <- paste(" and no default value is given for this mRNA per cell,", 
                    "so we cannot estimate the cellFractions, only", 
                    "the mRNA proportions")
    }
    else {
      tStr <- paste(" - using the default value of", 
                    mRNA_cell[defaultInd], "for these but this might bias the true cell proportions from", 
                    "all cell types.")
    }
    warning("mRNA_cell value unknown for some cell types: ", 
            paste(colnames(mRNAProportions)[is.na(tInds)], collapse = ", "), 
            tStr)
    tInds[is.na(tInds)] <- defaultInd
  }
  cellFractions <- t(t(mRNAProportions)/mRNA_cell[tInds])
  cellFractions <- cellFractions/rowSums(cellFractions, na.rm = FALSE)
  return(list(mRNAProportions = mRNAProportions, cellFractions = cellFractions, 
              fit.gof = fit.gof))
}

merge_duplicates <- function (mat, warn = TRUE, in_type = NULL) 
{
    dupl <- duplicated(rownames(mat))
    if (sum(dupl) > 0) {
        dupl_genes <- unique(rownames(mat)[dupl])
        if (warn) {
            warning("There are ", length(dupl_genes), " duplicated gene names", 
                ifelse(!is.null(in_type), paste(" in the", in_type), 
                  ""), ". We'll use the median value for each of these cases.")
        }
        mat_dupl <- mat[rownames(mat) %in% dupl_genes, , drop = F]
        mat_dupl_names <- rownames(mat_dupl)
        mat <- mat[!dupl, , drop = F]
        mat[dupl_genes, ] <- t(sapply(dupl_genes, FUN = function(cgene) apply(mat_dupl[mat_dupl_names == 
            cgene, , drop = F], MARGIN = 2, FUN = median)))
    }
    return(mat)
}


scaleCounts <- function (counts, sigGenes = NULL, renormGenes = NULL, normFact = NULL) 
{
    if (is.null(sigGenes)) 
        sigGenes <- 1:nrow(counts)
    if (is.null(normFact)) {
        if (is.null(renormGenes)) 
            renormGenes <- 1:nrow(counts)
        normFact <- colSums(counts[renormGenes, , drop = FALSE], 
            na.rm = TRUE)
    }
    counts <- t(t(counts[sigGenes, , drop = FALSE])/normFact) * 
        1e+06
    return(list(counts = counts, normFact = normFact))
}
