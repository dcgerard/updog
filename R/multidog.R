################
## Multiple SNP genotyping
################




#' Fit \code{\link{flexdog}} to multiple SNPs.
#'
#' This is a convenience function that will run \code{\link{flexdog}} over many SNPs.
#' Support is provided for parallel computing through the doParallel package.
#' This function has not been extensively tested. Please report any bugs to
#' \url{http://github.com/dcgerard/updog/issues}.
#'
#' You should format your reference counts and total read counts in two
#' separate matrices. The rows should index the markers (SNPs) and the
#' columns should index the individuals. Row names are how we ID the SNPs
#' and column names are how we ID the individuals, and so they are required
#' attributes.
#'
#' If your data are in VCF files, I would recommend importing them using the
#' VariantAnnotation package from Bioconductor
#' \url{https://bioconductor.org/packages/VariantAnnotation/}. It's a great
#' VCF parser.
#'
#' See the details of \code{\link{flexdog}} for the possible values of
#' \code{model}.
#'
#' If \code{model = "f1"}, \code{model = "s1"}, \code{model = "f1pp"}
#' or \code{model = "s1pp"} then the user may
#' provide the individual ID for parent(s) via the \code{p1_id}
#' and \code{p2_id} arguments.
#'
#' The output of the function is written in multiple tables parallely. This
#' can cause corruption of some lines in the writing process. These lines
#' will be eliminated, and thus also the markers that they contain. However,
#' the proportion of corrupted lines is very small (our tests indicate 0.17%)
#' and it occurs rarely in smaller datasets.
#'
#' SNPs that contain 0 reads (or all missing data) are entirely removed.
#'
#' @section Parallel Computation:
#'
#' The \code{multidog()} function supports parallel computing. It does
#' so through the \href{https://cran.r-project.org/package=future}{future}
#' package.
#'
#' If you are just running \code{multidog()} on a local machine, then you
#' can use the \code{nc} argument to specify the parallelization. Any value
#' of \code{nc} greater than 1 will result in multiple background R sessions to
#' genotype all of the SNPs. The maximum value of \code{nc} you should
#' try can be found by running \code{future::availableCores()}. Running
#' \code{multidog()} using \code{nc} is equivalent to setting the future
#' plan with \code{future::plan(future::multisession, workers = nc)}.
#'
#' Using the future package means that different evaluation strategies
#' are possible. In particular, if you are using a high performance machine,
#' you can explore using the
#' \href{https://cran.r-project.org/package=future.batchtools}{future.batchtools}
#' package to evaluate \code{multidog()} using schedulers like Slurm
#' or TORQUE/PBS.
#'
#' To use a different strategy, set \code{nc = NA} and then
#'  run \code{future::plan()} prior to
#' running \code{multidog()}. For example, to set up forked R processes
#' on your current machine (instead of using background R sessions), you would
#' run (will not work on Windows):
#' \code{future::plan(future::multicore)}, followed by
#' running \code{multidog()} with \code{nc = NA}. See the examples below.
#'
#' @inheritParams flexdog
#' @param refmat A matrix of reference read counts. The columns index
#'     the individuals and the rows index the markers (SNPs). This matrix must have
#'     rownames (for the names of the markers) and column names (for the names
#'     of the individuals). These names must match the names in \code{sizemat}.
#' @param sizemat A matrix of total read counts. The columns index
#'     the individuals and the rows index the markers (SNPs). This matrix must have
#'     rownames (for the names of the markers) and column names (for the names
#'     of the individuals). These names must match the names in \code{refmat}.
#' @param nc The number of computing cores to use when doing parallelization
#'     on your local machine. See the section "Parallel Computation" for how
#'     to implement more complicated evaluation strategies using the
#'     \code{future} package.
#'
#'     When you are specifying other evaluation strategies using the
#'     \code{future} package, you should also set \code{nc = NA}.
#'
#'     The value of \code{nc} should never be
#'     more than the number of cores available in your computing environment.
#'     You can determine the maximum number of available cores by running
#'     \code{future::availableCores()} in R.
#' @param p1_id The ID of the first parent. This should be a character of
#'     length 1. This should correspond to a single column name in \code{refmat}
#'     and \code{sizemat}.
#' @param p2_id The ID of the second parent. This should be a character of
#'     length 1. This should correspond to a single column name in \code{refmat}
#'     and \code{sizemat}.
#'
#'#' @param out Location where the output files will be written. A prefix for the file
#'     might also be added. Defaults to "updog".
#' @param outpars Character vector indicating which outputs from \code{flexdog}
#'     should be written out.
#'     Can take either of the following values:
#'\describe{
#'    \item{\code{"model"}:}{All snp genotyping parameters and model parameters, such as
#'    bias, overdispersion, sequencing error, log-likelihood, etc.}
#'    \item{\code{"all"}:}{All possible outputs}
#'    \item{Any of the outputs in \code{flexdog()}:}{bias, seq, od, num_iter, llike, prop_mis, gene_dist,
#'    par, geno, maxpostprob, postmean, postmat, genologlike}
#'    }. For a full description of what each parameter means see \code{\link{flexdog}()}.
#'    
#' @return A list-like object of two data frames.
#' \describe{
#' \item{\code{snpdf}}{A data frame containing properties of the SNPs (markers).
#'     The rows index the SNPs. The variables include:
#'     \describe{
#'     \item{\code{snp}}{The name of the SNP (marker).}
#'     \item{\code{bias}}{The estimated allele bias of the SNP.}
#'     \item{\code{seq}}{The estimated sequencing error rate of the SNP.}
#'     \item{\code{od}}{The estimated overdispersion parameter of the SNP.}
#'     \item{\code{prop_mis}}{The estimated proportion of individuals
#'         misclassified in the SNP.}
#'     \item{\code{num_iter}}{The number of iterations performed during
#'         the EM algorithm for that SNP.}
#'     \item{\code{llike}}{The maximum marginal likelihood of the SNP.}
#'     \item{\code{ploidy}}{The provided ploidy of the species.}
#'     \item{\code{model}}{The provided model for the prior genotype
#'         distribution.}
#'     \item{\code{p1ref}}{The user-provided reference read counts of parent 1.}
#'     \item{\code{p1size}}{The user-provided total read counts of parent 1.}
#'     \item{\code{p2ref}}{The user-provided reference read counts of parent 2.}
#'     \item{\code{p2size}}{The user-provided total read counts of parent 2.}
#'     \item{\code{Pr_k}}{The estimated frequency of individuals with genotype
#'         k, where k can be any integer between 0 and the ploidy level.}
#'     \item{Model specific parameter estimates}{See the return value of
#'         \code{par} in the help page of \code{\link{flexdog}}.}
#'     }}
#' \item{\code{inddf}}{A data frame containing the properties of the
#'     individuals at each SNP. The variables include:
#'     \describe{
#'     \item{\code{snp}}{The name of the SNP (marker).}
#'     \item{\code{ind}}{The name of the individual.}
#'     \item{\code{ref}}{The provided reference counts for that individual at
#'          that SNP.}
#'     \item{\code{size}}{The provided total counts for that individual at
#'          that SNP.}
#'     \item{\code{geno}}{The posterior mode genotype for that individual at
#'          that SNP. This is the estimated reference allele dosage for a
#'          given individual at a given SNP.}
#'     \item{\code{postmean}}{The posterior mean genotype for that individual
#'          at that SNP. This is a continuous genotype estimate of the
#'          reference allele dosage for a given individual at a given SNP.}
#'     \item{\code{maxpostprob}}{The maximum posterior probability. This
#'          is the posterior probability that the individual was genotyped
#'          correctly.}
#'     \item{\code{Pr_k}}{The posterior probability that a given individual
#'          at a given SNP has genotype k, where k can vary from 0 to the
#'          ploidy level of the species.}
#'     \item{\code{logL_k}}{The genotype \emph{log}-likelihoods for dosage
#'          k for a given individual at a given SNP, where k can vary f
#'          rom 0 to the ploidy level of the species.}
#'     }}
#' }
#'
#' @author David Gerard
#'
#' @seealso
#' \itemize{
#'   \item{\code{\link{flexdog}()}:}{For the underlying genotyping function.}
#'   \item{\code{\link{format_multidog}()}:}{For converting the output
#'       of \code{multidog()} to a matrix.}
#'   \item{\code{\link{filter_snp}()}:}{For filtering SNPs using the
#'       output of \code{multidog()}.}
#' }
#'
#' @examples
#' \dontrun{
#' data("uitdewilligen")
#'
#' ## Run multiple R sessions using the `nc` variable.
#' mout <- multidog(refmat = t(uitdewilligen$refmat),
#'                  sizemat = t(uitdewilligen$sizemat),
#'                  ploidy = uitdewilligen$ploidy,
#'                  nc = 2)
#' mout$inddf
#' mout$snpdf
#'
#' ## Run multiple external R sessions on the local machine.
#' ## Note that we set `nc = NA`.
#' cl <- parallel::makeCluster(2, timeout = 60)
#' future::plan(future::cluster, workers = cl)
#' mout <- multidog(refmat = t(uitdewilligen$refmat),
#'                  sizemat = t(uitdewilligen$sizemat),
#'                  ploidy = uitdewilligen$ploidy,
#'                  nc = NA)
#' mout$inddf
#' mout$snpdf
#'
#' ## Close cluster and reset future to current R process
#' parallel::stopCluster(cl)
#' future::plan(future::sequential)
#' }
#'
#' @export
multidog <- function(refmat,
                     sizemat,
                     ploidy,
                     model = c("norm",
                               "hw",
                               "bb",
                               "s1",
                               "s1pp",
                               "f1",
                               "f1pp",
                               "flex",
                               "uniform",
                               "custom"),
                     nc = 1,
                     p1_id = NULL,
                     p2_id = NULL,
                     bias_init = exp(c(-1, -0.5, 0, 0.5, 1)),
                     prior_vec = NULL,
                     ...) {

  cat(paste0(  "    |                                   *.#,%    ",
             "\n   |||                                 *******/  ",
             "\n |||||||    (**..#**.                  */   **/  ",
             "\n|||||||||    */****************************/*%   ",
             "\n   |||    &****..,*.************************/    ",
             "\n   |||     (....,,,*,...****%********/(******    ",
             "\n   |||                ,,****%////,,,,./.****/    ",
             "\n   |||                  /**//         .*///....  ",
             "\n   |||                  .*/*/%#         .,/   ., ",
             "\n   |||               , **/   #%         .*    .. ",
             "\n   |||                               ,,,*        ",
             "\n\nWorking on it..."))

  ## Check input --------------------------------------------------------------
  assertthat::assert_that(is.matrix(refmat))
  assertthat::assert_that(is.matrix(sizemat))
  assertthat::assert_that(is.numeric(refmat))
  assertthat::assert_that(is.numeric(sizemat))
  assertthat::assert_that(!is.null(colnames(refmat)))
  assertthat::assert_that(!is.null(colnames(sizemat)))
  assertthat::assert_that(!is.null(rownames(refmat)))
  assertthat::assert_that(!is.null(rownames(sizemat)))
  if (length(setdiff(colnames(refmat), colnames(sizemat))) != 0) {
    stop(paste0("multidog: refmat and sizemat must have the same column names\n",
                "The following column names are in the set difference:\n",
                setdiff(colnames(refmat), colnames(sizemat))))
  }
  if (length(setdiff(rownames(refmat), rownames(sizemat))) != 0) {
    stop(paste0("multidog: refmat and sizemat must have the same row names\n",
                "The following row names are in the set difference:\n",
                setdiff(rownames(refmat), rownames(sizemat))))
  }
  model <- match.arg(model)
  assertthat::assert_that(length(nc) == 1)
  if (!is.na(nc)) {
    assertthat::assert_that(is.numeric(nc))
    nc <- round(nc)
    assertthat::assert_that(nc >= 1)
  }
  if (!is.null(p1_id)) {
    stopifnot(is.character(p1_id))
    stopifnot(length(p1_id) == 1)
    stopifnot(is.element(el = p1_id, set = colnames(refmat)))
    stopifnot(model == "f1" | model == "s1" | model == "f1pp" | model == "s1pp")
  }
  if (!is.null(p2_id)) {
    stopifnot(is.character(p2_id))
    stopifnot(length(p2_id) == 1)
    stopifnot(is.element(el = p2_id, set = colnames(refmat)))
    stopifnot(model == "f1" | model == "f1pp")
  }
  if (!is.null(p2_id) & is.null(p1_id)) {
    warning("setting p1_id to be p2_id and setting p2_id to be NULL.")
    p1_id <- p2_id
    p2_id <- NULL
  }

  ## Get list of individuals ---------------------------------------------------
  indlist <- colnames(refmat)

  if (!is.null(p1_id)) {
    indlist <- indlist[indlist != p1_id]
  }
  if (!is.null(p2_id)) {
    indlist <- indlist[indlist != p2_id]
  }

  ## Remove NA SNPs ------------------------------------------------------------
  which_bad_size <- apply(X = (sizemat[, indlist, drop = FALSE] == 0) | is.na(sizemat[, indlist, drop = FALSE]),
                          MARGIN = 1,
                          FUN = all)
  which_bad_ref <- apply(X = is.na(refmat[, indlist, drop = FALSE]),
                         MARGIN = 1,
                         FUN = all)

  bad_snps <- unique(c(rownames(sizemat)[which_bad_size], rownames(refmat)[which_bad_ref]))

  if (length(bad_snps) > 0) {
    if (length(bad_snps) == nrow(sizemat)) {
      stop("multidog: All SNPs are missing.")
    }
    sizemat <- sizemat[!(rownames(sizemat) %in% bad_snps), , drop = FALSE]
    refmat  <- refmat[!(rownames(refmat) %in% bad_snps), , drop = FALSE]
  }

  ## Get list of SNPs ---------------------------------------------------------
  snplist <- rownames(refmat)

  ## Initialize write files -----------------------------------------------
  #I make write groups of 100 markers, but other numbers are possible
  #Larger numbers should reduce the chance of corruption
  #But the number should always be larger than snps/ncores
  write_groups <- split(snplist,ceiling(seq_along(1:length(snplist))/100))
  
  #We need to initialize the outfiles by writing the column names
  snp_cols <- c(snp_pars[!snp_pars %in% c("par","gene_dist")])
  if("gene_dist" %in% snp_pars){
    snp_pars <- c(snp_pars[!snp_pars %in% "gene_dist"],"gene_dist")
    snp_cols <- c(snp_cols,paste("gene_dist",0:ploidy,sep = "_"))
  }
  if("par" %in% snp_pars){
    snp_pars <- c(snp_pars[!snp_pars %in% "par"],"par")
    snp_cols <- c(snp_cols,par_in_model(model))
  }
  names_snp_ind_geno <- as.vector(sapply(snp_ind_geno_pars,paste,0:ploidy,sep = "_")  )
  outfiles <- paste0(out,c(snp_ind_pars,names_snp_ind_geno))
  names(outfiles) <- c(snp_ind_pars,names_snp_ind_geno)
  if(length(snp_pars) > 1){
    outfiles <- c(paste0(out,"model"),outfiles)
    names(outfiles)[1] <- "model"
  }
  
  for(n in names(outfiles)){
    if(n == "model"){
      cat("marker",snp_cols,"\n",file = outfiles[n])
    }else{
      cat("marker",indlist,"\n",file = outfiles[n])
    }
  }
  
  ## Register workers ----------------------------------------------------------
  if (nc == 1) {
    foreach::registerDoSEQ()
  } else {
    cl = parallel::makeCluster(nc)
    doParallel::registerDoParallel(cl = cl)
    if (foreach::getDoParWorkers() == 1) {
      stop("multidog: nc > 1 but only one core registered from foreach::getDoParWorkers().")
    }
  }
  
  ## Fit flexdog on all SNPs --------------------------------------------------
  #The result object N is an empty object, had to add it to prevent printing it
  N <- foreach::foreach(w = write_groups,
                        .export = "flexdog",
                        .combine = c) %dopar% {
                          
                          fout <- lapply(w,function(current_snp){
                            #here we calculate genotypes for groups of 100 markers
                            #m stands for marker
                            refvec <- refmat[current_snp, indlist, drop = TRUE]
                            sizevec <- sizemat[current_snp, indlist, drop = TRUE]
                            
                            if (!is.null(p1_id)) {
                              p1_ref <- refmat[current_snp, p1_id, drop = TRUE]
                              p1_size <- sizemat[current_snp, p1_id, drop = TRUE]
                            } else {
                              p1_ref <- NULL
                              p1_size <- NULL
                            }
                            
                            if (!is.null(p2_id)) {
                              p2_ref <- refmat[current_snp, p2_id, drop = TRUE]
                              p2_size <- sizemat[current_snp, p2_id, drop = TRUE]
                            } else {
                              p2_ref <- NULL
                              p2_size <- NULL
                            }
                            
                            fout <- flexdog(refvec    = refvec,
                                            sizevec   = sizevec,
                                            ploidy    = ploidy,
                                            model     = model,
                                            p1ref     = p1_ref,
                                            p1size    = p1_size,
                                            p2ref     = p2_ref,
                                            p2size    = p2_size,
                                            snpname   = current_snp,
                                            bias_init = bias_init,
                                            verbose   = FALSE,
                                            prior_vec = prior_vec)
                            
                            #We only return that which will be printed out. This saves a lot of memory
                            return(fout[c(snp_pars,snp_ind_pars,snp_ind_geno_pars)])
                          })
                          
                          #We process the results to obtain printable lines
                          #First the snp-lines. All these parameters go into one matrix
                          if(length(snp_pars) > 0){
                            snp <- sapply(snp_pars,function(s) sapply(fout,'[[',s))
                            snp <- t(do.call(rbind,snp))
                            snp <- round(matrix(unlist(snp),ncol = ncol(snp)),4)
                            snp <- cbind(w,snp)
                            write(t(snp),outfiles["model"],append = T,ncolumns = ncol(snp))
                          }
                          
                          #Then the snp x ind matrices
                          if(length(snp_ind_pars) > 0){
                            snp_ind <- lapply(snp_ind_pars,function(s){
                              output_matrix <- sapply(fout,'[[',s)
                              rbind(w,output_matrix)
                            })
                            names(snp_ind) <- snp_ind_pars
                            for(out in names(snp_ind)){
                              write(snp_ind[[out]],outfiles[out],
                                    append = T,
                                    ncol = length(indlist) + 1)
                            }
                          }
                          
                          #Then the snp x ind x geno which we will split into snp x ind matrices
                          if(length(snp_ind_geno_pars) > 0){
                            snp_ind_geno <- lapply(snp_ind_geno_pars,function(s){
                              output_matlist <- lapply(fout,'[[',s)
                              lapply(1:(ploidy+1),function(i){
                                output_1_col <- t(sapply(output_matlist,'[',,i))
                                cbind(w,round(output_1_col,4))
                              })
                            })
                            snp_ind_geno <- unlist(snp_ind_geno,recursive = F)
                            names(snp_ind_geno) <- names_snp_ind_geno
                            for(out in names(snp_ind_geno)){
                              write(t(snp_ind_geno[[out]]),
                                    file = outfiles[out],
                                    append = T,
                                    ncolumns = length(indlist) + 1)
                            }
                          }
                        }
  
  if (nc > 1) {
    parallel::stopCluster(cl)
  }
  
  ## Correct and sort outfiles ----------------------------------
  corrupt_lines <- sapply(outfiles,function(f){
    out_table <- corrupt_read(f)
    new_order <- match(snplist[snplist %in% out_table$marker],out_table$marker)
    out_table <- out_table[new_order,]
    write.table(out_table,f,quote = F, row.names = F)
    return(attr(out_table,"corrupt"))
  })
  
  cat("done!")
  
  return(list(files = outfiles,
              corrupt = corrupt_lines))
}


par_in_model <- function(model){
  parlist <- list(norm=c("mu","alpha"),
                  hw = "alpha",
                  bb=c("alpha","tau"),
                  s1 = c("p1geno","alpha"),
                  f1 = c("p1geno","p2geno","alpha"),
                  s1pp = c("ell1","tau1","gamma1","alpha"),
                  f1pp = c("ell1","ell2","tau1","tau2","gamma1","gamma2","alpha"),
                  flex = NULL,
                  uniform = NULL,
                  custom = NULL)
  return(parlist[[model]])
}

clean_table <- function(file, out = NULL){
  if(is.null(out)) out = paste0(file,"_clean")
  con <- file(file,"r")
  corrupt = c(NA)
  for( i in 1:R.utils::countLines(file)){
    line <- readLines(con,n = 1)
    if(i == 1){
      fields <- length(strsplit(line," ")[[1]])
      cat(paste0(line,"\n"), file = out)
    }else{
      n <- length(strsplit(line," ")[[1]])
      if(n == fields){
        cat(paste0(line,"\n"),file = out,append = T)
      }else{
        corrupt = c(corrupt,i)
      }
    }
  }
  close(con)
  if(length(corrupt) > 1){ corrupt <- corrupt[2:length(corrupt)]}
  attr(corrupt,"msg") <- paste("There were",sum(!is.na(corrupt)),"corrupted lines in",file)
  return(corrupt)
}

corrupt_read <- function(file){
  tab <- tryCatch(as.data.frame(data.table::fread(file)),
                  warning = function(w){
                    corrupt <- clean_table(file, out = paste0(file,"_clean"))
                    system(paste("mv",paste0(file,"_clean"),file))
                    res <- as.data.frame(data.table::fread(file))
                    attr(res,"corrupt") <- corrupt
                    return(res)
                  })
  if(!is.null(attr(tab,"corrupt"))){
    markers_lost <- length(attr(tab,"corrupt"))
  }else{
    markers_lost <- 0
  }
  attr(tab,"corrupt") <- markers_lost
  return(tab)
}



#' Filter SNPs based on the output of \code{\link{multidog}()}.
#'
#' Filter based on provided logical predicates in terms of the variable
#' names in \code{x$snpdf}. This function filters both \code{x$snpdf}
#' and \code{x$inddf}.
#'
#' @param x The output of \code{multidog}.
#' @param expr Logical predicate expression defined in terms of the variables
#'     in \code{x$snpdf}. Only SNPs where the condition evaluates to
#'     \code{TRUE} are kept.
#'
#' @examples
#' \dontrun{
#' data("uitdewilligen")
#' mout <- multidog(refmat = t(uitdewilligen$refmat),
#'                  sizemat = t(uitdewilligen$sizemat),
#'                  ploidy = uitdewilligen$ploidy,
#'                  nc = 2)
#'
#' ## The following filters are for educational purposes only and should
#' ## not be taken as a default filter:
#' mout2 <- filter_snp(mout, bias < 0.8 & od < 0.003)
#' }
#'
#' @seealso
#' \itemize{
#'   \item{\code{\link{multidog}()}:}{For the variables in \code{x$snpdf}
#'       which you can filter by.}
#' }
#'
#' @author David Gerard
#'
#' @export
filter_snp <- function(x, expr) {
  assertthat::assert_that(is.multidog(x))
  cond <- eval(expr = substitute(expr), envir = x$snpdf)
  x$snpdf <- x$snpdf[cond, , drop = FALSE]
  goodsnps <- x$snpdf$snp
  x$inddf <- x$inddf[x$inddf$snp %in% goodsnps, , drop = FALSE]
  return(x)
}




#' Save the output of \code{\link{multidog}()} to a VCF file.
#'
#' This is an experimental function to save the output of
#' \code{\link{multidog}()} to a VCF file.
#'
#' This function uses the Bioconductor packages VariantAnnotation,
#' GenomicRanges, S4Vectors, and IRanges. You can install these using the
#' BiocManager package via:
#'
#' \code{install.packages("BiocManager")}
#'
#' \code{BiocManager::install(c("VariantAnnotation", "GenomicRanges", "S4Vectors", "IRanges"))}
#'
#' To read more about the VCF format, see the official documentation
#' on the Samtools website: \url{http://samtools.github.io/hts-specs/}.
#'
#' @param obj An object of class \code{\link{multidog}()}.
#' @param filename A string, the path to save the VCF file.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(VariantAnnotation)
#' data("snpdat")
#' refmat <- reshape2::acast(data = snpdat,
#'                           formula = snp ~ id,
#'                           value.var = "counts")
#' sizemat <- reshape2::acast(data = snpdat,
#'                            formula = snp ~ id,
#'                            value.var = "size")
#' mout <- multidog(refmat = refmat,
#'                  sizemat = sizemat,
#'                  ploidy = 6,
#'                  model = "s1",
#'                  p1_id = "Xushu18")
#' export_vcf(obj = mout, filename = "./sweet_potato.vcf")
#' spvcf <- readVcf("./sweet_potato.vcf")
#' }
#'
#'
export_vcf <- function(obj, filename) {
  if (requireNamespace("VariantAnnotation", quietly = TRUE) &&
      requireNamespace("GenomicRanges", quietly = TRUE) &&
      requireNamespace("S4Vectors", quietly = TRUE) &&
      requireNamespace("IRanges", quietly = TRUE)) {

    ploidy <- unique(obj$snpdf$ploidy)
    stopifnot(length(ploidy) == 1)
    obj$inddf$alt <- obj$inddf$size - obj$inddf$ref
    geno <- S4Vectors::SimpleList(
      AD = format_multidog(x = obj, varname = c("alt", "ref")),
      DP = format_multidog(x = obj, varname = "size"),
      DS = format_multidog(x = obj, varname = "postmean"),
      GP = format_multidog(x = obj, varname = paste0("Pr_", 0:ploidy)),
      GL = format_multidog(x = obj, varname = paste0("logL_", 0:ploidy)) / log(10)
    )
    for (i in seq_along(geno)) {
      dimnames(geno[[i]]) <- NULL
    }

    nind <- ncol(geno$AD)
    nsnp <- nrow(obj$snpdf)

    infodf <- S4Vectors::DataFrame(
      row.names = c("snp",
                    "bias",
                    "seq",
                    "od",
                    "prop_mis",
                    "llike",
                    "ploidy",
                    "model",
                    paste0("Pr_", 0:ploidy)),
      Description = c("SNP name",
                      "Allele bias",
                      "Sequencing error rate",
                      "Overdispersion",
                      "Proportion of individuals genotyped incorrectly",
                      "Maximized marginal log-likelihood",
                      "Ploidy",
                      "Model used",
                      paste0("Prior probability of dosage ", 0:ploidy)),
      Type = c("String", "Float", "Float", "Float",
               "Float", "Float", "Integer", "String",
               rep("Float", ploidy + 1)),
      Number = 1)

    vcfobj <- VariantAnnotation::VCF(
      rowRanges = GenomicRanges::GRanges(
        seqnames=paste0("snp:", seq_len(nsnp)),
        ranges=NULL,
        strand=NULL,
        seqinfo = NULL,
        names = obj$snpdf$snp),
      colData = as(matrix(nrow = nind, ncol = 0), "DataFrame"),
      exptData = list(
        header = VariantAnnotation::VCFHeader(
          reference = character(),
          samples = character(),
          header = IRanges::DataFrameList(
            fileformat = S4Vectors::DataFrame(row.names = "fileformat", Value = "VCFv4.3"),
            fileDate = S4Vectors::DataFrame(row.names = "fileDate", Value = gsub("-", "", Sys.Date())),
            source = S4Vectors::DataFrame(row.names = "source", Value = paste0("updogv", utils::packageVersion("updog"))),
            FORMAT = S4Vectors::DataFrame(row.names = c("AD", "DP", "DS", "GP", "GL"),
                                          Number = c("R", "1", "1", "G", "G"),
                                          Type = c("Integer", "Integer", "Float", "Float", "Float"),
                                          Description = c("Read depth for each allele",
                                                          "Read depth",
                                                          "Posterior mean genotype",
                                                          "Genotype posterior probabilities",
                                                          "Genotype likelihoods")),
            INFO = infodf
            )
          )
        ),
      fixed = as(matrix(nrow = nsnp, ncol = 0), "DataFrame"),
      geno = geno,
      info = as(obj$snpdf[, row.names(infodf)], "DataFrame"),
      collapsed = TRUE,
      verbose = TRUE
      )

    VariantAnnotation::writeVcf(obj = vcfobj, filename = filename)
  } else {
    message(paste0("Need to have VariantAnnotation, S4Vectors,\n",
                   "GenomicRanges, and IRanges installed to run",
                   "\nexport_vcf()\n",
                   "To install, run in R:\n\n",
                   "install.packages(\"BiocManager\")\n",
                   "BiocManager::install(c(\"VariantAnnotation\",",
                   " \"GenomicRanges\", \"S4Vectors\", \"IRanges\"))\n"))
  }
}

