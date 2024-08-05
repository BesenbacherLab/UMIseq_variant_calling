"
\033[1;4;32m
Shearwater calls 
\033[0m

Usage:
  shearwater.R [-h] --vcfall=VCFALL --tempvcfall=TEMPVCFALL --pon=PON --sample=SP --region=BED [--model=MODEL] <pileup>
  
Options:
-h --help                 Show this screen.
--vcfall=VCFALL           Vcf OUTPUT of all positions called for this sample.
--tempvcfall=TEMPVCFALL   VCF all output after filtering but before intersection.
--model=MODEL             The shearwater model to use {AND,OR} [default: AND].
--pon=PON                 A PON RDS file.
--sample=SP               Sample for variant calling.
--region=BED              Interested regions from the bed file.


<pileup>         Pileup file.


Examples:
shearwater.R     --vcfall=out.all.vcf --tempvcfall=out.all.temp.vcf --model=AND --pon='pon.RDS' --sample='B00001' --region='.bed' 'fq_data2/B00001/B00001_pileup'
Requires:
docopt, abind, tidyverse, deepSNV

" -> doc
suppressWarnings(suppressMessages(require(docopt)))
suppressWarnings(suppressMessages(require(abind)))
suppressWarnings(suppressMessages(require(tidyverse)))
suppressWarnings(suppressMessages(require(deepSNV)))

opt <- docopt(doc)

cat("[shearwater.R] Running with arguments:\n")
print(opt)
cat("\n")

model <- opt$model
pon <- readRDS(opt$pon)
pileup <- opt$pileup
vcf_all <-opt$vcfall
tempvcf_all <-opt$tempvcfall
sample <- opt$sample
targets <-read.table(opt$region)

#library("deepSNV")
#library("abind")
#library("tidyverse")

#targets <-read.table("NEW_METHOD_hg38_08feb2016_capture_targets.bed")


new_loadAllData <-function(files, regions, ..., mc.cores = 1){
    #if (class(regions) == "GRanges") {
    #    regions = as.data.frame(regions)[, 1:3]
    #    colnames(regions) = c("chr", "start", "stop")
    #}
    nucleotides = c("A", "T", "C", "G", "DEL", "INS", "N", "a", "t", "c", "g", "del", "ins", "n")
    lengths = regions$stop - regions$start + 1
    rows = sum(lengths)
    beg = cumsum(c(1, lengths[-length(lengths)]))
    end = cumsum(lengths)
    c = mclapply(files, function(f) {
        test.matrix <- matrix(0, ncol = length(nucleotides), 
            nrow = rows)
        for (j in 1:nrow(regions)) {
            test.matrix[beg[j]:end[j], ] = bam2R(f, regions$chr[j], 
                regions$start[j], regions$stop[j], ...)[, nucleotides] ##bam2R: A    T C    G - N INS DEL HEAD TAIL   QUAL...
        }
        mode(test.matrix) <- "integer"
        test.matrix
    }, mc.cores = mc.cores)
    counts = array(0, c(length(files), rows, length(nucleotides)))
    mode(counts) <- "integer"
    for (i in 1:length(files)) counts[i, , ] = c[[i]]
    return(counts)
}

# piles_to_counts ---------------------------------------------------------
piles_to_counts <- function(file, regions, extended = F) {
  # v0.2.1 Use D and d for deletions (not - and _). Added extended (FALSE) which,
  # if TRUE, insertions (I, i) and N (N, n) are also returned.
  # 
  # v0.2.0 Removed dependency on package abind
  #INPUT:
  #files is a list of pileup files generated and named as by pysamstat.
  #regions is a data.frame with three first columns giving chrom, start, end
  #OUTPUT:
  # Array with dim pts x positions x bases where bases are A,T,C,G,D,a,t,c,g,d or
  # (if extended = T) A,T,C,G,D,I,N,a,t,c,g,d,i,n
  #
  # Requires pileup_reader, subset_ranges, dplyr

  #require(dplyr)

  #Make skeleton of positions using the regions argument.
  owarn <- getOption("warn")
  options(warn=-1)
  skeleton <- apply(regions, 1, function(x){
      data.frame('chrom' = x[1], 'pos' = seq.int(as.numeric(x[2]), as.numeric(x[3])), stringsAsFactors = F)
    }) %>%
    bind_rows() %>%
    unique()

  options(warn=owarn)
  pileup_reader <-function(x) {
    data <-read.table(x,header = T, stringsAsFactors = F)
    return (data)
    }

  nucleotide_freq <-
    function(x) {
      if(extended) {
        tmp <-
          pileup_reader(x) %>%
          dplyr::select(chrom, pos, ref,
                        A = A_pp_fwd, `T` = T_pp_fwd, C = C_pp_fwd, G = G_pp_fwd, D = deletions_pp_fwd, I = insertions_pp_fwd, N = N_pp_fwd,
                        a = A_pp_rev, t = T_pp_rev, c = C_pp_rev, g = G_pp_rev, d = deletions_pp_rev,  i = insertions_pp_rev, n = N_pp_rev)
        
      } else {
        tmp <-
          pileup_reader(x) %>%
          dplyr::select(chrom, pos, ref,
                        A = A_pp_fwd, `T` = T_pp_fwd, C = C_pp_fwd, G = G_pp_fwd, D = deletions_pp_fwd,
                        a = A_pp_rev, t = T_pp_rev, c = C_pp_rev, g = G_pp_rev, d = deletions_pp_rev)
      }
      
      # Note subset_ranges fails when nrow(tmp) == 0
      res <-
        tmp %>%
        dplyr::left_join(skeleton, ., by = c('chrom' = 'chrom', 'pos' = 'pos')) %>%
        dplyr::select(-c(1,2,3)) %>% #Just get the counts, delete columns of chrom,pos,ref
        as.matrix()

      #Region positions not sequenced (NA) is set to 0 for all bases
      res[ is.na(res) ] <- 0
      
      return(res)
    }
  res <-nucleotide_freq(file)
  # abind::abind(res, along = 0) #old abind dependent way
  #res <- aperm(array(unlist(freq), c(dim(freq[[1]]), length(freq) ) ), c(3,1,2) )
  res3d <-array(res, dim=c(1, dim(res)[1], dim(res)[2]))
  if(extended) {
    dimnames(res3d) <- list(file, NULL, c('A','T','C','G','D', 'I', 'N', 
                                        'a','t','c','g','d', 'i', 'n') )
    return(res3d)  
  }
  dimnames(res3d) <- list(file, NULL, c('A','T','C','G','-','a','t','c','g','_') )
  return(res3d)

}

#new_bf2Vcf <-function (BF, counts, coordinates, samples = 1:nrow(counts), err = NULL, 
#    mu = NULL, cutoff = 0.05, prior = 0.5){
#    #coordinates <- regions2Coordinates(regions) ##18094*2(chr&pos)
#
#    prior = array(rep(prior, each = length(BF)/length(prior)), 
#        dim = dim(BF)) ##47*18094*7, all 0.5
#    odds = prior/(1 - prior) ##47*18094*7, all 1
#    posterior = BF/(BF + odds) ##47*18094*7, A, T, C, G, D, I, N
#    w = which(posterior < cutoff, arr.ind = TRUE) ##10447*3(dim1, dim2, dim3)--posterior[dim1, dim2, dim3]<0.05
#    w = w[order(w[, 2], w[, 3], w[, 1]), , drop = FALSE] #first sort by dim2(pos), then sort by dim3(nulceotides), finally sort by dim1(sample)
#    if (dim(w)[1] == 0) {
#        w = matrix(1, ncol = 3)
#        isNull = TRUE
#    }
#    else {
#        isNull = FALSE
#    }
#    totCounts = counts[, , 1:7] + counts[, , 8:14] ##change from 5 to 7, from 47*18094*14 to 47*18094*7
#    if (is.null(err)) {
#        err = colSums(totCounts) ##18094*7, sum for the 47 samples
#        ##err[1,] A      T      C      G      D      I      N 
#                ##17     46     11 361410      0      0      4 
#        err = err/(rowSums(err) + .Machine$double.eps) #rowSums(err)--length=18094, [1]=G: 361488
#    }
#    if (is.null(mu)) {
#        mu = (totCounts)/rep(rowSums(counts, dims = 2) + .Machine$double.eps, 7) ##rowSums(counts, dims = 2): 47*18094
#        ##A    T    C    G    D    I    N    a    t    c    g    d    i    n 
#        ##0    0    0 2177    0    0    3    0    0    0 1437    0    0    1 --rowSums=3618
#    }        
#    #ref = factor(apply(err, 1, which.max), levels = 1:7, labels = c("A", "T", "C", "G", "D", "I", "N"))
#
#    u = !duplicated(w[, -1, drop = FALSE])
#    wu = w[u, , drop = FALSE]
#    pp = mapply(function(i, j) {posterior[, i, j]}, wu[, 2], wu[, 3])
#    geno = SimpleList(GT = t(pp < cutoff) + 0, GQ = t(round(-10 * log10(pp))), 
#    BF = signif(t(mapply(function(i, j) {BF[, i, j]}, wu[, 2], wu[, 3])), 3), 
#    VF = signif(t(mapply(function(i, j) {mu[, i, j]}, wu[, 2], wu[, 3])), 3), 
#    FW = t(mapply(function(i, j) {counts[, i, j]}, wu[, 2], wu[, 3])), 
#    BW = t(mapply(function(i, j) {counts[, i, j + 7]}, wu[, 2], wu[, 3])), 
#    FD = t(mapply(function(i, j) {rowSums(counts[, i, 1:7])}, wu[, 2], wu[, 3])), 
#    BD = t(mapply(function(i, j) {rowSums(counts[, i, 8:14])}, wu[, 2], wu[, 3])))
#    v = VCF(rowRanges = GRanges(coordinates$chr[wu[, 2]], 
#                                IRanges(coordinates$pos[wu[, 2]] - (wu[, 3] == 7), 
#                                width = 1 + (wu[, 3] == 7)), ), 
#            fixed = DataFrame(#REF = DNAStringSet(paste(ifelse(wu[, 3] == 7, as.character(ref[wu[, 2] - 1]), ""), ref[wu[, 2]], sep = "")), 
#                              REF = coordinates$ref,
#                              #ALT = do.call(DNAStringSetList, as.list(paste(ifelse(wu[, 3] == 7, as.character(ref[wu[, 2] - 1]), ""), c("A", "T", "C", "G", "D", "I", "N")[wu[, 3]], sep = ""))), QUAL = 1, FILTER = "PASS"), 
#                              ALT = c('A', 'T', 'C', 'G', 'D', 'I', 'N')[ w[,3] ]),
#                              info = DataFrame(ER = select(wu[, -1], err), 
#                              PI = select(wu, prior), 
#                              AF = rowMeans(geno$GT), LEN = 1), 
#                              geno = geno, 
#                              exptData = list(header = scanVcfHeader(system.file("extdata", "shearwater2.vcf", package = "deepSNV"))), 
#            colData = DataFrame(samples = 1:length(samples), row.names = samples), collapsed = TRUE)
#        colnames(v) = samples
#    }
#    metadata(v)$header@samples <- as.character(samples)
#    meta(metadata(v)$header)[[1]]["date", 1] <- paste(Sys.time())
#    if (isNull) 
#        v <- v[NULL]
#    return(v)
#}


# estimateRho_ ------------------------------------------------------------
estimateRho_ = function(counts, truncate = 0.1, rho.min = 1e-4, rho.max = 0.1,  pseudo.rho=.Machine$double.eps){
  #counts is the original array with dims pts x pos x 10 (whereas this is pts x pos x 5 in deepSNV::estimateRho)
  #mu is calculated directly from counts the + and - stranded summed array with dims pts x pos x 5
  #Output is a pos x 5 matrix of dispersions.

  ncol = dim(counts)[3]/2
  x.fw = counts[, , 1:ncol, drop = FALSE] #dim(x.fw)
  x.bw = counts[, , 1:ncol + ncol, drop = FALSE]
  x <- x.fw + x.bw #dim(x) #counts , + and - strands summed dim(x)

  #counts per pos per pt on fw and bw and the sum of the two
  n.fw = rep(rowSums(x.fw, dims = 2), dim(x.fw)[3]) #n.fw dim(rowSums(x.fw, dims = 2))
  n.bw = rep(rowSums(x.bw, dims = 2), dim(x.bw)[3])
  n = array(n.fw + n.bw, dim = dim(x)[1:2]) #Depth on all sites on all pts, aross# dim(n); n[1,1020]

  mu = (x + pseudo.rho)/(rep(n + ncol * pseudo.rho, dim(x)[3])) # dim(mu)

  ix = (mu < truncate) #truncate = 0.1; dim(ix) Boolean array dim(ix)

  n = array(rep(rowSums(x, dims=2), dim(x)[3]), dim=dim(x)[1:2]) #dim(n); dim(x) #Same as above actually
  Xix = colSums(x * ix, dims=1) #dim(Xix) ; Xix[1230,]
  nix = array(rep(n, dim(x)[3]), dim=dim(x)) #dim(nix)
  nix[!ix | nix == 0] = NA
  N = colSums(ix) #dim(N)
  nu = (Xix + pseudo.rho) /  (rep(rowSums(colSums(x, dims=1)) + dim(x)[3]* pseudo.rho, dim(x)[3])) #dim(nu) dim(mu)
  s2 =   N * colSums(nix  * (mu - rep(nu, each = nrow(x)))^2, na.rm=TRUE) / ( (N-1) * colSums(nix, na.rm=TRUE)) #str(rep(nu, each = nrow(x)))
  rho = (N * (s2/nu/(1-nu)) - colSums(1/nix, na.rm=TRUE)) / (N -  colSums(1/nix, na.rm=TRUE))
  rho = pmin(pmax(rho, 0), 1)

  rho = pmin( pmax(rho, rho.min), rho.max)

  rho[is.na(rho)] = rho.min

  return(rho)
}







# new_bbb---------------------------------------------------------------
new_bbb <-function (counts, counts_pon, rho, alternative = "greater", truncate = 0.1, 
    rho.min = 1e-04, rho.max = 0.1, pseudo = .Machine$double.eps, 
    return.value = c("BF", "P0", "err"), model = c("OR", "AND", "adaptive"), 
    min.cov = NULL, max.odds = 10, mu.min = 1e-06, 
    mu.max = 1 - mu.min) 
{
  logbb <- function(x, n, mu, disp) {
      lbeta(x + mu, n - x - mu + disp) - lbeta(mu, disp - mu)
  }
  pseudo.rho = .Machine$double.eps
  model = match.arg(model)
  return.value = match.arg(return.value)
  ncol = dim(counts)[3]/2
  x.fw = counts[, , 1:ncol, drop = FALSE]
  x.bw = counts[, , 1:ncol + ncol, drop = FALSE]
  n.fw = rep(rowSums(x.fw, dims = 2), dim(x.fw)[3])
  n.bw = rep(rowSums(x.bw, dims = 2), dim(x.bw)[3])
  x <- x.fw + x.bw
  n = array(n.fw + n.bw, dim = dim(x)[1:2])
  mu = (x + pseudo.rho)/(rep(n + ncol * pseudo.rho, dim(x)[3]))
  ix = (mu < truncate)
  X = colSums(x, dims = 1)
  bound = function(x, xmin, xmax) {
      x = pmax(x, xmin)
      x = pmin(x, xmax)
      return(x)
  }
  disp = (1 - rho)/rho
  rdisp <- rep(disp, each = nrow(counts))
  mu = (x + pseudo)/(rep(n + ncol * pseudo, dim(x)[3]))
  mu = bound(mu, mu.min, mu.max) * rdisp
  tr.fw = x.fw * ix
  X.fw = rep(colSums(tr.fw, dims = 1), each = nrow(counts)) - tr.fw
  N.fw = rep(colSums(n.fw * ix), each = nrow(counts)) - n.fw * ix
  nu0.fw <- (X.fw + x.fw + pseudo)/(N.fw + n.fw + ncol * pseudo)
  nu0.fw <- bound(nu0.fw, mu.min, mu.max) * rdisp
  mu0.bw <- (x.bw + pseudo)/(n.bw + ncol * pseudo)
  mu0.bw <- bound(mu0.bw, mu.min, mu.max) * rdisp
  nu.fw <- (X.fw + pseudo)/(N.fw + ncol * pseudo)
  nu.fw <- bound(nu.fw, mu.min, mu.max) * rdisp
  tr.bw = x.bw * ix
  X.bw = rep(colSums(tr.bw, dims = 1), each = nrow(counts)) - tr.bw
  N.bw = rep(colSums(n.bw * ix), each = nrow(counts)) - n.bw * ix
  nu0.bw <- (X.bw + x.bw + pseudo)/(N.bw + n.bw + ncol * pseudo)
  nu0.bw <- bound(nu0.bw, mu.min, mu.max) * rdisp
  mu0.fw <- (x.fw + pseudo)/(n.fw + ncol * pseudo)
  mu0.fw <- bound(mu0.fw, mu.min, mu.max) * rdisp
  nu.bw <- (X.bw + pseudo)/(N.bw + ncol * pseudo)
  nu.bw <- bound(nu.bw, mu.min, mu.max) * rdisp
  if (return.value == "err") {
      nu0 <- (X.bw + tr.fw + X.fw + tr.bw + pseudo)/(N.bw + 
          n.bw + N.fw + n.fw + ncol * pseudo)
      nu0 <- bound(nu0, mu.min, mu.max)
      return(list(nu = nu0[1, , ], nu.fw = (nu0.fw/rdisp)[1, , ], nu.bw = (nu0.bw/rdisp)[1, , ], rho = rho))
  }
  rm(tr.fw)
  rm(tr.bw)
  mu = pmax(mu, nu0.fw)
  mu = pmax(mu, nu0.bw)
  mu0.fw = pmax(mu0.fw, nu0.fw)
  mu0.bw = pmax(mu0.bw, nu0.bw)
  if (model %in% c("OR", "adaptive")) {
      Bf.fw <- logbb(x.fw, n.fw, nu0.fw, rdisp) + logbb(x.bw, 
          n.bw, mu0.bw, rdisp) + logbb(X.fw, N.fw, nu0.fw, 
          rdisp) - logbb(x.fw, n.fw, mu, rdisp) - logbb(x.bw, 
          n.bw, mu, rdisp) - logbb(X.fw, N.fw, nu.fw, rdisp)
      Bf.fw = exp(Bf.fw)
      Bf.both = logbb(x.fw, n.fw, nu0.fw, rdisp) + logbb(X.fw, 
          N.fw, nu0.fw, rdisp) - logbb(x.fw, n.fw, mu, rdisp) - 
          logbb(X.fw, N.fw, nu.fw, rdisp)
      rm(X.fw, N.fw, mu0.bw, nu.fw)
      Bf.bw <- logbb(x.fw, n.fw, mu0.fw, rdisp) + logbb(x.bw, 
          n.bw, nu0.bw, rdisp) + logbb(X.bw, N.bw, nu0.bw, 
          rdisp) - logbb(x.fw, n.fw, mu, rdisp) - logbb(x.bw, 
          n.bw, mu, rdisp) - logbb(X.bw, N.bw, nu.bw, rdisp)
      Bf.bw = exp(Bf.bw)
      Bf.both = Bf.both + logbb(x.bw, n.bw, nu0.bw, rdisp) + 
          logbb(X.bw, N.bw, nu0.bw, rdisp) - logbb(x.bw, n.bw, 
          mu, rdisp) - logbb(X.bw, N.bw, nu.bw, rdisp)
      Bf.both = exp(Bf.both)
      rm(X.bw, N.bw, mu0.fw, nu.bw)
      rm(mu, nu0.fw, nu0.bw)
      Bf = Bf.fw + Bf.bw - Bf.both + .Machine$double.xmin
  }
  else {
      Bf.both = logbb(x.fw, n.fw, nu0.fw, rdisp) + logbb(X.fw, 
          N.fw, nu0.fw, rdisp) - logbb(x.fw, n.fw, mu, rdisp) - 
          logbb(X.fw, N.fw, nu.fw, rdisp) + logbb(x.bw, n.bw, 
          nu0.bw, rdisp) + logbb(X.bw, N.bw, nu0.bw, rdisp) - 
          logbb(x.bw, n.bw, mu, rdisp) - logbb(X.bw, N.bw, 
          nu.bw, rdisp)
      Bf = exp(Bf.both)
  }
  if (model == "adaptive") {
      if (!is.null(min.cov)) 
          ix <- n.fw < min.cov | n.bw < min.cov
      else ix <- na.omit(abs(log10(n.fw/n.bw)) > log10(max.odds))
      Bf[ix] <- Bf.both[ix]
  }
  ###cons = apply(X, 1, which.max)
  ###for (i in 1:ncol(Bf)) Bf[, i, cons[i]] = NA
  if (return.value == "P0") {
      return(Bf/(1 + Bf))
  }
  else {
      return(Bf)
  }
}



# bbb_to_df ---------------------------------------------------------------
bbb_to_df <- function(BF, counts, prior=0.5, coordinates, sample, cutoff = 0.05, rho, pseudo.rho=.Machine$double.eps) {
  #Input:
  #Output:
  #Returns NULL when no call is found below cutoff

  prior_a = array(rep(prior, each = length(BF)/length(prior)), dim=dim(BF)) ##47*18094*7, all 0.5
  odds = prior_a/(1-prior_a) ##47*18094*7, all 1
  posterior = BF / (BF + odds) ##47*18094*7, A, T, C, G, D, I, N
  w = which(posterior < cutoff, arr.ind=TRUE) ##10447*3(dim1, dim2, dim3)--posterior[dim1, dim2, dim3]<0.05

  if(length(w)==0) { return(NULL) }

  #We want to return a data.frame with the following columns: 1) sample_name, 2-3)chr:pos,
  #4) ref, 5) mut, 6)maf, 7) depth, 8) posterior, 9)
  res <- data.frame('sample_name' = sample[ w[,1] ],  #sample_names[ w[,1] ],
                    'chr' = coordinates$chr[ w[,2] ],
                    'pos' = coordinates$pos[ w[,2] ],
                    'ref' = coordinates$ref[ w[,2] ],
                    'mut' = NA,
                    'maf' = NA,
                    'depth' = NA,
                    'supportive_counts' = NA,
                    'prior' = NA,
                    'bf' = NA,
                    'posterior' = NA,
                    'rho' = NA,
                    stringsAsFactors = F)


  #Make depth array x from counts
  ncol = dim(counts)[3]/2 ##=7
  x.fw = counts[, , 1:ncol, drop = FALSE] #dim(x.fw)=47*18094*7
  x.bw = counts[, , 1:ncol + ncol, drop = FALSE]
  x <- x.fw + x.bw #dim(x) #counts , + and - strands summed dim(x) 47*18094*7

  res$depth <-
    apply(X = w, 1, function(z){
      sum(x[ as.numeric(z[1]), as.numeric(z[2]),])
    })

  #Make a maf array mu
  #n.fw = rep(rowSums(x.fw, dims = 2), dim(x.fw)[3])
  #n.bw = rep(rowSums(x.bw, dims = 2), dim(x.bw)[3])
  #n = array(n.fw + n.bw, dim = dim(x)[1:2]) #Depth on all sites on all pts; dim(n) pts x pos 47*18094
  n.fw = rowSums(x.fw, dims = 2)
  n.bw = rowSums(x.bw, dims = 2)
  n = n.fw + n.bw

  mu = (x + pseudo.rho)/(rep(n + ncol * pseudo.rho, dim(x)[3])) # dim(mu)
  mu[ mu < pseudo.rho ] = 0 #We want to report 0 and for those with no counts

  res$supportive_counts <- 
    apply(X = w, 1, function(y){
      x[ as.numeric(y[1]), as.numeric(y[2]), as.numeric(y[3]) ]
    })


  res$maf <-
    apply(X = w, 1, function(x){
      mu[ as.numeric(x[1]), as.numeric(x[2]), as.numeric(x[3]) ]
    })

  res$bf <-
    apply(X = w, 1, function(x){
      BF[ as.numeric(x[1]), as.numeric(x[2]), as.numeric(x[3]) ]
    })

  res$posterior <-
    apply(X = w, 1, function(x){
      posterior[ as.numeric(x[1]), as.numeric(x[2]), as.numeric(x[3]) ]
    })

  res$rho <-
    apply(X = w, 1, function(x){
      rho[ as.numeric(x[2]), as.numeric(x[3]) ]
    })

  res$prior <-
    apply(X = w, 1, function(x){
      prior_a[ as.numeric(x[1]), as.numeric(x[2]), as.numeric(x[3]) ]
    })

  
  #u = !duplicated(w[, -1, drop = FALSE]) #regardless of dim1(sample pileup)
  #wu = w[u, , drop = FALSE]
  
  res$mut <-c('A', 'T', 'C', 'G', 'D')[ w[,3] ]
  filtered_res <-res %>% filter(mut!=ref & sample_name==sample) %>% mutate(site_mut=paste0(chr, ":", pos, "_", ref, "/", mut))
  call_res <-merge(filtered_res, mut_list, by.x = "site_mut", by.y = "sitemut")
  return(call_res)
}

colnames(targets) <-c("chr", "start", "stop", "Gene")
regions <-makeGRangesFromDataFrame(targets, seqnames.field = names(targets)[1],
                                    start.field = names(targets)[2],
                                    end.field = names(targets)[3], 
                                    keep.extra.column=TRUE)
count_sample <-piles_to_counts(pileup, targets) ##1*18094*14
#pon <-readRDS('210617_v0.1.2-hg38-novaseq-pon.RDS')


##count_pon <-abind(pon$pon[,,1:5], pon$pon[,,8:12],along=3) ##46*18094*14
##coordinates <-pon$coordinates
counts <-abind(count_sample, pon, along=1) ##47*18094*14, A    T    C    G    D    I    N    a    t    c    g    d    i    n
rho0 <-estimateRho_(pon)
#BF_or <-new_bbb(counts=counts, model="OR", rho=rho0) ##bbb function: if rho is null, would estimate rho by itself; 47*18094*7, A T C G D I N
BF <-new_bbb(counts=counts, model=model, rho=rho0) ##bbb: set the highest BF among 7 nucleotides in each position for each sample to NA--47*18094NA ??????
#or_df <-bbb_to_df(BF=BF_or, counts=counts, coordinates=coordinates, rho=rho0, sample_names='B00001', cutoff=0.05, prior=0.5) ##posterior < cutoff to call variants
###df_file <-bbb_to_df(BF=BF, counts=counts, coordinates=coordinates, rho=rho0, sample=sample, cutoff=0.05, prior=0.5)
#write.table(or_df, "or.df", row.names=FALSE, sep=" ")
###write.table(df_file, df, row.names=FALSE, sep=" ")

#vcf_res1 <-bf2Vcf(BF=BF, counts=counts, regions=regions, cutoff=0.05, samples=dimnames(counts)[[1]])
vcf_res2 <-bf2Vcf(BF=BF, counts=counts, regions=regions, cutoff=Inf, samples=dimnames(counts)[[1]])
#writeVcf(vcf_res1, vcf)
writeVcf(vcf_res2, vcf_all)

##modified shearwater output
######cutoff=0.05
#vcf_df <-read.table(vcf)
#var <-paste0("*", sample, "*")
bed <-targets %>% 
  rowwise() %>%
  mutate(POS = list(start:stop)) %>%
  unnest(cols = c("POS")) %>% 
  mutate(pos=paste0(chr,':',POS))

#filtered_vcf_df <-vcf_df[ grepl(var, vcf_df$V3),1:10] %>% filter(V4 !=V5) %>% filter(nchar(V4)==1) %>% mutate(site_mut=paste0(V1,":", V2, "_",V4, "/", V5 )) %>% filter(paste0(V1,":", V2) %in% bed$pos)
#write.table(filtered_vcf_df, tempvcf, row.names=FALSE, sep=" ", quote = FALSE)
#call_vcf_df <-merge(filtered_vcf_df, mut_list, by.x = "site_mut", by.y = "sitemut")
#colnames(call_vcf_df) <-c('site_mut','chr', 'pos', 'sample', 'ref', 'mut', 'qual', 'filter', 'info', 'format1', 'format2', 'pt_id', 'maf', 'vaf')
#call_vcf_df <-call_vcf_df[!duplicated(call_vcf_df$site_mut),] #%>% mutate(vf = str_split_fixed(format2, ":", 8) [,4])
#write.table(call_vcf_df, fvcf, row.names=FALSE, sep=" ", quote = FALSE)

######cutoff=1
ref_pileup <-read.table(pileup, sep="\t", header=TRUE) %>% mutate(place=paste0(chrom, ':', pos))
vcf_all_df <-read.table(vcf_all) %>% mutate(place=paste0(V1,':',V2)) %>% filter(nchar(V4)==1)
vcf_all_df_new <-merge(vcf_all_df, ref_pileup, by.x=c('place'), by.y=c('place')) %>% dplyr::select(V1, V2, V3, V4=ref, V5, V6, V7, V8, V9, V10)
filtered_vcf_all <-vcf_all_df_new %>% filter(V4!=V5, V1!='chrY') %>% mutate(site_mut=paste0(V1,":", V2, "_",V4, "/", V5 )) %>% filter(paste0(V1,":", V2) %in% bed$pos)#%>% mutate(vf = str_split_fixed(V10, ":", 8) [,4])
write.table(filtered_vcf_all, tempvcf_all, row.names=FALSE, sep=" ", quote = FALSE)

