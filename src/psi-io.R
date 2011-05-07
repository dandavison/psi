structure.read.q.files.2 <- function(dir, outdir.regexp="outdir", q.file="q", labels=NULL) {
    q <- list()
    dirs <- list.files(path=dir, pattern="K[0-9][0-9]?.*")
    K <- as.integer(sub("K([0-9][0-9]?).*", "\\1"))
    cat("found K values:", K, "\n")
    for(K in K) {
        Kdir <- file.path(dir,paste("K", K, sep=""))
        outdirs <- list.files(path=Kdir, pattern=outdir.regexp, full.names=TRUE)
        q[[K]] <- lapply(outdirs, function(d) as.matrix(read.table(file.path(d, q.file), row.names=labels)))
    }
    q
}

structure.read.q.files.1 <- function(dir, outdir.regexp="outdir", q.file="q", labels=NULL) {
    q <- list()
    K <- as.integer(sub("K", "", list.files(path=dir, pattern="K[0-9][0-9]?")))
    for(K in K) {
        Kdir <- file.path(dir,paste("K", K, sep=""))
        outdirs <- list.files(path=Kdir, pattern=outdir.regexp, full.names=TRUE)
        q[[K]] <- lapply(outdirs, function(d) as.matrix(read.table(file.path(d, q.file), row.names=labels)))
    }
    q
}

get.loglike <- function(file)
    scan(pipe(paste("~/bin/loglike", "<", file)))

read.genotypes <- function(file, na.label=0) {
    stop("same as above and what's deal with dimensions and na.label default")
    x <- strsplit(scan(file, what=""), split="")
    L <- length(x)
    n <- sapply(x, length)
    stopifnot(all(n == n[1]))
    n <- n[1]
    x <- structure(as.integer(unlist(x)), dim=c(n,L))
    x[x == 0] <- NA
    x - 1
}

allocation <- function(p, L, nprocs) {
    (L + nprocs - p - 1) %/% nprocs
}

read.hgdp <- function(file="/usr/popgen/jonathan_pritchard/conrad_et_al_data/phased_HGDP_regions1to36") {
    x <- read.table(file, skip=4, as.is=TRUE)
    labels <- x[[5]]
    x <- x[,-(1:7)]
    twon <- nrow(x)
    stopifnot(twon %% 2 == 0)
    L <- ncol(x)
    
    x[x == "?"] <- NA
    x <- lapply(x, function(snp) unclass(factor(snp)) - 1)
    x <- array(unlist(x), dim=c(twon,L))
    x <- t(x)
    colnames(x) <- labels
    x
}

