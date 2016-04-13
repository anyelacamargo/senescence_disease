alleleLine <- function(a, n)
{
  v <- a / sum(a);
  return(sprintf("allele %s %s", n, paste(as.character(v), collapse = " ")));
}


alleleBlock <- function(s, markerName)
{
  a0 <- as.numeric(s <= 1);
  a1 <- as.numeric(s >= 1);
  if (sum(a0) == 0)
  {
    b <- sprintf("# marker %s: no founder in state 0", markerName);
  }
  else if (sum(a1) == 0)
  {
    b <- sprintf("# marker %s: no founder in state 1", markerName);
  }
  else
  {
    b <- sprintf("marker %s", markerName);
    b <- c(b, alleleLine(rep(1, length(s)), "NA"));
    b <- c(b, alleleLine(a0, "0"));
    b <- c(b, alleleLine(a1, "1"));
  }
  return(b);
}


g2a <- function(gfname, aBasename)
{
  gAll <- read.table(gfname, header = TRUE, stringsAsFactors = FALSE, sep=',');
  chrNameList <- unique(gAll$pos);
  founderNames <- colnames(gAll)[9:ncol(gAll)];
  aList <- list();
  for (chrName in chrNameList)
  {
    g <- gAll[gAll$pos == chrName, ];
    a <- sprintf("markers %d strains %d", nrow(g), length(founderNames));
    a <- c(a, sprintf("strain_names %s", paste(founderNames, collapse = " ")));
    for (i in 1:nrow(g))
    {
      markerName <- g[i, "X"];
      s <- as.numeric(g[i, founderNames]);
      names(s) <- founderNames;
      a <- c(a, alleleBlock(s, markerName));
    }
    afname <- sprintf("chr%s%s", chrName, aBasename);
    write(a, afname);
    aList[[chrName]] <- a;
  }
  return(invisible(aList));
}

a <- g2a("founders.geno.csv", ".alleles");
