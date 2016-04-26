alleleLine <- function(a, n)
{
  v <- round(a / sum(a), 4);
  return(sprintf("allele %s %s", n, paste(as.character(v), collapse = " ")));
}


alleleBlock1 <- function(s, markerName, cM)
{
  b= c();
  alleles <- unique(s);
  if(length(alleles) == 1) 
  {
    b <- sprintf("# marker %s: no polymorphic", markerName);
    return(b);
  }
  else
  {
    b <- sprintf("marker %s %d %f", markerName, length(alleles), cM);
    for (m in alleles)
    {
      
      an = as.numeric(s == m);
      b <- c(b, alleleLine(an, m));
      
    }
  }
  return(b);
}

alleleBlock2 <- function(s, markerName, cM, pos)
{
  alleles <- unique(s);
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
    b <- sprintf("marker %s %d %s %f", markerName, length(alleles)+1, pos, cM);
   
    b <- c(b, alleleLine(a0, "0"));
    b <- c(b, alleleLine(a1, "1"));
    b <- c(b, alleleLine(rep(1, length(s)), "NA"));
    
  }
  return(b);
}


g2a <- function(gfname, aBasename)
{
  
  #gAll <- read.table(gfname, header = TRUE, stringsAsFactors = FALSE, sep=',');
  gAll = gfname;
  chrNameList <- unique(gAll$pos);
  
  founderNames <- colnames(gAll)[11:ncol(gAll)];
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
      cM = as.numeric(g[i, 'cM'])
      pos = g[i, 'pos']
      names(s) <- founderNames;
      a <- c(a, alleleBlock2(s, markerName, cM, pos));
    }
    afname <- sprintf("chr%s%s", chrName, aBasename);
    write(a, afname);
    aList[[chrName]] <- a;
  }
  return(invisible(aList));
}


