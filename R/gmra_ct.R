# GMRA Using Cover Trees


#' Orthogonal projection
#'
#' TODO: Usage of eigen is expensive; fast option needed
#' @param A Symmetric matrix
#' @param dim Number of eigenvectors to retain
#'
#' @export
projector = function(A,dim) {
  # Get the PCA projector up to dim of A
  x = eigen(A,symmetric=TRUE);
  return(crossprod(t(x$vectors[,1:dim])));
}

#' GMRA from cover trees
#'
#' TODO: Usage of eigen is expensive; fast option needed
#' @param ct The cover tree used to build the GMRA
#' @param sc The scale for determing the cover tree partition
#' @param dim The rank of the local projections
#'
#' @export
covertree_GMRA = function(ct,sc,dim) {

  e = new.env();
  e$dim = dim;
  e$ct = ct;
  e$sc = sc;

  v = c(1);
  for (k in ct$csc:(sc-1)) {
    v = c(v,ct$ch$get(v,k));
  }

  e$v = v;
  e$mu = new.env();
  e$P = new.env();
  e$counts = new.env();

  for (i in v) {
    j = as.character(i);
    e$mu[[j]] = matrix(0,1,dim(ct$nodes)[2]);
    e$P[[j]] = matrix(0,dim(ct$nodes)[2],dim(ct$nodes)[2]);
    e$counts[[j]] = 0;
  }

  # first, extract all the NN indices for the nodes in e$ct
  for (i in 1:dim(ct$nodes)[1]) {
    j = as.character(ct$nn(ct$nodes[i,],sc));
    e$mu[[j]] = e$mu[[j]] + ct$nodes[i,];
    e$P[[j]] = e$P[[j]] + crossprod(t(ct$nodes[i,]));
    e$counts[[j]] = e$counts[[j]] + 1;
  }

  # loop through and perform replacements
  for (i in v) {
    j = as.character(i);
    # might consider a zero check here...
    #print("Means and covariances:");
    #print(e$mu[[j]]);
    #print(e$P[[j]]);
    e$P[[j]] = e$P[[j]] - crossprod(e$mu[[j]]);
    e$mu[[j]] = e$mu[[j]]/e$counts[[j]];
    e$P[[j]] = projector(e$P[[j]],dim);
    x = eigen(e$P[[j]]);
    #print(x$values);
    #print(c("Projection Loop in GMRA...",i));
  }

  e$approx = function(x) {
    j = as.character(e$ct$nn(x,e$sc));
    #print("In GMRA approximation: ");
    #print(e$P[[j]]);
    return( ( (x-e$mu[[j]]) %*% e$P[[j]] ) + e$mu[[j]] );
  }

  e$getMu = function(i) {
    j = as.character(i);
    return(e$mu[[j]]);
  }

  e$getP = function(i) {
    j = as.character(i);
    return(e$P[[j]]);
  }

  return(e);
}
