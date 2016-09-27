# GMRA Using Cover Trees


#' Orthogonal projection
#'
#' TODO: Implement randomized low-rank SVD approach
#' @param A matrix with data in rows
#' @param The shift value to be subtracted from each row of X
#' @param The rank of the projection to be extracted
#'
#' @export
projector = function(X,mu,dim) {
  x = svd(t(apply(X,1,function(x){x-mu})),nu=0,nv=dim)
  return(x$v);
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
  e$ambient_dim = dim(ct$nodes)[2];
  e$dim = dim;
  e$ct = ct;
  e$sc = sc;

  e$v = ct$getIndices(sc); # Indices of the cover at scale sc; induces the Voronoi cells
  e$mu = new.env(); # Local means
  e$Q = new.env(); # Local projections for encoding
  e$P = new.env(); # Local orthogonal projections for approximations
  e$counts = new.env(); # Counts of data in each cell

  # Initialize information for each partition cell
  print('Initializing local projection information...')
  for (i in e$v) {
    j = as.character(i);
    e$mu[[j]] = matrix(0,1,e$ambient_dim);
    e$P[[j]] = matrix(0,e$ambient_dim,e$ambient_dim);
    e$nnIndices[[j]] = vector();
    e$counts[[j]] = 0;
  }

  # first, extract all the NN indices for the nodes in e$ct
  print('Determining indices of cell members...')
  for (i in 1:dim(ct$nodes)[1]) {
    j = as.character(ct$nn(ct$nodes[i,],sc));
    e$nnIndices[[j]] = c(e$nnIndices[[j]],i);
  }

  # loop through and perform replacements
  print('Performing local svds...')
  for (i in e$v) {
    j = as.character(i);
    e$counts[[j]] = length(e$nnIndices[[j]])
    e$mu[[j]] = colSums(ct$nodes[e$nnIndices[[j]],])/e$counts[[j]] # This needs a zero check
    e$Q[[j]] = projector(ct$nodes[e$nnIndices[[j]],],e$mu[[j]],dim)
    e$P[[j]] = e$Q[[j]] %*% t(e$Q[[j]]);
  }

  # Encode function returns index of cell and coordinates of the projection
  e$encode = function(x) {
    i = e$ct$nn(x,e$sc)
    j = as.character(i);
    return(list(j,(x - e$mu[[j]]) %*% e$Q[[j]]))
  }

  # Approximation function projects onto the local affine approximation
  e$approx = function(x) {
    j = as.character(e$ct$nn(x,e$sc));
    return( ( (x-e$mu[[j]]) %*% e$P[[j]] ) + e$mu[[j]] );
  }

  # For regression, emit a closure given responses Y (as rows) from the ct$nodes values
  # The regression function is linear in each cell
  e$regression = function(Y) {
    betas = new.env()
    for (i in e$v) {
      j = as.character(i);
      coded_data = t(apply(e$ct$nodes[e$nnIndices[[j]],],1,function(x) {return(x-e$mu[[j]])})) %*% e$Q[[j]]
      design = cbind(matrix(1,nrow=e$counts[[j]],ncol=1),coded_data)
      betas[[j]] = solve(t(design) %*% design, t(design) %*% Y[e$nnIndices[[j]],])
    }
    f = function(x) {
      p = e$encode(x)
      return(cbind(1,p[[2]]) %*% betas[[p[[1]]]] )
    }
    return(f)
  }

  e$getMu = function(i) {
    j = as.character(i);
    return(e$mu[[j]]);
  }

  e$getP = function(i) {
    j = as.character(i);
    return(e$P[[j]]);
  }

  print('GMRA construction complete.')

  return(e);
}
