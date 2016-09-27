# Examples

#' Uniform random samples from unit sphere
#'
#' @param n Number of samples
#' @param dim Dimension of the Euclidean space
#'
runifsphere = function(n,dim=2) {
  X = matrix(0,n,dim);
  for (i in 1:n) {
    X[i,] = rnorm(dim);
    X[i,] = X[i,] / norm(as.matrix(X[i,]),'f');
  }
  return(X);
}

#' Testing function for cover trees
#'
#' TODO: Usage of eigen is expensive; fast option needed
#'
#'
#' @export
test_covertree = function() {

  print("Testing cover trees...");

  testdata1 = cbind(c(1,-1/32,1/12, 1/16,1/8));
  ct1 = covertree(testdata1);
  print("Test 1 cover tree structure:");
  ct1$ch$print(ct1$csc,ct1$fsc+1);
  print("Nearest neighbor test 1:");
  print(c(ct1$nn(1/2,ct1$fsc+1),ct1$nn(1/4,ct1$fsc+1)));

  testdata2 = rbind(c(1,1),c(0,1/2),c(1/4,3/4),c(0,1/4),c(0,1/8));
  ct2 = covertree(testdata2)
  print("Test 2 cover tree structure:");
  ct2$ch$print(ct2$csc,ct2$fsc);
  print("Nearest neighbor test 2:");
  print(c(ct2$nn(c(0.9,0.9),ct1$fsc+1),ct2$nn(c(0,1/4),ct1$fsc+1)));



  scalefactor=3 # Test different scaling factors
  numTest3 = 512;
  testdata3 = cbind(runif(numTest3),runif(numTest3));
  # plot the test data coming in
  emptyplot(c(0,1),frame.plot=TRUE,main="Samples");
  points(testdata3[,1],testdata3[,2],pch=20);

  ct3 = covertree(testdata3,scalefactor = scalefactor);
  print("Test 3 cover tree structure:");
  #ct3$ch$print(ct3$csc,ct3$csc+2);

  print("Plotting...");
  # plot each level with separating and covering circles, and the level below it
  scale = 3;
  sc = ct3$csc;
  cover = c(1);
  for (sc in ct3$csc:scale) {
    emptyplot(c(0,1),frame.plot = TRUE,main=sc);
    newCover = cover;
    #print(c("Cover:",cover,"Scale:",sc));

    # Plot each point in the cover, it's circle, and it's children, and add the newCover
    for (i in 1:length(cover)) {
      pt = ct3$nodes[cover[i],];
      plotcircle(mid=pt,r=scalefactor^(-sc));
      points(pt[1],pt[2],pch=20);
      ch = ct3$ch$get(cover[i],sc);
      #print(c("Children of ",cover[i],":",ch));
      for (j in 1:length(ch)) {
        #print("The trigger:");
        #print(ch[[j]])
        pt = ct3$nodes[ch[j],];
        points(pt[1],pt[2]);
      }
      newCover = c(newCover, ch);
    }

    cover = newCover;
  }

  # Get the metric graph at scale 3 and plot
  G = ct3$extractMetricGraph(scale);
  indices = ct3$getIndices(scale);

  # Plot this graph
  emptyplot(c(0,1),frame.plot = TRUE,main="Metric Graph");

  # First, plot the edges
  for (a in G$E) {
    index = scan(text=a,sep=",");
    i = indices[index[1]];
    j = indices[index[2]];
    lines(c(testdata3[i,1],testdata3[j,1]),c(testdata3[i,2],testdata3[j,2]));
  }

  # and the vertices
  for (i in G$V) {
    pt = ct3$nodes[indices[i],];
    points(pt[1],pt[2],pch=20);
  }

  # Now get the shortest path tree emanating from the root and display
  parents = shortest_path_tree(G,G$getMaxStationary());
  #parents$print()

  # Plot this graph
  emptyplot(c(0,1),frame.plot = TRUE,main="Shortest Path Tree");

  # First, plot the edges
  for (a in 1:length(G$V)) {
    i = indices[a];
    j = indices[parents$get(a)];
    lines(c(testdata3[i,1],testdata3[j,1]),c(testdata3[i,2],testdata3[j,2]));
  }

  # and the vertices
  for (i in G$V) {
    pt = ct3$nodes[indices[i],];
    points(pt[1],pt[2],pch=20);
  }

}

#' Testing function for cover trees in 3D
#'
#' @export
test_covertree3D = function() {

  N = 1024;
  data = runifsphere(N,3);
  df = data.frame(x=data[,1],y=data[,2],z=data[,3]);

  ct = covertree(data);

  sc = 0;
  ind = ct$getIndices(sc);
  ch = ct$ch$get(ind,sc);


  print("Extracting metric graph...");
  G = ct$extractMetricGraph(sc);

  print("Finding a graph center...");
  u = G$getMaxStationary();

  print("Extracting shortest path tree...")
  sptParents = shortest_path_tree(G,u);

  print("Building edge matrix for MG plot...");
  # Build the edge matrix for the metric graph
  edgeDataMG = matrix(0,2*length(G$E),3);
  cnt = 0;
  for (a in G$E) {
    index = scan(text=a,sep=","); # Prints out warnings about how many items are scanned
    i = index[1];
    j = index[2];
    edgeDataMG[2*cnt+1,] = data[ind[i],];
    edgeDataMG[2*cnt+2,] = data[ind[j],];
    cnt=cnt+1;
  }
  dfMG = data.frame(x=edgeDataMG[,1],y=edgeDataMG[,2],z=edgeDataMG[,3])

  print("Building edge matrix for SPT plot...");
  # And the shortest path tree
  edgeDataSPT = matrix(0,2*length(sptParents)-2,3);
  cnt = 0;
  for (index in 1:length(sptParents)) {
    if (!is.null(sptParents$get(index))) {
      i = ind[index];
      j = ind[sptParents$get(index)];
      edgeDataSPT[2*cnt+1,] = data[i,];
      edgeDataSPT[2*cnt+2,] = data[j,];
      cnt=cnt+1;
    }
  }
  dfSPT = data.frame(x=edgeDataSPT[,1],y=edgeDataSPT[,2],z=edgeDataSPT[,3])

  print("Constructing the GMRA...");
  gmra = covertree_GMRA(ct,sc,2);

  print("Projecting the data");
  pData = matrix(0,N,3);
  for (i in 1:N) {
    pData[i,] = gmra$approx(data[i,]);
  }
  dfGMRA = data.frame(x=pData[,1],y=pData[,2],z=pData[,3]);

  print("Plotting the scale of interest...");
  open3d();
  par3d("windowRect" = c(0,0,1000,1000));
  plot3d(df[ind,],type="s",radius=0.03,col='blue');
  plot3d(df[ind,],type="s",radius=2^(-sc),alpha=0.1,col='gray',add=TRUE);
  plot3d(df[ch,],type="s",radius=0.02,col="gray",add=TRUE);
  rgl.pop('lights');
  light3d(theta=45,phi=45,ambient='white',diffuse='white');

  print("Plotting the metric graph...");
  open3d();
  par3d("windowRect" = c(0,0,1000,1000));
  plot3d(df[ind,],type="s",radius=0.03,col='blue');
  plot3d(df[ind[u],],type="s",radius=0.05,col="red",add=TRUE);
  segments3d(dfMG,add=TRUE);
  rgl.pop('lights');
  light3d(theta=45,phi=45,ambient='white',diffuse='white');

  print("Plotting the shortest path tree...");
  open3d();
  par3d("windowRect" = c(0,0,1000,1000));
  plot3d(df[ind,],type="s",radius=0.03,col='blue');
  plot3d(df[ind[u],],type="s",radius=0.05,col="red",add=TRUE);
  segments3d(dfSPT,add=TRUE);
  rgl.pop('lights');
  light3d(theta=45,phi=45,ambient='white',diffuse='white');

  cols = rainbow(256)

  print("Plotting the GMRA approximations...");
  open3d();
  par3d("windowRect" = c(0,0,1000,1000));
  #plot3d(dfGMRA,type="s",radius=0.05,col='gray');
  plot3d(dfGMRA,type="s",radius=0.05,col=cols[cut(256*dfGMRA$x,256)]);
  rgl.pop('lights');
  light3d(theta=45,phi=45,ambient='white',diffuse='white');

  f = gmra$regression(matrix(dfGMRA$x,nrow=length(dfGMRA$x),ncol=1))
  regVals = apply(dfGMRA,1,f)

  print("Plotting the GMRA approximations with linear regression test...");
  open3d();
  par3d("windowRect" = c(0,0,1000,1000));
  #plot3d(dfGMRA,type="s",radius=0.05,col='gray');
  plot3d(dfGMRA,type="s",radius=0.05,col=cols[cut(256*regVals,256)]);
  rgl.pop('lights');
  light3d(theta=45,phi=45,ambient='white',diffuse='white');


  print(c('Sum of Squares Regression Error:', sum((dfGMRA$x-regVals)^2)))
}
