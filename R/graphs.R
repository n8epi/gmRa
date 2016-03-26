# Graphs

#' Weighted Graph
#'
#' @export
wgraph = function() {

  # simple weighted graph on vertices 1:length(e$V)

  e = new.env();
  e$V = c();
  e$E = c();

  e$addEdge = function(i,j,w) {
    key = toString(c(i,j));
    e$E = c(e$E,key);
    e[[key]] = w;
    e[[toString(c(j,i))]] = w;
    e[[as.character(i)]] = c(e[[as.character(i)]],j);
    e[[as.character(j)]] = c(e[[as.character(j)]],i);
  }

  e$getWeight = function(i,j) {
    return(e[[toString(c(i,j))]]);
  }

  e$getDiffusionMatrix = function() {
    Z = matrix(0,length(e$V),length(e$V));
    for (i in e$V) {
      s = 0;
      for (j in e[[as.character(i)]]) {
        Z[i,j] = 1/(sqrt(e$getWeight(i,j)));
        s = s + Z[i,j];
      }
      for (j in e[[as.character(i)]]) {
        Z[i,j] = Z[i,j]/s;
      }
    }
    return(Z);
  }

  e$getMaxStationary = function(maxIter=100) {
    Z = t(e$getDiffusionMatrix());
    u = matrix(1,length(e$V),1);
    up = Z %*% u;
    r = crossprod(u-up);
    iter = 0;
    while (r > 1e-18 && iter<maxIter) {
      u = up;
      up = Z %*% u;
      r = crossprod(u-up);
      iter = iter+1;
    }

    #print("Approx Stationary:");
    #print(u);
    #print(sum(u));

    return(which.max(u));

  }

  e$getVertices = function() {
    verts = c();
    for (v in e$V) {
      verts = c(verts,as.integer(v));
    }
    return(verts);
  }

  e$getEdges = function() {
    edges = new.env();
    edges$v0 = c();
    edges$v1 = c();
    for (a in e$E) {
      index = scan(text=a,sep=",");
      i = index[1];
      j = index[2];
      edges$v0 = c(edges$v0,i);
      edges$v1 = c(edges$v1,j);
    }
    return(edges);
  }

  return(e);

}


#' Shortest path tree
#'
#' @export
shortest_path_tree = function(G,v) {
  # Get the shortest path tree with root v in the weighted symmetric graph G
  d = c();
  parents = xlist();
  pq = priorityqueue();

  for(i in G$V) {
    d = c(d,Inf);
    parents$add(NULL);
    pq$push(i,Inf);
  }

  d[v] = 0;
  pq$decrease(v,0);

  iter = 0;

  while (pq$size > 0) {

    #print(toString(c("At iteration ",iter)));
    #print("Before pop:");
    #pq$printQueue();
    u = pq$pop();
    #print("After pop: ");
    #pq$printQueue();

    for (i in G[[as.character(u)]]) {
      q = d[u] + G$getWeight(u,i);
      #print(c(u,i,q,d[i],length(d)));
      #print(G$V);
      if (q < d[i]) {
        d[i] = q;
        parents$set(i,u);
        pq$decrease(i,q);
      }
    }

    #print("After decrease in the loop: ");
    #pq$printQueue();
    iter = iter+1;
  }

  print("Shortest Path Parents from the xlist:");
  print(parents$asvec());

  return(parents);

}

###############
# Hash map with keys having two values
###############
hmap2 = function() {
  e = new.env();
  e$keys = xlist();

  e$put = function(index,scale,n) {
    i = as.character(index);
    s = as.character(scale);
    if(is.null(e[[i]])) {
      e[[i]] = hmap();
    }

    if(is.null(e[[i]][[s]])) {
      x = xlist();
      x$add(n);
      e[[i]]$put(scale,x);
    } else {
      e[[i]][[s]]$add(n);
    }
  }

  e$get = function(indices,scale) {
    #print("Attempting get...");
    g = c();
    s = as.character(scale);
    for (i in 1:length(indices)) {
      index = as.character(indices[i]);
      #print("Before asvec...");
      if(!is.null(e[[index]][[s]])) {
        g = c(g,e[[index]][[s]]$asvec());
      }
      #print("After asvec...");
    }

    #print(c("Inside get:",g));
    return(g);
  }

  e$print = function(start,end) {
    cover = c(1);
    for(sc in start:end) {
      print(c("Scale:",sc));
      for(i in cover) {
        print(c(i," ",e$get(i,sc)));
      }
      cover = c(cover,e$get(cover,sc));
    }

  }

  return(e);
}

#' Test for shortest path
#'
#'
test_dijkstra = function() {

  G1 = wgraph();
  G1$V = 1:8;
  G1$addEdge(1,2,1);
  G1$addEdge(2,3,1);
  G1$addEdge(3,4,1);
  G1$addEdge(4,5,1);
  G1$addEdge(5,6,1);
  G1$addEdge(6,7,1);
  G1$addEdge(7,8,1);
  G1$addEdge(8,1,1);

  p1 = shortest_path_tree(G1,1);

  # Print the parents
  print("Parents from the first test:");
  print(p1$asvec());

  G2 = wgraph();

  # 1 2 3
  # 4 5 6
  # 7 8 9

  G2$V = 1:9;

  # the perimeter
  G2$addEdge(1,2,1);
  G2$addEdge(2,3,1);
  G2$addEdge(3,6,1);
  G2$addEdge(6,9,1);
  G2$addEdge(9,8,1);
  G2$addEdge(8,7,1);
  G2$addEdge(7,4,1);
  G2$addEdge(4,1,1);

  # inside cross
  G2$addEdge(2,5,1);
  G2$addEdge(4,5,1);
  G2$addEdge(8,5,1);
  G2$addEdge(6,5,1);

  p2 = shortest_path_tree(G2,1);

  # Print the parents
  print("Parents from the second test:");
  print(p2$asvec());

}

