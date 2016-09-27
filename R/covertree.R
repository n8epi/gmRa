

#' Euclidean distance
#'
#' @export
dist2 = function(a,b) { sqrt(crossprod(a-b)) }


#' Index conversion for the dist function
#'
#' @export
distIndex = function(i,j,n) {
  return((n-(i/2))*(i-1) + (j-i));
}

#' Cover Tree
#'
#' TODO: Implement insert, remove, and tree balancing
#' TODO: Accommodate data request model in order to deal with large datasets
#' @param data Matrix of data in rows
#' @param scalefactor Zooming factor for the cover tree
#' @param metric Function used to compute distances
#'
#' @export
covertree = function(data,scalefactor=2,metric=dist2) {
  # The main components of the covertree are the coarsest scale, the finest
  # scale, the root, and the children set indices
  ct = new.env();
  ct$nodes = data; # Doubles the memory footprint
  ct$root = data[1,]; # Root of the cover tree
  ct$ch = hmap2(); # Hash map for children associations
  ct$sf = scalefactor; # The scaling, or "zoom" factor for the cover tree
  ct$dist = metric; # Note that a slight tweak can avoid square roots for Euclidean distances
  # TODO: Implement faster evaluations for p-norm distances
  ct$csc = Inf; # Coursest scale
  ct$fsc = -Inf; # Finest scale

  # Loop through points, inserting into the cover tree
  # TODO: Make this robust to redundant data -- either it is not inserted or loops infinitely
  for (i in 2:dim(data)[1]) {

    print(c("Processing node",i));
    dat = data[i,];
    di = ct$dist(ct$root,dat);

    if (di > 0) {
      sc = floor(-log(di,ct$sf)); # This breaks if d2 = 0 -- redundant points

      # Must start at a slightly lower scale to make sure conflict doesn't
      # occur for insertions at a slightly coarser scale
      if (sc < ct$csc-2) {
        ct$csc = sc;
        ct$ch$put(1,sc,i);
        #print(c("Add node",i,"root",sc));
      } else {
        Qk = c(1);
        Qk1 = c(1);
        sc = ct$csc-2;
        scFac = ct$sf^(-sc);
        parent = 1;
        parentSc = sc;

        duplicate_flag = FALSE

        while (length(Qk1) > 0 && !duplicate_flag) {
          # Loop through the cover to identify any potential parent
          # Note that this is wasteful
          # print('In the Qk1 loop...')
          for (j in 1:length(Qk)) {
            d2j=ct$dist(dat,data[Qk[j],])
            if ( d2j < scFac ) {
              parent = Qk[j];
              parentSc = sc;
            }
          }

          # Update Qk and the scale
          Qk = Qk1;
          sc = sc + 1;
          scFac = scFac / ct$sf;

          # The new Qk1 consists of Qk and children within 2^-sc
          Qk1 = c();
          for (j in 1:length(Qk)) {
            # Check this point of the cover
            d2j = ct$dist(dat,data[Qk[j],]);
            if (d2j < scFac) {
              Qk1 = c(Qk1,Qk[j]);
            }
            if (d2j == 0) {
              #print(c(i,' is duplicate of ',Qk[j],' with index ', j))
              duplicate_flag = TRUE
            }
          }

          # Check the children of the cover
          #print(c("Right before traceback...",Qk,sc));
          children = ct$ch$get(Qk,sc);
          #print(c("Children from get: ",children));
          if(length(children)>0) {
            for (j in 1:length(children)) {
              #print("In the check children loop:");
              #print(children[l]);
              d2j = ct$dist(dat,data[children[j],]);
              if (d2j < scFac) {
                Qk1 = c(Qk1,children[j]);
              }
              if (d2j == 0) {
                #print(c(i,' is duplicate of ',children[j],' with index ',j))
                duplicate_flag = TRUE
              }
            }
          }

        }

        # Only add if we have not duplicated another point...
        if (!duplicate_flag) {
          # Update the children index
          ct$ch$put(parent,parentSc,i);

          # Update the coarsest and finest scale if necessary
          if (parentSc < ct$csc) {
            ct$csc = parentSc;
          }

          if (ct$fsc < parentSc) {
            ct$fsc = parentSc;
          }
        }


        #print(c("Added node",i,parent,parentSc));
      }
    }
  }

  # Nearest neighbor
  ct$nn = function(point,scale) {

    minI = 1;
    minD = ct$dist(ct$root,point);

    if (minD > 0) {

      d = hmap();
      d$put(1,minD);
      cover = c(1);
      factor = ct$sf^(-ct$csc);
      sc = ct$csc;

      while (sc < scale && sc <= ct$fsc && !is.null(cover) ) {
        thresh = minD+factor;
        #print(c("Cover and children in NN scale ",sc));
        #print(cover);
        children = ct$ch$get(cover,sc);
        #print(children);
        newCover = c();

        for (i in cover) {
          if (d$get(i) <= thresh) {
            newCover = c(newCover,i);
          } else {
            d$remove(i);
          }
        }

        for (i in children) {
          dCh = ct$dist(ct$nodes[i,],point);
          if (dCh <= thresh) {
            newCover = c(newCover,i);
            d$put(i,dCh);
            if (dCh < minD) {
              minI = i;
              minD = dCh;
            }
          }
        }

        cover = newCover;
        sc = sc + 1;
        factor = factor/ct$sf;
      }


    }

    return(minI);

  }

  # Construct a metric graph on the cover points based on proximities
  ct$extractMetricGraph = function(sc,fac=2) {
    indices = ct$getIndices(sc);
    #print("Indices found for Metric Graph:");
    #print(indices);
    D = dist(ct$nodes[indices,]);
    #print(D);
    G = wgraph();
    G$V = 1:length(indices);
    thresh = fac*(ct$sf^(-sc));
    for (i in 1:(length(indices)-1)) {
      for (j in (i+1):length(indices)) {
        w = D[distIndex(i,j,length(indices))];
        if (w < thresh) {
          G$addEdge(i,j,w);
        }
      }
    }

    return(G);
  }

  # Get cover indices at a specified scale
  ct$getIndices = function(sc) {
    indices = c(1);
    for (i in ct$csc:(sc-1)) {
      indices = c(indices,ct$ch$get(indices,i));
    }
    return(indices);
  }

  return(ct);

}

