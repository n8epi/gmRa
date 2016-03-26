# Data Structures

#' List data structure
#'
#' @export
xlist = function() {

  e = new.env(); # Initialize as environment
  e$length = 0;

  e$add = function(item) {
    e$length = e$length + 1;
    n = e$length;
    e[[as.character(n)]] = item;
  }

  e$get = function(n) {
    return(e[[as.character(n)]]);
  }

  e$set = function(n,value) {
    if (n <= e$length) {
      e[[as.character(n)]] = value;
    } else {
      print("ERROR: Set value outside xlist range.")
    }
  }

  e$swap = function(n1,n2) {
    if (n1 <= e$length && n2 <= e$length) {
      v = e[[as.character(n1)]];
      e[[as.character(n1)]] = e[[as.character(n2)]];
      e[[as.character(n2)]] = v;
    } else {
      print("ERROR: Swapping values outside xlist range.")
    }
  }

  e$pop = function() {
    # Bubble the value up and remove
    index = 1;
    while (index < e$length) {
      e$swap(index,index+1);
    }
    index = as.character(e$length);
    x = e[[index]];
    rm(list=index,envir=e);
    e$length = e$length-1;
    return(x);
  }

  e$asvec = function() {
    data = c();
    if (e$length > 0) {
      for(i in 1:e$length) {
        if(is.null(e[[as.character(i)]])) {
          data = c(data,"NULL");
        } else {
          data = c(data,e[[as.character(i)]]);
        }

      }
    }
    return(data);
  }

  e$print = function() {
    print(e$asvec());
  }

  e$map = function(fcn) {
    for(i in 1:e$length) {
      e[[as.character(i)]] = fcn(e[[as.character(i)]]);
    }
  }

  e$fold = function(init,fcn) {
    temp = init;
    for(i in 1:e$length) {
      temp = fcn(temp,e[[as.character(i)]]);
    }
  }

  e$remove = function(val) {
    # Remove particular values from the list by bubbling
    index = 1;
    while (index <= e$length) {
      scanlist = e$asvec();
      if (val == scanlist[index]) {
        # Bubble up the point and then remove
        subIndex = index;
        while (subIndex < e$length) {
          e$swap(subIndex,subIndex+1);
          subIndex = subIndex+1;
        }
        subIndex = as.character(e$length);
        x = e[[subIndex]];
        rm(list=subIndex,envir=e);
        e$length = e$length-1;
      }
      index = index + 1;
      #print("Xlist removal info: val, index, scanlist");
      #print(c(val,index,scanlist));
    }
  }

  return(e);
}

#' Hash Map data structure
#'
#'
#' @export
hmap = function() {

  e = new.env();
  e$keys = xlist();

  e$put = function(key,val) {
    if(is.null(e[[as.character(key)]])) {
      e$keys$add(key);
    }
    e[[as.character(key)]] = val;
  }

  e$apply = function(key,fcn) {
    e[[as.character(key)]] = fcn(e[[as.character(key)]]);
  }

  e$get = function(key) {
    return(e[[as.character(key)]]);
  }

  e$getKeys = function() {
    return(e$keys$asvec());
  }

  e$getVals = function() {
    kyz = e$getKeys();
    vals = c();
    for (i in kyz) {
      vals = c(vals,e$get(i));
    }
    return(vals);
  }

  e$remove = function(key) {
    rm(list=as.character(key),envir=e);
    e$keys$remove(key);
  }

  return(e);

}

#' Priority Queue data structure
#'
#' @export
priorityqueue = function() {

  e = new.env();
  e$size = 0;
  e$values = new.env(); # ordered according to the min heap structure
  e$keys = new.env(); # map from ordered index to the associated key
  e$index = new.env(); # map from keys to the ordered index

  e$getVal = function(i) {
    return(e$values[[as.character(i)]]); # Should do a null check..
  }

  e$getKey = function(i) {
    return(e$keys[[as.character(i)]]);
  }

  e$getIndex = function(key) {
    return(e$index[[as.character(key)]]);
  }

  # Set methods have thin veneer of write protection through null checks

  e$setVal = function(i,val) {
    if (!is.null(e$values[[as.character(i)]])) {
      e$values[[as.character(i)]] = val;
    }
  }

  e$setKey = function(i,key) {
    if (!is.null(e$keys[[as.character(i)]])) {
      e$keys[[as.character(i)]] = key;
    }
  }

  e$setIndex = function(key,i) {
    if (!is.null(e$index[[as.character(key)]])) {
      e$index[[as.character(key)]] = i;
    }
  }

  e$parent = function(n) {
    return(floor(n/2));
  }

  e$left = function(n) {
    return(2*n);
  }

  e$right = function(n) {
    return(2*n+1);
  }

  e$add = function(key,val) {
    e$size = e$size + 1;
    e$index[[as.character(key)]] = e$size;
    e$values[[as.character(e$size)]] = val;
    e$keys[[as.character(e$size)]] = key;
  }

  e$swap = function(i1,i2) {

    k1 = e$getKey(i1);
    k2 = e$getKey(i2);

    e$setIndex(k1,i2);
    e$setIndex(k2,i1);

    e$setKey(i1,k2);
    e$setKey(i2,k1);

    v  = e$getVal(i1);
    e$setVal(i1,e$getVal(i2));
    e$setVal(i2,v);

  }

  e$push = function(key,val) {
    e$add(key,val);
    n=e$size;
    isheap=FALSE;

    while(n>1 && !isheap) {
      pn = e$parent(n);
      if(e$getVal(pn)>val) {
        e$swap(n,pn);
        n = pn;
      } else {
        isheap = TRUE;
      }
    }

  }

  e$pop = function() {

    key = e$getKey(1);

    # Empty check
    if(e$size > 0) {

      # Swap data at the first index with data at the last index
      e$swap(1,e$size);

      # Delete the information at the last index
      k = as.character(e$getKey(e$size));
      s = as.character(e$size);
      rm(list=s,envir=e$values);
      rm(list=s,envir=e$keys);
      rm(list=k,envir=e$index);
      e$size = e$size - 1;

      if (e$size > 0) {
        # Bubble the point downwards
        n=1;
        v = e$getVal(1);
        isheap = FALSE;

        #print("Beginning bubble down...");

        while(!isheap) {
          minn = n;
          minv = v;

          ln = e$left(n);
          lv = e$getVal(ln);
          if(!is.null(lv)) {
            if(lv < minv) {
              minn = ln;
              minv = lv;
            }
          }

          rn = e$right(n);
          rv = e$getVal(rn);
          if(!is.null(rv)) {
            if(rv < minv) {
              minn = rn;
              minv = rv;
            }
          }


          # print(toString(c("Pop data: ",n,v,minn,minv)));
          if (minv < v) {
            e$swap(n,minn);
            n = minn;
            v = minv;
          } else {
            isheap = TRUE;
          }

        }
      }

    }

    return(key);

  }

  e$decrease = function(key,val) {
    n = e$getIndex(key);

    if (is.null(n)) {
      e$push(key,val);
    } else {
      v = e$getVal(n);

      if(val < v) {
        v = val;
        e$setVal(n,v);
        isheap = FALSE;

        while(n>1 && !isheap) {
          pn = e$parent(n);
          if(e$getVal(pn)>v) {
            e$swap(n,pn);
            n = pn;
          } else {
            isheap = TRUE;
          }
        }

      } else {
        print("Key value is not lower.");
      }

    }

  }

  e$printQueue = function() {
    k = c();
    v = c();
    for(i in 1:e$size) {
      k = c(k,e$getKey(i));
      v = c(v,e$getVal(i));
    }
    print("Priority Queue Contents:");
    print(k);
    print(v);
  }

  return(e);
}
