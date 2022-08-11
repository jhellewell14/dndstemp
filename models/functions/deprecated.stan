


// Replicates mydouble type arithmetic from Danny's code
real sillyplus(real x, real y) {
  real a;
  if(x == 0){
    a = y;
  }else if(y == 0){
    a = x;
  }else{
    real diff = x - y;
    if(diff == 0){
      a = log(2.0) + x;
    }else if(diff < 0){
      a = y + log(1.0 + exp(diff));
    }else{
      a = x + log(1.0 + exp(-diff));
    }
  }
  return a;
}


real sillymult(real x, real y) {
  real a;
  // if(y == 0 || x == 0){
    // a = 0;
    // }else{
      a = x + y;
      // }
      return a;
}