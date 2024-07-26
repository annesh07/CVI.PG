f0 = function(x){cumsum(x)[length(x)]-cumsum(x)}
f1 = function(x){cumsum(x)*(cumsum(x)[length(x)]-cumsum(x))}
