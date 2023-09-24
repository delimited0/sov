r <- matrix(c(4, 3, 2, 1, 
              3, 5, -1, 1,
              2, -1, 4, 2,
              1, 1, 2, 5),
            nrow = 4)
a <- -rep(Inf, 4)
b <- 1:4

sov::mvnxpb(r, a, b)
mvtnorm::pmvnorm(a, b, sigma = r)

d <- 100
r <- .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
a <- rep(0, d)
b <- rep(Inf, d)

sov::mvnxpb(r, a, b)

sov::bvnu(-3, -1, .35)

bvnu(-3, Inf, .35)

bvnu(-Inf, -1, .35)

sov::bvnu(1, 1, .1)



bvnu(-Inf, -Inf, 0)



foo <- Rcpp::cppFunction(
  "
  arma::uvec foo() {
    arma::uvec a = arma::regspace<arma::uvec>(4, 1, 3);  
    return(a);
  }
  ", depends = "RcppArmadillo")

baz <- Rcpp::cppFunction(
  "
  arma::vec foo() {
    arma::mat A = arma::eye(2, 2);
    return A.submat(1, 0, 1, -1);
  }
  ", depends = "RcppArmadillo"
)

# lattice rule comparison ----
lat <- mined::Lattice(10000, p = 2)
