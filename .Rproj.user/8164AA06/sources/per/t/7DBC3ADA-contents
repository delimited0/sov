#include "RcppArmadillo.h"

// singularity tolerance
const double ep = 1e-10;
// [[Rcpp::depends(RcppArmadillo)]]

// Rcpp::List qscmvnv(int m, arma::mat r, arma::vec a, arma::mat cn, arma::vec b) {
//   
//   Rcpp::List chlrsrt_result = chlsrt(r, a, cn, b);
//   
//   arma::vec as = chlsrt_chlsrt;
//   arma::vec bs = chlsrt_result;
//   
//   double ci = R::pnorm(as(0));
//   double dci = R::pnorm(bs(0)) - ci;
//   double p = 0.0;
//   double e = 0.0;
//   
//   int ns = 12;
//   int nv = std::trunc( std::max( (double)m / ns, 1.0) );
//   
// }


arma::vec norminv(const arma::vec & w) {
  
  arma::vec quantiles(w.n_elem);
  for (int i = 0; i < quantiles.n_elem; i++) {
    quantiles(i) = R::qnorm(w(i), 0.0, 1.0, true, false);
  }
  return quantiles;
}

//' Transformed integrand for computation of MVN probabilities
// [[Rcpp::export]]
arma::vec mvndnv(int n, arma::vec a, arma::mat ch, arma::vec b, arma::uvec clg,
                 double ci, double dci, arma::mat x, int nv) {
  
  arma::mat y = arma::zeros(n-1, nv);
  arma::vec on = arma::ones(nv);
  // arma::vec c = ci * on;
  // arma::vec dc = dci * on;
  double c = ci;
  double dc = dci;
  arma::vec p = dci * on;
  int li = 1;
  int lf = 0;
  arma::vec s;
  double ai;
  double bi;
  
  for (int i = 1; i < n; i++) {
    y.row(i-1) = norminv(c + x.row(i-1) * dc);
    lf += clg(i);
    
    if (lf < li) {
      c = 0.0;
      dc = 1.0;
    }
    else {
      s = ch.submat(li, 0, lf, (i-1)) * y.rows(0, i-1);
      ai = std::max( arma::max(a.rows(li, lf) * on - s), -36.);
      bi = std::max(ai, std::min( arma::min(b.rows(li, lf) * on - s), 9.) );
      c = arma::normcdf(ai);
      dc = arma::normcdf(bi) - c;
      p = p * dc;
    }
    
    li += clg(i);
  }
  
  return p;
}

//' r covariance matrix n x n
//' a lower constraint m x 1 
//' b upper constraint m x 1
//' cn polytope constraint m x n 
// [[Rcpp::export]]
Rcpp::List chlsrt(const arma::mat & r, const arma::vec & a, const arma::mat & cn,
                  const arma::vec & b) {
  
  int n = r.n_rows;
  arma::vec y = arma::zeros(n);
  arma::rowvec clg = arma::zeros(n);
  arma::vec ap = a;
  arma::vec bp = b;
  
  int m = cn.n_rows;
  arma::mat ch = cn;
  arma::mat c = r;
  arma::vec d = arma::sqrt(arma::max( arma::diagvec(c), arma::zeros(n) ));
  
  for (int i = 0; i < n; i++) {
    double di = d(i);
    if (di > 0) {
      c.col(i) /= di;
      c.row(i) /= di;
      ch.col(i) *= di;
    }
  }
  
  // determine r factors and form revised constraint matrix ch
  arma::vec D;
  arma::mat V;
  arma::eig_sym(D, V, c);
  arma::uvec pm = arma::sort_index(D, "descend");
  d = arma::sqrt( arma::max( D(pm), arma::zeros(pm.n_elem) ) );
  
  arma::uvec pos_idx = arma::find(d > 0.0);
  int np = pos_idx.n_elem;
  c = V.cols(pm.rows(0, np)) % ( arma::ones(n) * d.rows(0, np).t() );
  ch *= c;
  
  // use right reflectors to reduce ch to lower triangular
  double mn;
  double vr;
  arma::rowvec v;
  double ss;
  for (arma::uword i = 0; i < std::min(np-1, m); i++) {
    
    double epi = ep * i;
    double vm = 1;
    arma::uword lm = i;
    
    for (int l = i; i < m; l++) {
      
      v = ch.row(l).cols(0, np);
      double s = arma::as_scalar(v.cols(0, i-1) * y.rows(0, i-1));
      ss = std::max( std::sqrt( arma::sum( arma::square(v.cols(i, np)) ) ), epi);
      double al = (ap(l) - s) / ss;
      double bl = (bp(l) - s) / ss;
      double dna = 0.0;
      double dsa = 0.0;
      double dnb = 0.0;
      double dsb = 1.0;
      
      if (al > -9.) {
        dna = arma::normpdf(al);
        dsa = arma::normcdf(al);
      }
      if (bl < 9.) {
        dnb = arma::normpdf(bl);
        dsb = arma::normcdf(bl);
      }
      if ((dsb - dsa) > epi) {
        if (al <= -9) {
          mn = -dnb;
          vr = -bl * dnb;
        } 
        else if (bl >= 9) {
          mn = dna;
          vr = al * dna;
        }
        else {
          mn = dna - dnb;
          vr = al * dna - bl * dnb;
        }
        
        mn /= (dsb - dsa);
        vr = 1.0 + vr / (dsb - dsa) - std::pow(mn, 2.0);
      } 
      else {
        mn = (al + bl) / 2.;
        vr = 0.;
        if (al <= -9) {
          mn = bl;
        }
        else if (bl >= 9) {
          mn = al;
        }
      }
      
      if (vr <= vm) {
        lm = l;
        vm = vr;
        y(i) = mn;
      }
    }
    
    v = ch.row(lm).cols(0, np);
    
    if (lm > i) {
      arma::uvec ilm = {i, lm};
      arma::uvec lmi = {lm, i};
      arma::uvec col_select = arma::regspace<arma::uvec>(0, 1, np);
      ch.submat(ilm, col_select) = ch.submat(lmi, col_select);
      ap(ilm) = ap(lmi);
      bp(ilm) = bp(lmi);
    }
    
    ch.row(i).cols(i+1, np) = 0.0;
    ss = arma::sum( arma::square(v.rows(i+1, np)) );
    
  }
}