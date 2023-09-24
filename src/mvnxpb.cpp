#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include <pbv.h>

std::vector<int> pmone = {-1, 1};

#include <cmath>

const double eps = 1e-10;

//' compute P(x > dh, y > dk) for x,y ~ bivariate normal with correlation r
// [[Rcpp::export]]
double bvnu(double dh, double dk, double r) {
  
  double p;  // the bivariate normal probability
  
  int ng, lg;  // Gauss Legendre point and weight matrix dimensions
  
  arma::mat w; // Gauss Legendre weights
  arma::mat x; // Gauss Legendre points
  
  double h, k, hk, bvn, hs, asr, sn, sp, as, a, b, bs, c, d, xs, rs, ep, L;
  
  if (dh == INFINITY || dk == INFINITY) {
    p = 0.0;
  }
  else if (dh == - INFINITY) {
    if (dk == - INFINITY) {
      p = 1.0;
    }
    else {
      p = arma::normcdf(-dk);
    }
  }
  else if (dk == -INFINITY) {
    p = arma::normcdf(-dh);
  }
  else {
    if (std::abs(r) < 0.3) {
      ng = 1;
      lg = 3;
      w = arma::mat(lg, ng);
      x = arma::mat(lg, ng);
      w.rows(0, lg-1) = arma::vec({0.1713244923791705, 0.3607615730481384, 0.4679139345726904});
      x.rows(0, lg-1) = arma::vec({0.9324695142031522, 0.6612093864662647, 0.2386191860831970});
    }
    else if (std::abs(r) < 0.75) {
      ng = 2;
      lg = 6;
      w = arma::mat(lg, ng);
      x = arma::mat(lg, ng);
      w.submat(0, 1, 2, 1) = arma::vec({.04717533638651177, 0.1069393259953183, 0.1600783285433464});
      w.submat(3, 1, 5, 1) = arma::vec({0.2031674267230659, 0.2334925365383547, 0.2491470458134029});
      x.submat(0, 1, 2, 1) = arma::vec({0.9815606342467191, 0.9041172563704750, 0.7699026741943050});
      x.submat(3, 1, 5, 1) = arma::vec({0.5873179542866171, 0.3678314989981802, 0.1252334085114692});
    }
    else {
      ng = 3;
      lg = 10;
      w = arma::mat(lg, ng);
      x = arma::mat(lg, ng);
      w.submat(0, 2, 2, 2) = arma::vec({.01761400713915212, .04060142980038694, .06267204833410906});
      w.submat(3, 2, 5, 2) = arma::vec({.08327674157670475, 0.1019301198172404, 0.1181945319615184});
      w.submat(6, 2, 8, 2) = arma::vec({0.1316886384491766, 0.1420961093183821, 0.1491729864726037});
      w(9, 2) = 0.1527533871307259;
      x.submat(0, 2, 2, 2) = arma::vec({0.9931285991850949, 0.9639719272779138, 0.9122344282513259});
      x.submat(3, 2, 5, 2) = arma::vec({0.8391169718222188, 0.7463319064601508, 0.6360536807265150});
      x.submat(6, 2, 8, 2) = arma::vec({0.5108670019508271, 0.3737060887154196, 0.2277858511416451});
      x(9, 2) = 0.07652652113349733;
    }
    
    h = dh;
    k = dk;
    hk = h * k;
    bvn = 0.0;
    
    if (std::abs(r) < 0.925) {
      hs = (h * h + k * k) / 2;
      asr = std::asin(r);
      
      for (int i = 0; i < lg; i++) {
        sn = std::sin( asr * (1 - x(i, ng-1) ) / 2 );
        bvn += w(i, ng-1) * std::exp( (sn*hk - hs) / (1 - sn*sn) );
        sn = std::sin( asr * (1 - x(i, ng-1) ) / 2 );
        bvn += w(i, ng-1) * std::exp( (sn*hk - hs) / (1 - sn*sn) );
      }
      
      bvn *= asr / (4 * arma::datum::pi);
      bvn += arma::normcdf(-h) * arma::normcdf(-k);
    }
    else {
      if (r < 0.0) {
        k = -k;
        hk = -k;
      }
      
      if (std::abs(r) < 1.0) {
        as = (1-r) * (1+r);
        a = std::sqrt(as);
        bs = std::pow(h-k, 2.0);
        c = (4.0 - hk) / 8.0;
        d = (12.0 - hk) / 16.0;
        asr = -((bs / as) + hk) / 2.;
        
        if (asr > -100.) {
          bvn = a*std::exp(asr) * (1 - c*(bs - as) * (1 - d*bs/5.)/3. + c*d*as*as/5.);
        }
        if (hk > -100.) {
          b = std::sqrt(bs);
          sp = std::sqrt(2.*arma::datum::pi) * arma::normcdf(-b/a);
          bvn -= std::exp(-hk/2.) * sp * b * (1 - c*bs*(1 - d*bs/5.) / 3.);
        }
        
        a /= 2.;
        for (int i = 0; i < lg; i++) {
          for (int is : pmone) {
            xs = std::pow(a + a*is*x(i, ng-1), 2.);
            rs = std::sqrt(1 - xs);
            asr = -(bs/xs + hk) / 2;
            if (asr > -100) {
              sp = (1 + c*xs*(1 + d*xs));
              ep = std::exp( -hk*xs / (2*std::pow(1 + rs, 2.)) ) / rs;
              bvn += a*w(i, ng-1)*std::exp(asr)*( ep - sp);
            }
          }
        }
        
        bvn /= -2 * arma::datum::pi;
      }
      
      if (r > 0.) {
        bvn += arma::normcdf(-std::max(h, k));
      }
      else if (h >= k) {
        bvn *= -1;
      }
      else {
        if (h < 0.) {
          L = arma::normcdf(k) - arma::normcdf(h);
        }
        else {
          L = arma::normcdf(-h) - arma::normcdf(-k);
        }
        bvn = L - bvn;
      }
    }
    
    p = std::max(0.0, std::min(1.0, bvn));
  }
  return p;
}

std::tuple<double, double, double> bvnmom(double xl, double xu, 
                                          double yl, double yu,
                                          double r) {
  
  double rs = 1 / std::sqrt(1 - std::pow(r, 2.));
  double p = bvnu(xl, yl, r) - bvnu(xu, yl, r) - bvnu(xl, yu, r) + bvnu(xu, yu, r);
  // double p =
  //   pbv::pbv_rcpp_pbvnorm0(xu, yu, r) - pbv::pbv_rcpp_pbvnorm0(xu, yl, r) -
  //   pbv::pbv_rcpp_pbvnorm0(xl, yu, r) + pbv::pbv_rcpp_pbvnorm0(xl, yl, r);
  
  double Ex, Ey;
  
  double Phiyxu, Phixyu, Phixyl, Phiyxl;
  double pPhixy, pPhiyx;
  
  if (xl == -INFINITY && yl == -INFINITY && xu == INFINITY && yu == INFINITY) {
    Ey = 0.0, Ex = 0.0;
  }
  else if (xl == -INFINITY && xu == INFINITY && yl == -INFINITY) {
    Ey = -arma::normpdf(yu);
    Ex = 0.0;
  }
  else if (xl == -INFINITY && xu == INFINITY && yu == INFINITY) {
    Ey = arma::normpdf(yl);
    Ex = 0.0;
  }
  else if (xl == -INFINITY && yl == -INFINITY && yu == INFINITY) {
    Ex = -arma::normpdf(xu);
    Ey = 0.0;
  }
  else if (xu == INFINITY && yl == -INFINITY && yu == INFINITY) {
    Ex = arma::normpdf(xl);
    Ey = 0.0;
  }
  else if (xl == -INFINITY && xu == INFINITY) {
    Ey = arma::normpdf(yl) - arma::normpdf(yu);
    Ex = 0.0;
  }
  else if (yl == -INFINITY && yu == INFINITY) {
    Ex = arma::normpdf(xl) - arma::normpdf(xu);
    Ey = 0.0;
  }
  else {
    if (xl == -INFINITY && yl == -INFINITY) {
      Phiyxu = arma::normcdf( (yu - r*xu)*rs );
      pPhixy = -arma::normpdf(xu) * Phiyxu;
      Phixyu = arma::normcdf( (xu - r*yu)*rs );
      pPhiyx = -arma::normpdf(yu) * Phixyu;
    }
    else if (xu == INFINITY && yu == INFINITY) {
      Phiyxl = arma::normcdf(-(yl - r*xl)*rs);
      pPhixy = arma::normpdf(xl) * Phiyxl;
      Phixyl = arma::normcdf(-(xl - r*yl)*rs);
      pPhiyx = arma::normpdf(yl) * Phixyl;
    }
    else if (xl == -INFINITY && yu == INFINITY) {
      Phiyxu = arma::normcdf(-(yl - r*xu)*rs);
      pPhixy = -arma::normpdf(xu) * Phiyxu;
      Phixyl = arma::normcdf( (xu - r*yl)*rs);
      pPhiyx = arma::normpdf(yl) * Phixyl;
    }
    else if (xu == INFINITY && yl == -INFINITY) {
      Phiyxl = arma::normcdf( (yu - r*xl)*rs);
      pPhixy = arma::normpdf(xl) * Phiyxl;
      Phixyu = arma::normcdf(-(xl - r*yu)*rs);
      pPhiyx = -arma::normpdf(yu) * Phixyu;
    }
    else if (xl == -INFINITY) {
      Phiyxu = arma::normcdf((yu - r*xu)*rs) - arma::normcdf((yl - r*xu)*rs);
      pPhixy = -arma::normpdf(xu) * Phiyxu;
      Phixyl = arma::normcdf((xu - r*yl)*rs);
      Phixyu = arma::normcdf((xu - r*yu)*rs);
      pPhiyx = arma::normpdf(yl) * Phixyl - arma::normpdf(yu) * Phixyu;
    }
    else if (xu == INFINITY) {
      Phiyxl = arma::normcdf((yu - r*xl)*rs) - arma::normcdf((yl - r*xl)*rs);
      pPhixy = arma::normpdf(xl) * Phiyxl;
      Phixyl = arma::normcdf(-(xl - r*yl)*rs);
      Phixyu = arma::normcdf(-(xl - r*yu)*rs);
      pPhiyx = arma::normpdf(yl) * Phixyl - arma::normpdf(yu) * Phixyu;
    }
    else if (yl == -INFINITY) {
      Phiyxl = arma::normcdf((yu - r*xl)*rs);
      Phiyxu = arma::normcdf((yu - r*xu)*rs);
      pPhixy = arma::normpdf(xl) * Phiyxl - arma::normpdf(xu) * Phiyxu;
      Phixyu = arma::normcdf((xu - r*yu)*rs) - arma::normcdf((xl - r*yu)*rs);
      pPhiyx = -arma::normpdf(yu) * Phixyu;
    }
    else if (yu == INFINITY) {
      Phiyxl = arma::normcdf(-(yl - r*xl)*rs);
      Phiyxu = arma::normcdf(-(yl - r*xu)*rs);
      pPhixy = arma::normpdf(xl) * Phiyxl - arma::normpdf(xu) * Phiyxu;
      Phixyl = arma::normcdf((xu - r*yl)*rs) - arma::normcdf((xl - r*yl)*rs);
      pPhiyx = arma::normpdf(yl) * Phixyl;
    }
    else {
      Phiyxl = arma::normcdf((yu - r*xl)*rs) - arma::normcdf((yl - r*xl)*rs);
      Phiyxu = arma::normcdf((yu - r*xu)*rs) - arma::normcdf((yl - r*xu)*rs);
      pPhixy = arma::normpdf(xl) * Phiyxl - arma::normpdf(xu) * Phiyxu;
      Phixyl = arma::normcdf((xu - r*yl)*rs) - arma::normcdf((xl - r*yl)*rs);
      Phixyu = arma::normcdf((xu - r*yu)*rs) - arma::normcdf((xl - r*yu)*rs);
      pPhiyx = arma::normpdf(yl) * Phixyl - arma::normpdf(yu) * Phixyu;
    }
    
    Ex = pPhixy + r*pPhiyx;
    Ey = r*pPhixy + pPhiyx;
  }
  
  Ex /= p;
  Ey /= p;
  
  std::tuple<double, double, double> result = std::make_tuple(Ex, Ey, p);
  return result;
}

std::tuple<double, double, double> bvnmmg(arma::vec a, arma::vec b, arma::mat sg) {
  double cx = std::sqrt(sg(0, 0));
  double cy = std::sqrt(sg(1, 1));
  double r = sg(1, 0) / (cx * cy);
  double xl = a(0) / cx;
  double xu = b(0) / cx;
  double yl = a(1) / cy;
  double yu = b(1) / cy;

  std::tuple<double, double, double> mom_result = bvnmom(xl, xu, yl, yu, r);
  double Ex = std::get<0>(mom_result);
  double Ey = std::get<1>(mom_result);
  double p = std::get<2>(mom_result);
  
  return std::make_tuple(Ex * cx, Ey * cy, p);
}

// [[Rcpp::export]]
double mvnxpb(arma::mat r, arma::vec a, arma::vec b) {
  
  int n = r.n_rows;
  arma::mat c = r;
  arma::vec ap = a;
  arma::vec bp = b;
  arma::vec d = arma::sqrt(arma::max(arma::diagvec(c), arma::zeros(c.n_rows)));
  
  // scale c to make it a correlation matrix
  for (int i = 0; i < n; i++) {
    if (d(i) > 0) {
      ap(i) /= d(i);
      bp(i) /= d(i);
      c.col(i) /= d(i);
      c.row(i) /= d(i);
    }
  }
  
  // dynamically determine Cholesky factor
  // permute variables to minimize outer integrals
  arma::vec y = arma::zeros(n);
  double p = 1.0;
  double pb = 1.0;
  arma::mat D = arma::eye(n, n);
  
  double dem;
  double cii, ai, bi, de, ckk, am, bm;
  arma::uvec km, ip, ki;
  
  // Rcpp::Rcout << "c: " << c << std::endl;
  // Rcpp::Rcout << "y: " << y << std::endl;
  // Rcpp::Rcout << "ap: " << ap << std::endl;
  // Rcpp::Rcout << "bp: " << bp << std::endl;
  
  for (arma::uword k = 0; k < n; k++) {
    
    Rcpp::Rcout << "k: " << k << std::endl;
    
    arma::uword im = k;
    double ckk = 0.0;
    dem = 1.0;
    double s = 0.0;
    
    for (int i = k; i < n; i++) {
      
      Rcpp::Rcout << "i: " << i << std::endl;
      
      if (c(i, i) > eps) {
        cii = std::sqrt( std::max(c(i, i), 0.0) );
        // Rcpp::Rcout << "cii: " << cii << std::endl;
        if (i > 0) {
          if (k == 0) {
            s = 0.0;
          }
          else {
            s = arma::as_scalar(c.submat(i, 0, i, k-1) * y.subvec(0, k-1));
          }
        }
        ai = ( ap(i) - s ) / cii;
        bi = ( bp(i) - s ) / cii;
        de = arma::normcdf(bi) - arma::normcdf(ai);
        if (de <= dem) {
          ckk = cii;
          dem = de;
          am = ai;
          bm = bi;
          im = i;
        }
      }
    }
    
    // Rcpp::Rcout << "im: " << im << std::endl;
    // Rcpp::Rcout << "dem: " << dem << std::endl;
    
    arma::uvec kvec = arma::uvec({k});
    arma::uvec imvec = arma::uvec({im});
    if (im > k) {
      
      Rcpp::Rcout << "im > k" << std::endl;
      
      km = arma::regspace<arma::uvec>(0, 1, k-1);
      ip = arma::regspace<arma::uvec>(im+1, 1, n-1);
      ki = arma::regspace<arma::uvec>(k+1, 1, im-1);
      
      // Rcpp::Rcout << "km: " << km << std::endl;
      // Rcpp::Rcout << "ip: " << ip << std::endl;
      // Rcpp::Rcout << "ki: " << ki << std::endl;
       
      arma::uvec imk = arma::uvec({im, k});
      arma::uvec kim = arma::uvec({k, im});
      
      ap(imk) = ap(kim);
      bp(imk) = bp(kim);
      c.submat(imk, km) = c.submat(kim, km);
      c.submat(ip, imk) = c.submat(ip, kim);
      
      arma::vec t = c.submat(ki, kvec);
      c.submat(ki, kvec) = c.submat(imvec, ki).t();
      c.submat(imvec, ki) = t.t();
      c(im, im) = c(k, k);
    }
    
    // Rcpp::Rcout << "c after im > k: " << std::endl << c << std::endl;
    
    // Rcpp::Rcout << "should we enter ckk: " << ckk << ", eps*k: " << eps*(k+1) << std::endl;
    
    if (ckk > eps * (k+1)) {
      
      // Rcpp::Rcout << "------ ----- " << std::endl;
      // Rcpp::Rcout << "ckk > eps * k" << std::endl;
      
      c(k, k) = ckk;
      
      if (k+1 < n)
        c.submat(k, k+1, k, n-1).fill(0.0);
      
      // Rcpp::Rcout << "ckk: " << ckk << std::endl;
      
      for (int i = k+1; i < n; i++) {
        c(i, k) = c(i, k) / ckk;
        c.submat(i, k+1, i, i) = 
          c.submat(i, k+1, i, i) - c(i, k) * c.submat(k+1, k, i, k).t();
        // Rcpp::Rcout << "c: " << c << std::endl;
      }
      
      if (std::abs(dem) > eps) {
        y(k) = ( arma::normpdf(am) - arma::normpdf(bm) ) / dem;
      }
      else {
        y(k) = ( am + bm) / 2;
        if (am < -9.) {
          y(k) = bm;
        }
        else if (bm > 9.) {
          y(k) = am;
        }
      }
    }
    else {
      c.submat(k, k, n-1, k).fill(0.0);
      y(k) = ( ap(k) + bp(k) ) / 2;
    }
    
    // Rcpp::Rcout << "c after ckk > eps: " << std::endl << c << std::endl;
    
    p *= dem;
    // compute bivariate product
    if (k % 2 == 1) {
      
      Rcpp::Rcout << "doing a prob update" << std::endl;
      
      double u = c(k-1, k-1);
      double v = c(k, k);
      double w = c(k, k-1);
      
      arma::mat uvw_mat = arma::mat(2, 2);
      uvw_mat(0, 0) = 1 / u; uvw_mat(0, 1) = 0.0;
      uvw_mat(1, 0) = -w / (u * v);
      uvw_mat(1, 1) = 1 / v;
      c.submat(k-1, k-1, n-1, k) =  c.submat(k-1, k-1, n-1, k) * uvw_mat;
      
      arma::vec ab = ap.subvec(k-1, k);
      arma::vec bb = bp.subvec(k-1, k);
      arma::vec cb = arma::zeros(ab.n_elem);
      
      if (k > 1) {
        cb = c.submat(k-1, 0, k, k-2) * y.subvec(0, k-2);
      }
      
      arma::mat sg = arma::mat(2, 2);
      sg(0, 0) = u * u; sg(0, 1) = u * w;
      sg(1, 0) = u * w; sg(1, 1) = w*w + v*v;
      D.submat(k-1, k-1, k, k) = sg;
      
      Rcpp::Rcout << "ab: " << ab << std::endl;
      Rcpp::Rcout << "bb: " << bb << std::endl;
      Rcpp::Rcout << "cb: " << cb << std::endl;
      Rcpp::Rcout << "sg: " << std::endl << sg << std::endl;
      
      std::tuple<double, double, double> moments = bvnmmg(ab-cb, bb-cb, sg);
      double Ex = std::get<0>(moments);
      double Ey = std::get<1>(moments);
      double bv = std::get<2>(moments);
      pb *= bv;
      Rcpp::Rcout << "called interior functions" << std::endl;
      y.subvec(k-1, k) = arma::vec({Ex, Ey});
      
      Rcpp::Rcout << "Ex: " << Ex << std::endl;
      Rcpp::Rcout << "Ey: " << Ey << std::endl;
      Rcpp::Rcout << "bv: " << bv << std::endl;
    }
    
    // Rcpp::Rcout << "c after prob update: " << std::endl << c << std::endl;
  }
  
  if (n % 2 == 1) {
    pb *= dem;
  } 
  
  return(pb);
}
