//barrier_quant_revised.cpp
// ------------------------------------------------------------------
// 🌟 Fast Binomial Summation Formulas for computing some barrier option prices created by Dr. Burak YILDIZ, Ankara/Türkiye
//     Email: tr.burak.yildiz@gmail.com
//     Google Scholar: https://scholar.google.com/citations?user=DFV7KBQAAAAJ&hl=tr
//     Example How to compile:
//     g++ barrier_quant_revised.cpp -o barrier_quant_revised -std=c++17 -I/opt/homebrew/include -L/opt/homebrew/lib -lQuantLib -lgmp -lm
//     Example run:
//     ./barrier_quant_revised 10000
// ------------------------------------------------------------------
#define _USE_MATH_DEFINES

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <cmath>

// ---------------- QuantLib includes ----------------
#include <ql/quantlib.hpp>
using namespace QuantLib;

using namespace std;

// ---------------- Global variables ----------------

int H = 0;


void downoutputprice(int n,const double& Barrier, const double& S, const double& K, const double& r, const double& v, const double& T){
  int j;
  double u=exp(v*sqrt(T/n));
  double d=exp(-v*sqrt(T/n));
  double p=(exp(r*T/n)-d)/(u-d);
  double mydiscount=exp(-r*T);
    
  mpf_t priceput,pj,qnj,uj,dnj,SS,KK,nchoosek,UU,VV,WW,MM;
    mpz_t binomial;
    mpz_t binomial2;
    mpz_init(binomial);
    mpz_init(binomial2);
    mpf_init(nchoosek);
    mpf_init(priceput);
    mpf_init(pj);
    mpf_init(qnj);
    mpf_init(UU);
    mpf_init(VV);
    mpf_init(WW);
    mpf_init(MM);
    mpf_init(uj);
    mpf_init(dnj);
    mpf_init(SS);
    mpf_init(KK);
    mpf_set_d(priceput,0.0);
    mpf_set_d(SS,S);
    mpf_set_d(KK,K);

     H=-(floor(log(Barrier/S)/log(d))+1);
  for (j=0; j <= n ; j++){
      mpz_bin_uiui(binomial,n,j);
      if(j-H<=n){
          mpz_bin_uiui(binomial2,n,j-H);
          if(mpz_cmp(binomial,binomial2)>0){
              mpz_sub(binomial,binomial,binomial2);
          }else{
              continue;
          }
      }
      mpf_set_z(nchoosek,binomial);
      mpf_set_d(pj,p);
      mpf_set_d(qnj,1-p);
      mpf_set_d(uj,u);
      mpf_set_d(dnj,d);
      mpf_pow_ui(pj,pj,j);
      mpf_pow_ui(uj,uj,j);
      mpf_pow_ui(qnj,qnj,n-j);
      mpf_pow_ui(dnj,dnj,n-j);
      mpf_mul(UU,nchoosek,pj);
      mpf_mul(UU,UU,qnj);
      mpf_mul(VV,SS,uj);
      mpf_mul(VV,VV,dnj);
      if(mpf_cmp_d(VV,K)<0){
          mpf_sub(WW,KK,VV);
          mpf_mul(MM,UU,WW);
          mpf_add(priceput,priceput,MM);
      }
  }
  printf ("%f", mpf_get_d(priceput)*mydiscount);
    
}
void downoutcallprice(int n,const double& Barrier, const double& S, const double& K, const double& r, const double& v, const double& T){
  int j;
  double u=exp(v*sqrt(T/n));
  double d=exp(-v*sqrt(T/n));
  double p=(exp(r*T/n)-d)/(u-d);
  double mydiscount=exp(-r*T);
    
  mpf_t pricecall,pj,qnj,uj,dnj,SS,KK,nchoosek,UU,VV,WW,MM;
    mpz_t binomial;
    mpz_t binomial2;
    mpz_init(binomial);
    mpz_init(binomial2);
    mpf_init(nchoosek);
    mpf_init(pricecall);
    mpf_init(pj);
    mpf_init(qnj);
    mpf_init(UU);
    mpf_init(VV);
    mpf_init(WW);
    mpf_init(MM);
    mpf_init(uj);
    mpf_init(dnj);
    mpf_init(SS);
    mpf_init(KK);
    mpf_set_d(pricecall,0.0);
    mpf_set_d(SS,S);
    mpf_set_d(KK,K);
     H=-(floor(log(Barrier/S)/log(d))+1);

  for (j=0; j <= n ; j++){
      mpz_bin_uiui(binomial,n,j);
      if(j-H<=n){
          mpz_bin_uiui(binomial2,n,j-H);
          if(mpz_cmp(binomial,binomial2)>0){
              mpz_sub(binomial,binomial,binomial2);
          }else{
              continue;
          }
      }
      mpf_set_z(nchoosek,binomial);
      mpf_set_d(pj,p);
      mpf_set_d(qnj,1-p);
      mpf_set_d(uj,u);
      mpf_set_d(dnj,d);
      mpf_pow_ui(pj,pj,j);
      mpf_pow_ui(uj,uj,j);
      mpf_pow_ui(qnj,qnj,n-j);
      mpf_pow_ui(dnj,dnj,n-j);
      mpf_mul(UU,nchoosek,pj);
      mpf_mul(UU,UU,qnj);
      mpf_mul(VV,SS,uj);
      mpf_mul(VV,VV,dnj);
      if(mpf_cmp_d(VV,K)>0){
          mpf_sub(WW,VV,KK);
          mpf_mul(MM,UU,WW);
          mpf_add(pricecall,pricecall,MM);
      }
  }
  printf ("%f", mpf_get_d(pricecall)*mydiscount);
    
}

void upoutcallprice(int n,const double& Barrier, const double& S, const double& K, const double& r, const double& v, const double& T){
  int j;
  double u=exp(v*sqrt(T/n));
  double d=exp(-v*sqrt(T/n));
  double p=(exp(r*T/n)-d)/(u-d);
  double mydiscount=exp(-r*T);
    
  mpf_t pricecall,pj,qnj,uj,dnj,SS,KK,nchoosek,UU,VV,WW,MM;
    mpz_t binomial;
    mpz_t binomial2;
    mpz_init(binomial);
    mpz_init(binomial2);
    mpf_init(nchoosek);
    mpf_init(pricecall);
    mpf_init(pj);
    mpf_init(qnj);
    mpf_init(UU);
    mpf_init(VV);
    mpf_init(WW);
    mpf_init(MM);
    mpf_init(uj);
    mpf_init(dnj);
    mpf_init(SS);
    mpf_init(KK);
    mpf_set_d(pricecall,0.0);
    mpf_set_d(SS,S);
    mpf_set_d(KK,K);
     H=floor(log(Barrier/S)/log(u))+1;
  for (j=0; j<= n ; j++){
      mpz_bin_uiui(binomial,n,j);
      if(H<=j){
          mpz_bin_uiui(binomial2,n,j-H);
          if(mpz_cmp(binomial,binomial2)>0){
              mpz_sub(binomial,binomial,binomial2);
          }else{
              continue;
          }
      }
      mpf_set_z(nchoosek,binomial);
      mpf_set_d(pj,p);
      mpf_set_d(qnj,1-p);
      mpf_set_d(uj,u);
      mpf_set_d(dnj,d);
      mpf_pow_ui(pj,pj,j);
      mpf_pow_ui(uj,uj,j);
      mpf_pow_ui(qnj,qnj,n-j);
      mpf_pow_ui(dnj,dnj,n-j);
      mpf_mul(UU,nchoosek,pj);
      mpf_mul(UU,UU,qnj);
      mpf_mul(VV,SS,uj);
      mpf_mul(VV,VV,dnj);
      if(mpf_cmp_d(VV,K)>0){
          mpf_sub(WW,VV,KK);
          mpf_mul(MM,UU,WW);
          mpf_add(pricecall,pricecall,MM);
      }
  }
  printf ("%f", mpf_get_d(pricecall)*mydiscount);
    
}
void upoutputprice(int n,const double& Barrier, const double& S, const double& K, const double& r, const double& v, const double& T){
  int j;
  double u=exp(v*sqrt(T/n));
  double d=exp(-v*sqrt(T/n));
  double p=(exp(r*T/n)-d)/(u-d);
  double mydiscount=exp(-r*T);
    
  mpf_t priceput,pj,qnj,uj,dnj,SS,KK,nchoosek,UU,VV,WW,MM;
    mpz_t binomial;
    mpz_t binomial2;
    mpz_init(binomial);
    mpz_init(binomial2);
    mpf_init(nchoosek);
    mpf_init(priceput);
    mpf_init(pj);
    mpf_init(qnj);
    mpf_init(UU);
    mpf_init(VV);
    mpf_init(WW);
    mpf_init(MM);
    mpf_init(uj);
    mpf_init(dnj);
    mpf_init(SS);
    mpf_init(KK);
    mpf_set_d(priceput,0.0);
    mpf_set_d(SS,S);
    mpf_set_d(KK,K);

 H=floor(log(Barrier/S)/log(u))+1;
  for (j=0; j<= n ; j++){
      mpz_bin_uiui(binomial,n,j);
      if(H<=j){
          mpz_bin_uiui(binomial2,n,j-H);
          if(mpz_cmp(binomial,binomial2)>0){
              mpz_sub(binomial,binomial,binomial2);
          }else{
              continue;
          }
      }
      mpf_set_z(nchoosek,binomial);
      mpf_set_d(pj,p);
      mpf_set_d(qnj,1-p);
      mpf_set_d(uj,u);
      mpf_set_d(dnj,d);
      mpf_pow_ui(pj,pj,j);
      mpf_pow_ui(uj,uj,j);
      mpf_pow_ui(qnj,qnj,n-j);
      mpf_pow_ui(dnj,dnj,n-j);
      mpf_mul(UU,nchoosek,pj);
      mpf_mul(UU,UU,qnj);
      mpf_mul(VV,SS,uj);
      mpf_mul(VV,VV,dnj);
      if(mpf_cmp_d(VV,K)<0){
          mpf_sub(WW,KK,VV);
          mpf_mul(MM,UU,WW);
          mpf_add(priceput,priceput,MM);
      }
  }
  printf ("%f", mpf_get_d(priceput)*mydiscount);
    
}

// ---------------- Main ----------------
int main(int argc, char * argv[]){
    int n;
    double S = 130.0;     // Spot
    double K = 140.0;    // Strike
    double Bdown = 100.0; // Barrier down
    double Bup = 170.0;  // Barrier up
    double r = 0.04;     // Risk-free
    double sigma = 0.45; // Volatility
    double T = 0.1;      // Maturity
    double q = 0.0;      // Dividend

  if (argc <= 1){
    printf ("Usage: %s <number> \n", argv[0]);
    return 2;
  }

  n = atoi(argv[1]);
  assert( n >= 0);

    // ---------------- Analytic prices using QuantLib ----------------
Calendar calendar = TARGET();
Date today = Date::todaysDate();
Settings::instance().evaluationDate() = today;

Handle<Quote> spot(boost::make_shared<SimpleQuote>(S));
Handle<YieldTermStructure> rTS(boost::make_shared<FlatForward>(today, r, Actual365Fixed()));
Handle<YieldTermStructure> qTS(boost::make_shared<FlatForward>(today, q, Actual365Fixed()));
Handle<BlackVolTermStructure> volTS(boost::make_shared<BlackConstantVol>(today, calendar, sigma, Actual365Fixed()));

boost::shared_ptr<GeneralizedBlackScholesProcess> process(
    boost::make_shared<GeneralizedBlackScholesProcess>(spot, qTS, rTS, volTS)
);

// Barrier Options
boost::shared_ptr<StrikedTypePayoff> payoff_call(boost::make_shared<PlainVanillaPayoff>(Option::Call,K));
boost::shared_ptr<StrikedTypePayoff> payoff_put(boost::make_shared<PlainVanillaPayoff>(Option::Put,K));
boost::shared_ptr<Exercise> exercise(boost::make_shared<EuropeanExercise>(today + Period(int(T*365), Days)));

BarrierOption dl_call(Barrier::DownOut, Bdown, 0.0, payoff_call, exercise);
dl_call.setPricingEngine(boost::make_shared<AnalyticBarrierEngine>(process));
double analytic_down_out_call = dl_call.NPV();

BarrierOption dl_put(Barrier::DownOut, Bdown, 0.0, payoff_put, exercise);
dl_put.setPricingEngine(boost::make_shared<AnalyticBarrierEngine>(process));
double analytic_down_out_put = dl_put.NPV();

BarrierOption ul_call(Barrier::UpOut, Bup, 0.0, payoff_call, exercise);
ul_call.setPricingEngine(boost::make_shared<AnalyticBarrierEngine>(process));
double analytic_up_out_call = ul_call.NPV();

BarrierOption ul_put(Barrier::UpOut, Bup, 0.0, payoff_put, exercise);
ul_put.setPricingEngine(boost::make_shared<AnalyticBarrierEngine>(process));
double analytic_up_out_put = ul_put.NPV();

    // ---------------- Print parameters and prices ----------------
    std::cout << "-----------Parameters------------------------------" << std::endl;
    std::cout << "Underlying:      " << S << std::endl;
    std::cout << "Strike:          " << K << std::endl;
    std::cout << "Barrier Up:      " << Bup << std::endl;
    std::cout << "Barrier down:    " << Bdown << std::endl;
    std::cout << "Risk-Free Rate:  " << r << std::endl;
    std::cout << "Volatility:      " << sigma << std::endl;
    std::cout << "Maturity:        " << T << std::endl;

    std::cout << "------------Analytic Solutions Begin --------------" << std::endl;
    std::cout << "QuantLib Down-Out Call: " << analytic_down_out_call << std::endl;
    std::cout << "QuantLib Down-Out Put:  " << analytic_down_out_put << std::endl;
    std::cout << "QuantLib Up-Out Call:   " << analytic_up_out_call << std::endl;
    std::cout << "QuantLib Up-Out Put:    " << analytic_up_out_put << std::endl;
    std::cout << "------------Binomial Scheme Begin -----------------" << std::endl;
std::cout << "------ Cupout ----- Pupout --- Cdownout --- Pdownout" << std::endl;
    int i=1000;
    while(i<=n){
        std::cout <<"n="<<i<<" [";
        upoutcallprice(i, Bup, S, K, r, sigma, T);
        std::cout <<" , ";
        upoutputprice(i, Bup, S, K, r, sigma, T);
        std::cout <<" , ";
        downoutcallprice(i, Bdown, S, K, r, sigma, T);
        std::cout <<" , ";
        downoutputprice(i, Bdown, S, K, r, sigma, T);
        std::cout <<"]"<< std::endl;
        i=i+1000;
    }
    std::cout << "------------Binomial Scheme End ---------------------" << std::endl;

    return 0;
}
