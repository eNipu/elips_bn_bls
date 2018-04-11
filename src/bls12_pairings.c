//
//  bls12_pairings.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/6/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#include <ELiPS_bn_bls/bls12_pairings.h>

void BLS12_tate(Fp12 *ANS,EFp12 *P, EFp12 *Q){
    //miller
    gettimeofday(&t0,NULL);
    BLS12_Miller_algo_for_tate(ANS,P,Q);
    gettimeofday(&t1,NULL);
    BLS12_MILLER_TATE=timedifference_msec(t0,t1);
    
    //final exp
    gettimeofday(&t0,NULL);
    BLS12_Final_exp_optimal(ANS,ANS);
    gettimeofday(&t1,NULL);
    BLS12_FINALEXP_OPT=timedifference_msec(t0,t1);
}

void BLS12_plain_ate(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    //miller
    gettimeofday(&t0,NULL);
    BLS12_Miller_algo_for_plain_ate(ANS,P,Q);
    gettimeofday(&t1,NULL);
    BLS12_MILLER_PLAINATE=timedifference_msec(t0,t1);
    
    //final exp
    gettimeofday(&t0,NULL);
    BLS12_Final_exp_optimal(ANS,ANS);
    gettimeofday(&t1,NULL);
    BLS12_FINALEXP_OPT=timedifference_msec(t0,t1);
}

void BLS12_opt_ate(Fp12 *ANS,EFp12 *P, EFp12 *Q){
    //miller
    gettimeofday(&t0,NULL);
    BLS12_Miller_algo_for_opt_ate(ANS,P,Q);
    gettimeofday(&t1,NULL);
    BLS12_MILLER_OPTATE=timedifference_msec(t0,t1);
    
    //final exp
    gettimeofday(&t0,NULL);
    BLS12_Final_exp_optimal(ANS,ANS);
    gettimeofday(&t1,NULL);
    BLS12_FINALEXP_OPT=timedifference_msec(t0,t1);
}

