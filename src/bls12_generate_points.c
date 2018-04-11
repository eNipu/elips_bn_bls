//
//  bls12_generate_points.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/6/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#include <ELiPS_bn_bls/bls12_generate_points.h>
#include <ELiPS_bn_bls/curve_settings.h>

//rational points
void BLS12_generate_G1_point(EFp12 *P){
    EFp Tmp_P;
    EFp_init(&Tmp_P);
    mpz_t scalar;
    mpz_init(scalar);
    
    EFp_rational_point_BLS12(&Tmp_P);
    BLS12_EFp_to_EFp12(P,&Tmp_P);
    mpz_tdiv_q(scalar,curve_parameters.EFp_total,curve_parameters.order);
    EFp12_SCM(P,P,scalar);
    P->infinity=Tmp_P.infinity;
    
    EFp_clear(&Tmp_P);
    mpz_clear(scalar);
}

void BLS12_generate_G2_point(EFp12 *Q){
    EFp12 random_P,P,frobenius_P;
    EFp12_init(&random_P);
    EFp12_init(&P);
    EFp12_init(&frobenius_P);
    mpz_t exp;
    mpz_init(exp);
    
    EFp12_rational_point_BLS12(&random_P);
    mpz_pow_ui(exp,curve_parameters.order,2);
    mpz_tdiv_q(exp,curve_parameters.EFpd_total,exp);
    EFp12_SCM(&P,&random_P,exp);
    BLS12_Fp12_frobenius_map_p1(&frobenius_P.x,&P.x);
    BLS12_Fp12_frobenius_map_p1(&frobenius_P.y,&P.y);
    frobenius_P.infinity=P.infinity;
    EFp12_set_neg(&P,&P);
    EFp12_ECA(Q,&P,&frobenius_P);
    
    mpz_clear(exp);
    EFp12_clear(&random_P);
    EFp12_clear(&P);
    EFp12_clear(&frobenius_P);
}

void BLS12_generate_random_point(EFp12 *R){
    EFp12_init(R);
    EFp12_rational_point_BLS12(R);
}
