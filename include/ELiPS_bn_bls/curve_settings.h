//
//  curve_settings.h
//
//  Created by Khandaker Md. Al-Amin on 1/25/18.
//  Copyright © 2018 Khandaker Md. Al-Amin. All rights reserved.
//

#ifndef bn_settings_h
#define bn_settings_h

#include <ELiPS_bn_bls/Commont_headers.h>
//#include "bn_fp.h"
//#include "bn_fp2.h"

extern int isCleared; // 1 is cleared. 0 is not cleared

struct curve_params {
    mpz_t prime;
    mpz_t X;
    int sign;
    mpz_t trace_t;
    mpz_t order;
    mpz_t EFp_total;
    mpz_t EFp2_total;
    mpz_t EFp6_total;
    mpz_t EFp12_total;
    mpz_t EFpd_total;
    mpz_t curve_a;
    mpz_t curve_b;
};
/**
 * @brief BN curve's systematically obtained parameters.
 *
 * It's a global variable that give access to BN public parameters.
 */
extern struct curve_params curve_parameters;

#define BN_X_length 114

extern gmp_randstate_t state;

extern char X_binary[BN_X_length+1];
extern char X_binary_opt[BN_X_length+3];

//bls12
extern int bls12_X_length;
extern mpz_t bls12_X;
extern int bls12_X_binary[78];


//init
extern void init_bn_settings(void);
extern void init_bn_parameters( void);
extern void generate_bn_mother_parameter(void);
extern int  generate_bn_prime(void);
extern int  generate_bn_order(void);
extern void generate_bn_trace(void);
extern void set_bn_curve_parameter(void);
extern void weil(void);
extern void print_curve_parameters(void);

extern void init_bls12_settings(void);
extern void init_bls12_parameters( void);
extern void bls12_generate_X(void);
extern int  bls12_generate_prime(void);
extern int  bls12_generate_order(void);
extern void bls12_generate_trace(void);
extern void bls12_set_curve_parameter(void);
extern void bls12_weil(void);
extern void bls12_print_parameters(void);


//KSS-16
#define TRUE 1
#define FALSE 0
#define KSS16_X_length 35

extern char X_bit_binary_kss16[KSS16_X_length+1];

extern mpz_t C1_INV; //c=2 its iverse value
//mpz_t PRIME_P,order_r,trace_t, order_EFp, a_x;
//mpz_t tmp_a;

extern void init_kss16_settings(void);
extern void generate_kss16_motherparam(void);
extern void generate_kss16_parameters(void);
extern void init_kss16_parameters(void);
#endif /* bn_settings_h */
