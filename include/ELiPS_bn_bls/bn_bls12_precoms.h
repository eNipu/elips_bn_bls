//
//  bn_bls12_precoms.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 1/31/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_precoms_h
#define bn_precoms_h

#include <ELiPS_bn_bls/bn_fp6.h>
#include <ELiPS_bn_bls/curve_settings.h>

#define d12 12
#define d24 24
//#define BN_X_length 114

enum state{
    f_p1,f_p2,f_p3,f_p4,f_p5,f_p6,f_p7,f_p8,f_p9,f_p10,f_p11,f_p12
};

extern struct Fp Fp_basis;
extern struct Fp2 Fp2_basis;
extern struct Fp2 Fp2_basis_inv;
extern struct Fp6 Fp6_basis;

extern mpz_t epsilon1,epsilon2;

extern  Fp2 d12_frobenius_constant[d12][6];
extern  Fp2 d12_skew_frobenius_constant[d12][2];

extern void init_precoms(int curvetype);
extern void get_epsilon(void);
extern void set_basis(void);
extern void set_frobenius_constant(void);
#endif /* bn_precoms_h */
