//
//  bls12_timeprint.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/5/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bls12_timeprint_h
#define bls12_timeprint_h

#include <stdio.h>
//#include <sys/time.h>

extern double BLS12_MILLER_TATE,BLS12_MILLER_PLAINATE,BLS12_MILLER_OPTATE;
extern double BLS12_FINALEXP_PLAIN,BLS12_FINALEXP_OPT;
extern double BLS12_G1SCM_PLAIN,BLS12_G1SCM_2SPLIT;
extern double BLS12_G2SCM_PLAIN,BLS12_G2SCM_2SPLIT,BLS12_G2SCM_4SPLIT;
extern double BLS12_G3EXP_PLAIN,BLS12_G3EXP_2SPLIT,BLS12_G3EXP_4SPLIT;

//print pairing
extern void BLS12_print_parameters(void);
extern void BLS12_print_G1_point(void);
extern void BLS12_print_G2_point(void);
extern void BLS12_print_tate_time(void);
extern void BLS12_print_plain_ate_time(void);
extern void BLS12_print_opt_ate_time(void);
extern void BLS12_print_final_exp_plain_time(void);
extern void BLS12_print_final_exp_optimal_time(void);
extern void BLS12_print_plain_G1_scm_time(void);
extern void BLS12_print_2split_G1_scm_time(void);
extern void BLS12_print_plain_G2_scm_time(void);
extern void BLS12_print_2split_G2_scm_time(void);
extern void BLS12_print_4split_G2_scm_time(void);
extern void BLS12_print_plain_G3_exp_time(void);
extern void BLS12_print_2split_G3_exp_time(void);
extern void BLS12_print_4split_G3_exp_time(void);

#endif /* bls12_timeprint_h */
