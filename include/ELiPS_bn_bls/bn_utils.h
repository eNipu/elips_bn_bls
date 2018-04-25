//
//  bn_utils.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/1/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_utils_h
#define bn_utils_h

#include <stdio.h>
#include <sys/time.h>

extern struct timeval t0,t1;

extern double MILLER_TATE,MILLER_PLAINATE,MILLER_OPTATE,MILLER_XATE;
extern double FINALEXP_PLAIN,FINALEXP_OPT;
extern double G1SCM_PLAIN,G1SCM_2SPLIT;
extern double G2SCM_PLAIN,G2SCM_2SPLIT,G2SCM_4SPLIT;
extern double G3SCM_PLAIN,G3SCM_2SPLIT,G3SCM_4SPLIT;

extern double timedifference_msec(struct timeval t0, struct timeval t1);
extern double timedifference_usec(struct timeval t0, struct timeval t1);

//print pairing
extern void BN12_print_parameters(void);
extern void BN12_print_G1_point(void);
extern void BN12_print_G2_point(void);
extern void BN12_print_tate_time(void);
extern void BN12_print_plain_ate_time(void);
extern void BN12_print_opt_ate_time(void);
extern void BN12_print_x_ate_time(void);
extern void BN12_print_final_exp_plain_time(void);
extern void BN12_print_final_exp_optimal_time(void);
extern void BN12_print_G1_plain_time(void);
extern void BN12_print_G1_2split_time(void);
extern void BN12_print_G2_plain_time(void);
extern void BN12_print_G2_2split_time(void);
extern void BN12_print_G2_4split_time(void);
extern void BN12_print_G3_plain_time(void);
extern void BN12_print_G3_2split_time(void);
extern void BN12_print_G3_4split_time(void);

#endif /* bn_utils_h */
