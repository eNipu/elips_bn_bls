//
//  bls12_timeprint.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/5/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#include <ELiPS_bn_bls/bls12_timeprint.h>

 double BLS12_MILLER_TATE,BLS12_MILLER_PLAINATE,BLS12_MILLER_OPTATE;
 double BLS12_FINALEXP_PLAIN,BLS12_FINALEXP_OPT;
 double BLS12_G1SCM_PLAIN,BLS12_G1SCM_2SPLIT;
 double BLS12_G2SCM_PLAIN,BLS12_G2SCM_2SPLIT,BLS12_G2SCM_4SPLIT;
 double BLS12_G3EXP_PLAIN,BLS12_G3EXP_2SPLIT,BLS12_G3EXP_4SPLIT;

void BLS12_print_tate_time(){
    printf("miller tate time:%.2f[ms]\n",BLS12_MILLER_TATE);
}

void BLS12_print_plain_ate_time(){
    printf("miller plain-ate time:%.2f[ms]\n",BLS12_MILLER_PLAINATE);
}

void BLS12_print_opt_ate_time(){
    printf("miller opt-ate time:%.2f[ms]\n",BLS12_MILLER_OPTATE);
}

void BLS12_print_final_exp_plain_time(){
    printf("plain-final exp time:%.2f[ms]\n",BLS12_FINALEXP_PLAIN);
}

void BLS12_print_final_exp_optimal_time(){
    printf("optimal-final exp time:%.2f[ms]\n",BLS12_FINALEXP_OPT);
}

void BLS12_print_plain_G1_scm_time(){
    printf("plain-G1-SCM time:%.2f[ms]\n",BLS12_G1SCM_PLAIN);
}

void BLS12_print_2split_G1_scm_time(){
    printf("2split-G1-SCM time:%.2f[ms]\n",BLS12_G1SCM_2SPLIT);
}

void BLS12_print_plain_G2_scm_time(){
    printf("plain-G2-SCM time:%.2f[ms]\n",BLS12_G2SCM_PLAIN);
}

void BLS12_print_2split_G2_scm_time(){
    printf("2split-G2-SCM time:%.2f[ms]\n",BLS12_G2SCM_2SPLIT);
}

void BLS12_print_4split_G2_scm_time(){
    printf("4split-G2-SCM time:%.2f[ms]\n",BLS12_G2SCM_4SPLIT);
}

void BLS12_print_plain_G3_exp_time(){
    printf("plain-G3-SCM time:%.2f[ms]\n",BLS12_G3EXP_PLAIN);
}

void BLS12_print_2split_G3_exp_time(){
    printf("2split-G3-SCM time:%.2f[ms]\n",BLS12_G3EXP_2SPLIT);
}

void BLS12_print_4split_G3_exp_time(){
    printf("4split-G3-SCM time:%.2f[ms]\n",BLS12_G3EXP_4SPLIT);
}
