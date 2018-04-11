//
//  bn_utils.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/1/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#include <ELiPS_bn_bls/bn_utils.h>

struct timeval t0,t1;

float MILLER_TATE,MILLER_PLAINATE,MILLER_OPTATE,MILLER_XATE;
float FINALEXP_PLAIN,FINALEXP_OPT;
float G1SCM_PLAIN,G1SCM_2SPLIT;
float G2SCM_PLAIN,G2SCM_2SPLIT,G2SCM_4SPLIT;
float G3SCM_PLAIN,G3SCM_2SPLIT,G3SCM_4SPLIT;

float timedifference_msec(struct timeval t0, struct timeval t1){
    return (t1.tv_sec - t0.tv_sec) * 1000.0f + (t1.tv_usec - t0.tv_usec) / 1000.0f;
}

float timedifference_usec(struct timeval t0, struct timeval t1){
    return (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec);
}

void BN12_print_tate_time(){
    printf("miller tate time:%.2f[ms]\n",MILLER_TATE);
}

void BN12_print_plain_ate_time(){
    printf("miller plain-ate time:%.2f[ms]\n",MILLER_PLAINATE);
}

void BN12_print_opt_ate_time(){
    printf("miller opt-ate time:%.2f[ms]\n",MILLER_OPTATE);
}

void BN12_print_x_ate_time(){
    printf("miller x-ate time:%.2f[ms]\n",MILLER_XATE);
}

void BN12_print_final_exp_plain_time(){
    printf("plain-final exp time:%.2f[ms]\n",FINALEXP_PLAIN);
}

void BN12_print_final_exp_optimal_time(){
    printf("optimal-final exp time:%.2f[ms]\n",FINALEXP_OPT);
}

void BN12_print_G1_plain_time(){
    printf("plain-G1-SCM time:%.2f[ms]\n",G1SCM_PLAIN);
}

void BN12_print_G1_2split_time(){
    printf("2split-G1-SCM time:%.2f[ms]\n",G1SCM_2SPLIT);
}

void BN12_print_G2_plain_time(){
    printf("plain-G2-SCM time:%.2f[ms]\n",G2SCM_PLAIN);
}

void BN12_print_G2_2split_time(){
    printf("2split-G2-SCM time:%.2f[ms]\n",G2SCM_2SPLIT);
}

void BN12_print_G2_4split_time(){
    printf("4split-G2-SCM time:%.2f[ms]\n",G2SCM_4SPLIT);
}

void BN12_print_G3_plain_time(){
    printf("plain-G3-SCM time:%.2f[ms]\n",G3SCM_PLAIN);
}

void BN12_print_G3_2split_time(){
    printf("2split-G3-SCM time:%.2f[ms]\n",G3SCM_2SPLIT);
}

void BN12_print_G3_4split_time(){
    printf("4split-G3-SCM time:%.2f[ms]\n",G3SCM_4SPLIT);
}


