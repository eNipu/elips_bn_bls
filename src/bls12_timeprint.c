//
//  bls12_timeprint.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/5/18.

/*
 * ELiPS is an Efficient Library for Pairing-based Systems
 * Copyright (C) 2008-2018 ELiPS Authors. Please refer to the COPYRIGHT file
 * for contact information.
 *
 * This file is part of ELiPS Project. ELiPS is legal property of its
 * developers, whose names are not listed here. 
 *
 * ELiPS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * ELiPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ELiPS. If not, see <http://www.gnu.org/licenses/>.
 */

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
