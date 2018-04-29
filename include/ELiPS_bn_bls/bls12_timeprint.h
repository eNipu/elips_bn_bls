//
//  bls12_timeprint.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/5/18.

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
