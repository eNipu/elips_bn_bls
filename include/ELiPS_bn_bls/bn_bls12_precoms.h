//
//  bn_bls12_precoms.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 1/31/18.

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
