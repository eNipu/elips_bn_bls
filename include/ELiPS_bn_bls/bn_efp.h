//
//  bn_efp.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/1/18.

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

#ifndef bn_efp_h
#define bn_efp_h

#include <ELiPS_bn_bls/curve_dtypes.h>

extern void EFp_init(EFp *P);
extern void EFp_clear(EFp *P);
extern void EFp_printf(EFp *P,char *str);
extern void EFp_set(EFp *ANS,EFp *P);
extern void EFp_set_ui(EFp *ANS,unsigned long int UI);
extern void EFp_set_mpz(EFp *ANS,mpz_t A);
extern void EFp_set_neg(EFp *ANS,EFp *P);
extern void EFp_rational_point_bn(EFp *P);
extern void EFp_rational_point_bls12(EFp *P);
extern void EFp_ECD(EFp *ANS,EFp *P);
extern void EFp_ECA(EFp *ANS,EFp *P1,EFp *P2);
extern void EFp_SCM(EFp *ANS,EFp *P,mpz_t scalar);

#endif /* bn_efp_h */
