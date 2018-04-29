//
//  bn_efp2.h
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

#ifndef bn_efp2_h
#define bn_efp2_h

#include <ELiPS_bn_bls/bn_efp.h>

extern void EFp2_init(EFp2 *P);
extern void EFp2_clear(EFp2 *P);
extern void EFp2_printf(EFp2 *P,char *str);
extern void EFp2_set(EFp2 *ANS,EFp2 *P);
extern void EFp2_set_ui(EFp2 *ANS,unsigned long int UI);
extern void EFp2_set_mpz(EFp2 *ANS,mpz_t A);
extern void EFp2_set_neg(EFp2 *ANS,EFp2 *P);
extern void EFp2_rational_point(EFp2 *P);
extern void EFp2_ECD(EFp2 *ANS,EFp2 *P);
extern void EFp2_ECA(EFp2 *ANS,EFp2 *P1,EFp2 *P2);
extern void EFp2_SCM(EFp2 *ANS,EFp2 *P,mpz_t scalar);

#endif /* bn_efp2_h */
