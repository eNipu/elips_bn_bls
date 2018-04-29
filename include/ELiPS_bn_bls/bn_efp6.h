//
//  bn_efp6.h
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

#ifndef bn_efp6_h
#define bn_efp6_h

#include <ELiPS_bn_bls/bn_efp2.h>


extern void EFp6_init(EFp6 *P);
extern void EFp6_clear(EFp6 *P);
extern void EFp6_printf(EFp6 *P,char *str);
extern void EFp6_set(EFp6 *ANS,EFp6 *P);
extern void EFp6_set_ui(EFp6 *ANS,unsigned long int UI);
extern void EFp6_set_mpz(EFp6 *ANS,mpz_t A);
extern void EFp6_set_neg(EFp6 *ANS,EFp6 *P);
extern void EFp6_rational_point(EFp6 *P);
extern void EFp6_ECD(EFp6 *ANS,EFp6 *P);
extern void EFp6_ECA(EFp6 *ANS,EFp6 *P1,EFp6 *P2);
extern void EFp6_SCM(EFp6 *ANS,EFp6 *P,mpz_t scalar);

#endif /* bn_efp6_h */
