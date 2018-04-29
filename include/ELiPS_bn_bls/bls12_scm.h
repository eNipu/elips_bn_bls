//
//  bls12_scm.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/6/18.

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

#ifndef bls12_scm_h
#define bls12_scm_h

#include <ELiPS_bn_bls/bn_efp12.h>
#include <ELiPS_bn_bls/bn_utils.h>
#include <ELiPS_bn_bls/bls12_timeprint.h>
#include <ELiPS_bn_bls/bls12_twist.h>
#include <ELiPS_bn_bls/bls12_skew_frobenius.h>



//G1
extern void BLS12_plain_G1_scm(EFp12 *ANS,EFp12 *P,mpz_t scalar);
extern void BLS12_2split_G1_scm(EFp12 *ANS,EFp12 *P,mpz_t scalar);
//G2
extern void BLS12_plain_G2_scm(EFp12 *ANS,EFp12 *Q,mpz_t scalar);
extern void BLS12_2split_G2_scm(EFp12 *ANS,EFp12 *Q,mpz_t scalar);
extern void BLS12_4split_G2_scm(EFp12 *ANS,EFp12 *Q,mpz_t scalar);

#endif /* bls12_scm_h */
