//
//  bls12_G3_exp.h
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

#ifndef bls12_G3_exp_h
#define bls12_G3_exp_h

#include <ELiPS_bn_bls/bn_utils.h>
#include <ELiPS_bn_bls/bls12_timeprint.h>
#include <ELiPS_bn_bls/bls12_frobenius.h>

//G3
extern void bls12_plain_G3_exp(Fp12 *ANS,Fp12 *A,mpz_t scalar);
extern void bls12_2split_G3_exp(Fp12 *ANS,Fp12 *A,mpz_t scalar);
extern void bls12_4split_G3_exp( Fp12 *ANS,Fp12 *A,mpz_t scalar);

#endif /* bls12_G3_exp_h */
