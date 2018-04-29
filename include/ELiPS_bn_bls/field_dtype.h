//
//  field_dtype.h
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

#ifndef field_dtype_h
#define field_dtype_h

#include <ELiPS_bn_bls/Commont_headers.h>


typedef struct Fp Fp;
struct Fp{
    mpz_t x0;
};

typedef struct Fp2 Fp2;
struct Fp2{
    struct Fp x0, x1;
};

typedef struct Fp4 Fp4;
struct Fp4{
    struct Fp2 x0,x1;
};

typedef struct Fp6 Fp6;
struct Fp6{
    Fp2 x0,x1,x2;
};

typedef struct Fp8 Fp8;
struct Fp8{
    struct Fp4 x0,x1;
};

typedef struct Fp12 Fp12;
struct Fp12{
    Fp6 x0,x1;
};

typedef struct Fp16 Fp16;
struct Fp16{
    struct Fp8 x0,x1;
};


#endif /* field_dtype_h */
