//
//  curve_dtypes.h
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

#ifndef curve_dtypes_h
#define curve_dtypes_h

#include <ELiPS_bn_bls/bn_fp12.h>
#include <ELiPS_bn_bls/field_dtype.h>

typedef struct EFp EFp;
struct EFp{
    Fp x,y;
    int infinity;
};

typedef struct EFp2 EFp2;
struct EFp2{
    Fp2 x,y;
    int infinity;
};

typedef struct EFp6 EFp6;
struct EFp6{
    Fp6 x,y;
    int infinity;
};

typedef struct EFp12 EFp12;
struct EFp12{
    Fp12 x,y;
    int infinity;
};

typedef struct EFp4 EFp4;
struct EFp4{
    struct Fp4 x,y;
    int infinity;
};

typedef struct EFp8 EFp8;
struct EFp8{
    struct Fp8 x,y;
    int infinity;
};

typedef struct EFp16 EFp16;
struct EFp16{
    struct Fp16 x,y;
    int infinity;
};


#endif /* curve_dtypes_h */
