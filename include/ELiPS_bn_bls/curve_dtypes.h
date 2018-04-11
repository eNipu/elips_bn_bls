//
//  curve_dtypes.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/1/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

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
