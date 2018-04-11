//
//  field_dtype.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 1/31/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

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
