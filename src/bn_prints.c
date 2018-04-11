//
//  bn_prints.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/1/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#include <ELiPS_bn_bls/bn_prints.h>

void BN12_print_G1_point(){
    EFp12_printf(&G1_P,"G1_P\n");
    printf("\n");
}

void BN12_print_G2_point(){
    EFp12_printf(&G2_Q,"G2_Q\n");
    printf("\n");
}
