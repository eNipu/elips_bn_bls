//
//  bn_inits.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/1/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#include <ELiPS_bn_bls/bn_inits.h>

void init_bn(void){
    init_bn_settings();
    init_precoms(1);
}
