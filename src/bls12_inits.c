//
//  bls12_inits.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/5/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#include <ELiPS_bn_bls/bls12_inits.h>

void bls12_inits(){
    init_bls12_settings();
    init_precoms(2);
}
