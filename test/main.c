//
//  main.c
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 1/25/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#include <ELiPS_bn_bls/bn_inits.h>
#include <ELiPS_bn_bls/bn_clears.h>
#include <ELiPS_bn_bls/bn_pairing_test.h>

#include <ELiPS_bn_bls/bls12_inits.h>
#include <ELiPS_bn_bls/bls12_pairings.h>
#include <ELiPS_bn_bls/bls12_test_pairings.h>


int main(int argc, const char * argv[]) {
    
    bls12_inits();
    BLS12_print_parameters();
    BLS12_test_opt_ate_pairing();
    clear_parameters();
    
//    init_bn();
//    print_curve_parameters();
//    BN12_test_x_ate_pairing();
//    BN12_test_plain_ate_pairing();
//    BN12_test_opt_ate_pairing();
//    clear_parameters();
 
    
    return 0;
}
