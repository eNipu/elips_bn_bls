//
//  bn_pairing_test.h
//  elips_refactoring
//
//  Created by Khandaker Md. Al-Amin on 2/2/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#ifndef bn_pairing_test_h
#define bn_pairing_test_h

#include <ELiPS_bn_bls/bn_pairings.h>
#include <ELiPS_bn_bls/bn_generate_points.h>

extern void BN12_test_tate_pairing(void);
extern void BN12_test_plain_ate_pairing(void);
extern void BN12_test_opt_ate_pairing(void);
extern void BN12_test_x_ate_pairing(void);

#endif /* bn_pairing_test_h */
