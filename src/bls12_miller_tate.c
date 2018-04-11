//
//  bls12_miller_tate.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/6/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#include <ELiPS_bn_bls/bls12_miller_tate.h>


void BLS12_Miller_algo_for_tate(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    EFp Tmp_P,T;
    EFp_init(&Tmp_P);
    EFp_init(&T);
    Fp12 f;
    Fp12_init(&f);
    int i,length;
    length=(int)mpz_sizeinbase(curve_parameters.order,2);
    char binary[length];
    mpz_get_str(binary,2,curve_parameters.order);
    
    //set
    BLS12_EFp12_to_EFp(&Tmp_P,P);
    EFp_set(&T,&Tmp_P);
    Fp12_set_ui(&f,0);
    Fp_set_ui(&f.x0.x0.x0,1);
    
    //miller
    for(i=1; binary[i]!='\0'; i++){
        BLS12_ff_ltt_vtt_for_tate(&f,&T,Q);
        if(binary[i]=='1'){
            BLS12_f_ltp_vtp_for_tate(&f,&T,&Tmp_P,Q);
        }
    }
    Fp12_set(ANS,&f);
    
    Fp12_clear(&f);
    EFp_clear(&T);
    EFp_clear(&Tmp_P);
}

