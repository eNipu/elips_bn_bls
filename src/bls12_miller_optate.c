//
//  bls12_miller_optate.c
//  elips_refactoring
//
//  Created by Y.Nanjo and M.A.A Khandaker on 2/6/18.
//  Copyright © 2018 ISec Lab. All rights reserved.
//

#include <ELiPS_bn_bls/bls12_miller_optate.h>

void BLS12_Miller_algo_for_opt_ate(Fp12 *ANS,EFp12 *P,EFp12 *Q){
    EFp2 T;
    EFp2_init(&T);
    EFp2 mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2_neg;
    EFp2_init(&mapped_Q);
    EFp2_init(&mapped_Q_neg);
    EFp2_init(&mapped_Q1);
    EFp2_init(&mapped_Q2_neg);
    EFp mapped_P;
    EFp_init(&mapped_P);
    Fp12 f;
    Fp12_init(&f);
    Fp L;
    Fp_init(&L);
    int i;
    
    //set
    BLS12_EFp12_to_EFp(&mapped_P,P);  //set P
    BLS12_EFp12_to_EFp2(&mapped_Q,Q);//set mapped_Q
    BLS12_Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
    
    EFp2_set_neg(&mapped_Q_neg,&mapped_Q);//set mapped_Q_neg
    
    EFp2_set(&T,&mapped_Q_neg);        //set T
    Fp12_set_ui(&f,0);            //set f
    Fp_set_ui(&f.x0.x0.x0,1);
    //miller
    for(i=BLS12_X_length-1; i>=0; i--){
        switch(BLS12_X_binary[i]){
            case 0:
                BLS12_ff_ltt(&f,&T,&mapped_P,&L);
                break;
            case 1:
                BLS12_ff_ltt(&f,&T,&mapped_P,&L);
                BLS12_f_ltq(&f,&T,&mapped_Q,&mapped_P,&L);
                break;
            case -1:
                BLS12_ff_ltt(&f,&T,&mapped_P,&L);
                BLS12_f_ltq(&f,&T,&mapped_Q_neg,&mapped_P,&L);
                break;
            default:
                break;
        }
        
    }
    
    Fp12_set(ANS,&f);
    
    Fp12_clear(&f);
    EFp2_clear(&T);
    EFp2_clear(&mapped_Q);
    EFp2_clear(&mapped_Q_neg);
    EFp2_clear(&mapped_Q1);
    EFp2_clear(&mapped_Q2_neg);
    EFp_clear(&mapped_P);
    Fp_clear(&L);
}
