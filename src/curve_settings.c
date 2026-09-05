//
//  bn_settings.c
//  bn_Single_File_Elips
//
//  Created by Y.Nanjo and M.A.A Khandaker on 1/25/18.
//  Copyright © 2018 Khandaker Md. Al-Amin. All rights reserved.
//

#include <ELiPS_bn_bls/curve_settings.h>
//#include "bn_bls12_precoms.h"


struct curve_params curve_parameters;
char X_binary[bn_X_length+1];
char X_binary_opt[bn_X_length+3];
//---------------------------------------------------------------------
int bls12_X_length;
mpz_t bls12_X;
int bls12_X_binary[78];

int init_bn_settings(){
    init_bn_parameters();
    
    generate_bn_mother_parameter();
    /* A6: these can fail (a candidate may not be prime) and their results were
     * discarded, so setup carried on with an unusable curve. */
    if(generate_bn_prime()!=1) return 0;
    if(generate_bn_order()!=1) return 0;
    generate_bn_trace();
    
    weil();
    set_bn_curve_parameter();
    return 1;
}

int init_bls12_settings(){
    init_bls12_parameters();
    
    bls12_generate_X();
    if(bls12_generate_prime()!=1) return 0;   /* A6 */
    if(bls12_generate_order()!=1) return 0;
    bls12_generate_trace();
    
    bls12_weil();
    bls12_set_curve_parameter();
    return 1;
}


/*============================================================================*/
/*  init                                                                      */
/*============================================================================*/
void init_bn_parameters(){
    //parameters
    mpz_init(curve_parameters.prime);
    mpz_init(curve_parameters.X);
    mpz_init(curve_parameters.trace_t);
    mpz_init(curve_parameters.order);
    mpz_init(curve_parameters.EFp_total);
    mpz_init(curve_parameters.EFp2_total);
    mpz_init(curve_parameters.EFp6_total);
    mpz_init(curve_parameters.EFp12_total);
    /* M2: EFpd_total holds #E(F_p^12); weil() writes it and both G2
     * generators read it, but nothing ever initialised it. */
    mpz_init(curve_parameters.EFpd_total);
    /* M3: curve_a is read by print_curve_parameters() and released by
     * clear_parameters(), so it must exist even though BN does not use it. */
    mpz_init(curve_parameters.curve_a);
    mpz_init(curve_parameters.curve_b);
    
    int i;
    for(i=0; i<bn_X_length+1; i++){
        X_binary[i]=0;
    }
    for(i=0; i<bn_X_length+3; i++){
        X_binary_opt[i]=0;
    }
}

void init_bls12_parameters(){
    //parameters
    mpz_init(curve_parameters.prime);
    mpz_init(curve_parameters.X);
    mpz_init(curve_parameters.trace_t);
    mpz_init(curve_parameters.order);
    mpz_init(curve_parameters.EFp_total);
    mpz_init(curve_parameters.EFp2_total);
    mpz_init(curve_parameters.EFp6_total);
    mpz_init(curve_parameters.EFp12_total);
    mpz_init(curve_parameters.EFpd_total);   /* M2 */
    mpz_init(curve_parameters.curve_a);      /* M3 */
    mpz_init(curve_parameters.curve_b);
    
    int i;
    for(i=0; i<bn_X_length+1; i++){
        X_binary[i]=0;
    }
    for(i=0; i<bn_X_length+3; i++){
        X_binary_opt[i]=0;
    }
    
    bls12_X_length=77;
    for(i=0; i<bls12_X_length+1; i++){
        bls12_X_binary[i]=0;
    }
}




void generate_bn_mother_parameter(){
    int i;
    mpz_t buf;
    mpz_init(buf);
    
    //X_binary
    X_binary[114]=1;
    X_binary[101]=1;
    X_binary[14]=-1;
    X_binary[0]=-1;
    
    //X_binary_opt
    X_binary_opt[116]=1;
    X_binary_opt[115]=1;
    X_binary_opt[103]=1;
    X_binary_opt[102]=1;
    X_binary_opt[16]=-1;
    X_binary_opt[15]=-1;
    X_binary_opt[2]=-1;
    
    //X
    mpz_set_ui(curve_parameters.X,0);
    for(i=bn_X_length; i>=0; i--){
        if(X_binary[i]==1){
            mpz_ui_pow_ui(buf,2,i);
            mpz_add(curve_parameters.X,curve_parameters.X,buf);
        }else if(X_binary[i]==-1){
            mpz_ui_pow_ui(buf,2,i);
            mpz_sub(curve_parameters.X,curve_parameters.X,buf);
        }
    }
    
    mpz_clear(buf);
}

int generate_bn_prime(){
    mpz_t buf,result;
    mpz_init(buf);
    mpz_init(result);
    
    //prime
    mpz_pow_ui(buf,curve_parameters.X,4);
    mpz_mul_ui(buf,buf,36);
    mpz_set(result,buf);
    mpz_pow_ui(buf,curve_parameters.X,3);
    mpz_mul_ui(buf,buf,36);
    mpz_add(result,result,buf);
    mpz_pow_ui(buf,curve_parameters.X,2);
    mpz_mul_ui(buf,buf,24);
    mpz_add(result,result,buf);
    mpz_mul_ui(buf,curve_parameters.X,6);
    mpz_add(result,result,buf);
    mpz_add_ui(result,result,1);
    
    //isprime
    if(mpz_probab_prime_p(result,25)==0){
        mpz_clear(buf);
        mpz_clear(result);
        return 0;
    }else{
        mpz_set(curve_parameters.prime,result);
        mpz_clear(buf);
        mpz_clear(result);
        return 1;
    }
}

int  generate_bn_order(){
    mpz_t buf,result;
    mpz_init(buf);
    mpz_init(result);
    
    //prime
    mpz_pow_ui(buf,curve_parameters.X,4);
    mpz_mul_ui(buf,buf,36);
    mpz_set(result,buf);
    mpz_pow_ui(buf,curve_parameters.X,3);
    mpz_mul_ui(buf,buf,36);
    mpz_add(result,result,buf);
    mpz_pow_ui(buf,curve_parameters.X,2);
    mpz_mul_ui(buf,buf,18);
    mpz_add(result,result,buf);
    mpz_mul_ui(buf,curve_parameters.X,6);
    mpz_add(result,result,buf);
    mpz_add_ui(result,result,1);
    
    //isprime
    if(mpz_probab_prime_p(result,25)==0){
        mpz_clear(buf);
        mpz_clear(result);
        return 0;
    }else{
        mpz_set(curve_parameters.order,result);
        mpz_clear(buf);
        mpz_clear(result);
        return 1;
    }
}


void generate_bn_trace(){
    mpz_t buf;
    mpz_init(buf);
    
    mpz_pow_ui(buf,curve_parameters.X,2);
    mpz_mul_ui(buf,buf,6);
    mpz_add_ui(curve_parameters.trace_t,buf,1);
    
    mpz_clear(buf);
}

void set_bn_curve_parameter(){
    mpz_set_ui(curve_parameters.curve_b,4);
}

void weil(){
    mpz_t t2,t6,t12,p2,p6,buf;
    mpz_init(t2);
    mpz_init(t6);
    mpz_init(t12);
    mpz_init(p2);
    mpz_init(p6);
    mpz_init(buf);
    
    //EFp_total
    mpz_add_ui(buf,curve_parameters.prime,1);
    mpz_sub(curve_parameters.EFp_total,buf,curve_parameters.trace_t);
    
    //t2←α^2+β^2
    mpz_pow_ui(t2,curve_parameters.trace_t,2);
    mpz_mul_ui(buf,curve_parameters.prime,2);
    mpz_sub(t2,t2,buf);
    mpz_pow_ui(p2,curve_parameters.prime,2);
    
    //α^6+β^6
    mpz_pow_ui(t6,t2,3);
    mpz_mul(buf,t2,p2);
    mpz_mul_ui(buf,buf,3);
    mpz_sub(t6,t6,buf);
    mpz_pow_ui(p6,p2,3);
    
    //α^12+β^12
    mpz_pow_ui(t12,t6,2);
    mpz_mul_ui(buf,p6,2);
    mpz_sub(t12,t12,buf);
    
    //EFp12_total
    mpz_pow_ui(buf,p6,2);
    mpz_sub(buf,buf,t12);
    mpz_add_ui(curve_parameters.EFpd_total,buf,1);
    
    mpz_clear(t2);
    mpz_clear(t6);
    mpz_clear(t12);
    mpz_clear(p2);
    mpz_clear(p6);
    mpz_clear(buf);
}


static void print_bn_blscurve() {
    gmp_printf("\nelliptic curve\n");
    gmp_printf("E:y^2=x^3-%Zd\n",curve_parameters.curve_b);
}

void print_curve_parameters(){
    printf("====================================================================================\n");
    printf("bn12\n\n");
    gmp_printf("parameters\n");
    gmp_printf("X     (%dbit length) : %Zd \n",(int)mpz_sizeinbase(curve_parameters.X,2),curve_parameters.X);
    gmp_printf("prime (%dbit length) : %Zd \n",(int)mpz_sizeinbase(curve_parameters.prime,2),curve_parameters.prime);
    gmp_printf("order (%dbit length) : %Zd \n",(int)mpz_sizeinbase(curve_parameters.order,2),curve_parameters.order);
    gmp_printf("trace (%dbit length) : %Zd \n",(int)mpz_sizeinbase(curve_parameters.trace_t,2),curve_parameters.trace_t);
    
    if (mpz_cmp_ui(curve_parameters. curve_b, 0) > 0) {
         print_bn_blscurve();
    }
//    gmp_printf("\nmodulo polynomial\n");
//    gmp_printf("Fp2  : f(x) = x^2+1\n");
//    gmp_printf("Fp6  : f(x) = x^3-(alpha+1)\n");
//    gmp_printf("Fp12 : f(x) = x^2-beta\n");
    
//    gmp_printf("\nnumber of the total rational points\n");
//    gmp_printf("EFp total   : %Zd\n",curve_parameters.EFp_total);
//    gmp_printf("EFp12 total : %Zd\n",curve_parameters.EFpd_total);
    
//    gmp_printf("\ncubic root of 1\n");
//    gmp_printf("epsilon1 : %Zd\n",epsilon1);
//    gmp_printf("epsilon2 : %Zd\n",epsilon2);
    
}

void bls12_generate_X(){
    int i;
    mpz_t buf,set_2;
    mpz_init(buf);
    mpz_init(set_2);
    mpz_set_ui(set_2,2);
    
    //bls12_X_binary
    bls12_X_binary[77]=-1;
    bls12_X_binary[50]=1;
    bls12_X_binary[33]=1;
    
    //bls12_X
    mpz_init(bls12_X);
    mpz_set_ui(bls12_X,0);
    for(i=bls12_X_length; i>=0; i--){
        if(bls12_X_binary[i]==1){
            mpz_pow_ui(buf,set_2,i);
            mpz_add(bls12_X,bls12_X,buf);
        }else if(bls12_X_binary[i]==-1){
            mpz_pow_ui(buf,set_2,i);
            mpz_sub(bls12_X,bls12_X,buf);
        }
    }
    
    mpz_clear(buf);
    mpz_clear(set_2);
}

int  bls12_generate_prime(){
    mpz_t result,buf1,buf2,modtest;
    mpz_init(result);
    mpz_init(buf1);
    mpz_init(buf2);
    mpz_init(modtest);
    
    mpz_sub_ui(result,bls12_X,1);
    mpz_pow_ui(result,result,2);
    
    mpz_pow_ui(buf1,bls12_X,4);
    mpz_pow_ui(buf2,bls12_X,2);
    mpz_sub(buf1,buf1,buf2);
    mpz_add_ui(buf1,buf1,1);
    
    mpz_mul(result,result,buf1);
    
    //check div3
    mpz_mod_ui(modtest,result,3);
    if(mpz_cmp_ui(modtest,0)!=0){
        mpz_clear(result);
        mpz_clear(buf1);
        mpz_clear(buf2);
        mpz_clear(modtest);
        return 0;
    }
    
    mpz_tdiv_q_ui(result,result,3);
    mpz_add(result,result,bls12_X);
    
    //isprime
    if(mpz_probab_prime_p(result,25)==0){
        mpz_clear(result);
        mpz_clear(buf1);
        mpz_clear(buf2);
        mpz_clear(modtest);
        return 0;
    }
    
    mpz_set(curve_parameters.prime,result);
    
    mpz_clear(result);
    mpz_clear(buf1);
    mpz_clear(buf2);
    mpz_clear(modtest);
    return 1;
}

int  bls12_generate_order(){
    mpz_t buf1,buf2;
    mpz_init(buf1);
    mpz_init(buf2);
    
    mpz_pow_ui(buf1,bls12_X,4);
    mpz_pow_ui(buf2,bls12_X,2);
    mpz_sub(curve_parameters.order,buf1,buf2);
    mpz_add_ui(curve_parameters.order,curve_parameters.order,1);
    
    mpz_clear(buf1);
    mpz_clear(buf2);
    /* A6: 1 means success, matching generate_bn_prime and generate_bn_order.
     * This used to return 0, the opposite of every sibling. */
    return 1;
}

void bls12_generate_trace(){
    mpz_add_ui(curve_parameters.trace_t,bls12_X,1);
}

void bls12_set_curve_parameter(){
    mpz_set_ui(curve_parameters.curve_b,4);
}

void bls12_weil(){
    mpz_t t2,t6,t12,p2,p6,buf;
    mpz_init(t2);
    mpz_init(t6);
    mpz_init(t12);
    mpz_init(p2);
    mpz_init(p6);
    mpz_init(buf);
    
    //EFp_total
    mpz_add_ui(buf,curve_parameters.prime,1);
    mpz_sub(curve_parameters.EFp_total,buf,curve_parameters.trace_t);
    
    //t2←α^2+β^2
    mpz_pow_ui(t2,curve_parameters.trace_t,2);
    mpz_mul_ui(buf,curve_parameters.prime,2);
    mpz_sub(t2,t2,buf);
    mpz_pow_ui(p2,curve_parameters.prime,2);
    
    //α^6+β^6
    mpz_pow_ui(t6,t2,3);
    mpz_mul(buf,t2,p2);
    mpz_mul_ui(buf,buf,3);
    mpz_sub(t6,t6,buf);
    mpz_pow_ui(p6,p2,3);
    
    //α^12+β^12
    mpz_pow_ui(t12,t6,2);
    mpz_mul_ui(buf,p6,2);
    mpz_sub(t12,t12,buf);
    
    //EFp12_total
    mpz_pow_ui(buf,p6,2);
    mpz_sub(buf,buf,t12);
    mpz_add_ui(curve_parameters.EFpd_total,buf,1);
    
    mpz_clear(t2);
    mpz_clear(t6);
    mpz_clear(t12);
    mpz_clear(p2);
    mpz_clear(p6);
    mpz_clear(buf);
}

void bls12_print_parameters(){
    printf("====================================================================================\n");
    printf("bls12\n\n");
    gmp_printf("parameters\n");
    gmp_printf("X     (%dbit length) : %Zd \n",(int)mpz_sizeinbase(bls12_X,2),bls12_X);
    gmp_printf("prime (%dbit length) : %Zd \n",(int)mpz_sizeinbase(curve_parameters.prime,2),curve_parameters.prime);
    gmp_printf("order (%dbit length) : %Zd \n",(int)mpz_sizeinbase(curve_parameters.order,2),curve_parameters.order);
    gmp_printf("trace (%dbit length) : %Zd \n",(int)mpz_sizeinbase(curve_parameters.trace_t,2),curve_parameters.trace_t);
    
    gmp_printf("\nelliptic curve\n");
    gmp_printf("E:y^2=x^3+%Zd\n",curve_parameters.curve_b);   /* A7 */
    
    gmp_printf("\nmodulo polynomial\n");
    gmp_printf("Fp2  : f(x) = x^2+1\n");
    gmp_printf("Fp6  : f(x) = x^3-(alpha+1)\n");
    gmp_printf("Fp12 : f(x) = x^2-beta\n");
    
}
//---------------------------------------------------------------------




