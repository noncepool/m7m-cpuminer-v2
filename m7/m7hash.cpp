
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <gmp.h>
#include <limits>
#include "uint256.h"
#include "sph_sha2.h"
#include "sph_keccak.h"
#include "sph_haval.h"
#include "sph_tiger.h"
#include "sph_whirlpool.h"
#include "sph_ripemd.h"
#include "magimath.h"

static void mpz_set_uint256(mpz_t r, char *u)
{
    mpz_import(r, 32 / sizeof(unsigned long), -1, sizeof(unsigned long), -1, 0, u);
}

static void mpz_get_uint256(mpz_t r, char *u)
{
    u=0;
    mpz_export(u, 0, -1, sizeof(unsigned long), -1, 0, r);
}

static void mpz_set_uint512(mpz_t r, char *u)
{
    mpz_import(r, 64 / sizeof(unsigned long), -1, sizeof(unsigned long), -1, 0, u);
}

static void set_one_if_zero(char *hash512) {
    for (int i = 0; i < 32; i++) {
        if (hash512[i] != 0) {
            return;
        }
    }
    hash512[0] = 1;
}

#define BITS_PER_DIGIT 3.32192809488736234787
#define EPS (std::numeric_limits<double>::epsilon())

#define NM7M 5
#define SW_DIVS 5

template<typename T1>
inline uint256 hash_M7M_v2(const T1 pbegin, const T1 pend, const unsigned int nnNonce)
{
    sph_sha256_context       ctx_sha256;
    sph_sha512_context       ctx_sha512;
    sph_keccak512_context    ctx_keccak;
    sph_whirlpool_context    ctx_whirlpool;
    sph_haval256_5_context   ctx_haval;
    sph_tiger_context        ctx_tiger;
    sph_ripemd160_context    ctx_ripemd;
    static unsigned char pblank[1];
    int bytes, nnNonce2 = (int)(nnNonce/2);

    uint512 hash[7];
    uint256 finalhash;
    for(int i=0; i < 7; i++)
	    hash[i] = 0;

    const void* ptr = (pbegin == pend ? pblank : static_cast<const void*>(&pbegin[0]));
    size_t sz = (pend - pbegin) * sizeof(pbegin[0]);

    sph_sha256_init(&ctx_sha256);
    sph_sha256 (&ctx_sha256, ptr, sz);
    sph_sha256_close(&ctx_sha256, static_cast<void*>(&hash[0]));

    sph_sha512_init(&ctx_sha512);
    sph_sha512 (&ctx_sha512, ptr, sz);
    sph_sha512_close(&ctx_sha512, static_cast<void*>(&hash[1]));
    
    sph_keccak512_init(&ctx_keccak);
    sph_keccak512 (&ctx_keccak, ptr, sz);
    sph_keccak512_close(&ctx_keccak, static_cast<void*>(&hash[2]));

    sph_whirlpool_init(&ctx_whirlpool);
    sph_whirlpool (&ctx_whirlpool, ptr, sz);
    sph_whirlpool_close(&ctx_whirlpool, static_cast<void*>(&hash[3]));
    
    sph_haval256_5_init(&ctx_haval);
    sph_haval256_5 (&ctx_haval, ptr, sz);
    sph_haval256_5_close(&ctx_haval, static_cast<void*>(&hash[4]));

    sph_tiger_init(&ctx_tiger);
    sph_tiger (&ctx_tiger, ptr, sz);
    sph_tiger_close(&ctx_tiger, static_cast<void*>(&hash[5]));

    sph_ripemd160_init(&ctx_ripemd);
    sph_ripemd160 (&ctx_ripemd, ptr, sz);
    sph_ripemd160_close(&ctx_ripemd, static_cast<void*>(&hash[6]));

    mpz_t bns[8];

    for(int i=0; i < 7; i++){
	if(hash[i]==0)
	    hash[i] = 1;
	mpz_init(bns[i]);
	mpz_set_uint512(bns[i],hash[i]);
    }

    mpz_init(bns[7]);
    mpz_set_ui(bns[7],0);
    for(int i=0; i < 7; i++)
	mpz_add(bns[7], bns[7], bns[i]);

    mpz_t product;
    mpz_init(product);
    mpz_set_ui(product,1);
    for(int i=0; i < 8; i++){
	mpz_mul(product,product,bns[i]);
    }
    mpz_pow_ui(product, product, 2);

    bytes = mpz_sizeinbase(product, 256);
    char *data = (char*)malloc(bytes);
    mpz_export(data, NULL, -1, 1, 0, 0, product);

    sph_sha256_init(&ctx_sha256);
    sph_sha256 (&ctx_sha256, data, bytes);
    sph_sha256_close(&ctx_sha256, static_cast<void*>(&finalhash));
    free(data);

    int digits=(int)((sqrt((double)(nnNonce2))*(1.+EPS))/9000+75);
    int iterations=20;
    mpf_set_default_prec((long int)(digits*BITS_PER_DIGIT+16));

    mpz_t magipi;
    mpz_t magisw;
    mpf_t magifpi;
    mpf_t mpa1, mpb1, mpt1, mpp1;
    mpf_t mpa2, mpb2, mpt2, mpp2;
    mpf_t mpsft;

    mpz_init(magipi);
    mpz_init(magisw);
    mpf_init(magifpi);
    mpf_init(mpsft);

    mpf_init(mpa1);
    mpf_init(mpb1);
    mpf_init(mpt1);
    mpf_init(mpp1);

    mpf_init(mpa2);
    mpf_init(mpb2);
    mpf_init(mpt2);
    mpf_init(mpp2);

    uint32_t usw_;
    usw_ = sw_(nnNonce2, SW_DIVS);
    if (usw_ < 1) usw_ = 1;
    mpz_set_ui(magisw, usw_);
    uint32_t mpzscale=mpz_size(magisw);

    for(int i=0; i < NM7M; i++){
        if (mpzscale > 1000) {
            mpzscale = 1000;
        }
        else if (mpzscale < 1) {
            mpzscale = 1;
        }

        mpf_set_ui(mpa1, 1);
        mpf_set_ui(mpb1, 2);
        mpf_set_d(mpt1, 0.25*mpzscale);
        mpf_set_ui(mpp1, 1);
        mpf_sqrt(mpb1, mpb1);
        mpf_ui_div(mpb1, 1, mpb1);
        mpf_set_ui(mpsft, 10);

        for(int i=0; i <= iterations; i++){
	        mpf_add(mpa2, mpa1,  mpb1);
	        mpf_div_ui(mpa2, mpa2, 2);
	        mpf_mul(mpb2, mpa1, mpb1);
	        mpf_abs(mpb2, mpb2);
	        mpf_sqrt(mpb2, mpb2);
	        mpf_sub(mpt2, mpa1, mpa2);
	        mpf_abs(mpt2, mpt2);
	        mpf_sqrt(mpt2, mpt2);
	        mpf_mul(mpt2, mpt2, mpp1);
	        mpf_sub(mpt2, mpt1, mpt2);
	        mpf_mul_ui(mpp2, mpp1, 2);
	        mpf_swap(mpa1, mpa2);
	        mpf_swap(mpb1, mpb2);
	        mpf_swap(mpt1, mpt2);
	        mpf_swap(mpp1, mpp2);
        }

        mpf_add(magifpi, mpa1, mpb1);
        mpf_pow_ui(magifpi, magifpi, 2);
        mpf_div_ui(magifpi, magifpi, 4);
        mpf_abs(mpt1, mpt1);
        mpf_div(magifpi, magifpi, mpt1);

        mpf_pow_ui(mpsft, mpsft, digits/2);
        mpf_mul(magifpi, magifpi, mpsft);

        mpz_set_f(magipi, magifpi);

        mpz_add(product,product,magipi);
        mpz_add(product,product,magisw);
    
        if(finalhash==0) finalhash = 1;
        mpz_set_uint256(bns[0],finalhash);
        mpz_add(bns[7], bns[7], bns[0]);

        mpz_mul(product,product,bns[7]);
        mpz_cdiv_q (product, product, bns[0]);
        if (mpz_sgn(product) <= 0) mpz_set_ui(product,1);

        bytes = mpz_sizeinbase(product, 256);
        mpzscale=bytes;
        char *bdata = (char*)malloc(bytes);
        mpz_export(bdata, NULL, -1, 1, 0, 0, product);

        sph_sha256_init(&ctx_sha256);
        sph_sha256 (&ctx_sha256, bdata, bytes);
        sph_sha256_close(&ctx_sha256, static_cast<void*>(&finalhash));
        free(bdata);
    }

    //Free the memory
    for(int i=0; i < 8; i++){
	    mpz_clear(bns[i]);
    }

    mpz_clear(product);

    mpz_clear(magipi);
    mpz_clear(magisw);
    mpf_clear(magifpi);
    mpf_clear(mpsft);

    mpf_clear(mpa1);
    mpf_clear(mpb1);
    mpf_clear(mpt1);
    mpf_clear(mpp1);

    mpf_clear(mpa2);
    mpf_clear(mpb2);
    mpf_clear(mpt2);
    mpf_clear(mpp2);
    
    return finalhash;
}

extern "C" void M7MHash_v2(char *finalhash, const unsigned char *input, int len, const unsigned int nnNonce)
{
    uint256 hash = hash_M7M_v2(input, input + len, nnNonce);
    memcpy(finalhash, &hash, 32);
}

