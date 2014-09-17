
#include <cstring>
#include <cstdlib>
#include <gmp.h>
#include "uint256.h"
#include "sph_sha2.h"
#include "sph_keccak.h"
#include "sph_haval.h"
#include "sph_tiger.h"
#include "sph_whirlpool.h"
#include "sph_ripemd.h"

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

#define NM7M 5
template<typename T1>
inline uint256 Hash7(const T1 pbegin, const T1 pend)
{
    sph_sha256_context       ctx_sha256;
    sph_sha512_context       ctx_sha512;
    sph_keccak512_context    ctx_keccak;
    sph_whirlpool_context    ctx_whirlpool;
    sph_haval256_5_context   ctx_haval;
    sph_tiger_context        ctx_tiger;
    sph_ripemd160_context    ctx_ripemd;
    static unsigned char pblank[1];

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

    //Take care of zeros and load gmp
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

    int bytes = mpz_sizeinbase(product, 256);
    char *data = (char*)malloc(bytes);
    mpz_export(data, NULL, -1, 1, 0, 0, product);

    sph_sha256_init(&ctx_sha256);
    sph_sha256 (&ctx_sha256, data, bytes);
    sph_sha256_close(&ctx_sha256, static_cast<void*>(&finalhash));
    free(data);

    for(int i=0; i < NM7M; i++){
        if(finalhash==0) finalhash = 1;
        mpz_set_uint256(bns[0],finalhash);
        mpz_add(bns[7], bns[7], bns[0]);

        mpz_mul(product,product,bns[7]);
        mpz_cdiv_q (product, product, bns[0]);
        if (mpz_sgn(product) <= 0) mpz_set_ui(product,1);

        bytes = mpz_sizeinbase(product, 256);
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

    return finalhash;
}

extern "C" void M7MHash(char *finalhash, const unsigned char *input, int len)
{
    uint256 hash = Hash7(input, input + len);
    memcpy(finalhash, &hash, 32);
}

// unsigned char *hash, const unsigned char *data, int len
extern "C" void m7hash(const char *finalhash, const unsigned char *input, int len)
{
    sph_sha256_context       ctx_sha256;
    sph_sha512_context       ctx_sha512;
    sph_keccak512_context    ctx_keccak;
    sph_whirlpool_context    ctx_whirlpool;
    sph_haval256_5_context   ctx_haval;
    sph_tiger_context        ctx_tiger;
    sph_ripemd160_context    ctx_ripemd;

    char hash[7][64];
    memset(hash, 0, 7 * 64);

    const void* ptr = input;
    size_t sz = len;

    sph_sha256_init(&ctx_sha256);
    // ZSHA256;
    sph_sha256 (&ctx_sha256, ptr, sz);
    sph_sha256_close(&ctx_sha256, (void*)(hash[0]));
    
    sph_sha512_init(&ctx_sha512);
    // ZSHA512;
    sph_sha512 (&ctx_sha512, ptr, sz);
    sph_sha512_close(&ctx_sha512, (void*)(hash[1]));
    
    sph_keccak512_init(&ctx_keccak);
    // ZKECCAK;
    sph_keccak512 (&ctx_keccak, ptr, sz);
    sph_keccak512_close(&ctx_keccak, (void*)(hash[2]));

    sph_whirlpool_init(&ctx_whirlpool);
    // ZWHIRLPOOL;
    sph_whirlpool (&ctx_whirlpool, ptr, sz);
    sph_whirlpool_close(&ctx_whirlpool, (void*)(hash[3]));
    
    sph_haval256_5_init(&ctx_haval);
    // ZHAVAL;
    sph_haval256_5 (&ctx_haval, ptr, sz);
    sph_haval256_5_close(&ctx_haval, (void*)(hash[4]));

    sph_tiger_init(&ctx_tiger);
    // ZTIGER;
    sph_tiger (&ctx_tiger, ptr, sz);
    sph_tiger_close(&ctx_tiger, (void*)(hash[5]));

#if 1
    sph_ripemd160_init(&ctx_ripemd);
    // ZRIPEMD;
    sph_ripemd160 (&ctx_ripemd, ptr, sz);
    sph_ripemd160_close(&ctx_ripemd, (void*)(hash[6]));

    //printf("%s\n", hash[6].GetHex().c_str());
#endif
    mpz_t bns[7];

    //Take care of zeros and load gmp
    for(int i=0; i < 7; i++){
        set_one_if_zero(hash[i]);
        mpz_init(bns[i]);
        mpz_set_uint512(bns[i],hash[i]);
    }
 
    mpz_t product;
    mpz_init(product);
    mpz_set_ui(product,1);
    for(int i=0; i < 7; i++){
        mpz_mul(product,product,bns[i]);
    }

    int bytes = mpz_sizeinbase(product, 256);
    char *data = (char*)malloc(bytes);
    mpz_export((void *)data, NULL, -1, 1, 0, 0, product);

    //Free the memory
    for(int i=0; i < 7; i++){
        mpz_clear(bns[i]);
    }
    mpz_clear(product);

    sph_sha256_init(&ctx_sha256);
    // ZSHA256;
    sph_sha256 (&ctx_sha256, data,bytes);
    sph_sha256_close(&ctx_sha256, (void*)(finalhash));

    free(data);
}
