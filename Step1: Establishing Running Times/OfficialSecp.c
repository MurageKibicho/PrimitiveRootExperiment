#include <secp256k1.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
//Run: clear && gcc OfficialSecp.c -lm -lsecp256k1 -o m.o && ./m.o
// Benchmark is 100,000 private sequential keyss 
#define NUM_ITERATIONS 100000

void IntToPrivateKeyString(uint64_t n, unsigned char *privateKey)
{
	memset(privateKey, 0, 32);
	for(int i = 0; i < 8; i++)
	{
		privateKey[31 - i] = (n >> (8 * i)) & 0xFF;
	}
}

int main() 
{
	secp256k1_context *ctx;
	secp256k1_pubkey pubkey;
	clock_t start, end;
	double cpuTimeCount;
	unsigned char privateKey[32];
	int returnValue;

	// Create secp context
	ctx = secp256k1_context_create(SECP256K1_CONTEXT_SIGN | SECP256K1_CONTEXT_VERIFY);

	printf("Benchmarking Libsecp256k1 public key generation...\n");
	printf("Iterations: %d\n", NUM_ITERATIONS);
	//Warm up
	for(uint64_t i = 1; i <= 10000; i++) 
	{
		IntToPrivateKeyString(i, privateKey);
		returnValue = secp256k1_ec_pubkey_create(ctx, &pubkey, privateKey);
	}

	start = clock();

	for(uint64_t i = 1; i <= NUM_ITERATIONS; i++)
	{
		IntToPrivateKeyString(i, privateKey);
		returnValue = secp256k1_ec_pubkey_create(ctx, &pubkey, privateKey);
	}

	end = clock(); 

	cpuTimeCount = ((double)(end - start)) / CLOCKS_PER_SEC;
	double iterationsPerSecond = NUM_ITERATIONS / cpuTimeCount;

	printf("Total time: %.4f seconds\n", cpuTimeCount);
	printf("Operations per second: %.2f\n", iterationsPerSecond);
	printf("Operations per minute: %.2f\n", iterationsPerSecond * 60);
	printf("Operations per hour: %.2f\n", iterationsPerSecond * 3600);
	printf("Operations per day: %.2f\n", iterationsPerSecond * 86400);

	secp256k1_context_destroy(ctx);

	return 0;
}
