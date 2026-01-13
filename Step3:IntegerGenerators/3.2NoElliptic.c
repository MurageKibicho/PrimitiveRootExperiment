#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <time.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fq.h>
#include <flint/fq_poly.h>
#include <flint/fq_poly_factor.h>
#include "ff_asm_primes.h"
#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"
#define PRIME_COUNT 5259

//clear && gcc 3.2NoElliptic.c -lm -lgmp -lmpfr -lflint -o m.o && ./m.o

bool FindGX(fmpz_t temp, fmpz_t numerator, fmpz_t denominator, fmpz_t generator, fmpz_t t0,fmpz_t t1,fmpz_t h0,fmpz_t h1, fmpz_t pointOrder)
{
	bool found = false;
	fmpz_powm(temp, generator, t0, pointOrder);
	fmpz_powm(numerator, generator, h0, pointOrder);
	fmpz_add(numerator, numerator, temp);
	fmpz_sub(numerator, pointOrder, numerator);
	
	fmpz_powm(temp, generator, t1, pointOrder);
	fmpz_powm(denominator, generator, h1, pointOrder);
	fmpz_add(denominator, denominator, temp);
	
	if(fmpz_cmp_ui(denominator, 0) != 0)
	{
		fmpz_invmod(denominator, denominator, pointOrder);
		fmpz_mul(numerator, numerator, denominator);
		fmpz_mod(numerator, numerator, pointOrder);
		found = true;
	}
	return found;	
}

void RandomWalkUpdate(fmpz_t exponentAnimal0, fmpz_t exponentAnimal1, fmpz_t pointOrderMinusOne)
{
	int modulo = 3;
	int moduloResult = rand() % modulo;
	if(moduloResult == 0)
	{
		fmpz_add_ui(exponentAnimal0, exponentAnimal0, 1);
		fmpz_mod(exponentAnimal0, exponentAnimal0, pointOrderMinusOne);
	}
	else if(moduloResult == 1)
	{
		fmpz_add_ui(exponentAnimal1, exponentAnimal1, 1);
		fmpz_mod(exponentAnimal1, exponentAnimal1, pointOrderMinusOne);
	}
	else
	{
		int count = 100;
		fmpz_add_ui(exponentAnimal1, exponentAnimal1, count);
		fmpz_mod(exponentAnimal1, exponentAnimal1, pointOrderMinusOne);	
	}
}

void NoEllipticTest()
{
	char *pointOrderHexadecimal = "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141";	
	int mersenneExponent = 7;
	int primeNumberHolder = 43573;
	int pointOrderHolder = 2293;
	assert(primeNumberHolder % 3 == 1);
	assert(primeNumberHolder < first5239[PRIME_COUNT-1]);
	
	fmpz_t pointOrder, pointOrderMinusOne, primeNumber,mersenneGenerator,lhs,rhs;
	fmpz_t t0,t1,h0,h1, temp, numerator, denominator;
	fmpz_init(pointOrderMinusOne);fmpz_init(temp);fmpz_init(numerator);fmpz_init(denominator);
	fmpz_init(lhs);fmpz_init(rhs);
	fmpz_init(t0);fmpz_init(t1);
	fmpz_init(h0);fmpz_init(h1);
	fmpz_init(pointOrder);
	fmpz_init(primeNumber);
	fmpz_init(mersenneGenerator);
	fmpz_set_ui(primeNumber, primeNumberHolder);
	
	fmpz_one(mersenneGenerator);
	fmpz_mul_2exp(mersenneGenerator, mersenneGenerator, mersenneExponent); 
	fmpz_sub_ui(mersenneGenerator, mersenneGenerator, 1);
	
	fmpz_set_ui(pointOrder, pointOrderHolder);
	//fmpz_set_str(pointOrder, pointOrderHexadecimal, 16);
	fmpz_sub_ui(pointOrderMinusOne, pointOrder, 1);
	
	fmpz_set_ui(t0, 2);
	fmpz_set_ui(t1, 714);
	fmpz_set_ui(h0, 11);
	fmpz_set_ui(h1, 1619);
	
	fmpz_set_ui(t0, 0);
	fmpz_set_ui(t1, 1);
	fmpz_set_ui(h0, 1);
	fmpz_set_ui(h1, 100);
	srand(43246156);
	int keyCount = 100000;
	clock_t start, end;
	double cpuTimeCount;
	printf("Benchmarking MySecp public key generation...\n");
	printf("Iterations: %d\n", keyCount);
	start = clock();
	
	for(int step = 0; step < 2293; step++)
	{
		bool foundGX = FindGX(temp, numerator, denominator, mersenneGenerator, t0, t1, h0, h1, pointOrder);
		if(foundGX)
		{
			fmpz_powm(temp, mersenneGenerator, t1, pointOrder);
			fmpz_add(lhs, lhs, temp);
			fmpz_mul(lhs, lhs, numerator);
			fmpz_powm(temp, mersenneGenerator, t0, pointOrder);
			fmpz_add(lhs, lhs, temp);
			fmpz_mod(lhs, lhs, pointOrder);
			
			fmpz_powm(temp, mersenneGenerator, h1, pointOrder);
			fmpz_add(rhs, rhs, temp);
			fmpz_mul(rhs, rhs, numerator);
			fmpz_powm(temp, mersenneGenerator, h0, pointOrder);
			fmpz_add(rhs, rhs, temp);
			fmpz_mod(rhs, rhs, pointOrder);
			
			fmpz_add(temp, lhs, rhs);
			fmpz_mod(denominator,temp, pointOrder);
			//fmpz_mod(temp, temp, pointOrder);
			bool foundCollision = false;
			if(fmpz_cmp_ui(denominator, 0) == 0)
			{
				foundCollision = true;
			}
			printf("(%d : %s) : g^x=",step, foundCollision ? "true" : "false");fmpz_print(numerator);printf(" ");
			printf("lhs = ");fmpz_print(lhs);printf(" ");
			printf("rhs = ");fmpz_print(rhs);printf(" ");
			printf("sum = ");fmpz_print(temp);printf("\n");
			
			RandomWalkUpdate(t0, t1, pointOrderMinusOne);
			RandomWalkUpdate(h0, h1, pointOrderMinusOne);
			RandomWalkUpdate(h0, h1, pointOrderMinusOne);
		}
	}
	
	end = clock(); 

	cpuTimeCount = ((double)(end - start)) / CLOCKS_PER_SEC;
	double iterationsPerSecond = keyCount / cpuTimeCount;

	printf("Total time MSB: %.4f seconds\n", cpuTimeCount);
	printf("Operations per second: %.2f\n", iterationsPerSecond);

	fmpz_clear(pointOrderMinusOne);fmpz_clear(lhs);fmpz_clear(rhs);
	fmpz_clear(temp);fmpz_clear(numerator);fmpz_clear(denominator);
	fmpz_clear(t0);fmpz_clear(t1);
	fmpz_clear(h0);fmpz_clear(h1);
	fmpz_clear(pointOrder);
	fmpz_clear(primeNumber);
	fmpz_clear(mersenneGenerator);
}

int main()
{
	NoEllipticTest();
	flint_cleanup();
	return 0;
}
