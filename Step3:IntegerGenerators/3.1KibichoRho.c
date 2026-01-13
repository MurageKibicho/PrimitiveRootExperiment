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

//clear && gcc 3.1KibichoRho.c -lm -lgmp -lmpfr -lflint -o m.o && ./m.o

typedef struct elliptic_point_struct *Elliptic;
struct elliptic_point_struct
{
	fmpz_t x;
	fmpz_t y;
	int infinity;
};

Elliptic CreatePoint()
{
	Elliptic point = malloc(sizeof(struct elliptic_point_struct));
	fmpz_init(point->x);
	fmpz_init(point->y);
	point->infinity = 0;
	return point;
}

void CopyPoint(Elliptic source, Elliptic destination)
{
	fmpz_set(destination->x, source->x);
	fmpz_set(destination->y, source->y);
	destination->infinity = source->infinity;
}

void FindPointInverse(Elliptic source, Elliptic destination, fmpz_t primeNumber)
{
	fmpz_set(destination->x, source->x);
	fmpz_sub(destination->y, primeNumber,source->y);
	fmpz_mod(destination->y,destination->y,primeNumber);
}

int TestPointEquality(Elliptic a, Elliptic b)
{
	if(a->infinity && b->infinity) return 1;
	if(a->infinity || b->infinity) return 0;
	return fmpz_cmp(a->x, b->x) == 0 && fmpz_cmp(a->y, b->y) == 0;	
}

void PrintPoint(Elliptic point)
{
	fmpz_print(point->x);printf(" ");fmpz_print(point->y);printf("\n");
}

void PrintPointSpace(Elliptic point)
{
	printf("(");fmpz_print(point->x);printf(" ");fmpz_print(point->y);printf(") ");
}

void DestroyPoint(Elliptic point)
{
	if(point)
	{
		fmpz_clear(point->x);
		fmpz_clear(point->y);
		free(point);
	}
}

bool AddCurvePoints(Elliptic R, Elliptic P, Elliptic Q, fmpz_t primeNumber)
{
	//Case 0: Handle Points at infinity
	if(P->infinity != 0){CopyPoint(Q,R);return true;}
	if(Q->infinity != 0){CopyPoint(P,R);return true;}
	
	//Case 1: Handle P->x == Q->x
	if(fmpz_cmp(P->x, Q->x) == 0) 
	{
		//Case 1.1: X values similar but Y differ or 0
		if(fmpz_cmp(P->y, Q->y) != 0 || fmpz_cmp_ui(P->y, 0) == 0){R->infinity = 1;return true;}
		
		//Case 1.2: Point Doubling
		fmpz_t s, num, den, denominatorInverse, tmp;
		fmpz_init(s);
		fmpz_init(num);
		fmpz_init(den);
		fmpz_init(denominatorInverse);
		fmpz_init(tmp);

		//num = x^2
		fmpz_mul(num, P->x, P->x);
		//num = 3x^2   
		fmpz_mul_ui(num, num, 3);       

		//den = 2y
		fmpz_mul_ui(den, P->y, 2); 
		//Find denominator inverse     
		if(!fmpz_invmod(denominatorInverse, den, primeNumber))
		{
			R->infinity = 1;
			fmpz_clear(s);
			fmpz_clear(num);
			fmpz_clear(den);
			fmpz_clear(denominatorInverse);
			fmpz_clear(tmp);
			return true;
		}

		//s = (3x^2)/(2y)
		fmpz_mul(s, num, denominatorInverse);     
		fmpz_mod(s, s, primeNumber);

		//x3 = s^2 - 2x
		fmpz_mul(tmp, s, s);
		fmpz_sub(tmp, tmp, P->x);
		fmpz_sub(tmp, tmp, Q->x);
		fmpz_mod(R->x, tmp, primeNumber);

		//y3 = s*(x - x3) - y
		fmpz_sub(tmp, P->x, R->x);
		fmpz_mul(tmp, s, tmp);
		fmpz_sub(tmp, tmp, P->y);
		fmpz_mod(R->y, tmp, primeNumber);

		R->infinity = 0;

		fmpz_clear(s);
		fmpz_clear(num);
		fmpz_clear(den);
		fmpz_clear(denominatorInverse);
		fmpz_clear(tmp);
		return true;

	}
	//Case 2: Handle P != Q (Point Addition)
	else
	{
		fmpz_t s, num, den, denominatorInverse, tmp;
		fmpz_init(s);
		fmpz_init(num);
		fmpz_init(den);
		fmpz_init(denominatorInverse);
		fmpz_init(tmp);
		
		//num = y2 - y1
		fmpz_sub(num, Q->y, P->y);
		//den = x2 - x1   
		fmpz_sub(den, Q->x, P->x); 
		//Find inverse of denominator mod p    
		if(!fmpz_invmod(denominatorInverse, den, primeNumber))
		{
			R->infinity = 1;
			fmpz_clear(s);
			fmpz_clear(num);
			fmpz_clear(den);
			fmpz_clear(denominatorInverse);
			fmpz_clear(tmp);
			return true;
		}

		//s = (y2 - y1)/(x2 - x1) mod p
		fmpz_mul(s, num, denominatorInverse);fmpz_mod(s, s, primeNumber);

		//x3 = s^2 - x1 - x2 mod p
		fmpz_mul(tmp, s, s);fmpz_sub(tmp, tmp, P->x);fmpz_sub(tmp, tmp, Q->x);fmpz_mod(R->x, tmp, primeNumber);

		//y3 = s*(x1 - x3) - y1
		fmpz_sub(tmp, P->x, R->x);fmpz_mul(tmp, s, tmp);fmpz_sub(tmp, tmp, P->y);fmpz_mod(R->y, tmp, primeNumber);
		R->infinity = 0;
		//Free memory
		fmpz_clear(s);
		fmpz_clear(num);
		fmpz_clear(den);
		fmpz_clear(denominatorInverse);
		fmpz_clear(tmp);
		return true;
	}
}


void MSBScalarMultiplication(Elliptic resultant, Elliptic generator, fmpz_t privateKey, fmpz_t primeNumber)
{
	//Create pointAtInfinity, temp0
	Elliptic pointAtInfinity = CreatePoint();
	pointAtInfinity->infinity = 1;  
	
	Elliptic temp0 = CreatePoint();
	size_t binaryLength = fmpz_sizeinbase(privateKey, 2);
	
	//loop from MSB to LSB
	for(ssize_t i = binaryLength - 1; i >= 0; --i)
	{
		//temp0 = 2*pointAtInfinity
		AddCurvePoints(temp0, pointAtInfinity, pointAtInfinity, primeNumber);
		CopyPoint(temp0, pointAtInfinity);
		
		//Test the current bit's parity
		if(fmpz_tstbit(privateKey, i) != 0)
		{
			//temp0 = pointAtInfinity + generator
			AddCurvePoints(temp0, pointAtInfinity, generator, primeNumber);
			CopyPoint(temp0, pointAtInfinity);	
		}
	}
	//Save pointAtInfinity to the result variable
	CopyPoint(pointAtInfinity, resultant);

	//Free memory
	DestroyPoint(pointAtInfinity);
	DestroyPoint(temp0);
}

void LSBScalarMultiplication(Elliptic resultant, Elliptic generator, fmpz_t privateKey, fmpz_t primeNumber)
{
	//Create pointAtInfinity, temp0, temp1
	Elliptic pointAtInfinity = CreatePoint();
	pointAtInfinity->infinity = 1;  
	
	Elliptic temp0 = CreatePoint();
	//Save generator to temp0
	CopyPoint(generator, temp0); 
	Elliptic temp1 = CreatePoint();

	size_t binaryLength = fmpz_sizeinbase(privateKey, 2);
	//loop from LSB to MSB
	for(size_t i = 0; i < binaryLength; ++i)
	{
		int bit = fmpz_tstbit(privateKey, i);
		if(bit)
		{
			//temp1 = pointAtInfinity + currentGenerator
			AddCurvePoints(temp1, pointAtInfinity, temp0, primeNumber);
			CopyPoint(temp1, pointAtInfinity);
		}
		//Update temp0(currentGenerator) = 2*temp0(previousGenerator)
		AddCurvePoints(temp1, temp0, temp0, primeNumber);
		CopyPoint(temp1, temp0);
	}
	//Save pointAtInfinity to the result variable
	CopyPoint(pointAtInfinity, resultant);
	DestroyPoint(pointAtInfinity);
	DestroyPoint(temp0);
	DestroyPoint(temp1);
}

void BenchMarkMSB()
{
	int keyCount = 1000;
	char *primeNumberHexadecimal = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";	
	char *generatorXHexadecimal  = "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798";
	char *generatorYHexadecimal  = "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8";
	
	fmpz_t primeNumber, privateKey;
	fmpz_init(primeNumber);
	fmpz_init(privateKey);
	Elliptic generator = CreatePoint();
	Elliptic result = CreatePoint();
	fmpz_set_str(primeNumber, primeNumberHexadecimal, 16);
	fmpz_set_str(generator->x, generatorXHexadecimal, 16);
	fmpz_set_str(generator->y, generatorYHexadecimal, 16);
	
	clock_t start, end;
	double cpuTimeCount;
	printf("Benchmarking MySecp public key generation...\n");
	printf("Iterations: %d\n", keyCount);
	start = clock();
	for(int i = 1; i < keyCount; i++)
	{
		fmpz_set_ui(privateKey, i);
		MSBScalarMultiplication(result, generator, privateKey, primeNumber);	
		//PrintPoint(result);
	}
	end = clock(); 

	cpuTimeCount = ((double)(end - start)) / CLOCKS_PER_SEC;
	double iterationsPerSecond = keyCount / cpuTimeCount;

	printf("Total time MSB: %.4f seconds\n", cpuTimeCount);
	
	DestroyPoint(generator);
	DestroyPoint(result);
	fmpz_clear(primeNumber);
	fmpz_clear(privateKey);
}

Elliptic *CacheDoubles(Elliptic generator, fmpz_t primeNumber, int bitCount)
{
	Elliptic *doubles = malloc(bitCount * sizeof(Elliptic));
	assert(doubles != NULL);

	Elliptic temp = CreatePoint();
	Elliptic temp2 = CreatePoint();
	CopyPoint(generator, temp);

	for(int i = 0; i < bitCount; i++)
	{
		doubles[i] = CreatePoint();
		//Store current power of 2 * G
		CopyPoint(temp, doubles[i]); 

		//Compute next power of 2: temp = 2 * temp
	
		AddCurvePoints(temp2, temp, temp, primeNumber);
		CopyPoint(temp2, temp);
	}
	DestroyPoint(temp2);
	DestroyPoint(temp);
	return doubles;
}

void LSBCachedMultiplication(int bitCount, Elliptic *cachedDoubles, Elliptic resultant, Elliptic generator, fmpz_t privateKey, fmpz_t primeNumber)
{
	//Create pointAtInfinity, temp0, temp1
	Elliptic pointAtInfinity = CreatePoint();
	pointAtInfinity->infinity = 1;  
	
	Elliptic temp0 = CreatePoint();
	//Save generator to temp0
	CopyPoint(generator, temp0); 
	Elliptic temp1 = CreatePoint();

	size_t binaryLength = fmpz_sizeinbase(privateKey, 2);
	//loop from LSB to MSB
	for(size_t i = 0; i < binaryLength; ++i)
	{
		int bit = fmpz_tstbit(privateKey, i);
		if(bit)
		{
			//temp1 = pointAtInfinity + currentGenerator
			AddCurvePoints(temp1, pointAtInfinity, cachedDoubles[i], primeNumber);
			CopyPoint(temp1, pointAtInfinity);
		}
		
	}
	//Save pointAtInfinity to the result variable
	CopyPoint(pointAtInfinity, resultant);
	DestroyPoint(pointAtInfinity);
	DestroyPoint(temp0);
	DestroyPoint(temp1);
}

bool FindRandomSecpPoint(Elliptic result, fmpz_t primeNumber, flint_rand_t randomState)
{
	fmpz_t x, y, rhs, lhs, tmp;
	fmpz_init(x);
	fmpz_init(y);
	fmpz_init(rhs);
	fmpz_init(lhs);
	fmpz_init(tmp);

	bool found = false;
	int attempts = 0;
	int maxAttempts = 1000;

	while(!found && attempts < maxAttempts)
	{
		attempts++;
		//Generate random x coordinate
		fmpz_randm(x, randomState, primeNumber);

		// Compute right side: x^3 + 7 mod primeNumber
		fmpz_powm_ui(rhs, x, 3, primeNumber);      // rhs = x^3 mod p
		fmpz_add_ui(rhs, rhs, 7);        // rhs = x^3 + 7
		fmpz_mod(rhs, rhs, primeNumber);           // rhs = (x^3 + 7) mod p

		// Check if rhs is a quadratic residue mod p
		int legendre = fmpz_jacobi(rhs, primeNumber);

		if(legendre == 1) 
		{
			if(fmpz_sqrtmod(y, rhs, primeNumber))
			{
				fmpz_set(result->x, x);
				fmpz_set(result->y, y);
				found = true;
			}
		}
	}
	fmpz_clear(x);
	fmpz_clear(y);
	fmpz_clear(rhs);
	fmpz_clear(lhs);
	fmpz_clear(tmp);

	return found;
}

bool BruteforcePointOrder(Elliptic generator, fmpz_t primeNumber, fmpz_t pointOrder)
{
	int bitCount = 32;
	bool foundPointOrder = false;
	assert(generator->infinity != 1);
	assert(fmpz_sizeinbase(primeNumber, 2) < bitCount);
	Elliptic *cachedDoubles = CacheDoubles(generator, primeNumber, bitCount);
	fmpz_t privateKey;
	Elliptic result0 = CreatePoint();
	fmpz_init(privateKey);
	for(int i = 1; i < fmpz_get_ui(primeNumber); i++)
	{
		fmpz_set_ui(privateKey, i);
		LSBCachedMultiplication(bitCount, cachedDoubles, result0, generator, privateKey, primeNumber);
		if(result0->infinity == 1)
		{
			fmpz_set_ui(pointOrder, i);
			foundPointOrder = true;
			break;
		}
	}
	for(int i = 0; i < bitCount; i++){DestroyPoint(cachedDoubles[i]);}free(cachedDoubles);
	fmpz_clear(privateKey);
	DestroyPoint(result0);
	return foundPointOrder;
}

bool RaiseGeneratorTest(fmpz_t testPrime, fmpz_t pointOrder, fmpz_t pointOrderMinusOne, fmpz_factor_t factorization)
{
	bool result = false;
	fmpz_t testResult, exponent;
	fmpz_init(testResult);
	fmpz_init(exponent);
	for(slong i = 0; i < factorization->num; i++)
	{
		fmpz_mod(exponent, pointOrderMinusOne, &factorization->p[i]);
		assert(fmpz_cmp_ui(exponent, 0 ) == 0);
		fmpz_divexact(exponent, pointOrderMinusOne, &factorization->p[i]);
		fmpz_powm(testResult, testPrime, exponent, pointOrder);
		if(fmpz_cmp_ui(testResult, 1) == 0)
		{
			result = false;
			return result;
		}
	}
	
	fmpz_clear(exponent);
	fmpz_clear(testResult);
	result = true;
	return result;
}

void BenchMarkLSB()
{
	int keyCount = 100000;
	int bitCount = 20;
	char *primeNumberHexadecimal = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";	
	char *generatorXHexadecimal  = "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798";
	char *generatorYHexadecimal  = "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8";
	
	fmpz_t primeNumber, privateKey;
	fmpz_init(primeNumber);
	fmpz_init(privateKey);
	Elliptic generator = CreatePoint();
	Elliptic result0 = CreatePoint();
	Elliptic result1 = CreatePoint();
	fmpz_set_str(primeNumber, primeNumberHexadecimal, 16);
	fmpz_set_str(generator->x, generatorXHexadecimal, 16);
	fmpz_set_str(generator->y, generatorYHexadecimal, 16);
	
	Elliptic *cachedDoubles = CacheDoubles(generator, primeNumber, bitCount);
	clock_t start, end;
	double cpuTimeCount;
	printf("Benchmarking MySecp public key generation...\n");
	printf("Iterations: %d\n", keyCount);
	start = clock();
	for(int i = 1; i < keyCount; i++)
	{
		fmpz_set_ui(privateKey, i);
		LSBCachedMultiplication(bitCount, cachedDoubles, result1, generator, privateKey, primeNumber);
	}
	end = clock(); 

	cpuTimeCount = ((double)(end - start)) / CLOCKS_PER_SEC;
	double iterationsPerSecond = keyCount / cpuTimeCount;

	printf("Total time MSB: %.4f seconds\n", cpuTimeCount);
	
	DestroyPoint(generator);
	DestroyPoint(result0);
	DestroyPoint(result1);
	fmpz_clear(primeNumber);
	fmpz_clear(privateKey);
	for(int i = 0; i < bitCount; i++){DestroyPoint(cachedDoubles[i]);}free(cachedDoubles);
}


void FindPuzzleGenerator()
{
	int startUpperBound = 71;
	int mersenneExponents [] = {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933};
	int mersenneLength = sizeof(mersenneExponents) / sizeof(int);
	fmpz_factor_t factorization;
	fmpz_factor_init(factorization);
	fmpz_t pointOrder,pointOrderMinusOne, largestFactor, testPrime;
	fmpz_init(pointOrder);
	fmpz_init(pointOrderMinusOne);
	fmpz_init(largestFactor);
	fmpz_init(testPrime);
	for(int i = startUpperBound; i <= 160; i++)
	{
		//Avoid already solveD or bad minus one
		if(i == 75 || i == 80 ||  i == 83 || i == 85 || i == 88 || i == 89 || i == 90 || i == 93 || i == 95 || i == 97 || i == 98 || i == 100 || i == 101 || i == 103 || i == 105 || i == 107 || i == 109 || i == 115 || i == 118 || i == 119 || i == 120 || i == 121 || i == 122 || i == 123||  i == 125 || i == 127 || i == 129 || i == 130 || i == 131 || i == 133 || i == 134 || i == 137 || i == 138 || i == 139 || i == 141 || i == 143 || i == 145 || i == 146 || i == 147 || i == 148 || i == 149 || i == 151 || i == 153 || i == 154 || i == 155 || i == 157 || i == 158 || i == 159 || i == 160)
		{
			//https://btcpuzzle.info/puzzle
			//printf("Bad Factor Range (2^%3d to 2^%d)\n", i - 1, i);	
		}
		else
		{
			fmpz_one(pointOrder);
			fmpz_mul_2exp(pointOrder, pointOrder, i); 
			fmpz_sub_ui(pointOrderMinusOne, pointOrder, 1);
			fmpz_factor(factorization, pointOrderMinusOne);
			printf("Range (2^%3d to 2^%d)\n", i - 1, i);
			printf("Prime Factors of p-1:(%ld): ",factorization->num);
			if(factorization->num > 0)
			{
				fmpz_set(largestFactor, &factorization->p[factorization->num - 1]);
				for(slong i = 0; i < factorization->num; i++)
				{
					printf("(");fmpz_print(&factorization->p[i]);flint_printf("^%wd),", factorization->exp[i]);		
					if(fmpz_cmp(largestFactor, &factorization->p[i]) < 0)
					{
							fmpz_set(largestFactor, &factorization->p[i]);
					}
				}
				size_t factorSize = fmpz_sizeinbase(largestFactor, 2);
				printf("\nLargest factor : ");fmpz_print(largestFactor);
				printf("(%ld bits)\n", factorSize);	
				if(i == 135 || i == 140 || i == 145 || i == 150 || i == 155 || i == 160)
				{
					printf("Public Key Known\n");	
	 			}
	 			int mersenneFoundCount = 0;
	 			int *foundMersenne = calloc(mersenneLength, sizeof(int));
				for(int j = 0; j < mersenneLength; j++)
				{
					fmpz_one(testPrime);
					fmpz_mul_2exp(testPrime, testPrime, mersenneExponents[j]); 
					fmpz_sub_ui(testPrime, testPrime, 1);
					if(fmpz_cmp(testPrime, pointOrder) > 0){break;}
					bool testResult = RaiseGeneratorTest(testPrime, pointOrder, pointOrderMinusOne, factorization);
					foundMersenne[i] = 0;
					if(testResult == true)
					{
						foundMersenne[i] = 1;
						mersenneFoundCount += 1;
					}
					printf("(2^%d) - 1: ",mersenneExponents[j]);fmpz_print(testPrime);printf(" %s\n", testResult ? "True" : "False");
		}
				
				//We focus on mersenne primes since more than fermat
				if(mersenneFoundCount > 0)
				{
					printf("\nFound %d mersenne primes\n",mersenneFoundCount);
				}
				else
				{
					printf("\nNo good generator\n");
				}
				free(foundMersenne);
			}
			else
			{
				printf("No factors\n");
			}
			
			printf("\n");
		}
		
	}
	fmpz_factor_clear(factorization);
	fmpz_clear(largestFactor);
	fmpz_clear(pointOrder);
	fmpz_clear(pointOrderMinusOne);
	fmpz_clear(testPrime);
}

void PollardRho_EllipticUpdate(fmpz_t setElement, fmpz_t exponentAnimal0, fmpz_t exponentAnimal1, Elliptic temp0, Elliptic temp1, Elliptic desiredGenerator, Elliptic desiredResult, fmpz_t pointOrderMinusOne, fmpz_t pointOrder, fmpz_t primeNumber, fmpz_t temporary, int mersenneExponent)
{
	
	int modulo = 2;
	int moduloResult = fmpz_mod_ui(temporary, setElement, modulo);
	if(moduloResult == 0)
	{
		//Update Coordinates by recurrence
		CopyPoint(desiredGenerator, temp0);
		for(int p = 0; p < mersenneExponent; p++)
		{
			AddCurvePoints(temp1, temp0, temp0, primeNumber);
			CopyPoint(temp1, temp0);
		}
		//Find Previous point's inverse
		FindPointInverse(desiredGenerator, temp1, primeNumber);
		//Add doubled point to inverse
		AddCurvePoints(desiredGenerator, temp0, temp1, primeNumber);
		
		//Update exponentAnimal0 by 1 mod primeNumberMinOne
		fmpz_mul(setElement, setElement, desiredGenerator->x);
		fmpz_mod(setElement, setElement, primeNumber);
		fmpz_add_ui(exponentAnimal0, exponentAnimal0, 1);
		fmpz_mod(exponentAnimal0, exponentAnimal0, pointOrderMinusOne);
	}
	else if(moduloResult == 1)
	{
		//Update Coordinates by recurrence
		CopyPoint(desiredResult, temp0);
		for(int p = 0; p < mersenneExponent; p++)
		{
			AddCurvePoints(temp1, temp0, temp0, primeNumber);
			CopyPoint(temp1, temp0);
		}
		//Find Previous point's inverse
		FindPointInverse(desiredResult, temp1, primeNumber);
		//Add doubled point to inverse
		AddCurvePoints(desiredResult, temp0, temp1, primeNumber);
		//Update exponentAnimal1 by 1 mod primeNumberMinOne
		fmpz_mul(setElement, setElement, desiredResult->x);
		fmpz_mod(setElement, setElement, primeNumber);
		fmpz_add_ui(exponentAnimal1, exponentAnimal1, 1);
		fmpz_mod(exponentAnimal1, exponentAnimal1, pointOrderMinusOne);
	}
	else
	{
	
	}
	//Update set element

}

void PrintAnimalStep(int i, char *animalString, fmpz_t mersenneGenerator, fmpz_t exponent0, fmpz_t exponent1, fmpz_t setElement,Elliptic generator, Elliptic target)
{
	printf("%s: (", animalString);fmpz_print(mersenneGenerator);printf("^");fmpz_print(exponent0);printf(" : ");PrintPointSpace(generator);
	printf(" * ");
	fmpz_print(mersenneGenerator);printf("^(x+");fmpz_print(exponent1);printf(") : ");PrintPointSpace(target);
	printf(" = ");fmpz_print(setElement);
	printf(")\n");	
}
bool PollardRho(Elliptic generator, Elliptic target, fmpz_t pointOrder, fmpz_t primeNumber, fmpz_t mersenneGenerator, int mersenneExponent)
{
	bool result = false;
	fmpz_t pointOrderMinusOne,temporary;
	fmpz_t exponentTortoise0, exponentTortoise1,setElementTortoise;
	fmpz_t exponentHare0, exponentHare1, setElementHare;
	Elliptic temp0 = CreatePoint();
	Elliptic temp1 = CreatePoint();
	Elliptic tempGeneratorTortoise = CreatePoint();
	Elliptic tempTargetTortoise  = CreatePoint();
	Elliptic tempGeneratorHare = CreatePoint();
	Elliptic tempTargetHare = CreatePoint();
	CopyPoint(generator, tempGeneratorTortoise);
	CopyPoint(target, tempTargetTortoise);
	CopyPoint(generator, tempGeneratorHare);
	CopyPoint(target, tempTargetHare);

	fmpz_init(pointOrderMinusOne);fmpz_init(temporary);fmpz_init(exponentTortoise0);fmpz_init(exponentTortoise1);fmpz_init(setElementTortoise);fmpz_init(exponentHare0);fmpz_init(exponentHare1);fmpz_init(setElementHare);
	
	fmpz_sub_ui(pointOrderMinusOne, pointOrder, 1);
	//Initialize exponents to 0
	fmpz_set_ui(exponentTortoise0, 0);fmpz_set_ui(exponentHare0, 0);
	fmpz_set_ui(exponentTortoise1, 0);fmpz_set_ui(exponentHare1, 0);
	
	//Initialize our current set element to 1
	fmpz_set_ui(setElementTortoise, 1);fmpz_set_ui(setElementHare, 1);
	
	for(int step = 0; step < 100; step++)
	{
		//Update tortoise once
		PollardRho_EllipticUpdate(setElementTortoise, exponentTortoise0, exponentTortoise1, temp0, temp1, tempGeneratorTortoise, tempTargetTortoise, pointOrderMinusOne, pointOrder, primeNumber, temporary,mersenneExponent);
		//Update hare twice
		PollardRho_EllipticUpdate(setElementHare, exponentHare0, exponentHare1, temp0,temp1, tempGeneratorHare, tempTargetHare, pointOrderMinusOne, pointOrder, primeNumber, temporary,mersenneExponent);
		PollardRho_EllipticUpdate(setElementHare, exponentHare0, exponentHare1, temp0,temp1, tempGeneratorHare, tempTargetHare, pointOrderMinusOne, pointOrder, primeNumber, temporary,mersenneExponent);

		printf("%d\n",step);
		PrintAnimalStep(step, "Tortoise", mersenneGenerator, exponentTortoise0, exponentTortoise1, setElementTortoise,tempGeneratorTortoise,tempTargetTortoise);
		PrintAnimalStep(step, "Hare    ", mersenneGenerator, exponentHare0, exponentHare1, setElementHare,tempGeneratorHare,tempTargetHare);
		if(fmpz_cmp(setElementTortoise, setElementHare) == 0)
		{
			printf("Found Collision\n");
			break;
		}

	}
	
	DestroyPoint(temp0);
	DestroyPoint(temp1);
	DestroyPoint(tempGeneratorTortoise);
	DestroyPoint(tempTargetTortoise);
	DestroyPoint(tempGeneratorHare);
	DestroyPoint(tempTargetHare);
	fmpz_clear(pointOrderMinusOne);fmpz_clear(temporary);fmpz_clear(exponentTortoise0);fmpz_clear(exponentTortoise1);fmpz_clear(setElementTortoise);fmpz_clear(exponentHare0);fmpz_clear(exponentHare1);fmpz_clear(setElementHare);
	return result;
}

void PrintAnimalStep2(int i, char *animalString, fmpz_t mersenneGenerator, fmpz_t exponent0, fmpz_t exponent1, Elliptic setElement,Elliptic generator, Elliptic target)
{
	printf("%s: (", animalString);fmpz_print(mersenneGenerator);printf("^");fmpz_print(exponent0);printf(" : ");PrintPointSpace(generator);
	printf(" + ");
	fmpz_print(mersenneGenerator);printf("^(x+");fmpz_print(exponent1);printf(") : ");PrintPointSpace(target);
	printf(" = ");PrintPointSpace(setElement);
	printf(")\n");	
}
void PollardRho_EllipticUpdate2(Elliptic setElement, fmpz_t exponentAnimal0, fmpz_t exponentAnimal1, Elliptic temp0, Elliptic temp1, Elliptic desiredGenerator, Elliptic desiredResult, fmpz_t pointOrderMinusOne, fmpz_t pointOrder, fmpz_t primeNumber, fmpz_t temporary, int mersenneExponent)
{
	
	int modulo = 3;
	int moduloResult = rand() % 3;
	if(moduloResult == 0)
	{
		//Update Coordinates by recurrence
		CopyPoint(desiredGenerator, temp0);
		for(int p = 0; p < mersenneExponent; p++)
		{
			AddCurvePoints(temp1, temp0, temp0, primeNumber);
			CopyPoint(temp1, temp0);
		}
		//Find Previous point's inverse
		FindPointInverse(desiredGenerator, temp1, primeNumber);
		//Add doubled point to inverse
		AddCurvePoints(desiredGenerator, temp0, temp1, primeNumber);
		
		fmpz_add_ui(exponentAnimal0, exponentAnimal0, 1);
		fmpz_mod(exponentAnimal0, exponentAnimal0, pointOrderMinusOne);
	}
	else if(moduloResult == 1)
	{
		//Update Coordinates by recurrence
		CopyPoint(desiredResult, temp0);
		for(int p = 0; p < mersenneExponent; p++)
		{
			AddCurvePoints(temp1, temp0, temp0, primeNumber);
			CopyPoint(temp1, temp0);
		}
		//Find Previous point's inverse
		FindPointInverse(desiredResult, temp1, primeNumber);
		//Add doubled point to inverse
		AddCurvePoints(desiredResult, temp0, temp1, primeNumber);
		
		fmpz_add_ui(exponentAnimal1, exponentAnimal1, 1);
		fmpz_mod(exponentAnimal1, exponentAnimal1, pointOrderMinusOne);
	}
	else
	{
		int count = 100;
		for(int j = 0; j < count; j++)
		{
			//Update Coordinates by recurrence
			CopyPoint(desiredResult, temp0);
			for(int p = 0; p < mersenneExponent; p++)
			{
				AddCurvePoints(temp1, temp0, temp0, primeNumber);
				CopyPoint(temp1, temp0);
			}
			//Find Previous point's inverse
			FindPointInverse(desiredResult, temp1, primeNumber);
			//Add doubled point to inverse
			AddCurvePoints(desiredResult, temp0, temp1, primeNumber);
		}
		fmpz_add_ui(exponentAnimal1, exponentAnimal1, count);
		fmpz_mod(exponentAnimal1, exponentAnimal1, pointOrderMinusOne);
	}
	//Update set element
	AddCurvePoints(setElement, desiredGenerator, desiredResult, primeNumber);
}

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

bool PollardRho2(Elliptic generator, Elliptic target, fmpz_t pointOrder, fmpz_t primeNumber, fmpz_t mersenneGenerator, int mersenneExponent)
{
	bool result = false;
	fmpz_t pointOrderMinusOne,temporary,numerator, denominator;
	fmpz_t exponentTortoise0, exponentTortoise1;
	fmpz_t exponentHare0, exponentHare1;
	Elliptic temp0 = CreatePoint();
	Elliptic temp1 = CreatePoint();
	Elliptic tempGeneratorTortoise = CreatePoint();
	Elliptic tempTargetTortoise  = CreatePoint();
	Elliptic tempGeneratorHare = CreatePoint();
	Elliptic tempTargetHare = CreatePoint();
	Elliptic setElementHare = CreatePoint();
	Elliptic setElementTortoise = CreatePoint();
	CopyPoint(generator, tempGeneratorTortoise);
	CopyPoint(target, tempTargetTortoise);
	CopyPoint(generator, tempGeneratorHare);
	CopyPoint(target, tempTargetHare);

	fmpz_init(numerator);fmpz_init(denominator);fmpz_init(pointOrderMinusOne);fmpz_init(temporary);fmpz_init(exponentTortoise0);fmpz_init(exponentTortoise1);fmpz_init(exponentHare0);fmpz_init(exponentHare1);
	
	fmpz_sub_ui(pointOrderMinusOne, pointOrder, 1);
	//Initialize exponents to 0
	fmpz_set_ui(exponentTortoise0, 0);fmpz_set_ui(exponentHare0, 0);
	fmpz_set_ui(exponentTortoise1, 0);fmpz_set_ui(exponentHare1, 0);
	
	//Initialize our current set element to 1
	fmpz_set_ui(setElementTortoise->x, 1);fmpz_set_ui(setElementHare->x, 1);
	
	for(int step = 0; step < 10000; step++)
	{
		//Update tortoise once
		PollardRho_EllipticUpdate2(setElementTortoise, exponentTortoise0, exponentTortoise1, temp0, temp1, tempGeneratorTortoise, tempTargetTortoise, pointOrderMinusOne, pointOrder, primeNumber, temporary,mersenneExponent);
		//Update hare twice
		PollardRho_EllipticUpdate2(setElementHare, exponentHare0, exponentHare1, temp0,temp1, tempGeneratorHare, tempTargetHare, pointOrderMinusOne, pointOrder, primeNumber, temporary,mersenneExponent);
		PollardRho_EllipticUpdate2(setElementHare, exponentHare0, exponentHare1, temp0,temp1, tempGeneratorHare, tempTargetHare, pointOrderMinusOne, pointOrder, primeNumber, temporary,mersenneExponent);


		bool found = FindGX(temporary, numerator, denominator, mersenneGenerator, exponentTortoise0, exponentTortoise1,exponentHare0, exponentHare1, pointOrder);
		printf("%d (%s) g^x = ",step,found ? "true" : "false");fmpz_print(numerator);printf("\n");
		PrintAnimalStep2(step, "Tortoise", mersenneGenerator, exponentTortoise0, exponentTortoise1, setElementTortoise,tempGeneratorTortoise,tempTargetTortoise);
		PrintAnimalStep2(step, "Hare    ", mersenneGenerator, exponentHare0, exponentHare1, setElementHare,tempGeneratorHare,tempTargetHare);
		if(fmpz_cmp(setElementTortoise->x, setElementHare->x) == 0)// && fmpz_cmp(setElementTortoise->y, setElementHare->y) != 0)
		{
			found = FindGX(temporary, numerator, denominator, mersenneGenerator, exponentTortoise0, exponentTortoise1,exponentHare0, exponentHare1, pointOrder);
			printf("Found Collision: g^x = ");
			fmpz_print(numerator);
			printf("\n");
			break;
		}
		printf("\n");

	}
	 
	DestroyPoint(setElementHare);
	DestroyPoint(setElementTortoise);
	DestroyPoint(temp0);
	DestroyPoint(temp1);
	DestroyPoint(tempGeneratorTortoise);
	DestroyPoint(tempTargetTortoise);
	DestroyPoint(tempGeneratorHare);
	DestroyPoint(tempTargetHare);
	fmpz_clear(numerator);fmpz_clear(denominator);fmpz_clear(pointOrderMinusOne);fmpz_clear(temporary);fmpz_clear(exponentTortoise0);fmpz_clear(exponentTortoise1);fmpz_clear(exponentHare0);fmpz_clear(exponentHare1);
	return result;
}

void TestSmallCurve()
{
	int primeNumberHolder = 43573;
	assert(primeNumberHolder % 3 == 1);
	assert(primeNumberHolder < first5239[PRIME_COUNT-1]);
	
	int mersenneExponents [] = {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933};
	int fermatExponents[] = {1, 2, 4, 8, 16};
	int fermatLength = sizeof(fermatLength) / sizeof(int);
	int mersenneLength = sizeof(mersenneExponents) / sizeof(int);
	int *foundMersenne = calloc(mersenneLength, sizeof(int));
	int *foundFermat   = calloc(mersenneLength, sizeof(int));
	
	
	Elliptic generator = CreatePoint();
	Elliptic target = CreatePoint();
	Elliptic result0 = CreatePoint();
	flint_rand_t randomState;
	fmpz_factor_t factorization;
	fmpz_t pointOrder,pointOrderMinusOne, testQuotient, primeNumber, testPrime, targetKey, mersenneGenerator;
	
	fmpz_init(testQuotient);
	fmpz_init(pointOrderMinusOne);
	fmpz_init(pointOrder);
	fmpz_init(primeNumber);
	fmpz_init(testPrime);
	fmpz_init(targetKey);
	fmpz_init(mersenneGenerator);
	fmpz_factor_init(factorization);
	flint_randinit(randomState);
	flint_randseed(randomState, 7194586,126043);//flint_randseed(randomState, 79929,1223); mod 43573
	fmpz_set_ui(primeNumber, primeNumberHolder);
	srand(43246156);
	//Loop to find point with prime order
	bool foundPrimePointOrder = false;
	int attempts = 1000;
	while(foundPrimePointOrder == false)
	{
		bool foundSecpPoint = FindRandomSecpPoint(generator, primeNumber, randomState);
		assert(foundSecpPoint);
		bool foundPointOrder= BruteforcePointOrder(generator, primeNumber,pointOrder);
		assert(foundPointOrder);
		//printf("Found P(");fmpz_print(generator->x);printf(",");fmpz_print(generator->y);printf(") mod ");fmpz_print(primeNumber);printf("\n");
		//printf("Order is: ");fmpz_print(pointOrder);printf("\n");
		for(int i = 0; i < PRIME_COUNT; i++)
		{
			int currentPrime = first5239[i];
			if(fmpz_cmp_ui(pointOrder, currentPrime) == 0)
			{
				foundPrimePointOrder = true;
				break;
			}
			else if(fmpz_cmp_ui(pointOrder, currentPrime) < 0){break;}
		}
		attempts -= 1;
		if(attempts == 0){break;}
	}
	
	assert(foundPrimePointOrder == true);
	fmpz_sub_ui(pointOrderMinusOne, pointOrder, 1);
	fmpz_factor(factorization, pointOrderMinusOne);
	
	printf("Found P(");fmpz_print(generator->x);printf(",");fmpz_print(generator->y);printf(") mod ");fmpz_print(primeNumber);printf("\n");
	printf("Prime order is: ");fmpz_print(pointOrder);printf("\n");
	printf("Prime Factors of p-1:(%ld): ",factorization->num);
	for(slong i = 0; i < factorization->num; i++)
	{
		printf("(");fmpz_print(&factorization->p[i]);flint_printf("^%wd),", factorization->exp[i]);		
	}
	printf("\n");
	//Check Fermat Exponents
	for(int i = 0; i < fermatLength; i++)
	{
		fmpz_one(testPrime);
		fmpz_mul_2exp(testPrime, testPrime, fermatExponents[i]); 
		fmpz_add_ui(testPrime, testPrime, 1);
		if(fmpz_cmp(testPrime, primeNumber) > 0){break;}
		bool testResult = RaiseGeneratorTest(testPrime, pointOrder, pointOrderMinusOne, factorization);
		foundFermat[i] = 0;
		if(testResult == true)
		{
			foundFermat[i] = 1;
		}
		printf("(2^%d) + 1: ",fermatExponents[i]);fmpz_print(testPrime);printf(" %s\n", testResult ? "True" : "False");
	}
	//Check Mersenne Exponents
	int mersenneFoundCount = 0;
	int mersenneExponent = 0;
	for(int i = 0; i < mersenneLength; i++)
	{
		fmpz_one(testPrime);
		fmpz_mul_2exp(testPrime, testPrime, mersenneExponents[i]); 
		fmpz_sub_ui(testPrime, testPrime, 1);
		if(fmpz_cmp(testPrime, pointOrder) > 0){break;}
		bool testResult = RaiseGeneratorTest(testPrime, pointOrder, pointOrderMinusOne, factorization);
		foundMersenne[i] = 0;
		if(testResult == true)
		{
			foundMersenne[i] = 1;
			fmpz_set(mersenneGenerator, testPrime);
			mersenneExponent = mersenneExponents[i];
			mersenneFoundCount += 1;
		}
		printf("(2^%d) - 1: ",mersenneExponents[i]);fmpz_print(testPrime);printf(" %s\n", testResult ? "True" : "False");
	}
	
	//We focus on mersenne primes since more than fermat
	printf("\nFound %d mersenne primes\n",mersenneFoundCount);
	assert(mersenneFoundCount > 0);
	//Set to 
	
	//Set Target Key
	int bitCount = fmpz_sizeinbase(primeNumber, 2) + 5;
	fmpz_randm(targetKey, randomState, pointOrder);
	Elliptic *cachedDoubles = CacheDoubles(generator, primeNumber, bitCount);
	LSBCachedMultiplication(bitCount, cachedDoubles, target, generator, targetKey, primeNumber);

	PollardRho2(generator, target, pointOrder, primeNumber, mersenneGenerator, mersenneExponent);
	printf("Target key : "); fmpz_print(targetKey);
	printf("\nFound Q(");fmpz_print(target->x);printf(",");fmpz_print(target->y);printf(") ");printf("\n");
	
	

	for(int i = 0; i < bitCount; i++){DestroyPoint(cachedDoubles[i]);}free(cachedDoubles);
	free(foundMersenne);
	free(foundFermat);
	DestroyPoint(generator);
	DestroyPoint(target);
	DestroyPoint(result0);
	flint_randclear(randomState);
	fmpz_factor_clear(factorization);
	fmpz_clear(pointOrder);
	fmpz_clear(testQuotient);
	fmpz_clear(pointOrderMinusOne);
	fmpz_clear(primeNumber);
	fmpz_clear(testPrime);
	fmpz_clear(targetKey);
	fmpz_clear(mersenneGenerator);
}

int main()
{	
	//BenchMarkLSB();
	//FindPuzzleGenerator();
	TestSmallCurve();
	flint_cleanup();
	return 0;
}
