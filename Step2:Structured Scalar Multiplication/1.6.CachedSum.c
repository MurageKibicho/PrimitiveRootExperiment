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

//clear && gcc 1.6.CachedSum.c -lm -lgmp -lmpfr -lflint -o m.o && ./m.o

typedef struct elliptic_point_struct *Elliptic;
typedef struct delta_struct *Delta;
struct delta_struct
{
	fmpz_t x;
	fmpz_t y;
	fmpz_t D;
	fmpz_t deltaD;
	fmpz_t accumulatedD;
	int bit;
};
struct elliptic_point_struct
{
	fmpz_t x;
	fmpz_t y;
	fmpz_t D;
	int infinity;
};

Delta CreateDelta()
{
	Delta delta = malloc(sizeof(struct delta_struct));
	fmpz_init(delta->x);
	fmpz_init(delta->y);
	fmpz_init(delta->D);
	fmpz_init(delta->deltaD);
	fmpz_init(delta->accumulatedD);
	delta->bit = 0;
	return delta;
}

void DestroyDelta(Delta delta)
{
	if(delta)
	{
		fmpz_clear(delta->x);
		fmpz_clear(delta->y);
		fmpz_clear(delta->D);
		fmpz_clear(delta->deltaD);
		fmpz_clear(delta->accumulatedD);
		free(delta);
	}
}

Elliptic CreatePoint()
{
	Elliptic point = malloc(sizeof(struct elliptic_point_struct));
	fmpz_init(point->x);
	fmpz_init(point->y);
	fmpz_init(point->D);
	point->infinity = 0;
	return point;
}

void CopyPoint(Elliptic source, Elliptic destination)
{
	fmpz_set(destination->x, source->x);
	fmpz_set(destination->y, source->y);
	fmpz_set(destination->D, source->D);
	destination->infinity = source->infinity;
}

int TestPointEquality(Elliptic a, Elliptic b)
{
	if(a->infinity && b->infinity) return 1;
	if(a->infinity || b->infinity) return 0;
	return fmpz_cmp(a->x, b->x) == 0 && fmpz_cmp(a->y, b->y) == 0;	
}

void PrintPoint(Elliptic point)
{
	printf("X:");fmpz_print(point->x);printf("\nY:");fmpz_print(point->y);printf("\nD:");fmpz_print(point->D);printf("\n");
}

void PrintDelta(Delta delta)
{
	printf("X:");fmpz_print(delta->x);printf("\nY:");fmpz_print(delta->y);printf("\nD:");fmpz_print(delta->D);
	printf("\nΔ:");fmpz_print(delta->deltaD);printf("\nA:");fmpz_print(delta->accumulatedD);printf("\n");
}


void DestroyPoint(Elliptic point)
{
	if(point)
	{
		fmpz_clear(point->x);
		fmpz_clear(point->y);
		fmpz_clear(point->D);
		free(point);
	}
}

bool SolveD(Elliptic point, fmpz_t N)
{
	fmpz_t tmp1, tmp2;
	fmpz_init(tmp1);
	fmpz_init(tmp2);
	//tmp1 = x^2 mod N
	fmpz_powm_ui(tmp1, point->x, 2, N);
	//tmp2 = y^2 mod N
	fmpz_powm_ui(tmp2, point->y, 2, N);

	// Compute D = (x^2 - 1) * (y^2)^(-1) mod N
	fmpz_sub_ui(tmp1, tmp1, 1);    
	fmpz_mod(tmp1, tmp1, N); 

	if(fmpz_invmod(tmp2, tmp2, N) == 0)
	{
		// y^2 not invertible mod N
		fmpz_set_si(point->D, -1);
		fmpz_clear(tmp1);
		fmpz_clear(tmp2);
		return false;
	}

	fmpz_mul(point->D, tmp1, tmp2);
	fmpz_mod(point->D, point->D, N);

	fmpz_clear(tmp1);
	fmpz_clear(tmp2);
	return true;
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
		SolveD(doubles[i], primeNumber);
		//Compute next power of 2: temp = 2 * temp
		AddCurvePoints(temp2, temp, temp, primeNumber);
		CopyPoint(temp2, temp);
	}
	DestroyPoint(temp2);
	DestroyPoint(temp);
	return doubles;
}

void PrintDoubles(int bitCount, Elliptic *doubles)
{
	for(int i = 0; i < bitCount; i++)
	{
		printf("2^(%3d) G:\n", i+1);
		PrintPoint(doubles[i]);
	}

}

void LSBCachedMultiplication(int bitCount, Elliptic *cachedDoubles, Elliptic resultant, Elliptic generator, fmpz_t privateKey, fmpz_t primeNumber)
{
	Elliptic pointAtInfinity = CreatePoint();
	pointAtInfinity->infinity = 1;  
	Elliptic temp1 = CreatePoint();

	size_t binaryLength = fmpz_sizeinbase(privateKey, 2);
	Delta *deltas = malloc(binaryLength * sizeof(Delta));
	fmpz_t DPrevious, deltaD,accumulatedD;
	fmpz_init(DPrevious);
	fmpz_init(deltaD);
	fmpz_init(accumulatedD);
	fmpz_set(DPrevious, generator->D);
	fmpz_set(accumulatedD, generator->D);
	for(size_t i = 0; i < binaryLength; ++i)
	{
		int bit = fmpz_tstbit(privateKey, i);
		deltas[i] = CreateDelta();
		deltas[i]->bit = bit;
		if(bit)
		{
			AddCurvePoints(temp1, pointAtInfinity, cachedDoubles[i], primeNumber);
			CopyPoint(temp1, pointAtInfinity);
			SolveD(pointAtInfinity, primeNumber);
			//Copy
			fmpz_set(deltas[i]->x, pointAtInfinity->x);
			fmpz_set(deltas[i]->y, pointAtInfinity->y);
			fmpz_set(deltas[i]->D, pointAtInfinity->D);
			
			fmpz_sub(deltaD, pointAtInfinity->x, DPrevious);
			fmpz_mod(deltaD,deltaD,primeNumber);
			fmpz_set(DPrevious, pointAtInfinity->x);	
			fmpz_set(deltas[i]->deltaD, deltaD);	
			
			fmpz_add(accumulatedD, accumulatedD, deltas[i]->deltaD);
			//fmpz_mod(accumulatedD, accumulatedD, primeNumber);
			
			fmpz_set(deltas[i]->accumulatedD, accumulatedD);
		}
	}
	//Save pointAtInfinity to the result variable
	CopyPoint(pointAtInfinity, resultant);
	
	//Print Target
	for(size_t i = 0; i < binaryLength; ++i)
	{
		printf("(%d)\n",deltas[i]->bit);
		if(deltas[i]->bit == 1)
		{
			PrintDelta(deltas[i]);
		}
	}
	printf("\n");
	SolveD(resultant,primeNumber);
	PrintPoint(resultant);
	
	DestroyPoint(pointAtInfinity);
	DestroyPoint(temp1);
	fmpz_clear(DPrevious);
	fmpz_clear(deltaD);
	fmpz_clear(accumulatedD);
	for(size_t i = 0; i < binaryLength; ++i){DestroyDelta(deltas[i]);}free(deltas);
}


void BenchMarkLSB()
{
	int keyCount = 97566;
	int bitCount = 131;
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
	//PrintDoubles(bitCount, cachedDoubles);
	clock_t start, end;
	double cpuTimeCount;
	printf("Benchmarking MySecp public key generation...\n");
	printf("Test: %d|\n", keyCount-1);
	start = clock();
	for(int i = keyCount-1; i < keyCount; i++)
	{
		fmpz_set_ui(privateKey, i);
		LSBCachedMultiplication(bitCount, cachedDoubles, result1, generator, privateKey, primeNumber);
		//PrintPoint(result1);
		//printf("\n");
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

int main()
{	
	BenchMarkLSB();
	flint_cleanup();
	return 0;
}
