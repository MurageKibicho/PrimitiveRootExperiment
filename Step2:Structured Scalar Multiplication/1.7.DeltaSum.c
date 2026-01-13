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
#include <string.h> 
#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"
//clear && gcc 1.6.CachedSum.c -lm -lgmp -lmpfr -lflint -o m.o && ./m.o

typedef struct elliptic_point_struct *Elliptic;
typedef struct delta_struct *Delta;
struct delta_struct
{
	fmpz_t x;
	fmpz_t y;
	fmpz_t deltaX;
	fmpz_t deltaY;
	fmpz_t accumulatedX;
	fmpz_t accumulatedY;
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
	fmpz_init(delta->deltaX);
	fmpz_init(delta->deltaY);
	fmpz_init(delta->accumulatedX);
	fmpz_init(delta->accumulatedY);
	delta->bit = 0;
	return delta;
}

void DestroyDelta(Delta delta)
{
	if(delta)
	{
		fmpz_clear(delta->x);
		fmpz_clear(delta->y);
		fmpz_clear(delta->deltaX);
		fmpz_clear(delta->deltaY);
		fmpz_clear(delta->accumulatedX);
		fmpz_clear(delta->accumulatedY);
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
	printf("X:");fmpz_print(point->x);printf("\nY:");fmpz_print(point->y);printf("\n");
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

void PrintDelta(Delta delta)
{
	printf("X :");fmpz_print(delta->x);printf("\nY :");fmpz_print(delta->y);
	printf("\nΔX:");fmpz_print(delta->deltaX);printf("\nΔY:");fmpz_print(delta->deltaY);
	printf("\nTX:");fmpz_print(delta->accumulatedX);printf("\nTY:");fmpz_print(delta->accumulatedY);
	printf("\n");
}

void PrintDeltaArray(size_t length, Delta *deltas)
{
	for(int i = length-1; i >= 0; i--)
	{
		printf("(%d)\n",deltas[i]->bit);
		if(deltas[i]->bit == 1)
		{
			PrintDelta(deltas[i]);
		}
	}
	printf("\n");
}


void CachedMultiplication(int bitCount, Elliptic *cachedDoubles, Elliptic resultant, Elliptic generator, fmpz_t privateKey, fmpz_t primeNumber)
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

Delta *CopyDeltas(Delta* original)
{
	Delta* copy = NULL;
	for(size_t i = 0; i < arrlen(original); i++)
	{
		Delta newDelta = CreateDelta();
		fmpz_set(newDelta->x, original[i]->x);
		fmpz_set(newDelta->y, original[i]->y);
		fmpz_set(newDelta->deltaX, original[i]->deltaX);
		fmpz_set(newDelta->deltaY, original[i]->deltaY);
		fmpz_set(newDelta->accumulatedX, original[i]->accumulatedX);
		fmpz_set(newDelta->accumulatedY, original[i]->accumulatedY);
		newDelta->bit = original[i]->bit;
		arrput(copy, newDelta);
	}
	return copy;
}


void DeltaForwardPass(size_t targetLength, int bitCount, Elliptic lossLSB,Elliptic targetLSB, Elliptic *cachedDoubles, Delta *deltas,fmpz_t primeNumber)
{
	Elliptic pointAtInfinity = CreatePoint();
	pointAtInfinity->infinity = 1;  
	Elliptic temp1 = CreatePoint();

	fmpz_t xPrevious, yPrevious, delta,accumulatedX,accumulatedY;
	fmpz_init(xPrevious);fmpz_init(yPrevious);
	fmpz_init(delta);
	fmpz_init(accumulatedX);fmpz_init(accumulatedY);
	for(int i = targetLength - 1, j = 0; i >= 0; i--, j++)
	{
		int bit = deltas[i]->bit;
		if(bit)
		{
			AddCurvePoints(temp1, pointAtInfinity, cachedDoubles[j], primeNumber);
			CopyPoint(temp1, pointAtInfinity);
			//Copy
			fmpz_set(deltas[i]->x, pointAtInfinity->x);
			fmpz_set(deltas[i]->y, pointAtInfinity->y);
			//Find Delta x
			fmpz_sub(delta, pointAtInfinity->x, xPrevious);
			fmpz_mod(delta,delta,primeNumber);
			fmpz_set(xPrevious, pointAtInfinity->x);	
			fmpz_set(deltas[i]->deltaX, delta);	
			fmpz_add(accumulatedX, accumulatedX, deltas[i]->deltaX);
			fmpz_mod(accumulatedX, accumulatedX, primeNumber);
			fmpz_set(deltas[i]->accumulatedX, accumulatedX);
			//Find Delta y
			fmpz_sub(delta, pointAtInfinity->y, yPrevious);
			fmpz_mod(delta,delta,primeNumber);
			fmpz_set(yPrevious, pointAtInfinity->y);	
			fmpz_set(deltas[i]->deltaY, delta);	
			fmpz_add(accumulatedY, accumulatedY, deltas[i]->deltaY);
			fmpz_mod(accumulatedY, accumulatedY, primeNumber);
			fmpz_set(deltas[i]->accumulatedY, accumulatedY);
		}
	}
	//Find Loss
	CopyPoint(pointAtInfinity, lossLSB);
	DestroyPoint(pointAtInfinity);
	DestroyPoint(temp1);
	fmpz_clear(xPrevious);fmpz_clear(yPrevious);
	fmpz_clear(delta);
	fmpz_clear(accumulatedX);fmpz_clear(accumulatedY);
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

void DeltaRandomPass(size_t targetLength, int bitCount, Elliptic lossLSB,Elliptic targetLSB, Elliptic *cachedDoubles, Delta *deltas,fmpz_t primeNumber)
{
	Delta *copies = CopyDeltas(deltas);
	Elliptic tempLoss = CreatePoint();
	int bestIndex = -1;
	int bitFlipCount = 1;
	int *visited = calloc(targetLength, sizeof(int));
	int *flipped = calloc(targetLength, sizeof(int));
	int *bestFlip= calloc(targetLength, sizeof(int));
	fmpz_t temp0, totalLoss, currentLoss, bestLoss;
	fmpz_init(temp0);
	fmpz_init(totalLoss);
	fmpz_init(currentLoss);
	fmpz_init(bestLoss);

	SolveD(lossLSB, primeNumber);
	fmpz_set(bestLoss, lossLSB->D);
	printf("Initial Loss:%ld bits\n", fmpz_sizeinbase(totalLoss, 2));
	int stagnationCount = 0;
	for(int tests = 0; tests < 100000; tests++)
	{
		bestIndex = -1;
		for(int i = 0; i < targetLength; i++)
		{
			if(bitFlipCount == 1)
			{
				copies[i]->bit = !copies[i]->bit;
			}
			else
			{
				memset(flipped,0,targetLength * sizeof(int));
				for(int j = 0; j < bitFlipCount; j++)
				{
					int bitIndex = rand() % targetLength;
					copies[bitIndex]->bit = !copies[bitIndex]->bit;	
					flipped[bitIndex] = 1;
				}
				copies[0]->bit = 1;
			}
			DeltaForwardPass(targetLength, bitCount, tempLoss, targetLSB, cachedDoubles, copies, primeNumber);
			SolveD(tempLoss, primeNumber);
			fmpz_set(currentLoss, tempLoss->D);
			//printf("New Loss:%ld bits Best: %d\n", fmpz_sizeinbase(currentLoss, 2),bestIndex);
			if(fmpz_cmp(bestLoss, currentLoss) > 0 && visited[i] == 0)
			{
				bestIndex = i;
				fmpz_set(bestLoss, currentLoss);
				if(bitFlipCount > 1)
				{
					memcpy(bestFlip, flipped, targetLength * sizeof(int));	
					printf("{");
					for(int j = 0; j < targetLength; j++)
					{
						int bit = copies[j]->bit;
						
						printf("%d",bit);
					}
					printf("}\n");
				}
			}
			
			if(bitFlipCount == 1)
			{
				copies[i]->bit = !copies[i]->bit;
			}
			else
			{
				for(int j = 0; j < targetLength; j++)
				{
					if(flipped[j] == 1)
					{
						int bitIndex = j;
						copies[bitIndex]->bit = !copies[bitIndex]->bit;	
					}
				}
				copies[0]->bit = 1;
			}
		}
		if(bestIndex > -1)
		{
			stagnationCount = 0;
			if(bitFlipCount == 1)
			{
				copies[bestIndex]->bit = !copies[bestIndex]->bit;
				visited[bestIndex] = 1;	
			}
			else
			{
				for(int j = 0; j < targetLength; j++)
				{
					if(bestFlip[j] == 1)
					{
						int bitIndex = j;
						copies[bitIndex]->bit = !copies[bitIndex]->bit;	
					}
				}
			}

			printf("Best Loss:%ld bits:", fmpz_sizeinbase(bestLoss, 2));
			for(int j = arrlen(copies)-1; j >= 0; j--){printf("%d", copies[j]->bit);}printf("\n");
			if(fmpz_cmp_ui(bestLoss, 0) == 0)
			{
				printf("Found at iteration (%d)",tests);
				break;
			}
		}
		else
		{
			//printf("No improvements after(%d)\n", tests);
			bitFlipCount = rand() % targetLength;
			stagnationCount += 1;
			if(stagnationCount == 100000)
			{
				break;
			}
		}
	}
	fmpz_clear(temp0);
	fmpz_clear(totalLoss);
	fmpz_clear(currentLoss);
	fmpz_clear(bestLoss);
	free(visited);
	free(flipped);
	free(bestFlip);
	DestroyPoint(tempLoss);
	for(size_t i = 0; i < arrlen(copies); i++){DestroyDelta(copies[i]);}arrfree(copies);
}

void LearnPrivateKey(size_t targetLength, int bitCount, Elliptic *cachedDoubles, Elliptic publicKey, Elliptic generator, fmpz_t primeNumber)
{
	srand(787368678);
	printf("Target key:%ld\nTarget Point:\n",targetLength);
	PrintPoint(publicKey);
	printf("\n");
	Elliptic lossLSB = CreatePoint();
	Delta *deltas = NULL;
	//Initialize Deltas
	for(size_t i = 0; i < targetLength; i++)
	{
		Delta delta = CreateDelta();
		delta->bit = rand() % 2;
		arrput(deltas, delta);
	}
	deltas[0]->bit = 1;
	//ForwardPass
	DeltaForwardPass(targetLength, bitCount, lossLSB, publicKey, cachedDoubles, deltas, primeNumber);
	//PrintDeltaArray(targetLength, deltas);
	//PrintDelta(deltas[0]);
	//PrintPoint(lossLSB);
	DeltaRandomPass(targetLength, bitCount, lossLSB, publicKey, cachedDoubles, deltas, primeNumber);

	DestroyPoint(lossLSB);
	for(size_t i = 0; i < arrlen(deltas); i++){DestroyDelta(deltas[i]);}arrfree(deltas);
}


void BenchMarkLSB()
{
	int keyCount = 49274;
	int bitCount = 161;
	char *primeNumberHexadecimal = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F";	
	char *generatorXHexadecimal  = "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798";
	char *generatorYHexadecimal  = "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8";
	char *X160  = "e0a8b039282faf6fe0fd769cfbc4b6b4cf8758ba68220eac420e32b91ddfa673";
	char *Y160  = "c2d9690945dd98f6e0e45d4a1f760c9e85ed5ae5ffeeda74e121ee0d836a7c86";
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
	printf("Test: %d\n", keyCount-1);
	start = clock();
	for(int i = keyCount-1; i < keyCount; i++)
	{
		fmpz_set_ui(privateKey, 5);
		CachedMultiplication(bitCount, cachedDoubles, result0, generator, privateKey, primeNumber);
		//fmpz_set_str(result0->x, X160, 16);
		//fmpz_set_str(result0->y, Y160, 16);
		//size_t targetLength = 160;
		//LearnPrivateKey(targetLength, bitCount, cachedDoubles, result0, generator, primeNumber);
		PrintPoint(result0);
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
