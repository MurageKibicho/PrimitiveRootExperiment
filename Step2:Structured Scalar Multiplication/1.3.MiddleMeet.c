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

//clear && gcc 1.3.MiddleMeet.c -lm -lgmp -lmpfr -lflint -o m.o && ./m.o

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


void BothSideScalarMultiplicationMidpoint(Elliptic resultant, Elliptic generator, fmpz_t privateKey, fmpz_t primeNumber)
{
	size_t binaryLength = fmpz_sizeinbase(privateKey, 2);
	size_t midpoint = 5;
    
	Elliptic pointMSB = CreatePoint();
	Elliptic pointLSB = CreatePoint();
	Elliptic temp = CreatePoint();

	pointMSB->infinity = 1;
	pointLSB->infinity = 1;

	Elliptic baseLSB = CreatePoint();
	Elliptic baseMSB = CreatePoint();
	CopyPoint(generator, baseLSB);
	CopyPoint(generator, baseMSB);
    
	// Precompute 2^midpoint * G for MSB side
	for(size_t i = 0; i < midpoint; i++)
	{
		AddCurvePoints(temp, baseMSB, baseMSB, primeNumber);
		CopyPoint(temp, baseMSB);
	}
    
	for(size_t i = 0; i < binaryLength; i++)
	{
		if(i < midpoint)
		{
			//LSB processing
			if(fmpz_tstbit(privateKey, i))
			{
				AddCurvePoints(temp, pointLSB, baseLSB, primeNumber);
				CopyPoint(temp, pointLSB);
			}
			AddCurvePoints(temp, baseLSB, baseLSB, primeNumber);
			CopyPoint(temp, baseLSB);
		}

		if(i >= midpoint)
		{
			//MSB side processing
			size_t msbIndex = binaryLength - 1 - (i - midpoint);
			printf("%ld %ld\n", i, msbIndex);
			AddCurvePoints(temp, pointMSB, pointMSB, primeNumber);
			CopyPoint(temp, pointMSB);

			if(fmpz_tstbit(privateKey, msbIndex))
			{
				AddCurvePoints(temp, pointMSB, baseMSB, primeNumber);
				CopyPoint(temp, pointMSB);
			}
		}
	}

	// Combine results
	AddCurvePoints(temp, pointLSB, pointMSB, primeNumber);
	CopyPoint(temp, resultant);


	// Cleanup
	DestroyPoint(pointMSB);
	DestroyPoint(pointLSB);
	DestroyPoint(baseLSB);
	DestroyPoint(baseMSB);
	DestroyPoint(temp);
}

void MiddleMeetBenchmark()
{
	int keyCount = 100000;
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
	for(int i = keyCount - 1; i < keyCount; i++)
	{
		fmpz_set_ui(privateKey, i);
		BothSideScalarMultiplicationMidpoint(result, generator, privateKey, primeNumber);	
	}
	PrintPoint(result);
	end = clock(); 

	cpuTimeCount = ((double)(end - start)) / CLOCKS_PER_SEC;
	double iterationsPerSecond = keyCount / cpuTimeCount;

	printf("Total time MSB: %.4f seconds\n", cpuTimeCount);
	
	DestroyPoint(generator);
	DestroyPoint(result);
	fmpz_clear(primeNumber);
	fmpz_clear(privateKey);
}



int main()
{
	MiddleMeetBenchmark();
	
	flint_cleanup();
	return 0;
}
