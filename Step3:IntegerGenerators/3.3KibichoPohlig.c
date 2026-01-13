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

//clear && gcc 3.3KibichoPohlig.c -lm -lgmp -lmpfr -lflint -o m.o && ./m.o
typedef struct elliptic_point_struct *Elliptic;
struct elliptic_point_struct
{
	fmpz_t x;
	fmpz_t y;
	fmpz_t genX;
	fmpz_t genY;
	fmpz_t a;
	fmpz_t a0;
	fmpz_t a1;
	fmpz_t D;
	int legendre;
	int infinity;
};

Elliptic CreatePoint()
{
	Elliptic point = malloc(sizeof(struct elliptic_point_struct));
	fmpz_init(point->x);
	fmpz_init(point->y);
	fmpz_init(point->genX);
	fmpz_init(point->genY);
	fmpz_init(point->a);
	fmpz_init(point->a0);
	fmpz_init(point->a1);
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
	printf("(x:");fmpz_print(point->x);printf(" y:");fmpz_print(point->y);
	printf(")(g:");fmpz_print(point->genX);printf(",");fmpz_print(point->genY);
	printf(")A:");fmpz_print(point->a);
	printf(" A0:");fmpz_print(point->a0);
	printf(" A1:");fmpz_print(point->a1);
	printf(" D:");fmpz_print(point->D);
	printf("(%d)\n", point->legendre);
}

void DestroyPoint(Elliptic point)
{
	if(point)
	{
		fmpz_clear(point->x);
		fmpz_clear(point->y);
		fmpz_clear(point->genX);
		fmpz_clear(point->genY);
		fmpz_clear(point->a);
		fmpz_clear(point->a0);
		fmpz_clear(point->a1);
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

void ComputeSP()
{

}


int fmpz_mod_sqrt_tonelli(fmpz_t r, const fmpz_t a, const fmpz_t p)
{
    fmpz_t exp, legendre, Q, z, c, t, R, b, tmp, two_pow;
    ulong S = 0, M, i;

    // Initialize all fmpz variables individually
    fmpz_init(exp);
    fmpz_init(legendre);
    fmpz_init(Q);
    fmpz_init(z);
    fmpz_init(c);
    fmpz_init(t);
    fmpz_init(R);
    fmpz_init(b);
    fmpz_init(tmp);
    fmpz_init(two_pow);

    // Check if a is a quadratic residue mod p using Euler's criterion
    fmpz_sub_ui(exp, p, 1);                // exp = p - 1
    fmpz_fdiv_q_2exp(exp, exp, 1);         // exp = (p - 1) / 2
    fmpz_powm(legendre, a, exp, p);        // legendre = a^((p-1)/2) mod p

    if (!fmpz_is_one(legendre)) {
        // Not a quadratic residue
        fmpz_clear(exp); fmpz_clear(legendre); fmpz_clear(Q);
        fmpz_clear(z); fmpz_clear(c); fmpz_clear(t); fmpz_clear(R);
        fmpz_clear(b); fmpz_clear(tmp); fmpz_clear(two_pow);
        return 0;
    }

    // Shortcut for p ≡ 3 mod 4
    if (fmpz_fdiv_ui(p, 4) == 3) {
        fmpz_add_ui(exp, p, 1);
        fmpz_fdiv_q_2exp(exp, exp, 2);     // exp = (p + 1) / 4
        fmpz_powm(r, a, exp, p);
        fmpz_clear(exp); fmpz_clear(legendre); fmpz_clear(Q);
        fmpz_clear(z); fmpz_clear(c); fmpz_clear(t); fmpz_clear(R);
        fmpz_clear(b); fmpz_clear(tmp); fmpz_clear(two_pow);
        return 1;
    }

    // Factor p - 1 = Q * 2^S
    fmpz_sub_ui(Q, p, 1); // Q = p - 1
    while (fmpz_tstbit(Q, 0) == 0) { // While Q is even
        fmpz_fdiv_q_2exp(Q, Q, 1);
        S++;
    }

    // Find quadratic non-residue z
    fmpz_set_ui(z, 2);
    while (1) {
        fmpz_powm(legendre, z, exp, p);
        if (!fmpz_is_one(legendre))
            break;
        fmpz_add_ui(z, z, 1);
    }

    // Initialize variables
    fmpz_powm(c, z, Q, p);
    fmpz_powm(t, a, Q, p);

    fmpz_add_ui(exp, Q, 1);
    fmpz_fdiv_q_2exp(exp, exp, 1);  // (Q + 1) / 2
    fmpz_powm(R, a, exp, p);

    M = S;

    // Main loop
    while (!fmpz_is_one(t)) {
        fmpz_set(tmp, t);
        for (i = 1; i < M; i++) {
            fmpz_powm_ui(tmp, tmp, 2, p);
            if (fmpz_is_one(tmp))
                break;
        }

        if (i == M) {
            // No solution
            fmpz_clear(exp); fmpz_clear(legendre); fmpz_clear(Q);
            fmpz_clear(z); fmpz_clear(c); fmpz_clear(t); fmpz_clear(R);
            fmpz_clear(b); fmpz_clear(tmp); fmpz_clear(two_pow);
            return 0;
        }

        fmpz_set_ui(two_pow, 1UL << (M - i - 1));
        fmpz_powm(b, c, two_pow, p);     // b = c^(2^{M-i-1}) mod p

        fmpz_mul(R, R, b);
        fmpz_mod(R, R, p);

        fmpz_powm_ui(c, b, 2, p);        // c = b^2
        fmpz_mul(t, t, c);
        fmpz_mod(t, t, p);

        M = i;
    }

    fmpz_set(r, R);

    // Clean up
    fmpz_clear(exp); fmpz_clear(legendre); fmpz_clear(Q);
    fmpz_clear(z); fmpz_clear(c); fmpz_clear(t); fmpz_clear(R);
    fmpz_clear(b); fmpz_clear(tmp); fmpz_clear(two_pow);
    return 1;
}


void SolveA(Elliptic point, fmpz_t scalar, fmpz_t pointOrder)
{
	fmpz_t temp;
	fmpz_init(temp);
	fmpz_sub(point->a, scalar, point->x);
	fmpz_invmod(temp, point->y, pointOrder);
	fmpz_mul(point->a, point->a, temp);
	fmpz_mod(point->a, point->a, pointOrder);	
	fmpz_clear(temp);	
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
	point->legendre = fmpz_jacobi(point->D, N);
	if(point->legendre == 1)
	{
		fmpz_mod_sqrt_tonelli(point->a0, point->D, N);
		fmpz_sub(point->a1, N, point->a0);
	}
	else
	{
		fmpz_set_ui(point->a0, 0);
		fmpz_set_ui(point->a1, 0);
	}
	fmpz_clear(tmp1);
	fmpz_clear(tmp2);
	return true;
}

void Pohlig(Elliptic generator, Elliptic target, fmpz_t primeNumber, fmpz_t pointOrder, fmpz_t pointOrderMinusOne, fmpz_t mersenneGenerator, fmpz_factor_t factorization)
{
	fmpz_t temp0, exponentG, exponentT,privateKey;
	fmpz_init(temp0);
	fmpz_init(exponentG);
	fmpz_init(exponentT);
	fmpz_init(privateKey);
	Elliptic result = CreatePoint();
	int bitCount = 130;
	Elliptic *cachedDoubles = CacheDoubles(generator, primeNumber, bitCount);
	printf("Prime Factors of pointOrder-1:(%ld): ",factorization->num);
	for(slong i = 0; i < factorization->num; i++)
	{
		fmpz_powm(&factorization->p[i],&factorization->p[i],&factorization->exp[i],pointOrderMinusOne);
		fmpz_set_ui(&factorization->exp[i],1);
		printf("(");fmpz_print(&factorization->p[i]);flint_printf("^%wd),", factorization->exp[i]);		
	}
	printf("\n");
	int steps = 1;
	for(int exponent = 0; exponent < fmpz_get_ui(pointOrderMinusOne); exponent++)
	{
		fmpz_set_ui(exponentG, exponent);
		fmpz_powm(privateKey, mersenneGenerator, exponentG, pointOrder);
		LSBCachedMultiplication(bitCount, cachedDoubles, result, generator, privateKey, primeNumber);
		fmpz_mod(result->genX, result->x, pointOrder);fmpz_mod(result->genY, result->y, pointOrder);
		
		SolveA(result, privateKey, pointOrder);
		SolveD(result, pointOrder);
		//if(fmpz_cmp_ui(result->x, 22830) == 0)
		{
			fmpz_print(mersenneGenerator); printf(" ^ %d : (",exponent);
			fmpz_print(privateKey); printf(")");
			PrintPoint(result);
		}
	}
	for(int i = 0; i < bitCount; i++){DestroyPoint(cachedDoubles[i]);}free(cachedDoubles);
	fmpz_clear(privateKey);
	fmpz_clear(exponentG);
	fmpz_clear(exponentT);
	fmpz_clear(temp0);
	DestroyPoint(result);
}

void TestSmallCurve()
{
	int mersenneExponent = 7;
	int primeNumberHolder = 43573;
	int pointOrderHolder = 2293;
	assert(primeNumberHolder % 3 == 1);
	assert(primeNumberHolder < first5239[PRIME_COUNT-1]);
	
	Elliptic generator = CreatePoint();
	Elliptic target = CreatePoint();
	Elliptic result0 = CreatePoint();
	fmpz_t pointOrder,pointOrderMinusOne, primeNumber, targetKey, mersenneGenerator;
	fmpz_factor_t factorization;
	
	fmpz_factor_init(factorization);
	fmpz_init(pointOrder);
	fmpz_init(pointOrderMinusOne);
	fmpz_init(primeNumber);
	fmpz_init(targetKey);
	fmpz_init(mersenneGenerator);
	
	fmpz_one(mersenneGenerator);
	fmpz_mul_2exp(mersenneGenerator, mersenneGenerator, mersenneExponent); 
	fmpz_sub_ui(mersenneGenerator, mersenneGenerator, 1);
	
	fmpz_set_ui(primeNumber, primeNumberHolder);
	fmpz_set_ui(pointOrder, pointOrderHolder);
	fmpz_sub_ui(pointOrderMinusOne, pointOrder, 1);
	fmpz_factor(factorization, pointOrderMinusOne);
	
	//Set Generator
	fmpz_set_ui(generator->x, 9382);
	fmpz_set_ui(generator->y, 5436);
	//Set Target
	fmpz_set_ui(target->x, 1005);
	fmpz_set_ui(target->y, 1043);

	Pohlig(generator, target, primeNumber, pointOrder, pointOrderMinusOne, mersenneGenerator, factorization);
	
	fmpz_clear(pointOrder);
	fmpz_clear(pointOrderMinusOne);
	fmpz_clear(primeNumber);
	fmpz_clear(targetKey);
	fmpz_clear(mersenneGenerator);
	fmpz_factor_clear(factorization);
	DestroyPoint(generator);
	DestroyPoint(target);
	DestroyPoint(result0);
}

int main()
{
	TestSmallCurve();
	flint_cleanup();
	return 0;
}
