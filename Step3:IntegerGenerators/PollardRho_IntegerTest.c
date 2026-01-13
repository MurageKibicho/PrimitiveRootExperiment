#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <gmp.h>

//clear && gcc PollardRho_IntegerTest.c -lm -lgmp -o m.o && ./m.o

void PollardRho_IntegerUpdate(mpz_t setElement, mpz_t exponentAnimal0, mpz_t exponentAnimal1, mpz_t desiredGenerator, mpz_t desiredResult, mpz_t primeNumberMinOne, mpz_t primeNumber, mpz_t temporary)
{
	//Three update rules as defined in (Pollard, 1978)
	int modulo = 3;
	int moduloResult = mpz_mod_ui(temporary, setElement, modulo);
	if(moduloResult == 2)
	{
		//Multiply the set element by desiredResult mod primeNumber
		mpz_mul(setElement, setElement, desiredResult);
		mpz_mod(setElement, setElement, primeNumber);
	
		//Update exponentAnimal1 by 1 mod primeNumberMinOne
		mpz_add_ui(exponentAnimal1, exponentAnimal1, 1);
		mpz_mod(exponentAnimal1, exponentAnimal1, primeNumberMinOne);
	}
	else if(moduloResult == 1)
	{
		//Multiply the set element by desiredGenerator mod primeNumber
		mpz_mul(setElement, setElement, desiredGenerator);
		mpz_mod(setElement, setElement, primeNumber);

		//Update exponentAnimal0 by 1 mod primeNumberMinOne
		mpz_add_ui(exponentAnimal0, exponentAnimal0, 1);
		mpz_mod(exponentAnimal0, exponentAnimal0, primeNumberMinOne);
	}
	else
	{
		//Square the set element and mod primeNumber
		mpz_mul(setElement, setElement, setElement);
		mpz_mod(setElement, setElement, primeNumber);
		
		//Double both result and generator mod primeNumberMinOne
		mpz_mul_ui(exponentAnimal0, exponentAnimal0, 2);
		mpz_mul_ui(exponentAnimal1, exponentAnimal1, 2);
		mpz_mod(exponentAnimal0, exponentAnimal0, primeNumberMinOne);
		mpz_mod(exponentAnimal1, exponentAnimal1, primeNumberMinOne);
	}
}
bool PollardRho_IntegerDLP()
{
	bool result = false;
	mpz_t primeNumber, primeNumberMinOne, desiredGenerator, desiredResult;
	
	mpz_t exponentTortoise0, exponentTortoise1,setElementTortoise;
	mpz_t exponentHare0, exponentHare1, setElementHare,temporary;
	mpz_inits(exponentTortoise0, exponentTortoise1,setElementTortoise, exponentHare0, exponentHare1,setElementHare,temporary,primeNumber, primeNumberMinOne, desiredGenerator, desiredResult, NULL);
	
	//Set prime number and primeNumber Minus One and desires
	mpz_set_ui(primeNumber, 43573);
	mpz_sub_ui(primeNumberMinOne, primeNumber, 1);
	mpz_set_ui(desiredGenerator, 127);
	mpz_set_ui(desiredResult, 1984);
	
	//Set generator and desired result
	
	//Initialize generators, exponents and results to 0
	mpz_set_ui(exponentTortoise0, 0);mpz_set_ui(exponentHare0, 0);
	mpz_set_ui(exponentTortoise1, 0);mpz_set_ui(exponentHare1, 0);
	
	//Initialize our current set element to 1
	mpz_set_ui(setElementTortoise, 1);
	mpz_set_ui(setElementHare, 1);
	
	for(int i = 0; i < 100; i++)
	{
		//Update tortoise once
		PollardRho_IntegerUpdate(setElementTortoise, exponentTortoise0, exponentTortoise1, desiredGenerator, desiredResult, primeNumberMinOne, primeNumber, temporary);
		//Update hare twice
		PollardRho_IntegerUpdate(setElementHare, exponentHare0, exponentHare1, desiredGenerator, desiredResult, primeNumberMinOne, primeNumber, temporary);
		PollardRho_IntegerUpdate(setElementHare, exponentHare0, exponentHare1, desiredGenerator, desiredResult, primeNumberMinOne, primeNumber, temporary);
		
		//Print the current table
		gmp_printf("%d:\nTortoise: (%Zd^%3Zd * %Zd^%3Zd = %5Zd)\n", i,desiredGenerator, exponentTortoise0, desiredResult, exponentTortoise1, setElementTortoise);
		gmp_printf("Hare    : (%Zd^%3Zd * %Zd^%3Zd = %5Zd)\n", desiredGenerator, exponentHare0, desiredResult, exponentHare1, setElementHare);
		if(mpz_cmp(setElementTortoise, setElementHare) == 0)
		{
			printf("Found Collision\n");
			break;
		}
	}
	
	mpz_clears(exponentTortoise0, exponentTortoise1,setElementTortoise, exponentHare0, exponentHare1,setElementHare,temporary,primeNumber, primeNumberMinOne, desiredGenerator, desiredResult, NULL);
	return result;
}

int main()
{
	PollardRho_IntegerDLP();
	return 0;
}
