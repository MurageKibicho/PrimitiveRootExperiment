# Unofficial Repo for the Substack Post 'Using Group Theory to Speed up an Elliptic Curve Library from 352 million CPU years to 12 million CPU years'

We used all the group theory we know to speed up operations on the secp elliptic curve in this [Substack post](https://leetarxiv.substack.com/p/if-youre-smart-why-are-you-poor-elliptic).
Ultimately, our goal was to investigate the possibility of a p-1 logarithmic attack on Satoshi's wallet. 

## Conclusions
- This approach does not permit a sub-exponential p-1 attack because our final collision algorithm is additive, not multiplicative.
- We found an insanely fast algorithm for generating lots of valid, ordered (and unique) elliptic curve points.
  
## How to use this code
We assume you have FLINT and LibGMP installed. This library worked with FLINT version 3.0

It may not run with FLINT version 3.4 due to FLINT's breaking changes.

The folder numbers align with sections of the [original Substack article](https://leetarxiv.substack.com/p/if-youre-smart-why-are-you-poor-elliptic).

`Step1: Establishing Running Times`: This is the code displayed in [Section 1.0](https://leetarxiv.substack.com/i/175344062/10-introduction) where we establish it takes 352 million years.
`Step2:Structured Scalar Multiplication`: This aligns with [Section 2.0](https://leetarxiv.substack.com/i/175344062/20-rewriting-scalar-multiplication) where we show different ways to perform MSB, LSB and meet-in-the-middle scalar multiplication.
`Step3:IntegerGenerators`: This aligns with [Section 2.3](https://leetarxiv.substack.com/i/175344062/20-rewriting-scalar-multiplication) where we demonstrate how to find primitive roots and modify our scalar multiplication to find successive points super fast.

