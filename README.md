# lepore_factorization_nr_00
integer factorization

I have devised a factorization algorithm that:

Given M, a number to be factored,

From M, it generates an N of the order of magnitude of M^2.

It starts the complete factorization of N using trial division, ECM, and Pollard–Rho.

Each time it finds a new factor, it updates the current factorization.

For each updated factor, it checks whether a given value S is valid.

If S is valid, it is sent in real time to the main program via a callback.

The main program checks whether S produces a factorization of M; otherwise, it generates another N and restarts.

The source is written in C and uses the GMP and ECM libraries.

ECM uses B1 = 50000, 250000, 1000000, with 200 curves for each B1, for a total of 600 ECM curves, with standard GMP-ECM parameters and a timeout of 1800 seconds.

Optimal case:
N is factored in a few seconds
S is valid
S produces a factorization of M
