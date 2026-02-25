#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>
#include <ecm.h>


#define MAX_STREAM_S 10000

mpz_t STREAM_S[MAX_STREAM_S];
mpz_t STREAM_H[MAX_STREAM_S];
mpz_t STREAM_K[MAX_STREAM_S];
int STREAM_COUNT = 0;


void (*S_callback)(const mpz_t S, const mpz_t H, const mpz_t K) = NULL;

#define MAX_PRIME 100000

typedef struct {
    mpz_t p;
    unsigned long e;
} factor_t;

typedef struct {
    mpz_t H;
    mpz_t K;
    mpz_t S;
} solution_t;


int compute_S_from_current_factors(
    solution_t *sol,
    int max_solutions,
    factor_t *factors,
    int nf,
    const mpz_t A, const mpz_t B, const mpz_t C, const mpz_t N2,
    time_t start_time,
    int *timeout_flag
);

/* ---------- utilità base ---------- */

unsigned long v2(const mpz_t x) {
    if (mpz_cmp_ui(x, 0) == 0) return 0;
    return mpz_scan1(x, 0);
}

void add_factor(factor_t *factors, int *nf, int max_factors,
                const mpz_t p, unsigned long e,
                const mpz_t A, const mpz_t B, const mpz_t C, const mpz_t N2,
                time_t start_time, int *timeout_flag)
{
    for (int i = 0; i < *nf; i++) {
        if (mpz_cmp(factors[i].p, p) == 0) {
            factors[i].e += e;
            goto TRY_S;
        }
    }

    if (*nf >= max_factors) {
        fprintf(stderr, "Troppi fattori\n");
        exit(1);
    }

    mpz_init_set(factors[*nf].p, p);
    factors[*nf].e = e;
    (*nf)++;

TRY_S:
    if (*timeout_flag) return;

    solution_t sol[64];
    int count = compute_S_from_current_factors(
        sol, 64, factors, *nf, A, B, C, N2, start_time, timeout_flag
    );

    if (S_callback != NULL) {
        for (int i = 0; i < count; i++) {
            S_callback(sol[i].S, sol[i].H, sol[i].K);
        }
    }
}

/* ---------- generazione primi fino a 100000 ---------- */

int generate_primes(unsigned long **plist) {
    int size = MAX_PRIME + 1;
    char *sieve = calloc(size, 1);
    int count = 0;

    for (int i = 2; i * i <= MAX_PRIME; i++) {
        if (!sieve[i]) {
            for (int j = i * i; j <= MAX_PRIME; j += i)
                sieve[j] = 1;
        }
    }

    for (int i = 2; i <= MAX_PRIME; i++)
        if (!sieve[i]) count++;

    *plist = malloc(count * sizeof(unsigned long));
    int idx = 0;

    for (int i = 2; i <= MAX_PRIME; i++)
        if (!sieve[i]) (*plist)[idx++] = i;

    free(sieve);
    return count;
}

/* ---------- trial division fino a 100000 ---------- */

void trial_division(mpz_t n, factor_t *factors, int *nf, int max_factors,
                    const mpz_t A, const mpz_t B, const mpz_t C, const mpz_t N2,
                    time_t start_time, int *timeout_flag)
{
    unsigned long *primes;
    int pcount = generate_primes(&primes);

    mpz_t r, pz;
    mpz_inits(r, pz, NULL);

    for (int i = 0; i < pcount; i++) {
        if (*timeout_flag) break;
        if (time(NULL) - start_time >= 1800) {
            *timeout_flag = 1;
            break;
        }

        unsigned long p = primes[i];
        mpz_set_ui(pz, p);

        while (1) {
            mpz_mod_ui(r, n, p);
            if (mpz_cmp_ui(r, 0) != 0) break;

            mpz_divexact_ui(n, n, p);

            add_factor(factors, nf, max_factors, pz, 1,
                       A, B, C, N2, start_time, timeout_flag);

            if (*timeout_flag) break;
        }

        if (mpz_cmp_ui(n, 1) == 0) break;
    }

    mpz_clears(r, pz, NULL);
    free(primes);
}

/* ---------- Miller–Rabin ---------- */

int is_probable_prime(const mpz_t n, int reps) {
    return mpz_probab_prime_p(n, reps) > 0;
}

/* ---------- Pollard–Rho ---------- */

void pollard_rho(mpz_t factor, const mpz_t n,
                 time_t start_time, int *timeout_flag) {

    if (*timeout_flag) return;
    if (time(NULL) - start_time >= 1800) {
        *timeout_flag = 1;
        return;
    }

    if (mpz_even_p(n)) {
        mpz_set_ui(factor, 2);
        return;
    }

    gmp_randstate_t st;
    gmp_randinit_default(st);
    gmp_randseed_ui(st, (unsigned long)time(NULL));

    mpz_t x, y, c, d, tmp;
    mpz_inits(x, y, c, d, tmp, NULL);

    do {
        if (*timeout_flag) break;
        if (time(NULL) - start_time >= 1800) {
            *timeout_flag = 1;
            break;
        }

        mpz_urandomm(x, st, n);
        mpz_set(y, x);
        mpz_urandomm(c, st, n);
        if (mpz_cmp_ui(c, 0) == 0) mpz_set_ui(c, 1);
        mpz_set_ui(d, 1);

        while (mpz_cmp_ui(d, 1) == 0) {

            if (*timeout_flag) break;
            if (time(NULL) - start_time >= 1800) {
                *timeout_flag = 1;
                break;
            }

            mpz_mul(x, x, x);
            mpz_add(x, x, c);
            mpz_mod(x, x, n);

            mpz_mul(y, y, y);
            mpz_add(y, y, c);
            mpz_mod(y, y, n);

            mpz_mul(y, y, y);
            mpz_add(y, y, c);
            mpz_mod(y, y, n);

            mpz_sub(tmp, x, y);
            if (mpz_sgn(tmp) < 0) mpz_neg(tmp, tmp);
            mpz_gcd(d, tmp, n);

            if (mpz_cmp(d, n) == 0) break;
        }

    } while (mpz_cmp_ui(d, 1) == 0 || mpz_cmp(d, n) == 0);

    mpz_set(factor, d);

    mpz_clears(x, y, c, d, tmp, NULL);
    gmp_randclear(st);
}

/* ---------- fattorizzazione ibrida ECM + Pollard ---------- */

void factor_recursive(mpz_t n, factor_t *factors, int *nf, int max_factors,
                      const mpz_t A, const mpz_t B, const mpz_t C, const mpz_t N2,
                      time_t start_time, int *timeout_flag);

void factor_pollard(mpz_t n, factor_t *factors, int *nf, int max_factors,
                    const mpz_t A, const mpz_t B, const mpz_t C, const mpz_t N2,
                    time_t start_time, int *timeout_flag)
{
    if (*timeout_flag) return;
    if (time(NULL) - start_time >= 1800) {
        *timeout_flag = 1;
        return;
    }

    if (mpz_cmp_ui(n, 1) == 0) return;

    if (is_probable_prime(n, 40)) {
        add_factor(factors, nf, max_factors, n, 1,
                   A, B, C, N2, start_time, timeout_flag);
        return;
    }

    mpz_t f, q;
    mpz_inits(f, q, NULL);

    pollard_rho(f, n, start_time, timeout_flag);

    if (*timeout_flag) {
        mpz_clears(f, q, NULL);
        return;
    }

    if (mpz_cmp_ui(f, 0) == 0 ||
        mpz_cmp_ui(f, 1) == 0 ||
        mpz_cmp(f, n) == 0 ||
        !mpz_divisible_p(n, f))
    {
        add_factor(factors, nf, max_factors, n, 1,
                   A, B, C, N2, start_time, timeout_flag);
        mpz_clears(f, q, NULL);
        return;
    }

    mpz_divexact(q, n, f);

    factor_recursive(f, factors, nf, max_factors,
                     A, B, C, N2, start_time, timeout_flag);

    factor_recursive(q, factors, nf, max_factors,
                     A, B, C, N2, start_time, timeout_flag);

    mpz_clears(f, q, NULL);
}

void factor_ecm_single(mpz_t n, factor_t *factors, int *nf, int max_factors,
                       const mpz_t A, const mpz_t B, const mpz_t C, const mpz_t N2,
                       time_t start_time, int *timeout_flag)
{
    if (*timeout_flag) return;
    if (time(NULL) - start_time >= 1800) {
        *timeout_flag = 1;
        return;
    }

    if (mpz_cmp_ui(n, 1) == 0) return;

    if (is_probable_prime(n, 40)) {
        add_factor(factors, nf, max_factors, n, 1,
                   A, B, C, N2, start_time, timeout_flag);
        return;
    }

    ecm_params params;
    ecm_init(params);
    params->verbose = 0;

    mpz_t f, q;
    mpz_inits(f, q, NULL);

    int found = 0;
    double B1_list[] = { 50000.0, 250000.0, 1000000.0 };

    for (int bi = 0; bi < 3 && !found; bi++) {
        double B1 = B1_list[bi];

        for (int tries = 0; tries < 200 && !found; tries++) {

            if (*timeout_flag) break;
            if (time(NULL) - start_time >= 1800) {
                *timeout_flag = 1;
                break;
            }

            mpz_set_ui(f, 0);
            ecm_factor(f, n, B1, params);

            if (mpz_cmp_ui(f, 0) == 0) continue;
            if (mpz_cmp_ui(f, 1) == 0) continue;
            if (mpz_cmp(f, n) == 0) continue;
            if (!mpz_divisible_p(n, f)) continue;

            found = 1;
        }
    }

    if (!found) {
        mpz_clears(f, q, NULL);
        ecm_clear(params);
        factor_pollard(n, factors, nf, max_factors,
                       A, B, C, N2, start_time, timeout_flag);
        return;
    }

    mpz_divexact(q, n, f);

    factor_recursive(f, factors, nf, max_factors,
                     A, B, C, N2, start_time, timeout_flag);

    factor_recursive(q, factors, nf, max_factors,
                     A, B, C, N2, start_time, timeout_flag);

    mpz_clears(f, q, NULL);
    ecm_clear(params);
}

void factor_recursive(mpz_t n, factor_t *factors, int *nf, int max_factors,
                      const mpz_t A, const mpz_t B, const mpz_t C, const mpz_t N2,
                      time_t start_time, int *timeout_flag)
{
    if (*timeout_flag) return;
    if (time(NULL) - start_time >= 1800) {
        *timeout_flag = 1;
        return;
    }

    if (mpz_cmp_ui(n, 1) == 0) return;

    trial_division(n, factors, nf, max_factors,
                   A, B, C, N2, start_time, timeout_flag);

    if (*timeout_flag) return;
    if (mpz_cmp_ui(n, 1) == 0) return;

    if (is_probable_prime(n, 40)) {
        add_factor(factors, nf, max_factors, n, 1,
                   A, B, C, N2, start_time, timeout_flag);
        return;
    }

    factor_ecm_single(n, factors, nf, max_factors,
                      A, B, C, N2, start_time, timeout_flag);
}

/* ---------- test coppia (H,K) ---------- */

int test_pair(const mpz_t H, const mpz_t K,
              const mpz_t N2, const mpz_t A, const mpz_t BmodA,
              const mpz_t C2, unsigned long v2_N2) {
    (void)v2_N2;

    mpz_t S, tmp, halfS, diffHK, halfDiff;
    mpz_inits(S, tmp, halfS, diffHK, halfDiff, NULL);

    mpz_mul(tmp, H, K);
    if (mpz_cmp(tmp, N2) != 0) goto fail;

    if (mpz_cmp(H, K) == 0) goto fail;

    if (mpz_odd_p(H) != mpz_odd_p(K)) goto fail;

    mpz_add(S, H, K);
    if (mpz_odd_p(S)) goto fail;
    mpz_fdiv_q_2exp(halfS, S, 1);
    mpz_sub(tmp, halfS, BmodA);
    mpz_mod(tmp, tmp, A);
    if (mpz_cmp_ui(tmp, 0) != 0) goto fail;

    mpz_sub(diffHK, H, K);
    if (mpz_odd_p(diffHK)) goto fail;
    mpz_fdiv_q_2exp(halfDiff, diffHK, 1);
    mpz_mod(tmp, halfDiff, C2);
    if (mpz_cmp_ui(tmp, 0) != 0) goto fail;

    mpz_clears(S, tmp, halfS, diffHK, halfDiff, NULL);
    return 1;

fail:
    mpz_clears(S, tmp, halfS, diffHK, halfDiff, NULL);
    return 0;
}

/* ---------- ricorsione sui divisori di N^2 ---------- */

void search_divisors(
    int idx,
    int nf,
    factor_t *factors,
    const mpz_t N2,
    const mpz_t A,
    const mpz_t BmodA,
    const mpz_t C2,
    unsigned long v2_N2,
    mpz_t H,
    solution_t *sol,
    int *sol_count,
    int max_solutions,
    time_t start_time,
    int *timeout_flag
) {
    if (*timeout_flag) return;
    if (time(NULL) - start_time >= 1800) {
        *timeout_flag = 1;
        return;
    }

    if (idx == nf) {
        if (mpz_cmp_ui(H, 0) == 0) return;
        if (!mpz_divisible_p(N2, H)) return;

        mpz_t K;
        mpz_init(K);
        mpz_divexact(K, N2, H);

        if (test_pair(H, K, N2, A, BmodA, C2, v2_N2)) {
            if (*sol_count < max_solutions) {
                mpz_init_set(sol[*sol_count].H, H);
                mpz_init_set(sol[*sol_count].K, K);
                mpz_init(sol[*sol_count].S);
                mpz_add(sol[*sol_count].S, H, K);
                (*sol_count)++;
            }
        }

        mpz_clear(K);
        return;
    }

    mpz_t p_pow;
    mpz_init_set_ui(p_pow, 1);

    unsigned long max_exp = 2 * factors[idx].e;

    for (unsigned long e = 0; e <= max_exp; e++) {

        if (*timeout_flag) break;
        if (time(NULL) - start_time >= 1800) {
            *timeout_flag = 1;
            break;
        }

        mpz_t H_next;
        mpz_init(H_next);
        mpz_mul(H_next, H, p_pow);

        search_divisors(
            idx + 1, nf, factors,
            N2, A, BmodA, C2,
            v2_N2,
            H_next,
            sol, sol_count, max_solutions,
            start_time,
            timeout_flag
        );

        mpz_clear(H_next);
        mpz_mul(p_pow, p_pow, factors[idx].p);
    }

    mpz_clear(p_pow);
}

/* ---------- wrapper per S dai fattori correnti ---------- */

int compute_S_from_current_factors(
    solution_t *sol,
    int max_solutions,
    factor_t *factors,
    int nf,
    const mpz_t A, const mpz_t B, const mpz_t C, const mpz_t N2,
    time_t start_time,
    int *timeout_flag
) {
    mpz_t H0, C2, BmodA;
    mpz_inits(H0, C2, BmodA, NULL);

    mpz_set_ui(H0, 1);
    mpz_fdiv_q_2exp(C2, C, 1);
    mpz_mod(BmodA, B, A);

    unsigned long v2_N2 = v2(N2);

    int sol_count = 0;

    search_divisors(
        0, nf, factors,
        N2, A, BmodA, C2,
        v2_N2,
        H0,
        sol, &sol_count, max_solutions,
        start_time,
        timeout_flag
    );

    mpz_clears(H0, C2, BmodA, NULL);

    return sol_count;
}

/* ---------- wrapper principale ---------- */

int find_all_S(
    solution_t *sol,
    int max_solutions,
    const mpz_t A, const mpz_t B, const mpz_t C, const mpz_t N2
) {
    mpz_t N, tmp;
    mpz_inits(N, tmp, NULL);

    mpz_sqrt(N, N2);
    mpz_mul(tmp, N, N);
    if (mpz_cmp(tmp, N2) != 0) {
        fprintf(stderr, "N^2 non è un quadrato perfetto\n");
        mpz_clears(N, tmp, NULL);
        return 0;
    }

    const int MAX_FACTORS = 256;
    factor_t factors[MAX_FACTORS];
    int nf = 0;

    time_t start_time = time(NULL);
    int timeout_flag = 0;

    factor_recursive(N, factors, &nf, MAX_FACTORS,
                 A, B, C, N2, start_time, &timeout_flag);

    if (timeout_flag) {
        for (int i = 0; i < nf; i++) mpz_clear(factors[i].p);
        mpz_clears(N, tmp, NULL);
        return 0;
    }

    mpz_t H0, C2, BmodA;
    mpz_inits(H0, C2, BmodA, NULL);
    mpz_set_ui(H0, 1);
    mpz_fdiv_q_2exp(C2, C, 1);
    mpz_mod(BmodA, B, A);

    unsigned long v2_N2 = v2(N2);

    int sol_count = 0;

    search_divisors(
        0, nf, factors,
        N2, A, BmodA, C2,
        v2_N2,
        H0,
        sol, &sol_count, max_solutions,
        start_time,
        &timeout_flag
    );

    for (int i = 0; i < nf; i++) mpz_clear(factors[i].p);

    mpz_clears(N, tmp, H0, C2, BmodA, NULL);

    if (timeout_flag) return 0;

    return sol_count;
}

/* ---------- main di esempio ---------- */
/* Scrivi tu il main in base a come vuoi usare find_all_S */
/* ---------- main ---------- */


void on_new_S(const mpz_t S, const mpz_t H, const mpz_t K) {
    if (STREAM_COUNT >= MAX_STREAM_S) {
        printf("[MAIN] Buffer S pieno!\n");
        return;
    }

    mpz_set(STREAM_S[STREAM_COUNT], S);
    mpz_set(STREAM_H[STREAM_COUNT], H);
    mpz_set(STREAM_K[STREAM_COUNT], K);

    gmp_printf("[MAIN] S salvato: S=%Zd  (H=%Zd, K=%Zd)\n",
               STREAM_S[STREAM_COUNT],
               STREAM_H[STREAM_COUNT],
               STREAM_K[STREAM_COUNT]);

    STREAM_COUNT++;
}

int main() {
    S_callback = on_new_S;

    for (int i = 0; i < MAX_STREAM_S; i++) {
        mpz_init(STREAM_S[i]);
        mpz_init(STREAM_H[i]);
        mpz_init(STREAM_K[i]);
    }

    mpz_t A, B, C, N2, M, temp,temp2,temp3,p,P,a,M2,a2;
    mpz_inits(A, B, C, N2,M, temp,temp2,temp3,p,P,a,M2,a2, NULL);
int count =0;
  int i = 0;
  mpz_set_str(a, "3", 10);
    mpz_set_str(M, "187", 10);
    gmp_printf("RSA= %Zd\n", M);
    while (1) {
        STREAM_COUNT = 0;

        printf("\n--- Nuova chiamata a find_all_S ---\n");
        
  mpz_mul_ui(A,a,4);
      mpz_mul(M2,M,M);
      mpz_mul(a2,a,a);
      mpz_mul(temp,M2,a2);
    
      mpz_sub_ui(temp2,temp,1);/////
      mpz_div_ui(temp2,temp2,2);
      mpz_mul(N2,temp2,temp2);
      mpz_add_ui(temp2,temp,1);//////
      mpz_div_ui(temp2,temp2,2);
      mpz_add(B,temp2,a);
      mpz_sqrt(temp2,temp);
      mpz_mul_ui(C,temp2,2);
      gmp_printf("a: %Zd\n", a);
gmp_printf("A: %Zd\n", A);
 gmp_printf("B: %Zd\n", B);
 gmp_printf("C: %Zd\n", C);
 gmp_printf("N2: %Zd\n", N2);
        printf("\nCalcolo in corso...\n");
        printf("(Gli S verranno salvati in tempo reale nel buffer STREAM_S[])\n\n");

        solution_t sol[256];
        count = find_all_S(sol, 256, A, B, C, N2);

        printf("\n--- FINE CALCOLO ---\n");

        printf("Sono arrivati %d S in tempo reale.\n", STREAM_COUNT);

        for (int i = 0; i < STREAM_COUNT; i++) {
            gmp_printf("STREAM[%d]  S=%Zd  H=%Zd  K=%Zd\n",
                       i,
                       STREAM_S[i],
                       STREAM_H[i],
                       STREAM_K[i]);

mpz_div_ui(temp,STREAM_S[i],2);
	mpz_sub(temp,temp,B);
	mpz_div(temp,temp,A);gmp_printf("m trovato: %Zd\n", temp);
	mpz_mul_ui(temp,temp,4);
	mpz_add_ui(P,temp,1);
	mpz_mod(temp,P,M);
	if(mpz_cmp_ui(temp,0)==0){
	  mpz_div(P,P,M);
	}

	mpz_gcd(p,P,M);
	
	if(mpz_cmp_ui(p,1)==0 || mpz_cmp(p,M)==0){
	  printf("algoritmo non funzionante\n");
	}else{
	  gmp_printf("p= %Zd\n", p);
	  return 0;
	}

	    
        }

        if (count > 0) {
            printf("\nSoluzioni finali trovate:\n");
            for (i = 0; i < count; i++) {
                gmp_printf("S finale %d: %Zd\n", i+1, sol[i].S);
                mpz_clear(sol[i].H);
                mpz_clear(sol[i].K);
                mpz_clear(sol[i].S);
            }
        }
	i=0;
	count=0;
mpz_add_ui(a,a,4);
	
    }

    mpz_clears(A, B, C, N2, NULL);
    return 0;
}




