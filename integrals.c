/*
 * integrals.c -- analytic molecular integrals over Cartesian Gaussians.
 *
 * Compiled to a shared library (integrals.dll) and called from atoms.py
 * through ctypes. Everything is in atomic units (Bohr, Hartree).
 *
 * A primitive Cartesian Gaussian centred at A is
 *
 *     g(r) = N (x-Ax)^l (y-Ay)^m (z-Az)^n exp(-alpha |r-A|^2)
 *
 * and a basis function is a contraction (fixed linear combination) of
 * primitives. The integrals are evaluated with the Taketa-Huzinaga-O-ohata
 * scheme as presented in Cook, "Handbook of Computational Chemistry" (1998):
 *   - f_j(l, m, a, b): coefficient of x^j in (x+a)^l (x+b)^m   (f_j below)
 *   - 1-D overlap integrals                                     (S_integral)
 *   - A_{l,r,i} coefficients for the nuclear attraction         (calculate_A)
 *   - B_{l,l',r,r',i} coefficients for the two-electron
 *     repulsion integrals                                       (calculate_B)
 *   - Boys function F_n(x), which ties the three Cartesian
 *     directions together                                       (boys_recurrence)
 *
 * Notation used in the comments:
 *   A, B, C, D   centres of the (up to) four Gaussians
 *   P, Q         Gaussian product centres of the bra (A,B) and ket (C,D) pairs
 *   gamma1       alpha_a + alpha_b          gamma2   alpha_c + alpha_d
 *   delta        1/(4*gamma1) + 1/(4*gamma2)
 *
 * NOTE: repulsion_integral4 uses variable-length arrays (C99), so build with
 * gcc/clang (MSVC does not support VLAs).
 *
 * Build:  gcc -O3 -shared -o integrals.dll integrals.c      (Windows, MinGW-w64)
 *         gcc -O3 -shared -fPIC -o integrals.so integrals.c -lm (Linux)
 *         gcc -O3 -shared -fPIC -o integrals.dylib integrals.c -lm (macOS)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Row-major index of element (a,b,c,d) of a dim x dim x dim x dim tensor.
 * Must match the C-contiguous layout of the numpy array allocated in atoms.py. */
#define IDX4D(a, b, c, d, dim) ((a)*(dim)*(dim)*(dim) + (b)*(dim)*(dim) + (c)*(dim) + (d)) 
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Lookup table: FACT[n] = n!  (exact in double up to 20!). */
static const double FACT[] = {
    1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0,
    3628800.0, 39916800.0, 479001600.0, 6227020800.0, 87178291200.0,
    1307674368000.0, 20922789888000.0, 355687428096000.0,
    6402373705728000.0, 121645100408832000.0, 2432902008176640000.0
};

/* Lookup table: FACT2[n] = n!!  (double factorial), n = 0..20. */
static const double FACT2[] = {
    1.0, 1.0, 2.0, 3.0, 8.0, 15.0, 48.0, 105.0, 384.0, 945.0, 3840.0,
    10395.0, 46080.0, 135135.0, 645120.0, 2027025.0,
    10321920.0, 34459425.0, 185794560.0,
    654729075.0, 3715891200.0
};

/* Binomial coefficient C(n, k). Returns 0 when k is outside [0, n]
 * (the sums below rely on this) and also for n >= 21, which is beyond the
 * FACT table. */
static inline double comb(int n, int k) {
    if (k < 0 || k > n || n >= 21) return 0.0;
    return FACT[n] / (FACT[k] * FACT[n - k]);
}

/* Double factorial n!!. By convention (-1)!! = 0!! = 1, which is what the
 * normalisation and overlap formulas need (e.g. 2*l-1 with l = 0).
 * Uses the table for n <= 20 and a plain loop above that. */
static inline double fact2(int n) {
    if (n <= 1) return 1.0;
    if (n <= 20) return FACT2[n];
 
    double val = 1.0;
    for (int i = n; i > 0; i -= 2) {
        val *= i;
    }
    return val;
}

/* base^exp for an integer exponent (repeated multiplication, cheaper than pow).
 * Negative exponents return 1/base^|exp|. exp == 0 always gives 1, even for
 * base == 0 (0^0 = 1, needed by the polynomial expansions below); any other
 * power of a zero base returns 0. */
static inline double int_pow(double base, int exp) {
    if (exp == 0) return 1.0;
    if (base == 0.0) return 0.0;
    

    int abs_exp = (exp > 0) ? exp : -exp;
    double res = 1.0;
    
    for (int i = 0; i < abs_exp; i++) {
        res *= base;
    }
    

    return (exp > 0) ? res : (1.0 / res);
}

/* f_j(l, m, a, b): coefficient of x^j in the polynomial (x + a)^l (x + b)^m,
 *
 *     f_j = sum_k C(l,k) C(m,j-k) a^(l-k) b^(m+k-j)
 *
 * It comes from expanding the product of two Cartesian Gaussians' polynomial
 * parts around their product centre P, where a = P-A and b = P-B.
 * The sum runs only over k for which both binomials are non-zero. */
static inline double f_j(int l, int m, double a, double b, int j) {
    double value = 0.0;
    int start = (j - m > 0) ? j - m : 0;
    int end = (j < l) ? j : l;
    
    for (int k = start; k <= end; k++) {
        value += comb(l, k) * comb(m, j - k) * int_pow(a, l - k) * int_pow(b, m + k - j);
    }
    return value;
}

/* theta(j, l, m, a, b, r, gamma): one-centre expansion coefficient used by
 * the repulsion integrals,
 *
 *     theta = f_j(l, m, a, b) * j! * gamma^(r-j) / ( r! * (j-2r)! )
 *
 * with a = P-A, b = P-B (see calculate_B). Needs j >= 2r. */
static inline double theta(int j, int l, int m, double a, double b, int r, double gamma) {
    return f_j(l, m, a, b, j) * FACT[j] * int_pow(gamma, r - j) / (FACT[r] * FACT[j - 2 * r]);
}

/* One-dimensional overlap integral between two Cartesian Gaussians that have
 * angular momenta l and m along one axis (the 3-D overlap is the product of
 * three of these times the exponential prefactor, built in atoms.py):
 *
 *     S = sqrt(pi/g) * sum_j f_{2j}(l, m, Pa, Pb) (2j-1)!! / (2g)^j
 *
 * Pa = P-A and Pb = P-B along that axis, g = alpha_a + alpha_b.
 * Only even powers survive the integration, hence the 2j. */
double S_integral(int l, int m, double Pa, double Pb, double g) {
    double value = 0.0;
    int max_j = (l + m) / 2;
    double den = 1.0; 
    
    for (int j = 0; j <= max_j; j++) {
        double num = f_j(l, m, Pa, Pb, 2 * j) * fact2(2 * j - 1);
        value += num / den;
        den *= (2.0 * g);
    }
    
    value *= sqrt(M_PI / g);
    
    return value;
}

/* Normalisation constant of a single primitive Cartesian Gaussian,
 *
 *     N = (2 alpha / pi)^(3/4) * sqrt( (4 alpha)^(l+m+n) / ((2l-1)!! (2m-1)!! (2n-1)!!) )
 *
 * so that <g|g> = 1 for that primitive alone. The contracted function is
 * re-normalised afterwards in atoms.py (Matrix.__init__). */
double normalize_N(double alpha, int l, int m, int n) {
    double num = pow(4.0 * alpha, l + m + n);
    double den = fact2(2 * l - 1) * fact2(2 * m - 1) * fact2(2 * n - 1);
    
    double term1 = sqrt(num / den);
    double term2 = pow(2.0 * alpha / M_PI, 0.75);
    
    return term1 * term2;
}

/* Boys function F_n(x) = int_0^1 t^(2n) exp(-x t^2) dt, for all orders
 * 0..n at once. out_array must have room for n+1 doubles.
 *
 * Strategy (chosen for numerical stability):
 *   x ~ 0     : limit F_i(0) = 1/(2i+1).
 *   x <= 12   : power series for the highest order F_n, then DOWNWARD
 *               recursion F_{i-1} = (2x F_i + e^-x) / (2i-1), which is stable
 *               for small/moderate x.
 *   x  > 12   : F_0 from the error function, then UPWARD recursion
 *               F_{i+1} = ((2i+1) F_i - e^-x) / (2x), stable for large x.
 */
void boys_recurrence(int n, double x, double* out_array) {
    if (x < 1e-15) {
        for (int i = 0; i <= n; i++) {
            out_array[i] = 1.0 / (2.0 * i + 1.0);
        }
        return;
    }

    double exp_x = exp(-x);

    if (x <= 12.0) {
        /* Series: F_n(x) = sum_k (-x)^k / k! * 1/(2n+2k+1), accumulated
         * term by term using the ratio between consecutive terms. */
        double term = 1.0 / (2.0 * n + 1.0);
        double f_val = term;
        
        for (int k = 1; k < 500; k++) {
            term *= (-x) / k * (2.0 * n + 2.0 * k - 1.0) / (2.0 * n + 2.0 * k + 1.0);
            f_val += term;
            if (fabs(term) < 1e-17) break;
        }
        out_array[n] = f_val;
        
        for (int i = n; i > 0; i--) {
            out_array[i-1] = (2.0 * x * out_array[i] + exp_x) / (2.0 * i - 1.0);
        }
        
    } else {
        out_array[0] = sqrt(M_PI) * erf(sqrt(x)) / (2.0 * sqrt(x));
        
        for (int i = 0; i < n; i++) {
            out_array[i+1] = ((2.0 * i + 1.0) * out_array[i] - exp_x) / (2.0 * x);
        }
    }
}

/* Nuclear-attraction coefficients A_{l,r,i} along ONE Cartesian axis.
 *
 *   l1, l2   angular momenta of the two Gaussians on this axis
 *   A, B     their centres on this axis;  C = nucleus position on this axis
 *   P        Gaussian product centre on this axis
 *   eps      1/(4*gamma), gamma = alpha_a + alpha_b
 *
 * The result is accumulated (+=) into out_array, which must have l1+l2+1
 * entries and be zeroed by the caller. Entry u holds the sum of all terms
 * with u = index1 - 2*index2 - index3, i.e. the weight that multiplies the
 * Boys function F_(u + contributions from the other two axes). The 3-D
 * integral is  sum_{ux,uy,uz} A_x[ux] A_y[uy] A_z[uz] F_{ux+uy+uz}  (see
 * Matrix._nuclear_attraction in atoms.py).
 *
 *   A_{l,r,i} = (-1)^(l+i) f_l(l1,l2,P-A,P-B) l! (P-C)^(l-2r-2i) eps^(r+i)
 *               / ( r! i! (l-2r-2i)! )
 */
void calculate_A(int l1, int l2, double A, double B, double C, double P, double eps, double* out_array) {
    for (int index1 = 0; index1 <= l1 + l2; index1++) {              /* l */
        for (int index2 = 0; index2 <= index1 / 2; index2++) {       /* r */
            for (int index3 = 0; index3 <= (index1 - 2 * index2) / 2; index3++) {   /* i */
                
                double sign = ((index1 + index3) % 2 == 0) ? 1.0 : -1.0;
                
                double num = sign * f_j(l1, l2, P - A, P - B, index1) 
                             * FACT[index1] 
                             * int_pow(P - C, index1 - 2 * (index2 + index3)) 
                             * int_pow(eps, index2 + index3);
                             
                double den = FACT[index2] * FACT[index3] * FACT[index1 - 2 * (index2 + index3)];
                
                int target_idx = index1 - 2 * index2 - index3;
                out_array[target_idx] += num / den;
            }
        }
    }
}

/* Two-electron repulsion coefficients B_{l,l',r,r',i} along ONE Cartesian
 * axis, for the four Gaussians (A,B | C,D).
 *
 *   l1..l4         angular momenta on this axis (bra: l1,l2; ket: l3,l4)
 *   A,B,C,D        centres on this axis
 *   P, Q           product centres of the bra and ket pairs on this axis
 *   gamma1/gamma2  alpha_a+alpha_b  /  alpha_c+alpha_d
 *   delta          1/(4 gamma1) + 1/(4 gamma2)
 *
 * The bra and ket Gaussian products are expanded with theta (t1, t2) and
 * recombined through delta and (Q-P). As in calculate_A, the result is
 * accumulated (+=) into out_array, which must have l1+l2+l3+l4+1 entries and
 * be zeroed by the caller; entry u is the weight multiplying F_(u+...).
 *
 * Index map:  index1/index2 = (l, r) of the bra,  index4/index5 = (l', r')
 * of the ket,  index3 = i (the extra summation of the combined expansion).
 * The `continue` guarantees index1+index4-2(index3+index2+index5) >= 0,
 * i.e. that the factorial / power arguments below are never negative. */
void calculate_B(int l1, int l2, int l3, int l4, 
                 double A, double B, double C, double D, 
                 double P, double Q, 
                 double gamma1, double gamma2, double delta, 
                 double* out_array) {
                     
    for (int index1 = 0; index1 <= l1 + l2; index1++) {
        for (int index2 = 0; index2 <= index1 / 2; index2++) {
            for (int index3 = 0; index3 <= (index1 + l3 + l4 - 2 * index2) / 2; index3++) {
                for (int index4 = 0; index4 <= l3 + l4; index4++) {
                    for (int index5 = 0; index5 <= index4 / 2; index5++) {
                        
                        if (index3 > (index1 + index4 - 2 * (index2 + index5)) / 2) {
                            continue;
                        }
                        
                        double sign = ((index4 + index3) % 2 == 0) ? 1.0 : -1.0;
                        
                        double t1 = theta(index1, l1, l2, P - A, P - B, index2, gamma1);
                        double t2 = theta(index4, l3, l4, Q - C, Q - D, index5, gamma2);
                        
                        int fact_idx1 = index1 + index4 - 2 * (index2 + index5);
                        int fact_idx2 = index1 + index4 - 2 * (index3 + index2 + index5);
                        
                        double p1 = int_pow(2.0, 2 * (index2 + index5) - 2 * (index1 + index4));
                        double p2 = int_pow(delta, 2 * (index2 + index5) + index3 - (index1 + index4));
                        double p3 = int_pow(Q - P, fact_idx2);
                        
                        double term = sign * t1 * t2 * p1 * p2 * FACT[fact_idx1] * p3 
                                      / (FACT[index3] * FACT[fact_idx2]);
                        
                        int target_idx = fact_idx1 - index3;
                        out_array[target_idx] += term;
                    }
                }
            }
        }
    }
}

/* One contracted two-electron repulsion integral (mu nu | lambda sigma) in
 * chemists' notation, WITHOUT the contracted-function normalisation factors
 * (those are applied by repulsion_matrix).
 *
 *   lmn[12]        (l,m,n) of the four basis functions, flattened as
 *                  [l1 m1 n1  l2 m2 n2  l3 m3 n3  l4 m4 n4]
 *   A, B, C, D     centres, each a 3-vector (x, y, z) in Bohr
 *   num1..num4     number of primitives of each contracted function
 *   expK / coefK   exponents and coefficients of function K's primitives
 *                  (coefficients already include the primitive normalisation N)
 *
 * Loops over every combination of primitives (i,j,k,l) and, for each:
 *   1. builds the Gaussian product centres P (bra) and Q (ket);
 *   2. fills B_x, B_y, B_z with calculate_B;
 *   3. evaluates the Boys functions F_0..F_vmax at x = |P-Q|^2 / (4 delta);
 *   4. sums  B_x[ix] B_y[iy] B_z[iz] F_{ix+iy+iz}  over all (ix,iy,iz);
 *   5. multiplies by the prefactor omega
 *          omega = 2 pi^(5/2) / (gamma1 gamma2 sqrt(gamma1+gamma2))
 *                  * exp( -alpha_a alpha_b |AB|^2 / gamma1
 *                         -alpha_c alpha_d |CD|^2 / gamma2 )
 *      and by the four contraction coefficients.
 * Primitives with a zero coefficient are skipped. */
double repulsion_integral4(
    int* lmn, double* A, double* B, double* C, double* D,
    int num1, int num2, int num3, int num4,
    double* exp1, double* coef1, double* exp2, double* coef2,
    double* exp3, double* coef3, double* exp4, double* coef4
) {
    double total_value = 0.0;
    
    /* Squared A-B and C-D distances (enter the Gaussian product prefactor). */
    double ab2 = (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) + (A[2]-B[2])*(A[2]-B[2]);
    double cd2 = (C[0]-D[0])*(C[0]-D[0]) + (C[1]-D[1])*(C[1]-D[1]) + (C[2]-D[2])*(C[2]-D[2]);

    /* Total angular momentum along each axis, and overall (highest Boys order). */
    int lx_sum = lmn[0] + lmn[3] + lmn[6] + lmn[9];
    int ly_sum = lmn[1] + lmn[4] + lmn[7] + lmn[10];
    int lz_sum = lmn[2] + lmn[5] + lmn[8] + lmn[11];
    int vmax = lx_sum + ly_sum + lz_sum;

    /* Per-axis B coefficients and Boys values (stack VLAs, reused every loop). */
    double B_x[lx_sum + 1];
    double B_y[ly_sum + 1];
    double B_z[lz_sum + 1];
    double boys[vmax + 1];

    for(int i = 0; i < num1; i++) {
        if (coef1[i] == 0.0) continue; 
        
        for(int j = 0; j < num2; j++) {
            if (coef2[j] == 0.0) continue;
            
            /* Bra: Gaussian product of primitives i (on A) and j (on B). */
            double gamma1 = exp1[i] + exp2[j];
            double Px = (exp1[i]*A[0] + exp2[j]*B[0]) / gamma1;
            double Py = (exp1[i]*A[1] + exp2[j]*B[1]) / gamma1;
            double Pz = (exp1[i]*A[2] + exp2[j]*B[2]) / gamma1;

            for(int k = 0; k < num3; k++) {
                if (coef3[k] == 0.0) continue;
                
                for(int l = 0; l < num4; l++) {
                    if (coef4[l] == 0.0) continue;
                    
                    /* Ket: Gaussian product of primitives k (on C) and l (on D). */
                    double gamma2 = exp3[k] + exp4[l];
                    double delta = 1.0/(4.0*gamma1) + 1.0/(4.0*gamma2);

                    double Qx = (exp3[k]*C[0] + exp4[l]*D[0]) / gamma2;
                    double Qy = (exp3[k]*C[1] + exp4[l]*D[1]) / gamma2;
                    double Qz = (exp3[k]*C[2] + exp4[l]*D[2]) / gamma2;

                    /* |P-Q|^2, argument of the Boys function (divided by 4*delta below). */
                    double p2 = (Px-Qx)*(Px-Qx) + (Py-Qy)*(Py-Qy) + (Pz-Qz)*(Pz-Qz);

                    /* calculate_B accumulates, so the buffers must be cleared first. */
                    for(int idx = 0; idx <= lx_sum; idx++) B_x[idx] = 0.0;
                    for(int idx = 0; idx <= ly_sum; idx++) B_y[idx] = 0.0;
                    for(int idx = 0; idx <= lz_sum; idx++) B_z[idx] = 0.0;

                    calculate_B(lmn[0], lmn[3], lmn[6], lmn[9], A[0], B[0], C[0], D[0], Px, Qx, gamma1, gamma2, delta, B_x);
                    calculate_B(lmn[1], lmn[4], lmn[7], lmn[10], A[1], B[1], C[1], D[1], Py, Qy, gamma1, gamma2, delta, B_y);
                    calculate_B(lmn[2], lmn[5], lmn[8], lmn[11], A[2], B[2], C[2], D[2], Pz, Qz, gamma1, gamma2, delta, B_z);

                    boys_recurrence(vmax, p2/(4.0*delta), boys);

                    double omega = 2.0 * M_PI * M_PI * sqrt(M_PI/(gamma1 + gamma2)) * exp(-(exp1[i]*exp2[j]*ab2/gamma1 + exp3[k]*exp4[l]*cd2/gamma2)) / (gamma1*gamma2);

                    /* Combine the three axes: sum_{ix,iy,iz} Bx By Bz F_{ix+iy+iz}. */
                    double primitive_val = 0.0;
                    for (int ix = 0; ix <= lx_sum; ix++) {
                        for (int iy = 0; iy <= ly_sum; iy++) {
                            for (int iz = 0; iz <= lz_sum; iz++) {
                                primitive_val += B_x[ix] * B_y[iy] * B_z[iz] * boys[ix + iy + iz];
                            }
                        }
                    }
                    total_value += primitive_val * omega * coef1[i] * coef2[j] * coef3[k] * coef4[l];
                }
            }
        }
    }
    return total_value;
}

/* Fills the full dim^4 tensor (mu nu | lambda sigma) of two-electron repulsion
 * integrals over contracted basis functions.
 *
 *   dim      number of basis functions
 *   centers  dim x 3 array: centre (x,y,z) of each basis function
 *   lmns     dim x 3 array: (l,m,n) of each basis function
 *   lens     number of primitives of each basis function
 *   offsets  where each function's primitives start inside exps / coefs
 *   exps     all primitive exponents, concatenated function by function
 *   coefs    matching coefficients (contraction coefficient * primitive N)
 *   normalize  per-function factor 1/sqrt(<chi|chi>) of the contracted function
 *   matrix   output, dim^4 doubles, row-major (see IDX4D); must be zero-filled
 *
 * Only the symmetry-unique integrals are computed: i<=j, k<=l and the pair
 * (i,j) not greater than the pair (k,l) -- about 1/8 of dim^4. Each value is
 * then copied to the 8 equivalent positions:
 *   (ij|kl) = (ji|kl) = (ij|lk) = (ji|lk) = (kl|ij) = (lk|ij) = (kl|ji) = (lk|ji).
 *
 * This is the O(N^4) bottleneck of the whole program, hence its implementation
 * in C. (The return type is declared double but nothing is returned; atoms.py
 * sets restype = None, so the value is never read. `void` would be the clean
 * declaration.) */
double repulsion_matrix(
    int dim, double* centers, int* lmns, int* lens,         
    int* offsets, double* exps, double* coefs, 
    double* normalize, double* matrix  
) {
    for(int i = 0; i < dim; i++) {
        for(int j = i; j < dim; j++) {
            for(int k = 0; k < dim; k++) {
                for(int l = k; l < dim; l++) {
                    /* Skip pairs (k,l) that come before (i,j): the 8-fold
                     * symmetry already filled them when they were visited. */
                    if (k < i || (k == i && l < j)){
                        continue;
                    }

                    int lmn12[12] = {
                        lmns[i*3], lmns[i*3+1], lmns[i*3+2],
                        lmns[j*3], lmns[j*3+1], lmns[j*3+2],
                        lmns[k*3], lmns[k*3+1], lmns[k*3+2],
                        lmns[l*3], lmns[l*3+1], lmns[l*3+2]
                    };

                    double integral = repulsion_integral4(
                        lmn12,
                        &centers[i*3], &centers[j*3], &centers[k*3], &centers[l*3],
                        lens[i], lens[j], lens[k], lens[l],
                        &exps[offsets[i]], &coefs[offsets[i]],
                        &exps[offsets[j]], &coefs[offsets[j]],
                        &exps[offsets[k]], &coefs[offsets[k]],
                        &exps[offsets[l]], &coefs[offsets[l]]
                    );

                    /* Apply the contracted-function normalisation. */
                    double value = integral * normalize[i] * normalize[j] * normalize[k] * normalize[l];

                    /* Scatter to the 8 symmetry-equivalent positions. */
                    matrix[IDX4D(i, j, k, l, dim)] = value;
                    matrix[IDX4D(i, j, l, k, dim)] = value;
                    matrix[IDX4D(j, i, k, l, dim)] = value;
                    matrix[IDX4D(j, i, l, k, dim)] = value;
                    matrix[IDX4D(k, l, i, j, dim)] = value;
                    matrix[IDX4D(k, l, j, i, dim)] = value;
                    matrix[IDX4D(l, k, i, j, dim)] = value;
                    matrix[IDX4D(l, k, j, i, dim)] = value;
                }
            }
        }
    }
}