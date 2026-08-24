#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define IDX4D(a, b, c, d, dim) ((a)*(dim)*(dim)*(dim) + (b)*(dim)*(dim) + (c)*(dim) + (d)) 
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const double FACT[] = {
    1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0,
    3628800.0, 39916800.0, 479001600.0, 6227020800.0, 87178291200.0,
    1307674368000.0, 20922789888000.0, 355687428096000.0,
    6402373705728000.0, 121645100408832000.0, 2432902008176640000.0
};

static const double FACT2[] = {
    1.0, 1.0, 2.0, 3.0, 8.0, 15.0, 48.0, 105.0, 384.0, 945.0, 3840.0,
    10395.0, 46080.0, 135135.0, 645120.0, 2027025.0,
    10321920.0, 34459425.0, 185794560.0,
    654729075.0, 3715891200.0
};

static inline double comb(int n, int k) {
    if (k < 0 || k > n || n >= 21) return 0.0;
    return FACT[n] / (FACT[k] * FACT[n - k]);
}

static inline double fact2(int n) {
    if (n <= 1) return 1.0;
    if (n <= 20) return FACT2[n];
 
    double val = 1.0;
    for (int i = n; i > 0; i -= 2) {
        val *= i;
    }
    return val;
}

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

static inline double f_j(int l, int m, double a, double b, int j) {
    double value = 0.0;
    int start = (j - m > 0) ? j - m : 0;
    int end = (j < l) ? j : l;
    
    for (int k = start; k <= end; k++) {
        value += comb(l, k) * comb(m, j - k) * int_pow(a, l - k) * int_pow(b, m + k - j);
    }
    return value;
}

static inline double theta(int j, int l, int m, double a, double b, int r, double gamma) {
    return f_j(l, m, a, b, j) * FACT[j] * int_pow(gamma, r - j) / (FACT[r] * FACT[j - 2 * r]);
}

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

double normalize_N(double alpha, int l, int m, int n) {
    double num = pow(4.0 * alpha, l + m + n);
    double den = fact2(2 * l - 1) * fact2(2 * m - 1) * fact2(2 * n - 1);
    
    double term1 = sqrt(num / den);
    double term2 = pow(2.0 * alpha / M_PI, 0.75);
    
    return term1 * term2;
}

void boys_recurrence(int n, double x, double* out_array) {
    if (x < 1e-15) {
        for (int i = 0; i <= n; i++) {
            out_array[i] = 1.0 / (2.0 * i + 1.0);
        }
        return;
    }

    double exp_x = exp(-x);

    if (x <= 12.0) {
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

void calculate_A(int l1, int l2, double A, double B, double C, double P, double eps, double* out_array) {
    for (int index1 = 0; index1 <= l1 + l2; index1++) {
        for (int index2 = 0; index2 <= index1 / 2; index2++) {
            for (int index3 = 0; index3 <= (index1 - 2 * index2) / 2; index3++) {
                
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

double repulsion_integral4(
    int* lmn, double* A, double* B, double* C, double* D,
    int num1, int num2, int num3, int num4,
    double* exp1, double* coef1, double* exp2, double* coef2,
    double* exp3, double* coef3, double* exp4, double* coef4
) {
    double total_value = 0.0;
    
    double ab2 = (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) + (A[2]-B[2])*(A[2]-B[2]);
    double cd2 = (C[0]-D[0])*(C[0]-D[0]) + (C[1]-D[1])*(C[1]-D[1]) + (C[2]-D[2])*(C[2]-D[2]);

    int lx_sum = lmn[0] + lmn[3] + lmn[6] + lmn[9];
    int ly_sum = lmn[1] + lmn[4] + lmn[7] + lmn[10];
    int lz_sum = lmn[2] + lmn[5] + lmn[8] + lmn[11];
    int vmax = lx_sum + ly_sum + lz_sum;

    double B_x[lx_sum + 1];
    double B_y[ly_sum + 1];
    double B_z[lz_sum + 1];
    double boys[vmax + 1];

    for(int i = 0; i < num1; i++) {
        if (coef1[i] == 0.0) continue; 
        
        for(int j = 0; j < num2; j++) {
            if (coef2[j] == 0.0) continue;
            
            double gamma1 = exp1[i] + exp2[j];
            double Px = (exp1[i]*A[0] + exp2[j]*B[0]) / gamma1;
            double Py = (exp1[i]*A[1] + exp2[j]*B[1]) / gamma1;
            double Pz = (exp1[i]*A[2] + exp2[j]*B[2]) / gamma1;

            for(int k = 0; k < num3; k++) {
                if (coef3[k] == 0.0) continue;
                
                for(int l = 0; l < num4; l++) {
                    if (coef4[l] == 0.0) continue;
                    
                    double gamma2 = exp3[k] + exp4[l];
                    double delta = 1.0/(4.0*gamma1) + 1.0/(4.0*gamma2);

                    double Qx = (exp3[k]*C[0] + exp4[l]*D[0]) / gamma2;
                    double Qy = (exp3[k]*C[1] + exp4[l]*D[1]) / gamma2;
                    double Qz = (exp3[k]*C[2] + exp4[l]*D[2]) / gamma2;

                    double p2 = (Px-Qx)*(Px-Qx) + (Py-Qy)*(Py-Qy) + (Pz-Qz)*(Pz-Qz);

                    for(int idx = 0; idx <= lx_sum; idx++) B_x[idx] = 0.0;
                    for(int idx = 0; idx <= ly_sum; idx++) B_y[idx] = 0.0;
                    for(int idx = 0; idx <= lz_sum; idx++) B_z[idx] = 0.0;

                    calculate_B(lmn[0], lmn[3], lmn[6], lmn[9], A[0], B[0], C[0], D[0], Px, Qx, gamma1, gamma2, delta, B_x);
                    calculate_B(lmn[1], lmn[4], lmn[7], lmn[10], A[1], B[1], C[1], D[1], Py, Qy, gamma1, gamma2, delta, B_y);
                    calculate_B(lmn[2], lmn[5], lmn[8], lmn[11], A[2], B[2], C[2], D[2], Pz, Qz, gamma1, gamma2, delta, B_z);

                    boys_recurrence(vmax, p2/(4.0*delta), boys);

                    double omega = 2.0 * M_PI * M_PI * sqrt(M_PI/(gamma1 + gamma2)) * exp(-(exp1[i]*exp2[j]*ab2/gamma1 + exp3[k]*exp4[l]*cd2/gamma2)) / (gamma1*gamma2);

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

double repulsion_matrix(
    int dim, double* centers, int* lmns, int* lens,         
    int* offsets, double* exps, double* coefs, 
    double* normalize, double* matrix  
) {
    for(int i = 0; i < dim; i++) {
        for(int j = i; j < dim; j++) {
            for(int k = 0; k < dim; k++) {
                for(int l = k; l < dim; l++) {
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

                    double value = integral * normalize[i] * normalize[j] * normalize[k] * normalize[l];

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