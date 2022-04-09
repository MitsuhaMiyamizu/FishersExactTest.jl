/*#ifndef HTSLIB_KFUNC_H
#define HTSLIB_KFUNC_H 
*/

#include <cmath>
#include <stdlib.h>

#ifndef KFUNC_H_
#define KFUNC_H_

#include <cmath>
#include <stdlib.h>

#ifdef M_SQRT2
#undef M_SQRT2
#endif // M_SQRT2

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

double kf_lgamma(double z);
double kf_erfc(double x);
#define KF_GAMMA_EPS 1e-14
#define KF_TINY 1e-290

// regularized lower incomplete gamma function, by series expansion
double kf_betai(double a, double b, double x);
double kf_gammap(double s, double z);
double kf_gammaq(double s, double z);

static double kf_betai_aux(double a, double b, double x);
static double _kf_gammap(double s, double z);
static double _kf_gammaq(double s, double z);
static double lbinom(long long n, long long k);
static double hypergeo(long long n11, long long n1_, long long n_1, long long n);

typedef struct {
    long long n11, n1_, n_1, n;
    double p;
} hgacc_t;

// incremental version of hypergenometric distribution
double kt_fisher_exact(long long n11, long long n12, long long n21, long long n22, double *_left, double *_right, double *two);

static double hypergeo_acc(long long n11, long long n1_, long long n_1, long long n, hgacc_t *aux);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif /* KFUNC_H_ */
