/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   bmath.h
 * Author: veraalva
 *
 * Created on February 18, 2016, 2:49 PM
 */

#ifndef BMATH_H
#define BMATH_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI 0.918938533204672741780329736406 /* log(sqrt(2*pi)) == log(2*pi)/2 */
#endif

#ifndef M_LN2
#define M_LN2  0.693147180559945309417232121458 /* ln(2) */
#endif

#ifndef M_2PI
#define M_2PI  6.283185307179586476925286766559 /* 2*pi */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2 0.225791352644727432363097614947 /* log(sqrt(pi/2)) */
#endif

#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131E-16
#endif

#define give_log log_p
#define R_D__0 (log_p ? -INFINITY : 0.)  /* 0 */
#define R_D__1 (log_p ? 0. : 1.)   /* 1 */
#define R_DT_0 (lower_tail ? R_D__0 : R_D__1)  /* 0 */
#define R_DT_1 (lower_tail ? R_D__1 : R_D__0)  /* 1 */
#define R_D_nonint(x)    (fabs((x) - floor((x)+0.5)) > 1e-7)
#define R_D_negInonint(x) (x < 0. || R_D_nonint(x))
#define R_D_forceint(x)   floor((x) + 0.5)
#define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))
#define R_DT_Log(p) (lower_tail? (p) : R_Log1_Exp(p))
#define R_D_Lval(p) (lower_tail ? (p) : (0.5 - (p) + 0.5)) /*  p  */
#define R_D_exp(x) (log_p ?  (x)  : exp(x)) /* exp(x) */

#define ML_ERR_return_NAN { fprintf(stderr, "ML_ERR_return_NAN\n"); return NAN; }

    extern double stirlerr(double n);

    extern double lgammafn(double x);

    extern double gammafn(double x);

    extern double lgammacor(double x);

    extern double chebyshev_eval(double x, const double *a, const int n);

    extern double phyper(double x, double NR, double NB, double n,
            int lower_tail, int log_p);

#ifdef __cplusplus
}
#endif

#endif /* BMATH_H */

