#include <stdio.h>
#include <math.h>
#include "bmath.h"

double lgammafn_sign(double x, int *sgn) {
    double ans, y, sinpiy;

#ifdef NOMORE_FOR_THREADS
    static double xmax = 0.;
    static double dxrel = 0.;

    if (xmax == 0) {/* initialize machine dependent constants _ONCE_ */
        xmax = d1mach(2) / log(d1mach(2)); /* = 2.533 e305	 for IEEE double */
        dxrel = sqrt(d1mach(4)); /* sqrt(Eps) ~ 1.49 e-8  for IEEE double */
    }
#else
    /* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
       xmax  = DBL_MAX / log(DBL_MAX) = 2^1024 / (1024 * log(2)) = 2^1014 / log(2)
       dxrel = sqrt(DBL_EPSILON) = 2^-26 = 5^26 * 1e-26 (is *exact* below !)
     */
#define xmax  2.5327372760800758e+305
#define dxrel 1.490116119384765696e-8
#endif

    if (sgn != NULL) *sgn = 1;

#ifdef IEEE_754
    if (ISNAN(x)) return x;
#endif

    if (x < 0 && fmod(floor(-x), 2.) == 0)
        if (sgn != NULL) *sgn = -1;

    if (x <= 0 && x == trunc(x)) { /* Negative integer argument */
        return INFINITY; /* +Inf, since lgamma(x) = log|gamma(x)| */
    }

    y = fabs(x);

    if (y < 1e-306) return -log(x); // denormalized range, R change
    if (y <= 10) return log(fabs(gammafn(x)));
    /*
      ELSE  y = |x| > 10 ---------------------- */

    if (y > xmax) {
        return INFINITY;
    }

    if (x > 0) { /* i.e. y = x > 10 */
#ifdef IEEE_754
        if (x > 1e17)
            return (x * (log(x) - 1.));
        else if (x > 4934720.)
            return (M_LN_SQRT_2PI + (x - 0.5) * log(x) - x);
        else
#endif
            return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x + lgammacor(x);
    }
    /* else: x < -10; y = -x */
    sinpiy = fabs(sin(M_PI * y));

    if (sinpiy == 0) { /* Negative integer argument ===
			  Now UNNECESSARY: caught above */
        ML_ERR_return_NAN;
    }

    ans = M_LN_SQRT_PId2 + (x - 0.5) * log(y) - x - log(sinpiy) - lgammacor(y);

    if (fabs((x - trunc(x - 0.5)) * ans / x) < dxrel) {

        /* The answer is less than half precision because
         * the argument is too near a negative integer. */

        fprintf(stderr, "lgamma\n");
    }

    return ans;
}

double lgammafn(double x) {
    return lgammafn_sign(x, NULL);
}