#include <stdio.h>
#include "rtcore/raytrace.h"

extern double density_only(
    int iRay,
    double Param_I[],
    double x,
    double y,
    double z);

int main() {

    struct param prm;

    prm.r_chromo = 1.0129;
    prm.r_corona = 1.0158;

    double *Param_I = (double*)&prm;

    double Ne = density_only(
        0,
        Param_I,
        1.1,
        0.0,
        0.0);

    printf("Ne = %e\n", Ne);

    return 0;
}