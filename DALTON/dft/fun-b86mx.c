/*


!
!  Dalton, a molecular electronic structure program
!  Copyright (C) The Dalton Authors (see AUTHORS file for details).
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!

!

*/
/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* fun-b86mx.c:

   Automatically generated code implementing B86MX functional and
   its derivatives. It is generated by func-codegen.pl being a part of
   a "Automatic code generation framework for analytical functional
   derivative evaluation", Pawel Salek, 2005

    This functional is connected by making following changes:
    1. add "extern Functional b86mxFunctional;" to 'functionals.h'
    2. add "&b86mxFunctional," to 'functionals.c'
    3. add "fun-b86mx.c" to 'Makefile.am', 'Makefile.in' or 'Makefile'.

    This functional has been generated from following input:
    ------ cut here -------
Cx:   -3/4*(6/%PI)^(1/3);
bet:  0.00375;
lam:  0.007;

xa:   grada/(rhoa^(4/3));
xb:   gradb/(rhob^(4/3));

da:   (1+lam*xa^2)^(4/5);
db:   (1+lam*xb^2)^(4/5);

Exa:  (-rhoa^(4/3))*(Cx + bet*xa^2/da);
Exb:  (-rhob^(4/3))*(Cx + bet*xb^2/db);

K(rhoa,rhob,grada,gradb,gradab):=(Exa+Exb);


    ------ cut here -------
*/

 
/* strictly conform to XOPEN ANSI C standard */
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas and pow() prototype
 * is not specified. */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif
#include <math.h>
#include <stddef.h>
#include "general.h"

#define __CVERSION__
 
#include "functionals.h"
 
/* INTERFACE PART */
static integer b86mx_isgga(void) { return 1; } /* FIXME: detect! */
static integer b86mx_read(const char *conf_line);
static real b86mx_energy(const FunDensProp* dp);
static void b86mx_first(FunFirstFuncDrv *ds,   real factor,
                         const FunDensProp* dp);
static void b86mx_second(FunSecondFuncDrv *ds, real factor,
                          const FunDensProp* dp);
static void b86mx_third(FunThirdFuncDrv *ds,   real factor,
                         const FunDensProp* dp);
static void b86mx_fourth(FunFourthFuncDrv *ds,   real factor,
                          const FunDensProp* dp);
 
Functional B86mxFunctional = {
  "B86mx",       /* name */
  b86mx_isgga,   /* gga-corrected */
   1,
  b86mx_read,
  NULL,
  b86mx_energy,
  b86mx_first,
  b86mx_second,
  b86mx_third,
  b86mx_fourth
};
 
/* IMPLEMENTATION PART */
static integer
b86mx_read(const char *conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}

static real
b86mx_energy(const FunDensProp *dp)
{
    real res;
    real rhoa = dp->rhoa, rhob = dp->rhob;
    real grada = dp->grada, gradb = dp->gradb, gradab = dp->gradab;

    real t1, t2, t3, t4, t5;

    t1 = -0.75*pow(6.0,0.333333333333333)/pow(3.141592653589793,
        0.333333333333333);
    t2 = pow(grada,2.0);
    t3 = 1/pow(rhoa,2.666666666666667);
    t4 = pow(gradb,2.0);
    t5 = 1/pow(rhob,2.666666666666667);

   /* code */
    res = -1.0*(0.00375*t4*t5/pow(0.007*t4*t5+1.0,0.8)+t1)*
        pow(rhob,1.333333333333333)-1.0*(0.00375*t2*t3/pow(0.007*t2*
        t3+1.0,0.8)+t1)*pow(rhoa,1.333333333333333);

    return res;
}

static void
b86mx_first_helper(real rhoa, real grada, real *res)
{    real t1, t2, t3, t4, t5, t6;

    t1 = pow(grada,2.0);
    t2 = 1/pow(rhoa,2.666666666666667);
    t3 = 0.007*t1*t2+1.0;
    t4 = 1/pow(t3,0.8);
    t5 = 1/pow(t3,1.8);
    t6 = pow(rhoa,1.333333333333333);

   /* code */
    res[0] = -1.0*t6*(5.6E-5*t5*pow(grada,4.0)/pow(rhoa,6.333333333333333)-
        0.01*t1*t4/pow(rhoa,3.666666666666667))-1.333333333333333*
        (0.00375*t1*t4*t2-0.75*pow(6.0,0.333333333333333)/pow(3.141592653589793,
        0.333333333333333))*pow(rhoa,0.333333333333333);
    res[1] = -1.0*t6*(0.0075*grada*t4*t2-4.2E-5*t5*pow(grada,
        3.0)/pow(rhoa,5.333333333333333));
}

static void
b86mx_first(FunFirstFuncDrv *ds, real factor, const FunDensProp *dp)
{
    real res[2];

    b86mx_first_helper(dp->rhoa, dp->grada, res);
   /* Final assignment */
    ds->df1000 += factor*res[0];
    ds->df0010 += factor*res[1];


    if(fabs(dp->rhoa-dp->rhob)>1e-13 ||
       fabs(dp->grada-dp->gradb)>1e-13)
        b86mx_first_helper(dp->rhob, dp->gradb, res);
    ds->df0100 += factor*res[0];
    ds->df0001 += factor*res[1];

}

static void
b86mx_second_helper(real rhoa, real grada, real *res)
{
    real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    real t11, t12, t13, t14, t15, t16;

    t1 = pow(grada,2.0);
    t2 = 1/pow(rhoa,2.666666666666667);
    t3 = 0.007*t1*t2+1.0;
    t4 = 1/pow(t3,0.8);
    t5 = 0.00375*t1*t4*t2-0.75*pow(6.0,0.333333333333333)/
        pow(3.141592653589793,0.333333333333333);
    t6 = pow(rhoa,0.333333333333333);
    t7 = pow(grada,4.0);
    t8 = 1/pow(t3,1.8);
    t9 = 1/pow(rhoa,6.333333333333333);
    t10 = 1/pow(rhoa,3.666666666666667);
    t11 = 5.6E-5*t7*t8*t9-0.01*t1*t4*t10;
    t12 = pow(rhoa,1.333333333333333);
    t13 = pow(grada,3.0);
    t14 = 1/pow(rhoa,5.333333333333333);
    t15 = 0.0075*grada*t4*t2-4.2E-5*t13*t8*t14;
    t16 = 1/pow(t3,2.8);

   /* code */
    res[0] = -1.333333333333333*t5*t6-1.0*t11*t12;
    res[1] = -1.0*t12*t15;
    res[2] = -1.0*t12*(1.8816E-6*t16*pow(grada,6.0)/pow(rhoa,
        10.0)-5.04E-4*t7*t8/pow(rhoa,7.333333333333333)+0.036666666666667*
        t1*t4/pow(rhoa,4.666666666666667))-0.444444444444444*t5/pow(rhoa,
        0.666666666666667)-2.666666666666667*t11*t6;
    res[3] = -1.0*t12*(-1.4112E-6*t16*pow(grada,5.0)/pow(rhoa,
        9.0)+3.36E-4*t13*t8*t9-0.02*grada*t4*t10)-1.333333333333333*
        t15*t6;
    res[4] = -1.0*t12*(1.0584000000000001E-6*t16*t7/pow(rhoa,
        8.0)+0.0075*t4*t2-2.1000000000000004E-4*t1*t8*t14);

}

static void
b86mx_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real res[5];
 
    b86mx_second_helper(dp->rhoa, dp->grada, res);

    ds->df1000 += factor*res[0];
    ds->df0010 += factor*res[1];

    ds->df2000 += factor*res[2];
    ds->df1010 += factor*res[3];
    ds->df0020 += factor*res[4];


    if(fabs(dp->rhoa-dp->rhob)>1e-13 ||
       fabs(dp->grada-dp->gradb)>1e-13)
        b86mx_second_helper(dp->rhob, dp->gradb, res);
    ds->df0100 += factor*res[0];
    ds->df0001 += factor*res[1];

    ds->df0200 += factor*res[2];
    ds->df0101 += factor*res[3];
    ds->df0002 += factor*res[4];

}

static void
b86mx_third_helper(real rhoa, real grada, real *res)
{
    real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    real t11, t12, t13, t14, t15, t16, t17, t18;
    real t19, t20, t21, t22, t23, t24, t25, t26;
    real t27, t28;

    t1 = pow(grada,2.0);
    t2 = 1/pow(rhoa,2.666666666666667);
    t3 = 0.007*t1*t2+1.0;
    t4 = 1/pow(t3,0.8);
    t5 = 0.00375*t1*t4*t2-0.75*pow(6.0,0.333333333333333)/
        pow(3.141592653589793,0.333333333333333);
    t6 = pow(rhoa,0.333333333333333);
    t7 = pow(grada,4.0);
    t8 = 1/pow(t3,1.8);
    t9 = 1/pow(rhoa,6.333333333333333);
    t10 = 1/pow(rhoa,3.666666666666667);
    t11 = 5.6E-5*t7*t8*t9-0.01*t1*t4*t10;
    t12 = pow(rhoa,1.333333333333333);
    t13 = pow(grada,3.0);
    t14 = 1/pow(rhoa,5.333333333333333);
    t15 = 0.0075*grada*t4*t2-4.2E-5*t13*t8*t14;
    t16 = 1/pow(rhoa,0.666666666666667);
    t17 = pow(grada,6.0);
    t18 = 1/pow(t3,2.8);
    t19 = 1/pow(rhoa,10.0);
    t20 = 1/pow(rhoa,7.333333333333333);
    t21 = 1/pow(rhoa,4.666666666666667);
    t22 = 0.036666666666667*t1*t4*t21-5.04E-4*t7*t8*t20+1.8816E-6*
        t17*t18*t19;
    t23 = pow(grada,5.0);
    t24 = 1/pow(rhoa,9.0);
    t25 = -0.02*grada*t4*t10+3.36E-4*t13*t8*t9-1.4112E-6*
        t23*t18*t24;
    t26 = 1/pow(rhoa,8.0);
    t27 = 0.0075*t4*t2-2.1000000000000004E-4*t1*t8*t14+1.0584000000000001E-6*
        t7*t18*t26;
    t28 = 1/pow(t3,3.8);

   /* code */
    res[0] = -1.333333333333333*t5*t6-1.0*t11*t12;
    res[1] = -1.0*t12*t15;
    res[2] = -2.666666666666667*t11*t6-0.444444444444444*
        t16*t5-1.0*t12*t22;
    res[3] = -1.333333333333333*t15*t6-1.0*t12*t25;
    res[4] = -1.0*t12*t27;
    res[5] = -1.0*t12*(9.834495999999997E-8*t28*pow(grada,
        8.0)/pow(rhoa,13.66666666666667)-3.57504E-5*t17*t18/pow(rhoa,
        11.0)+0.004243555555556*t7*t8/pow(rhoa,8.333333333333334)-
        0.171111111111111*t1*t4/pow(rhoa,5.666666666666667))+0.296296296296296*
        t5/pow(rhoa,1.666666666666667)-4.0*t22*t6-1.333333333333333*
        t11*t16;
    res[6] = -1.0*t12*(-7.375872E-8*t28*pow(grada,7.0)/pow(rhoa,
        12.66666666666667)+0.073333333333333*grada*t4*t21-0.002426666666667*
        t13*t8*t20+2.39904E-5*t23*t18*t19)-2.666666666666667*t25*t6-
        0.444444444444444*t15*t16;
    res[7] = -1.0*t12*(5.531903999999999E-8*t17*t28/pow(rhoa,
        11.66666666666667)+0.001232*t1*t8*t9-1.55232E-5*t7*t18*t24-
        0.02*t4*t10)-1.333333333333333*t27*t6;
    res[8] = -1.0*t12*(-4.148928000000001E-8*t23*t28/pow(rhoa,
        10.66666666666667)+9.525600000000002E-6*t13*t18*t26-5.040000000000001E-4*
        grada*t8*t14);

}

static void
b86mx_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real res[9];
 
    b86mx_third_helper(dp->rhoa, dp->grada, res);

    ds->df1000 += factor*res[0];
    ds->df0010 += factor*res[1];

    ds->df2000 += factor*res[2];
    ds->df1010 += factor*res[3];
    ds->df0020 += factor*res[4];

    ds->df3000 += factor*res[5];
    ds->df2010 += factor*res[6];
    ds->df1020 += factor*res[7];
    ds->df0030 += factor*res[8];


    if(fabs(dp->rhoa-dp->rhob)>1e-13 ||
       fabs(dp->grada-dp->gradb)>1e-13)
        b86mx_third_helper(dp->rhob, dp->gradb, res);

    ds->df0100 += factor*res[0];
    ds->df0001 += factor*res[1];

    ds->df0200 += factor*res[2];
    ds->df0101 += factor*res[3];
    ds->df0002 += factor*res[4];

    ds->df0300 += factor*res[5];
    ds->df0201 += factor*res[6];
    ds->df0102 += factor*res[7];
    ds->df0003 += factor*res[8];

}

static void
b86mx_fourth_helper(real rhoa, real grada, real *res)
{
    real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    real t11, t12, t13, t14, t15, t16, t17, t18;
    real t19, t20, t21, t22, t23, t24, t25, t26;
    real t27, t28, t29, t30, t31, t32, t33, t34;
    real t35, t36, t37, t38, t39, t40, t41, t42;
    real t43;

    t1 = pow(grada,2.0);
    t2 = 1/pow(rhoa,2.666666666666667);
    t3 = 0.007*t1*t2+1.0;
    t4 = 1/pow(t3,0.8);
    t5 = 0.00375*t1*t4*t2-0.75*pow(6.0,0.333333333333333)/
        pow(3.141592653589793,0.333333333333333);
    t6 = pow(rhoa,0.333333333333333);
    t7 = pow(grada,4.0);
    t8 = 1/pow(t3,1.8);
    t9 = 1/pow(rhoa,6.333333333333333);
    t10 = 1/pow(rhoa,3.666666666666667);
    t11 = 5.6E-5*t7*t8*t9-0.01*t1*t4*t10;
    t12 = pow(rhoa,1.333333333333333);
    t13 = pow(grada,3.0);
    t14 = 1/pow(rhoa,5.333333333333333);
    t15 = 0.0075*grada*t4*t2-4.2E-5*t13*t8*t14;
    t16 = 1/pow(rhoa,0.666666666666667);
    t17 = pow(grada,6.0);
    t18 = 1/pow(t3,2.8);
    t19 = 1/pow(rhoa,10.0);
    t20 = 1/pow(rhoa,7.333333333333333);
    t21 = 1/pow(rhoa,4.666666666666667);
    t22 = 0.036666666666667*t1*t4*t21-5.04E-4*t7*t8*t20+1.8816E-6*
        t17*t18*t19;
    t23 = pow(grada,5.0);
    t24 = 1/pow(rhoa,9.0);
    t25 = -0.02*grada*t4*t10+3.36E-4*t13*t8*t9-1.4112E-6*
        t23*t18*t24;
    t26 = 1/pow(rhoa,8.0);
    t27 = 0.0075*t4*t2-2.1000000000000004E-4*t1*t8*t14+1.0584000000000001E-6*
        t7*t18*t26;
    t28 = 1/pow(rhoa,1.666666666666667);
    t29 = pow(grada,8.0);
    t30 = 1/pow(t3,3.8);
    t31 = 1/pow(rhoa,13.66666666666667);
    t32 = 1/pow(rhoa,11.0);
    t33 = 1/pow(rhoa,8.333333333333334);
    t34 = 1/pow(rhoa,5.666666666666667);
    t35 = -0.171111111111111*t1*t4*t34+0.004243555555556*
        t7*t8*t33-3.57504E-5*t17*t18*t32+9.834495999999997E-8*t29*
        t30*t31;
    t36 = pow(grada,7.0);
    t37 = 1/pow(rhoa,12.66666666666667);
    t38 = 0.073333333333333*grada*t4*t21-0.002426666666667*
        t13*t8*t20+2.39904E-5*t23*t18*t19-7.375872E-8*t36*t30*t37;
    t39 = 1/
        pow(rhoa,11.66666666666667);
    t40 = -0.02*t4*t10+0.001232*t1*t8*t9-1.55232E-5*t7*t18*
        t24+5.531903999999999E-8*t17*t30*t39;
    t41 = 1/pow(rhoa,10.66666666666667);
    t42 = -5.040000000000001E-4*grada*t8*t14+9.525600000000002E-6*
        t13*t18*t26-4.148928000000001E-8*t23*t30*t41;
    t43 = 1/pow(t3,4.8);

   /* code */
    res[0] = -1.333333333333333*t5*t6-1.0*t11*t12;
    res[1] = -1.0*t12*t15;
    res[2] = -2.666666666666667*t11*t6-0.444444444444444*
        t16*t5-1.0*t12*t22;
    res[3] = -1.333333333333333*t15*t6-1.0*t12*t25;
    res[4] = -1.0*t12*t27;
    res[5] = -4.0*t22*t6+0.296296296296296*t28*t5-1.0*t12*
        t35-1.333333333333333*t11*t16;
    res[6] = -2.666666666666667*t25*t6-1.0*t12*t38-0.444444444444444*
        t15*t16;
    res[7] = -1.333333333333333*t27*t6-1.0*t12*t40;
    res[8] = -1.0*t12*t42;
    res[9] = -1.0*t12*(6.97593582933333E-9*t43*pow(grada,
        10.0)/pow(rhoa,17.33333333333333)-3.212602026666666E-6*t29*
        t30/pow(rhoa,14.66666666666667)+5.358378666666667E-4*t17*t18/
        pow(rhoa,12.0)-0.037918222222222*t7*t8/pow(rhoa,9.333333333333334)+
        0.96962962962963*t1*t4/pow(rhoa,6.666666666666667))-5.333333333333333*
        t35*t6-0.493827160493827*t2*t5+1.185185185185185*t11*t28-2.666666666666667*
        t16*t22;
    res[10] = -1.0*t12*(-5.231951871999998E-9*t43*pow(grada,
        9.0)/pow(rhoa,16.33333333333333)-0.342222222222222*grada*t4*
        t34+0.018890666666667*t13*t8*t33-3.214399999999999E-4*t23*
        t18*t32+2.18817536E-6*t36*t30*t31)-4.0*t38*t6+0.296296296296296*
        t15*t28-1.333333333333333*t16*t25;
    res[11] = -1.0*t12*(3.923963904E-9*t29*t43/pow(rhoa,15.33333333333333)-
        1.45673472E-6*t17*t30*t37+0.073333333333333*t4*t21-0.008101333333333*
        t1*t8*t20+1.81104E-4*t7*t18*t19)-2.666666666666667*t40*t6-
        0.444444444444444*t16*t27;
    res[12] = -1.0*t12*(-2.942972928E-9*t36*t43/pow(rhoa,
        14.33333333333333)+0.002688*grada*t8*t9+9.4042368E-7*t23*t30*
        t39-9.31392E-5*t13*t18*t24)-1.333333333333333*t42*t6;
    res[13] = -1.0*t12*(2.207229696E-9*t17*t43/pow(rhoa,13.33333333333333)-
        5.8084992E-7*t7*t30*t41+4.127760000000001E-5*t1*t18*t26-5.040000000000001E-4*
        t8*t14);

}

static void
b86mx_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real res[14];
 
    b86mx_fourth_helper(dp->rhoa, dp->grada, res);

    ds->df1000 += factor*res[0];
    ds->df0010 += factor*res[1];

    ds->df2000 += factor*res[2];
    ds->df1010 += factor*res[3];
    ds->df0020 += factor*res[4];

    ds->df3000 += factor*res[5];
    ds->df2010 += factor*res[6];
    ds->df1020 += factor*res[7];
    ds->df0030 += factor*res[8];

    ds->df4000 += factor*res[9];
    ds->df3010 += factor*res[10];
    ds->df2020 += factor*res[11];
    ds->df1030 += factor*res[12];
    ds->df0040 += factor*res[13];


    if(fabs(dp->rhoa-dp->rhob)>1e-13 ||
       fabs(dp->grada-dp->gradb)>1e-13)
        b86mx_fourth_helper(dp->rhob, dp->gradb, res);

    ds->df0100 += factor*res[0];
    ds->df0001 += factor*res[1];

    ds->df0200 += factor*res[2];
    ds->df0101 += factor*res[3];
    ds->df0002 += factor*res[4];

    ds->df0300 += factor*res[5];
    ds->df0201 += factor*res[6];
    ds->df0102 += factor*res[7];
    ds->df0003 += factor*res[8];

    ds->df0400 += factor*res[9];
    ds->df0301 += factor*res[10];
    ds->df0202 += factor*res[11];
    ds->df0103 += factor*res[12];
    ds->df0004 += factor*res[13];

}
