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
/* fun-dk87x.c:

   Automatically generated code implementing DK87X functional and
   its derivatives. It is generated by func-codegen.pl being a part of
   a "Automatic code generation framework for analytical functional
   derivative evaluation", Pawel Salek, 2005

    This functional is connected by making following changes:
    1. add "extern Functional dk87xFunctional;" to 'functionals.h'
    2. add "&dk87xFunctional," to 'functionals.c'
    3. add "fun-dk87x.c" to 'Makefile.am', 'Makefile.in' or 'Makefile'.

    This functional has been generated from following input:
    ------ cut here -------
rho:  rhoa + rhob;
grad: sqrt(grada*grada + gradb*gradb + 2*gradab);
zeta: (rhoa-rhob)/(rhoa+rhob);

Bg:   0.00132327;
a1:   0.861213;
b1:   0.042076;
c:    0.98;

xa:  grada/(rhoa^(4/3));
xb:  gradb/(rhob^(4/3));

t1a: Bg*xa^2;
t1b: Bg*xb^2;

Exa: -rhoa^(4/3)*(t1a*(1+a1*xa^c)/(1+b1*xa^2));
Exb: -rhob^(4/3)*(t1b*(1+a1*xb^c)/(1+b1*xb^2));

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
static integer dk87x_isgga(void) { return 1; } /* FIXME: detect! */
static integer dk87x_read(const char *conf_line);
static real dk87x_energy(const FunDensProp* dp);
static void dk87x_first(FunFirstFuncDrv *ds,   real factor,
                         const FunDensProp* dp);
static void dk87x_second(FunSecondFuncDrv *ds, real factor,
                          const FunDensProp* dp);
static void dk87x_third(FunThirdFuncDrv *ds,   real factor,
                         const FunDensProp* dp);
static void dk87x_fourth(FunFourthFuncDrv *ds,   real factor,
                          const FunDensProp* dp);
 
Functional DK87xFunctional = {
  "DK87x",       /* name */
  dk87x_isgga,   /* gga-corrected */
   1,
  dk87x_read,
  NULL,
  dk87x_energy,
  dk87x_first,
  dk87x_second,
  dk87x_third,
  dk87x_fourth
};
 
/* IMPLEMENTATION PART */
static integer
dk87x_read(const char *conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}

static real
dk87x_energy(const FunDensProp *dp)
{
    real res;
    real rhoa = dp->rhoa, rhob = dp->rhob;
    real grada = dp->grada, gradb = dp->gradb, gradab = dp->gradab;

    real t1, t2, t3, t4;

    t1 = pow(grada,2.0);
    t2 = 1/pow(rhoa,1.333333333333333);
    t3 = pow(gradb,2.0);
    t4 = 1/pow(rhob,1.333333333333333);

   /* code */
    res = -0.00132327*t3*t4*(0.861213*pow(gradb*t4,0.98)+
        1.0)/(0.042076*t3/pow(rhob,2.666666666666667)+1.0)-0.00132327*
        t1*t2*(0.861213*pow(grada*t2,0.98)+1.0)/(0.042076*t1/pow(rhoa,
        2.666666666666667)+1.0);

    return res;
}

static void
dk87x_first_helper(real rhoa, real grada, real *res)
{    real t1, t2, t3, t4, t5, t6, t7, t8, t9;

    t1 = 1/pow(rhoa,1.333333333333333);
    t2 = grada*t1;
    t3 = 0.861213*pow(t2,0.98)+1.0;
    t4 = pow(grada,2.0);
    t5 = 1/pow(rhoa,2.666666666666667);
    t6 = 0.042076*t4*t5+1.0;
    t7 = 1/pow(t6,2.0);
    t8 = pow(grada,3.0);
    t9 = 1/t6;

   /* code */
    res[0] = -1.4847442272E-4*t3*t7*pow(grada,4.0)/pow(rhoa,
        5.0)+0.001489099973306*t8*t9/(pow(t2,0.02)*pow(rhoa,3.666666666666667))+
        0.00176436*t3*t4*t9/pow(rhoa,2.333333333333333);
    res[1] = -0.00111682497998*t5*t9*pow(grada,1.98)*pow(fabs(rhoa),
        0.026666666666667)+1.1135581704000001E-4*t3*t7*t8/pow(rhoa,
        4.0)-0.00264654*grada*t3*t9*t1;
}

static void
dk87x_first(FunFirstFuncDrv *ds, real factor, const FunDensProp *dp)
{
    real res[2];

    dk87x_first_helper(dp->rhoa, dp->grada, res);
   /* Final assignment */
    ds->df1000 += factor*res[0];
    ds->df0010 += factor*res[1];


    if(fabs(dp->rhoa-dp->rhob)>1e-13 ||
       fabs(dp->grada-dp->gradb)>1e-13)
        dk87x_first_helper(dp->rhob, dp->gradb, res);
    ds->df0100 += factor*res[0];
    ds->df0001 += factor*res[1];

}

static void
dk87x_second_helper(real rhoa, real grada, real *res)
{
    real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    real t11, t12, t13, t14, t15, t16, t17, t18;
    real t19;

    t1 = pow(grada,4.0);
    t2 = 1/pow(rhoa,1.333333333333333);
    t3 = grada*t2;
    t4 = 0.861213*pow(t3,0.98)+1.0;
    t5 = pow(grada,2.0);
    t6 = 1/pow(rhoa,2.666666666666667);
    t7 = 0.042076*t5*t6+1.0;
    t8 = 1/pow(t7,2.0);
    t9 = 1/pow(rhoa,5.0);
    t10 = pow(grada,3.0);
    t11 = 1/t7;
    t12 = 1/pow(rhoa,3.666666666666667);
    t13 = 1/pow(rhoa,2.333333333333333);
    t14 = 1/pow(rhoa,4.0);
    t15 = pow(grada,1.98);
    t16 = fabs(rhoa);
    t17 = pow(t16,0.026666666666667);
    t18 = 1/pow(t7,3.0);
    t19 = pow(grada,2.98);

   /* code */
    res[0] = -1.4847442272E-4*t1*t4*t8*t9+0.001489099973306*
        t10*t11*t12/pow(t3,0.02)+0.00176436*t5*t4*t11*t13;
    res[1] = -0.00111682497998*t15*t11*t6*t17-0.00264654*
        grada*t4*t11*t2+1.1135581704000001E-4*t10*t4*t8*t14;
    res[2] = -3.331845232195585E-5*t18*t4*pow(grada,6.0)/
        pow(rhoa,8.666666666666666)+3.341619758764805E-4*t17*t8*pow(grada,
        4.98)/pow(rhoa,7.333333333333333)+9.403380105600001E-4*t1*
        t4*t8/pow(rhoa,6.0)-0.007445499866532*t11*t17*t19/pow(rhoa,
        4.666666666666667)-0.00411684*t11*t4*t5/pow(rhoa,3.333333333333333)+
        3.970933262150404E-5*t11*t19*t6/pow(t16,1.973333333333333);
    res[3] = 2.4988839241466884E-5*t18*t4*pow(grada,5.0)/
        pow(rhoa,7.666666666666667)-2.506214819073604E-4*t17*t8*pow(grada,
        3.98)/pow(rhoa,6.333333333333333)-7.423721136000001E-4*t10*
        t4*t8*t9+0.005926617893759*t15*t11*t12*t17+0.00352872*grada*
        t4*t11*t13;
    res[4] = -1.8741629431100162E-5*t1*t18*t4/pow(rhoa,6.666666666666667)+
        1.8796611143052027E-4*t17*t19*t8/pow(rhoa,5.333333333333333)-
        0.00444496342032*t11*t17*t6*pow(grada,0.98)-0.00264654*t4*
        t11*t2+5.567790852E-4*t5*t4*t8*t14;

}

static void
dk87x_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real res[5];
 
    dk87x_second_helper(dp->rhoa, dp->grada, res);

    ds->df1000 += factor*res[0];
    ds->df0010 += factor*res[1];

    ds->df2000 += factor*res[2];
    ds->df1010 += factor*res[3];
    ds->df0020 += factor*res[4];


    if(fabs(dp->rhoa-dp->rhob)>1e-13 ||
       fabs(dp->grada-dp->gradb)>1e-13)
        dk87x_second_helper(dp->rhob, dp->gradb, res);
    ds->df0100 += factor*res[0];
    ds->df0001 += factor*res[1];

    ds->df0200 += factor*res[2];
    ds->df0101 += factor*res[3];
    ds->df0002 += factor*res[4];

}

static void
dk87x_third_helper(real rhoa, real grada, real *res)
{
    real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    real t11, t12, t13, t14, t15, t16, t17, t18;
    real t19, t20, t21, t22, t23, t24, t25, t26;
    real t27, t28, t29, t30, t31, t32, t33, t34;
    real t35;

    t1 = pow(grada,4.0);
    t2 = 1/pow(rhoa,1.333333333333333);
    t3 = grada*t2;
    t4 = 0.861213*pow(t3,0.98)+1.0;
    t5 = pow(grada,2.0);
    t6 = 1/pow(rhoa,2.666666666666667);
    t7 = 0.042076*t5*t6+1.0;
    t8 = 1/pow(t7,2.0);
    t9 = 1/pow(rhoa,5.0);
    t10 = pow(grada,3.0);
    t11 = 1/t7;
    t12 = 1/pow(rhoa,3.666666666666667);
    t13 = 1/pow(rhoa,2.333333333333333);
    t14 = 1/pow(rhoa,4.0);
    t15 = pow(grada,1.98);
    t16 = fabs(rhoa);
    t17 = pow(t16,0.026666666666667);
    t18 = pow(grada,6.0);
    t19 = 1/pow(t7,3.0);
    t20 = 1/pow(rhoa,8.666666666666666);
    t21 = 1/pow(rhoa,6.0);
    t22 = 1/pow(rhoa,3.333333333333333);
    t23 = pow(grada,2.98);
    t24 = 1/pow(t16,1.973333333333333);
    t25 = pow(grada,4.98);
    t26 = 1/pow(rhoa,7.333333333333333);
    t27 = 1/pow(rhoa,4.666666666666667);
    t28 = pow(grada,5.0);
    t29 = 1/pow(rhoa,7.666666666666667);
    t30 = pow(grada,3.98);
    t31 = 1/pow(rhoa,6.333333333333333);
    t32 = 1/pow(rhoa,6.666666666666667);
    t33 = 1/pow(rhoa,5.333333333333333);
    t34 = pow(grada,0.98);
    t35 = 1/pow(t7,4.0);

   /* code */
    res[0] = -1.4847442272E-4*t1*t4*t8*t9+0.001489099973306*
        t10*t11*t12/pow(t3,0.02)+0.00176436*t5*t4*t11*t13;
    res[1] = -0.00111682497998*t15*t11*t6*t17-0.00264654*
        grada*t4*t11*t2+1.1135581704000001E-4*t10*t4*t8*t14;
    res[2] = -0.007445499866532*t23*t11*t27*t17+3.341619758764805E-4*
        t25*t8*t26*t17+3.970933262150404E-5*t23*t11*t6*t24-0.00411684*
        t5*t4*t11*t22+9.403380105600001E-4*t1*t4*t8*t21-3.331845232195585E-5*
        t18*t4*t19*t20;
    res[3] = 0.005926617893759*t15*t11*t12*t17-2.506214819073604E-4*
        t30*t8*t31*t17+0.00352872*grada*t4*t11*t13-7.423721136000001E-4*
        t10*t4*t8*t9+2.4988839241466884E-5*t28*t4*t19*t29;
    res[4] = -0.00444496342032*t34*t11*t6*t17+1.8796611143052027E-4*
        t23*t8*t33*t17-0.00264654*t4*t11*t2+5.567790852E-4*t5*t4*t8*
        t14-1.8741629431100162E-5*t1*t4*t19*t32;
    res[5] = -1.1215257599188915E-5*t35*t4*pow(grada,8.0)/
        pow(rhoa,12.33333333333333)+1.1248159437583037E-4*t17*t19*
        pow(grada,6.98)/pow(rhoa,11.0)+4.997767848293376E-4*t18*t19*
        t4/pow(rhoa,9.666666666666666)-0.004344105686394*t17*t25*t8/
        pow(rhoa,8.333333333333334)-0.0061039484896*t1*t4*t8/pow(rhoa,
        7.0)+0.039378421516325*t11*t17*t23/pow(rhoa,5.666666666666667)+
        0.0137228*t11*t4*t5/pow(rhoa,4.333333333333333)-7.835974970643463E-5*
        t11*t23/(pow(t16,3.973333333333333)*pow(rhoa,1.666666666666667))+
        1.3366479035059231E-5*t25*t8*t31*t24-3.044382167648643E-4*
        t23*t11*t12*t24;
    res[6] = 8.411443199391686E-6*t35*t4*pow(grada,7.0)/pow(rhoa,
        11.33333333333333)-8.436119578187277E-5*t17*t19*pow(grada,
        5.98)/pow(rhoa,10.0)+0.001664126639865*t17*t26*t8*pow(grada,
        3.98)+1.1833381121208204E-4*t15*t11*t6*t24-3.341619758764808E-6*
        t30*t8*t33*t24-0.00823368*grada*t4*t11*t22+0.00410779236192*
        t10*t4*t8*t21-3.581733624610253E-4*t28*t4*t19*t20-0.025662156206647*
        t15*t11*t27*t17+0.001420188397475*t30*t8*t26*t17;
    res[7] = -6.308582399543764E-6*t18*t35*t4/pow(rhoa,10.33333333333333)+
        6.327089683640457E-5*t17*t19*t25/pow(rhoa,9.0)-0.00252406518624*
        t5*t4*t8*t9+2.498883924146688E-4*t1*t4*t19*t29-0.002122763951755*
        t23*t8*t31*t17+0.014712903376257*t34*t11*t12*t17+0.00352872*
        t4*t11*t13;
    res[8] = 4.731436799657823E-6*t28*t35*t4/pow(rhoa,9.333333333333334)-
        4.745317262730343E-5*t17*t19*t30/pow(rhoa,8.0)-0.006589714111873*
        t11*t17*t6/pow(grada,0.02)-1.6867466487990147E-4*t10*t4*t19*
        t32+0.001404106852386*t15*t8*t33*t17+0.00133626980448*grada*
        t4*t8*t14;

}

static void
dk87x_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real res[9];
 
    dk87x_third_helper(dp->rhoa, dp->grada, res);

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
        dk87x_third_helper(dp->rhob, dp->gradb, res);

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
dk87x_fourth_helper(real rhoa, real grada, real *res)
{
    real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    real t11, t12, t13, t14, t15, t16, t17, t18;
    real t19, t20, t21, t22, t23, t24, t25, t26;
    real t27, t28, t29, t30, t31, t32, t33, t34;
    real t35, t36, t37, t38, t39, t40, t41, t42;
    real t43, t44, t45, t46, t47, t48, t49, t50;
    real t51, t52, t53, t54, t55, t56, t57, t58;
    real t59;

    t1 = pow(grada,4.0);
    t2 = 1/pow(rhoa,1.333333333333333);
    t3 = grada*t2;
    t4 = 0.861213*pow(t3,0.98)+1.0;
    t5 = pow(grada,2.0);
    t6 = 1/pow(rhoa,2.666666666666667);
    t7 = 0.042076*t5*t6+1.0;
    t8 = 1/pow(t7,2.0);
    t9 = 1/pow(rhoa,5.0);
    t10 = pow(grada,3.0);
    t11 = 1/t7;
    t12 = 1/pow(rhoa,3.666666666666667);
    t13 = 1/pow(rhoa,2.333333333333333);
    t14 = 1/pow(rhoa,4.0);
    t15 = pow(grada,1.98);
    t16 = fabs(rhoa);
    t17 = pow(t16,0.026666666666667);
    t18 = pow(grada,6.0);
    t19 = 1/pow(t7,3.0);
    t20 = 1/pow(rhoa,8.666666666666666);
    t21 = 1/pow(rhoa,6.0);
    t22 = 1/pow(rhoa,3.333333333333333);
    t23 = pow(grada,2.98);
    t24 = 1/pow(t16,1.973333333333333);
    t25 = pow(grada,4.98);
    t26 = 1/pow(rhoa,7.333333333333333);
    t27 = 1/pow(rhoa,4.666666666666667);
    t28 = pow(grada,5.0);
    t29 = 1/pow(rhoa,7.666666666666667);
    t30 = pow(grada,3.98);
    t31 = 1/pow(rhoa,6.333333333333333);
    t32 = 1/pow(rhoa,6.666666666666667);
    t33 = 1/pow(rhoa,5.333333333333333);
    t34 = pow(grada,0.98);
    t35 = pow(grada,8.0);
    t36 = 1/pow(t7,4.0);
    t37 = 1/pow(rhoa,12.33333333333333);
    t38 = 1/pow(rhoa,9.666666666666666);
    t39 = 1/pow(rhoa,7.0);
    t40 = 1/pow(rhoa,4.333333333333333);
    t41 = 1/pow(rhoa,1.666666666666667);
    t42 = 1/pow(t16,3.973333333333333);
    t43 = pow(grada,6.98);
    t44 = 1/pow(rhoa,11.0);
    t45 = 1/pow(rhoa,8.333333333333334);
    t46 = 1/pow(rhoa,5.666666666666667);
    t47 = pow(grada,7.0);
    t48 = 1/pow(rhoa,11.33333333333333);
    t49 = pow(grada,5.98);
    t50 = 1/pow(rhoa,10.0);
    t51 = pow(grada,3.98);
    t52 = 1/pow(rhoa,10.33333333333333);
    t53 = 1/pow(rhoa,9.0);
    t54 = 1/pow(rhoa,9.333333333333334);
    t55 = 1/pow(rhoa,8.0);
    t56 = 1/pow(grada,0.02);
    t57 = 1/pow(t7,5.0);
    t58 = 1/pow(rhoa,10.66666666666667);
    t59 = 1/pow(rhoa,12.0);

   /* code */
    res[0] = -1.4847442272E-4*t1*t4*t8*t9+0.001489099973306*
        t10*t11*t12/pow(t3,0.02)+0.00176436*t5*t4*t11*t13;
    res[1] = -0.00111682497998*t15*t11*t6*t17-0.00264654*
        grada*t4*t11*t2+1.1135581704000001E-4*t10*t4*t8*t14;
    res[2] = -0.007445499866532*t23*t11*t27*t17+3.341619758764805E-4*
        t25*t8*t26*t17+3.970933262150404E-5*t23*t11*t6*t24-0.00411684*
        t5*t4*t11*t22+9.403380105600001E-4*t1*t4*t8*t21-3.331845232195585E-5*
        t18*t4*t19*t20;
    res[3] = 0.005926617893759*t15*t11*t12*t17-2.506214819073604E-4*
        t30*t8*t31*t17+0.00352872*grada*t4*t11*t13-7.423721136000001E-4*
        t10*t4*t8*t9+2.4988839241466884E-5*t28*t4*t19*t29;
    res[4] = -0.00444496342032*t34*t11*t6*t17+1.8796611143052027E-4*
        t23*t8*t33*t17-0.00264654*t4*t11*t2+5.567790852E-4*t5*t4*t8*
        t14-1.8741629431100162E-5*t1*t4*t19*t32;
    res[5] = 0.039378421516325*t23*t11*t46*t17-0.004344105686394*
        t25*t8*t45*t17+1.1248159437583037E-4*t43*t19*t44*t17-3.044382167648643E-4*
        t23*t11*t12*t24+1.3366479035059231E-5*t25*t8*t31*t24-7.835974970643463E-5*
        t23*t11*t41*t42+0.0137228*t5*t4*t11*t40-0.0061039484896*t1*
        t4*t8*t39+4.997767848293376E-4*t18*t4*t19*t38-1.1215257599188915E-5*
        t35*t4*t36*t37;
    res[6] = -0.025662156206647*t15*t11*t27*t17+0.001664126639865*
        t51*t8*t26*t17+0.001420188397475*t30*t8*t26*t17-8.436119578187277E-5*
        t49*t19*t50*t17+1.1833381121208204E-4*t15*t11*t6*t24-3.341619758764808E-6*
        t30*t8*t33*t24-0.00823368*grada*t4*t11*t22+0.00410779236192*
        t10*t4*t8*t21-3.581733624610253E-4*t28*t4*t19*t20+8.411443199391686E-6*
        t47*t4*t36*t48;
    res[7] = 0.014712903376257*t34*t11*t12*t17-0.002122763951755*
        t23*t8*t31*t17+6.327089683640457E-5*t25*t19*t53*t17+0.00352872*
        t4*t11*t13-0.00252406518624*t5*t4*t8*t9+2.498883924146688E-4*
        t1*t4*t19*t29-6.308582399543764E-6*t18*t4*t36*t52;
    res[8] = -0.006589714111873*t56*t11*t6*t17+0.001404106852386*
        t15*t8*t33*t17-4.745317262730343E-5*t30*t19*t55*t17+0.00133626980448*
        grada*t4*t8*t14-1.6867466487990147E-4*t10*t4*t19*t32+4.731436799657823E-6*
        t28*t4*t36*t54;
    res[9] = -5.033527239930376E-6*t4*t57*pow(grada,10.0)/
        pow(rhoa,16.0)+5.048293935954602E-5*t17*t36*pow(grada,8.98)/
        pow(rhoa,14.66666666666667)+3.065503743778303E-4*t35*t36*t4/
        pow(rhoa,13.33333333333333)+3.113494055002336E-4*t11*t23/(pow(t16,
        5.973333333333333)*pow(rhoa,0.666666666666667))-0.006200934182142*
        t18*t4*t19*t58+0.044267374181333*t1*t4*t8*t55+7.313576639267232E-4*
        t23*t11*t6*t42-3.516869150557807E-5*t25*t8*t33*t42-0.059465466666667*
        t5*t4*t11*t33+5.999018366710958E-6*t43*t19*t50*t24+0.002166364701906*
        t23*t11*t27*t24-2.3465596528215093E-4*t25*t8*t26*t24-0.002774545994604*
        t43*t19*t59*t17+0.047488129682891*t25*t8*t54*t17-0.238586906834203*
        t23*t11*t32*t17;
    res[10] = 3.775145429947782E-6*t4*t57*pow(grada,9.0)/
        pow(rhoa,15.0)-3.786220451965951E-5*t17*t36*pow(grada,7.98)/
        pow(rhoa,13.66666666666667)-2.335120541251752E-4*t15*t11*t41*
        t42+6.594129657295887E-6*t30*t8*t40*t42+0.0274456*grada*t4*
        t11*t40-0.025570595024*t10*t4*t8*t39+0.00402597965557*t28*
        t4*t19*t38-2.1589370878438658E-4*t47*t4*t36*t37-2.249631887516609E-6*
        t49*t19*t53*t24+6.656506559459498E-5*t51*t8*t31*t24+2.5619084817196862E-5*
        t30*t8*t31*t24-9.072258859592955E-4*t15*t11*t12*t24+0.12892958479992*
        t15*t11*t46*t17-0.021633646318243*t51*t8*t45*t17-0.008465436722204*
        t30*t8*t45*t17+0.001938057871096*t49*t19*t44*t17;
    res[11] = -2.831359072460837E-6*t35*t4*t57/pow(rhoa,14.0)+
        2.839665338974463E-5*t17*t36*t43/pow(rhoa,12.66666666666667)+
        0.006623224026662*t17*t26*t8*pow(grada,2.98)+1.493031167892024E-4*
        t18*t4*t36*t48+2.3430094619992243E-4*t34*t11*t6*t24+5.624079718791522E-7*
        t25*t19*t55*t24-2.3257673521003063E-5*t23*t8*t33*t24-0.00823368*
        t4*t11*t22+0.01301625772512*t5*t4*t8*t21-0.002482224697986*
        t1*t4*t19*t20-0.001325876793705*t25*t19*t50*t17-0.057760202497924*
        t34*t11*t27*t17+0.011278802090771*t23*t8*t26*t17;
    res[12] = 2.123519304345627E-6*t4*t47*t57/pow(rhoa,13.0)-
        2.1297490042308472E-5*t17*t36*t49/pow(rhoa,11.66666666666667)-
        0.00534507921792*grada*t4*t8*t9-1.0093731839270023E-4*t28*
        t4*t36*t52+0.001424363836764*t10*t4*t19*t29+3.150890662452948E-4*
        t51*t19*t53*t17+5.681726535909131E-4*t30*t19*t53*t17-0.009694239417362*
        t15*t8*t31*t17+0.017396845255344*t56*t11*t12*t17;
    res[13] = 1.3179428223745638E-4*t11*t17*t6/pow(grada,
        1.02)-1.5926394782592205E-6*t18*t4*t57*t59+6.624011519520952E-5*
        t1*t4*t36*t54-7.309235478129063E-4*t5*t4*t19*t32+1.5973117531731354E-5*
        t25*t36*t58*t17-5.675399446225491E-4*t23*t19*t55*t17+0.00446246585825*
        t34*t8*t33*t17+0.00133626980448*t4*t8*t14;

}

static void
dk87x_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real res[14];
 
    dk87x_fourth_helper(dp->rhoa, dp->grada, res);

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
        dk87x_fourth_helper(dp->rhob, dp->gradb, res);

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
