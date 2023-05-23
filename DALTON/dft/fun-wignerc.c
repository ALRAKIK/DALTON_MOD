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
/* fun-wignerc.c:

   Automatically generated code implementing WIGNERC functional and
   its derivatives. It is generated by func-codegen.pl being a part of
   a "Automatic code generation framework for analytical functional
   derivative evaluation", Pawel Salek, 2005

    This functional is connected by making following changes:
    1. add "extern Functional wignercFunctional;" to 'functionals.h'
    2. add "&wignercFunctional," to 'functionals.c'
    3. add "fun-wignerc.c" to 'Makefile.am', 'Makefile.in' or 'Makefile'.

    This functional has been generated from following input:
    ------ cut here -------
c:   0.04918;
d:   0.349;

rho:   rhoa + rhob;
rho13: rho^(1/3);
zeta:  (rhoa*rhob)/(rhoa+rhob);

aa: -4*c*zeta;
bb: (1+d/rho13);

K(rhoa,rhob,grada,gradb,gradab):=aa/bb;


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
static integer wignerc_isgga(void) { return 1; } /* FIXME: detect! */
static integer wignerc_read(const char *conf_line);
static real wignerc_energy(const FunDensProp* dp);
static void wignerc_first(FunFirstFuncDrv *ds,   real factor,
                         const FunDensProp* dp);
static void wignerc_second(FunSecondFuncDrv *ds, real factor,
                          const FunDensProp* dp);
static void wignerc_third(FunThirdFuncDrv *ds,   real factor,
                         const FunDensProp* dp);
static void wignerc_fourth(FunFourthFuncDrv *ds,   real factor,
                          const FunDensProp* dp);
 
Functional WignerFunctional = {
  "Wigner",       /* name */
  wignerc_isgga,   /* gga-corrected */
   1,
  wignerc_read,
  NULL,
  wignerc_energy,
  wignerc_first,
  wignerc_second,
  wignerc_third,
  wignerc_fourth
};
 
/* IMPLEMENTATION PART */
static integer
wignerc_read(const char *conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}

static real
wignerc_energy(const FunDensProp *dp)
{
    real res;
    real rhoa = dp->rhoa, rhob = dp->rhob;
    real grada = dp->grada, gradb = dp->gradb, gradab = dp->gradab;

    real t1;

    t1 = rhob+rhoa;

   /* code */
    res = -0.19672*rhoa*rhob/(t1*(0.349/pow(t1,0.333333333333333)+
        1.0));

    return res;
}

static void
wignerc_first(FunFirstFuncDrv *ds, real factor, const FunDensProp *dp)
{
    real dfdra, dfdrb, dfdga, dfdgb, dfdgab;
    real rhoa = dp->rhoa, rhob = dp->rhob;
    real grada = dp->grada, gradb = dp->gradb, gradab = dp->gradab;

    real t1, t2, t3, t4, t5, t6;

    t1 = rhob+rhoa;
    t2 = 0.349/pow(t1,0.333333333333333)+1.0;
    t3 = -0.022885093333333*rhoa*rhob/(pow(t1,2.333333333333333)*
        pow(t2,2.0));
    t4 = 1/t2;
    t5 = 0.19672*t4*rhoa*rhob/pow(t1,2.0);
    t6 = 1/t1;

   /* code */
    dfdra = -0.19672*rhob*t6*t4+t5+t3;
    dfdrb = -0.19672*rhoa*t6*t4+t5+t3;
    dfdga = 0.0;
    dfdgb = 0.0;
    dfdgab = 0.0;


    ds->df1000 += factor*dfdra;
    ds->df0100 += factor*dfdrb;
    ds->df0010 += factor*dfdga;
    ds->df0001 += factor*dfdgb;
    ds->df00001+= factor*dfdgab;
   
}

static void
wignerc_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real dfdra, dfdrb, dfdga, dfdgb, dfdgab;
    real d2fdrara, d2fdrarb, d2fdraga, d2fdragb, d2fdraab, d2fdrbrb,
        d2fdrbga, d2fdrbgb, d2fdrbgab, d2fdgaga, d2fdgagb, d2fdgagab,
        d2fdgbgb, d2fdgbgab, d2fdgabgab;
    real rhoa = dp->rhoa, rhob = dp->rhob;
    real grada = dp->grada, gradb = dp->gradb, gradab = dp->gradab;

    real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    real t11, t12;

    t1 = rhob+rhoa;
    t2 = 1/pow(t1,2.333333333333333);
    t3 = 0.349/pow(t1,0.333333333333333)+1.0;
    t4 = 1/pow(t3,2.0);
    t5 = -0.022885093333333*rhoa*rhob*t2*t4;
    t6 = 1/pow(t1,2.0);
    t7 = 1/t3;
    t8 = 0.19672*rhoa*rhob*t6*t7;
    t9 = 1/t1;
    t10 = -0.005324598382222*rhoa*rhob/(pow(t1,3.666666666666667)*
        pow(t3,3.0));
    t11 = 0.076283644444444*t4*rhoa*rhob/pow(t1,3.333333333333333);
    t12 = -
        0.39344*t7*rhoa*rhob/pow(t1,3.0);

   /* code */
    dfdra = -0.19672*rhob*t9*t7+t8+t5;
    dfdrb = -0.19672*rhoa*t9*t7+t8+t5;
    dfdga = 0.0;
    dfdgb = 0.0;
    dfdgab = 0.0;
    d2fdrara = 0.39344*rhob*t6*t7+t12-0.045770186666667*rhob*
        t2*t4+t11+t10;
    d2fdrarb = -0.19672*t9*t7+0.19672*rhob*t6*t7+0.19672*
        rhoa*t6*t7+t12-0.022885093333333*rhob*t2*t4-0.022885093333333*
        rhoa*t2*t4+t11+t10;
    d2fdraga = 0.0;
    d2fdragb = 0.0;
    d2fdraab = 0.0;
    d2fdrbrb = 0.39344*rhoa*t6*t7+t12-0.045770186666667*rhoa*
        t2*t4+t11+t10;
    d2fdrbga = 0.0;
    d2fdrbgb = 0.0;
    d2fdrbgab = 0.0;
    d2fdgaga = 0.0;
    d2fdgagb = 0.0;
    d2fdgagab = 0.0;
    d2fdgbgb = 0.0;
    d2fdgbgab = 0.0;
    d2fdgabgab = 0.0;


    ds->df1000 += factor*dfdra;
    ds->df0100 += factor*dfdrb;
    ds->df0010 += factor*dfdga;
    ds->df0001 += factor*dfdgb;
    ds->df00001+= factor*dfdgab;

    ds->df2000 += factor*d2fdrara;
    ds->df1100 += factor*d2fdrarb;
    ds->df1010 += factor*d2fdraga;
    ds->df1001 += factor*d2fdragb;
    ds->df10001+= factor*d2fdraab;
    ds->df0200 += factor*d2fdrbrb;
    ds->df0110 += factor*d2fdrbga;
    ds->df0101 += factor*d2fdrbgb;
    ds->df01001+= factor*d2fdrbgab;
    ds->df0020 += factor*d2fdgaga;
    ds->df0011 += factor*d2fdgagb;
    ds->df00101+= factor*d2fdgagab;
    ds->df0002 += factor*d2fdgbgb;
    ds->df00011+= factor*d2fdgbgab;
    ds->df00002+= factor*d2fdgabgab;

}

static void
wignerc_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real dfdra, dfdrb, dfdga, dfdgb, dfdgab;
    real d2fdrara, d2fdrarb, d2fdraga, d2fdragb, d2fdraab, d2fdrbrb,
        d2fdrbga, d2fdrbgb, d2fdrbgab, d2fdgaga, d2fdgagb, d2fdgagab,
        d2fdgbgb, d2fdgbgab, d2fdgabgab;
    real d3fdrarara, d3fdrararb, d3fdraraga, d3fdraragb, d3fdraraab,
         d3fdrarbrb, d3fdrarbga, d3fdrarbgb, d3fdrarbab, d3fdragaga,
         d3fdragagb, d3fdragaab, d3fdragbgb, d3fdragbab, d3fdraabab,
         d3fdrbrbrb, d3fdrbrbga, d3fdrbrbgb, d3fdrbrbab, d3fdrbgaga,
         d3fdrbgagb, d3fdrbgaab, d3fdrbgbgb, d3fdrbgbab, d3fdrbabab,
         d3fdgagaga, d3fdgagagb, d3fdgagaab, d3fdgagbgb, d3fdgagbab,
         d3fdgaabab, d3fdgbgbgb, d3fdgbgbab, d3fdgbabab, d3fdababab;
    real rhoa = dp->rhoa, rhob = dp->rhob;
    real grada = dp->grada, gradb = dp->gradb, gradab = dp->gradab;

    real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    real t11, t12, t13, t14, t15, t16, t17, t18;
    real t19, t20, t21, t22;

    t1 = rhob+rhoa;
    t2 = 1/pow(t1,2.333333333333333);
    t3 = 0.349/pow(t1,0.333333333333333)+1.0;
    t4 = 1/pow(t3,2.0);
    t5 = -0.022885093333333*rhoa*rhob*t2*t4;
    t6 = 1/pow(t1,2.0);
    t7 = 1/t3;
    t8 = 0.19672*rhoa*rhob*t6*t7;
    t9 = 1/t1;
    t10 = 1/pow(t1,3.666666666666667);
    t11 = 1/pow(t3,3.0);
    t12 = -0.005324598382222*rhoa*rhob*t10*t11;
    t13 = 1/pow(t1,3.333333333333333);
    t14 = 0.076283644444444*rhoa*rhob*t13*t4;
    t15 = 1/pow(t1,3.0);
    t16 = -0.39344*rhoa*rhob*t15*t7;
    t17 = -0.001858284835396*rhoa*rhob/(pow(t1,5.0)*pow(t3,
        4.0));
    t18 = 0.037272188675556*t11*rhoa*rhob/pow(t1,4.666666666666667);
    t19 = -
        0.300049001481481*t4*rhoa*rhob/pow(t1,4.333333333333333);
    t20 = 1.18032*
        t7*rhoa*rhob/pow(t1,4.0);
    t21 = -0.045770186666667*t2*t4;
    t22 = 0.39344*t6*t7;

   /* code */
    dfdra = -0.19672*rhob*t9*t7+t8+t5;
    dfdrb = -0.19672*rhoa*t9*t7+t8+t5;
    dfdga = 0.0;
    dfdgb = 0.0;
    dfdgab = 0.0;
    d2fdrara = 0.39344*rhob*t6*t7+t16-0.045770186666667*rhob*
        t2*t4+t14+t12;
    d2fdrarb = -0.19672*t9*t7+0.19672*rhob*t6*t7+0.19672*
        rhoa*t6*t7+t16-0.022885093333333*rhob*t2*t4-0.022885093333333*
        rhoa*t2*t4+t14+t12;
    d2fdraga = 0.0;
    d2fdragb = 0.0;
    d2fdraab = 0.0;
    d2fdrbrb = 0.39344*rhoa*t6*t7+t16-0.045770186666667*rhoa*
        t2*t4+t14+t12;
    d2fdrbga = 0.0;
    d2fdrbgb = 0.0;
    d2fdrbgab = 0.0;
    d2fdgaga = 0.0;
    d2fdgagb = 0.0;
    d2fdgagab = 0.0;
    d2fdgbgb = 0.0;
    d2fdgbgab = 0.0;
    d2fdgabgab = 0.0;
    d3fdrarara = -1.18032*rhob*t15*t7+t20+0.228850933333333*
        rhob*t13*t4+t19-0.015973795146667*rhob*t10*t11+t18+t17;
    d3fdrararb = t22-0.78688*rhob*t15*t7-0.39344*rhoa*t15*
        t7+t20+t21+0.152567288888889*rhob*t13*t4+0.076283644444444*
        rhoa*t13*t4+t19-0.010649196764444*rhob*t10*t11-0.005324598382222*
        rhoa*t10*t11+t18+t17;
    d3fdraraga = 0.0;
    d3fdraragb = 0.0;
    d3fdraraab = 0.0;
    d3fdrarbrb = t22-0.39344*rhob*t15*t7-0.78688*rhoa*t15*
        t7+t20+t21+0.076283644444444*rhob*t13*t4+0.152567288888889*
        rhoa*t13*t4+t19-0.005324598382222*rhob*t10*t11-0.010649196764444*
        rhoa*t10*t11+t18+t17;
    d3fdrarbga = 0.0;
    d3fdrarbgb = 0.0;
    d3fdrarbab = 0.0;
    d3fdragaga = 0.0;
    d3fdragagb = 0.0;
    d3fdragaab = 0.0;
    d3fdragbgb = 0.0;
    d3fdragbab = 0.0;
    d3fdraabab = 0.0;
    d3fdrbrbrb = -1.18032*rhoa*t15*t7+t20+0.228850933333333*
        rhoa*t13*t4+t19-0.015973795146667*rhoa*t10*t11+t18+t17;
    d3fdrbrbga = 0.0;
    d3fdrbrbgb = 0.0;
    d3fdrbrbab = 0.0;
    d3fdrbgaga = 0.0;
    d3fdrbgagb = 0.0;
    d3fdrbgaab = 0.0;
    d3fdrbgbgb = 0.0;
    d3fdrbgbab = 0.0;
    d3fdrbabab = 0.0;
    d3fdgagaga = 0.0;
    d3fdgagagb = 0.0;
    d3fdgagaab = 0.0;
    d3fdgagbgb = 0.0;
    d3fdgagbab = 0.0;
    d3fdgaabab = 0.0;
    d3fdgbgbgb = 0.0;
    d3fdgbgbab = 0.0;
    d3fdgbabab = 0.0;
    d3fdababab = 0.0;


    ds->df1000 += factor*dfdra;
    ds->df0100 += factor*dfdrb;
    ds->df0010 += factor*dfdga;
    ds->df0001 += factor*dfdgb;
    ds->df00001+= factor*dfdgab;

    ds->df2000 += factor*d2fdrara;
    ds->df1100 += factor*d2fdrarb;
    ds->df1010 += factor*d2fdraga;
    ds->df1001 += factor*d2fdragb;
    ds->df10001+= factor*d2fdraab;
    ds->df0200 += factor*d2fdrbrb;
    ds->df0110 += factor*d2fdrbga;
    ds->df0101 += factor*d2fdrbgb;
    ds->df01001+= factor*d2fdrbgab;
    ds->df0020 += factor*d2fdgaga;
    ds->df0011 += factor*d2fdgagb;
    ds->df00101+= factor*d2fdgagab;
    ds->df0002 += factor*d2fdgbgb;
    ds->df00011+= factor*d2fdgbgab;
    ds->df00002+= factor*d2fdgabgab;

    ds->df3000 += factor*d3fdrarara;
    ds->df2100  += factor*d3fdrararb;
    ds->df2010  += factor*d3fdraraga;
    ds->df2001  += factor*d3fdraragb;
    ds->df20001 += factor*d3fdraraab;
    ds->df1200  += factor*d3fdrarbrb;
    ds->df1110  += factor*d3fdrarbga;
    ds->df1101  += factor*d3fdrarbgb;
    ds->df11001 += factor*d3fdrarbab;
    ds->df1020  += factor*d3fdragaga;
    ds->df1011  += factor*d3fdragagb;
    ds->df10101 += factor*d3fdragaab;
    ds->df1002  += factor*d3fdragbgb;
    ds->df10011 += factor*d3fdragbab;
    ds->df10002 += factor*d3fdraabab;
    ds->df0300  += factor*d3fdrbrbrb;
    ds->df0210  += factor*d3fdrbrbga;
    ds->df0201  += factor*d3fdrbrbgb;
    ds->df02001 += factor*d3fdrbrbab;
    ds->df0120  += factor*d3fdrbgaga;
    ds->df0111  += factor*d3fdrbgagb;
    ds->df01101 += factor*d3fdrbgaab;
    ds->df0102  += factor*d3fdrbgbgb;
    ds->df01011 += factor*d3fdrbgbab;
    ds->df01002 += factor*d3fdrbabab;
    ds->df0030  += factor*d3fdgagaga;
    ds->df0021  += factor*d3fdgagagb;
    ds->df00201 += factor*d3fdgagaab;
    ds->df0012  += factor*d3fdgagbgb;
    ds->df00111 += factor*d3fdgagbab;
    ds->df00102 += factor*d3fdgaabab;
    ds->df0003  += factor*d3fdgbgbgb;
    ds->df00021 += factor*d3fdgbgbab;
    ds->df00012 += factor*d3fdgbabab;
    ds->df00003 += factor*d3fdababab;

}

static void
wignerc_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real dfdra, dfdrb, dfdga, dfdgb, dfdgab;
    real d2fdrara, d2fdrarb, d2fdraga, d2fdragb, d2fdraab, d2fdrbrb,
        d2fdrbga, d2fdrbgb, d2fdrbgab, d2fdgaga, d2fdgagb, d2fdgagab,
        d2fdgbgb, d2fdgbgab, d2fdgabgab;
    real d3fdrarara, d3fdrararb, d3fdraraga, d3fdraragb, d3fdraraab,
         d3fdrarbrb, d3fdrarbga, d3fdrarbgb, d3fdrarbab, d3fdragaga,
         d3fdragagb, d3fdragaab, d3fdragbgb, d3fdragbab, d3fdraabab,
         d3fdrbrbrb, d3fdrbrbga, d3fdrbrbgb, d3fdrbrbab, d3fdrbgaga,
         d3fdrbgagb, d3fdrbgaab, d3fdrbgbgb, d3fdrbgbab, d3fdrbabab,
         d3fdgagaga, d3fdgagagb, d3fdgagaab, d3fdgagbgb, d3fdgagbab,
         d3fdgaabab, d3fdgbgbgb, d3fdgbgbab, d3fdgbabab, d3fdababab;
    real d4fdrararara, d4fdrarararb, d4fdrararaga, d4fdrararagb,
         d4fdrararaab, d4fdrararbrb, d4fdrararbga, d4fdrararbgb, d4fdrararbab,
         d4fdraragaga, d4fdraragagb, d4fdraragaab, d4fdraragbgb, d4fdraragbab,
         d4fdraraabab, d4fdrarbrbrb, d4fdrarbrbga, d4fdrarbrbgb, d4fdrarbrbab,
         d4fdrarbgaga, d4fdrarbgagb, d4fdrarbgaab, d4fdrarbgbgb, d4fdrarbgbab,
         d4fdrarbabab, d4fdragagaga, d4fdragagagb, d4fdragagaab, d4fdragagbgb,
         d4fdragagbab, d4fdragaabab, d4fdragbgbgb, d4fdragbgbab, d4fdragbabab,
         d4fdraababab, d4fdrbrbrbrb, d4fdrbrbrbga, d4fdrbrbrbgb, d4fdrbrbrbab,
         d4fdrbrbgaga, d4fdrbrbgagb, d4fdrbrbgaab, d4fdrbrbgbgb, d4fdrbrbgbab,
         d4fdrbrbabab, d4fdrbgagaga, d4fdrbgagagb, d4fdrbgagaab, d4fdrbgagbgb,
         d4fdrbgagbab, d4fdrbgaabab, d4fdrbgbgbgb, d4fdrbgbgbab, d4fdrbgbabab,
         d4fdrbababab, d4fdgagagaga, d4fdgagagagb, d4fdgagagaab, d4fdgagagbgb,
         d4fdgagagbab, d4fdgagaabab, d4fdgagbgbgb, d4fdgagbgbab, d4fdgagbabab,
         d4fdgaababab, d4fdgbgbgbgb, d4fdgbgbgbab, d4fdgbgbabab, d4fdgbababab,
         d4fdabababab;
    real rhoa = dp->rhoa, rhob = dp->rhob;
    real grada = dp->grada, gradb = dp->gradb, gradab = dp->gradab;

    real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    real t11, t12, t13, t14, t15, t16, t17, t18;
    real t19, t20, t21, t22, t23, t24, t25, t26;
    real t27, t28, t29, t30, t31, t32, t33, t34;
    real t35;

    t1 = rhob+rhoa;
    t2 = 1/pow(t1,2.333333333333333);
    t3 = 0.349/pow(t1,0.333333333333333)+1.0;
    t4 = 1/pow(t3,2.0);
    t5 = -0.022885093333333*rhoa*rhob*t2*t4;
    t6 = 1/pow(t1,2.0);
    t7 = 1/t3;
    t8 = 0.19672*rhoa*rhob*t6*t7;
    t9 = 1/t1;
    t10 = 1/pow(t1,3.666666666666667);
    t11 = 1/pow(t3,3.0);
    t12 = -0.005324598382222*rhoa*rhob*t10*t11;
    t13 = 1/pow(t1,3.333333333333333);
    t14 = 0.076283644444444*rhoa*rhob*t13*t4;
    t15 = 1/pow(t1,3.0);
    t16 = -0.39344*rhoa*rhob*t15*t7;
    t17 = 1/pow(t1,5.0);
    t18 = 1/pow(t3,4.0);
    t19 = -0.001858284835396*rhoa*rhob*t17*t18;
    t20 = 1/pow(t1,4.666666666666667);
    t21 = 0.037272188675556*rhoa*rhob*t20*t11;
    t22 = 1/pow(t1,4.333333333333333);
    t23 = -0.300049001481481*rhoa*rhob*t22*t4;
    t24 = 1/pow(t1,4.0);
    t25 = 1.18032*rhoa*rhob*t24*t7;
    t26 = -0.045770186666667*t2*t4;
    t27 = 0.39344*t6*t7;
    t28 = -8.647218767373981E-4*rhoa*rhob/(pow(t1,6.333333333333333)*
        pow(t3,5.0));
    t29 = 0.022299418024747*t18*rhoa*rhob/pow(t1,6.0);
    t30 = -0.243748281497284*t11*rhoa*rhob/pow(t1,5.666666666666667);
    t31 = 1.437522899753086*
        t4*rhoa*rhob/pow(t1,5.333333333333333);
    t32 = -4.72128*rhoa*rhob*t17*t7;
    t33 = -0.015973795146667*t10*t11;
    t34 = 0.228850933333333*t13*t4;
    t35 = -1.18032*t15*t7;

   /* code */
    dfdra = -0.19672*rhob*t9*t7+t8+t5;
    dfdrb = -0.19672*rhoa*t9*t7+t8+t5;
    dfdga = 0.0;
    dfdgb = 0.0;
    dfdgab = 0.0;
    d2fdrara = 0.39344*rhob*t6*t7+t16-0.045770186666667*rhob*
        t2*t4+t14+t12;
    d2fdrarb = -0.19672*t9*t7+0.19672*rhob*t6*t7+0.19672*
        rhoa*t6*t7+t16-0.022885093333333*rhob*t2*t4-0.022885093333333*
        rhoa*t2*t4+t14+t12;
    d2fdraga = 0.0;
    d2fdragb = 0.0;
    d2fdraab = 0.0;
    d2fdrbrb = 0.39344*rhoa*t6*t7+t16-0.045770186666667*rhoa*
        t2*t4+t14+t12;
    d2fdrbga = 0.0;
    d2fdrbgb = 0.0;
    d2fdrbgab = 0.0;
    d2fdgaga = 0.0;
    d2fdgagb = 0.0;
    d2fdgagab = 0.0;
    d2fdgbgb = 0.0;
    d2fdgbgab = 0.0;
    d2fdgabgab = 0.0;
    d3fdrarara = -1.18032*rhob*t15*t7+t25+0.228850933333333*
        rhob*t13*t4+t23-0.015973795146667*rhob*t10*t11+t21+t19;
    d3fdrararb = t27-0.78688*rhob*t15*t7-0.39344*rhoa*t15*
        t7+t25+t26+0.152567288888889*rhob*t13*t4+0.076283644444444*
        rhoa*t13*t4+t23-0.010649196764444*rhob*t10*t11-0.005324598382222*
        rhoa*t10*t11+t21+t19;
    d3fdraraga = 0.0;
    d3fdraragb = 0.0;
    d3fdraraab = 0.0;
    d3fdrarbrb = t27-0.39344*rhob*t15*t7-0.78688*rhoa*t15*
        t7+t25+t26+0.076283644444444*rhob*t13*t4+0.152567288888889*
        rhoa*t13*t4+t23-0.005324598382222*rhob*t10*t11-0.010649196764444*
        rhoa*t10*t11+t21+t19;
    d3fdrarbga = 0.0;
    d3fdrarbgb = 0.0;
    d3fdrarbab = 0.0;
    d3fdragaga = 0.0;
    d3fdragagb = 0.0;
    d3fdragaab = 0.0;
    d3fdragbgb = 0.0;
    d3fdragbab = 0.0;
    d3fdraabab = 0.0;
    d3fdrbrbrb = -1.18032*rhoa*t15*t7+t25+0.228850933333333*
        rhoa*t13*t4+t23-0.015973795146667*rhoa*t10*t11+t21+t19;
    d3fdrbrbga = 0.0;
    d3fdrbrbgb = 0.0;
    d3fdrbrbab = 0.0;
    d3fdrbgaga = 0.0;
    d3fdrbgagb = 0.0;
    d3fdrbgaab = 0.0;
    d3fdrbgbgb = 0.0;
    d3fdrbgbab = 0.0;
    d3fdrbabab = 0.0;
    d3fdgagaga = 0.0;
    d3fdgagagb = 0.0;
    d3fdgagaab = 0.0;
    d3fdgagbgb = 0.0;
    d3fdgagbab = 0.0;
    d3fdgaabab = 0.0;
    d3fdgbgbgb = 0.0;
    d3fdgbgbab = 0.0;
    d3fdgbabab = 0.0;
    d3fdababab = 0.0;
    d4fdrararara = 4.72128*rhob*t24*t7+t32-1.200196005925926*
        rhob*t22*t4+t31+0.149088754702222*rhob*t20*t11+t30-0.007433139341582*
        rhob*t17*t18+t29+t28;
    d4fdrarararb = t35+3.54096*rhob*t24*t7+1.18032*rhoa*t24*
        t7+t32+t34-0.900147004444444*rhob*t22*t4-0.300049001481481*
        rhoa*t22*t4+t31+t33+0.111816566026667*rhob*t20*t11+0.037272188675556*
        rhoa*t20*t11+t30-0.005574854506187*rhob*t17*t18-0.001858284835396*
        rhoa*t17*t18+t29+t28;
    d4fdrararaga = 0.0;
    d4fdrararagb = 0.0;
    d4fdrararaab = 0.0;
    d4fdrararbrb = -1.57376*t15*t7+2.36064*rhob*t24*t7+2.36064*
        rhoa*t24*t7+t32+0.305134577777778*t13*t4-0.600098002962963*
        rhob*t22*t4-0.600098002962963*rhoa*t22*t4+t31-0.021298393528889*
        t10*t11+0.074544377351111*rhob*t20*t11+0.074544377351111*rhoa*
        t20*t11+t30-0.003716569670791*rhob*t17*t18-0.003716569670791*
        rhoa*t17*t18+t29+t28;
    d4fdrararbga = 0.0;
    d4fdrararbgb = 0.0;
    d4fdrararbab = 0.0;
    d4fdraragaga = 0.0;
    d4fdraragagb = 0.0;
    d4fdraragaab = 0.0;
    d4fdraragbgb = 0.0;
    d4fdraragbab = 0.0;
    d4fdraraabab = 0.0;
    d4fdrarbrbrb = t35+1.18032*rhob*t24*t7+3.54096*rhoa*t24*
        t7+t32+t34-0.300049001481481*rhob*t22*t4-0.900147004444444*
        rhoa*t22*t4+t31+t33+0.037272188675556*rhob*t20*t11+0.111816566026667*
        rhoa*t20*t11+t30-0.001858284835396*rhob*t17*t18-0.005574854506187*
        rhoa*t17*t18+t29+t28;
    d4fdrarbrbga = 0.0;
    d4fdrarbrbgb = 0.0;
    d4fdrarbrbab = 0.0;
    d4fdrarbgaga = 0.0;
    d4fdrarbgagb = 0.0;
    d4fdrarbgaab = 0.0;
    d4fdrarbgbgb = 0.0;
    d4fdrarbgbab = 0.0;
    d4fdrarbabab = 0.0;
    d4fdragagaga = 0.0;
    d4fdragagagb = 0.0;
    d4fdragagaab = 0.0;
    d4fdragagbgb = 0.0;
    d4fdragagbab = 0.0;
    d4fdragaabab = 0.0;
    d4fdragbgbgb = 0.0;
    d4fdragbgbab = 0.0;
    d4fdragbabab = 0.0;
    d4fdraababab = 0.0;
    d4fdrbrbrbrb = 4.72128*rhoa*t24*t7+t32-1.200196005925926*
        rhoa*t22*t4+t31+0.149088754702222*rhoa*t20*t11+t30-0.007433139341582*
        rhoa*t17*t18+t29+t28;
    d4fdrbrbrbga = 0.0;
    d4fdrbrbrbgb = 0.0;
    d4fdrbrbrbab = 0.0;
    d4fdrbrbgaga = 0.0;
    d4fdrbrbgagb = 0.0;
    d4fdrbrbgaab = 0.0;
    d4fdrbrbgbgb = 0.0;
    d4fdrbrbgbab = 0.0;
    d4fdrbrbabab = 0.0;
    d4fdrbgagaga = 0.0;
    d4fdrbgagagb = 0.0;
    d4fdrbgagaab = 0.0;
    d4fdrbgagbgb = 0.0;
    d4fdrbgagbab = 0.0;
    d4fdrbgaabab = 0.0;
    d4fdrbgbgbgb = 0.0;
    d4fdrbgbgbab = 0.0;
    d4fdrbgbabab = 0.0;
    d4fdrbababab = 0.0;
    d4fdgagagaga = 0.0;
    d4fdgagagagb = 0.0;
    d4fdgagagaab = 0.0;
    d4fdgagagbgb = 0.0;
    d4fdgagagbab = 0.0;
    d4fdgagaabab = 0.0;
    d4fdgagbgbgb = 0.0;
    d4fdgagbgbab = 0.0;
    d4fdgagbabab = 0.0;
    d4fdgaababab = 0.0;
    d4fdgbgbgbgb = 0.0;
    d4fdgbgbgbab = 0.0;
    d4fdgbgbabab = 0.0;
    d4fdgbababab = 0.0;
    d4fdabababab = 0.0;


    ds->df1000 += factor*dfdra;
    ds->df0100 += factor*dfdrb;
    ds->df0010 += factor*dfdga;
    ds->df0001 += factor*dfdgb;
    ds->df00001+= factor*dfdgab;

    ds->df2000 += factor*d2fdrara;
    ds->df1100 += factor*d2fdrarb;
    ds->df1010 += factor*d2fdraga;
    ds->df1001 += factor*d2fdragb;
    ds->df10001+= factor*d2fdraab;
    ds->df0200 += factor*d2fdrbrb;
    ds->df0110 += factor*d2fdrbga;
    ds->df0101 += factor*d2fdrbgb;
    ds->df01001+= factor*d2fdrbgab;
    ds->df0020 += factor*d2fdgaga;
    ds->df0011 += factor*d2fdgagb;
    ds->df00101+= factor*d2fdgagab;
    ds->df0002 += factor*d2fdgbgb;
    ds->df00011+= factor*d2fdgbgab;
    ds->df00002+= factor*d2fdgabgab;

    ds->df3000 += factor*d3fdrarara;
    ds->df2100  += factor*d3fdrararb;
    ds->df2010  += factor*d3fdraraga;
    ds->df2001  += factor*d3fdraragb;
    ds->df20001 += factor*d3fdraraab;
    ds->df1200  += factor*d3fdrarbrb;
    ds->df1110  += factor*d3fdrarbga;
    ds->df1101  += factor*d3fdrarbgb;
    ds->df11001 += factor*d3fdrarbab;
    ds->df1020  += factor*d3fdragaga;
    ds->df1011  += factor*d3fdragagb;
    ds->df10101 += factor*d3fdragaab;
    ds->df1002  += factor*d3fdragbgb;
    ds->df10011 += factor*d3fdragbab;
    ds->df10002 += factor*d3fdraabab;
    ds->df0300  += factor*d3fdrbrbrb;
    ds->df0210  += factor*d3fdrbrbga;
    ds->df0201  += factor*d3fdrbrbgb;
    ds->df02001 += factor*d3fdrbrbab;
    ds->df0120  += factor*d3fdrbgaga;
    ds->df0111  += factor*d3fdrbgagb;
    ds->df01101 += factor*d3fdrbgaab;
    ds->df0102  += factor*d3fdrbgbgb;
    ds->df01011 += factor*d3fdrbgbab;
    ds->df01002 += factor*d3fdrbabab;
    ds->df0030  += factor*d3fdgagaga;
    ds->df0021  += factor*d3fdgagagb;
    ds->df00201 += factor*d3fdgagaab;
    ds->df0012  += factor*d3fdgagbgb;
    ds->df00111 += factor*d3fdgagbab;
    ds->df00102 += factor*d3fdgaabab;
    ds->df0003  += factor*d3fdgbgbgb;
    ds->df00021 += factor*d3fdgbgbab;
    ds->df00012 += factor*d3fdgbabab;
    ds->df00003 += factor*d3fdababab;

    ds->df4000  += factor*d4fdrararara;
    ds->df3100  += factor*d4fdrarararb;
    ds->df3010  += factor*d4fdrararaga;
    ds->df3001  += factor*d4fdrararagb;
    ds->df30001 += factor*d4fdrararaab;
    ds->df2200  += factor*d4fdrararbrb;
    ds->df2110  += factor*d4fdrararbga;
    ds->df2101  += factor*d4fdrararbgb;
    ds->df21001 += factor*d4fdrararbab;
    ds->df2020  += factor*d4fdraragaga;
    ds->df2011  += factor*d4fdraragagb;
    ds->df20101 += factor*d4fdraragaab;
    ds->df2002  += factor*d4fdraragbgb;
    ds->df20011 += factor*d4fdraragbab;
    ds->df20002 += factor*d4fdraraabab;
    ds->df1300  += factor*d4fdrarbrbrb;
    ds->df1210  += factor*d4fdrarbrbga;
    ds->df1201  += factor*d4fdrarbrbgb;
    ds->df12001 += factor*d4fdrarbrbab;
    ds->df1120  += factor*d4fdrarbgaga;
    ds->df1111  += factor*d4fdrarbgagb;
    ds->df11101 += factor*d4fdrarbgaab;
    ds->df1102  += factor*d4fdrarbgbgb;
    ds->df11011 += factor*d4fdrarbgbab;
    ds->df11002 += factor*d4fdrarbabab;
    ds->df1030  += factor*d4fdragagaga;
    ds->df1021  += factor*d4fdragagagb;
    ds->df10201 += factor*d4fdragagaab;
    ds->df1012  += factor*d4fdragagbgb;
    ds->df10111 += factor*d4fdragagbab;
    ds->df10102 += factor*d4fdragaabab;
    ds->df1003  += factor*d4fdragbgbgb;
    ds->df10021 += factor*d4fdragbgbab;
    ds->df10012 += factor*d4fdragbabab;
    ds->df10003 += factor*d4fdraababab;
    ds->df0400  += factor*d4fdrbrbrbrb;
    ds->df0310  += factor*d4fdrbrbrbga;
    ds->df0301  += factor*d4fdrbrbrbgb;
    ds->df03001 += factor*d4fdrbrbrbab;
    ds->df0220  += factor*d4fdrbrbgaga;
    ds->df0211  += factor*d4fdrbrbgagb;
    ds->df02101 += factor*d4fdrbrbgaab;
    ds->df0202  += factor*d4fdrbrbgbgb;
    ds->df02011 += factor*d4fdrbrbgbab;
    ds->df02002 += factor*d4fdrbrbabab;
    ds->df0130  += factor*d4fdrbgagaga;
    ds->df0121  += factor*d4fdrbgagagb;
    ds->df01201 += factor*d4fdrbgagaab;
    ds->df0112  += factor*d4fdrbgagbgb;
    ds->df01111 += factor*d4fdrbgagbab;
    ds->df01102 += factor*d4fdrbgaabab;
    ds->df0103  += factor*d4fdrbgbgbgb;
    ds->df01021 += factor*d4fdrbgbgbab;
    ds->df01012 += factor*d4fdrbgbabab;
    ds->df01003 += factor*d4fdrbababab;
    ds->df0040  += factor*d4fdgagagaga;
    ds->df0031  += factor*d4fdgagagagb;
    ds->df00301 += factor*d4fdgagagaab;
    ds->df0022  += factor*d4fdgagagbgb;
    ds->df00211 += factor*d4fdgagagbab;
    ds->df00202 += factor*d4fdgagaabab;
    ds->df0013  += factor*d4fdgagbgbgb;
    ds->df00121 += factor*d4fdgagbgbab;
    ds->df00112 += factor*d4fdgagbabab;
    ds->df00103 += factor*d4fdgaababab;
    ds->df0004  += factor*d4fdgbgbgbgb;
    ds->df00031 += factor*d4fdgbgbgbab;
    ds->df00022 += factor*d4fdgbgbabab;
    ds->df00013 += factor*d4fdgbababab;
    ds->df00004 += factor*d4fdabababab;

}
