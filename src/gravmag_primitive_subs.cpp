/*
 * gravmag_primitive_subs.cpp -- primitive vertical prism gravmag subroutines
 *
 * These subroutines include:
 *   rect_prism_grav1()
 *   rect_prism_mag1()
 *   add_tangents()
 *
 *
 * CITATION
 * =======
 *
 * To cite this software in publications, please use:
 *
 *  Cronkite-Ratcliff, C., Phelps, G., Scheirer, D., 2022, gravmagsubs:
 *  Gravitational and magnetic attraction of 3-D vertical rectangular
 *  prisms, U.S. Geological Survey software release, version 1.0,
 *  https://doi.org/10.5066/P9HBSPRM
 *
 *
 */


#include "gravmag_primitive_subs.h"
#include <Rcpp.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cfloat>

using namespace Rcpp;
/******************************************************************************/
/************************* gravmag prism subroutines **************************/
/******************************************************************************/


/************************ rect_prism_grav1() ************************/
/*
 * Mimics the subroutine gpris() by Plouff (1975a).
 *
 * This subroutine calculates the vertical component of gravitational
 * attraction of a rectangular prism at any fieldpoint location (Plouff 1975a).
 *
 * Right-handed coordinate system, RHCS, utilized in this code is:
 *   +X=north, +Y=east,  +Z=down
 * Differs from left-handed coordinate system, LHCS, of calling program:
 *   +X=east,  +Y=north, +Z=up
 *
 * Input XYZ distances are in km; delta_rho is in g/cc.
 *
 * Returns vertical gravitational acceleration in mGal.
 * Return codes: 0=successful completion; <0=ERROR
 *
 * Utilizes subroutine: add_tangents()
 *
 * Calling variables:
 *      xstation          - double   - x-coordinate of station, in km, positive east
 *      ystation          - double   - y-coordinate of station, in km, positive north
 *      zstation          - double   - z-coordinate of station, in km, positive up
 *      xmin              - double   - x-coordinate minimum of prism, in km, positive east
 *      xmax              - double   - x-coordinate maximum of prism, in km, positive east
 *      ymin              - double   - y-coordinate minimum of prism, in km, positive north
 *      ymax              - double   - y-coordinate maximum of prism, in km, positive north
 *      z_deep            - double   - z-coordinate of bottom of prism, in km, positive up
 *      z_shallow         - double   - z-coordinate of top    of prism, in km, positive up
 *      delta_rho         - double   - density contrast, in g/cc
 *      *gravanom         - double * - gravity anomaly due to rectangular prism, in mGal
 *      subverbose        - int      - verbose/debug flag
 *
 */

int    rect_prism_grav1(double xstation, double ystation, double zstation,
                        double xmin, double xmax, double ymin, double ymax,
                        double z_deep, double z_shallow,
                        double delta_rho, double *gravanom, int subverbose)

{
    char   subname[NDIMSUBNAME];
    double z1, z2;
    double x1, y1, x2, y2;
    double r111, r121, r112, r122, r211, r221, r212, r222;
    double anom_add=0.;
    
    double tdw1=0., tdt1=0.;
    double tdw2=0., tdt2=0.;
    
    (void)strcpy(subname, "rect_prism_grav1\0");
    
    *gravanom = 0.;
    if (delta_rho == 0.) { /* no anomaly */
        if (subverbose > 1) Rprintf("%s: WARNING: delta_rho=%f; skipping calculations\n", subname, delta_rho);
        return 0;
    }
    if (fabs(z_deep - z_shallow) < eps) { /* no anomaly -- zero thickness */
        if (subverbose > 0) Rprintf("%s: WARNING: z_deep, %f, == z_shallow, %f; skipping calculations\n", subname, z_deep, z_shallow);
        return 0;
    }
    if (fabs(xmin - xmax) < eps) { /* no anomaly -- zero dimension */
        if (subverbose > 0) Rprintf("%s: WARNING: xmin, %f, == xmax, %f; skipping calculations\n", subname, xmin, xmax);
        return 0;
    }
    if (fabs(ymin - ymax) < eps) { /* no anomaly -- zero dimension */
        if (subverbose > 0) Rprintf("%s: WARNING: ymin, %f, == ymax, %f; skipping calculations\n", subname, ymin, ymax);
        return 0;
    }
    
    /* Convert LHCS to RHCS centered on observation point -- +Z=down */
    z1 = zstation - z_shallow;
    z2 = zstation - z_deep;
    
    if (fabs(z1 + z2) < eps) { /* no anomaly -- equal mass above and below */
        if (subverbose > 0) Rprintf("%s: WARNING: zstation, %f, is mid-way between z_deep, %f, and z_shallow, %f; skipping calculations\n", subname, zstation, z_deep, z_shallow);
        return 0;
    }
    
    if (subverbose > 1) {
        Rprintf("%s: DEBUG: zstation, z_shallow, z_deep are: %g %g %g\n",
                      subname, zstation, z_shallow, z_deep);
        Rprintf("%s: DEBUG: z1, z2 are: %g %g\n", subname, z1, z2);
    }
    
    /* Convert LHCS to RHCS centered on observation point -- +X=north, +Y=east */
    y1 = xmin - xstation;
    x1 = ymin - ystation;
    y2 = xmax - xstation;
    x2 = ymax - ystation;
    
    if (fabs(x1) < eps) x1 = 0.;
    if (fabs(x2) < eps) x2 = 0.;
    if (fabs(y1) < eps) y1 = 0.;
    if (fabs(y2) < eps) y2 = 0.;
    if (fabs(z1) < eps) z1 = 0.;
    if (fabs(z2) < eps) z2 = 0.;
    
    r111 = sqrt(x1*x1+y1*y1+z1*z1);
    r121 = sqrt(x1*x1+y2*y2+z1*z1);
    r112 = sqrt(x1*x1+y1*y1+z2*z2);
    r122 = sqrt(x1*x1+y2*y2+z2*z2);
    r211 = sqrt(x2*x2+y1*y1+z1*z1);
    r221 = sqrt(x2*x2+y2*y2+z1*z1);
    r212 = sqrt(x2*x2+y1*y1+z2*z2);
    r222 = sqrt(x2*x2+y2*y2+z2*z2);
    
    if (x1 != 0.) anom_add += x1 * log ( (y1+r111) * (y2+r122) / ( (y1+r112) * (y2+r121) ) );
    if (x2 != 0.) anom_add += x2 * log ( (y1+r212) * (y2+r221) / ( (y2+r222) * (y1+r211) ) );
    if (y1 != 0.) anom_add += y1 * log ( (x1+r111) * (x2+r212) / ( (x1+r112) * (x2+r211) ) );
    if (y2 != 0.) anom_add += y2 * log ( (x1+r122) * (x2+r221) / ( (x2+r222) * (x1+r121) ) );
    
    if (z1 != 0.) {
        if (x2*y2 != 0.) add_tangents(-x2*y2/(z1*r221), &tdw1, &tdt1);
        if (x1*y2 != 0.) add_tangents( x1*y2/(z1*r121), &tdw1, &tdt1);
        if (x2*y1 != 0.) add_tangents( x2*y1/(z1*r211), &tdw1, &tdt1);
        if (x1*y1 != 0.) add_tangents(-x1*y1/(z1*r111), &tdw1, &tdt1);
        tdt1 += atan(tdw1);
        anom_add += (z1 * tdt1);
    }
    
    if (z2 != 0.) {
        if (x2*y2 != 0.) add_tangents( x2*y2/(z2*r222), &tdw2, &tdt2);
        if (x1*y2 != 0.) add_tangents(-x1*y2/(z2*r122), &tdw2, &tdt2);
        if (x2*y1 != 0.) add_tangents(-x2*y1/(z2*r212), &tdw2, &tdt2);
        if (x1*y1 != 0.) add_tangents( x1*y1/(z2*r112), &tdw2, &tdt2);
        tdt2 += atan(tdw2);
        anom_add += (z2 * tdt2);
    }
    
    if (subverbose > 1) {
        Rprintf("%s: DEBUG: z1*tdt1, z2*tdt2 are: %g %g\n", subname, z1*tdt1, z2*tdt2);
    }
    
    if (subverbose > 1) {
        Rprintf("%s: DEBUG: xstation, ystation, zstation, anom_add are: %f %f %f %f\n",
                      subname, xstation, ystation, zstation, anom_add);
    }
    
    anom_add *= (BigGx1e11 * delta_rho);
    
    if (subverbose > 2) Rprintf("%s: DEBUG: anom_add is: %f\n", subname, anom_add);
    
    *gravanom = anom_add;
    if (!std::isfinite(*gravanom)) {
        (*gravanom) = NAN;
        return -10;
    }
    
    return 0 ;
} /* end of rect_prism_grav1() */
/************************ rect_prism_grav1() ************************/


/************************ rect_prism_mag1() *************************/
/*
 * Mimics the subroutine magw() by Plouff (1975b).
 *
 * Calculates the total magnetic anomaly generated by a vertical rectangular
 * prism for any observation point.
 *
 * Calling variables:
 *      xstation             - double   - x-coordinate of station, in km, positive east
 *      ystation             - double   - y-coordinate of station, in km, positive north
 *      zstation             - double   - z-coordinate of station, in km, positive up
 *      xmin                 - double   - x-coordinate minimum of prism, in km, positive east
 *      xmax                 - double   - x-coordinate maximum of prism, in km, positive east
 *      ymin                 - double   - y-coordinate minimum of prism, in km, positive north
 *      ymax                 - double   - y-coordinate maximum of prism, in km, positive north
 *      z_deep               - double   - z-coordinate of bottom of prism, in km, positive up
 *      z_shallow            - double   - z-coordinate of top    of prism, in km, positive up
 *      magnetization_intens - double   - magnetization of body, in A/m (both induced and remanent)
 *      magnetization_incl   - double   - magnetization inclination of body, in deg, positive below horizontal
 *      magnetization_decl   - double   - magnetization declination of body, in deg, positive east of true north
 *      field_total          - double   - Earth's field intensity, in nT
 *      field_incl           - double   - Earth's field inclination, in deg, positive below horizon
 *      field_decl           - double   - Earth's field declination, in deg, positive east of true north
 *      *maganom             - double * - magnetic anomaly due to rectangular prism, in nT
 *      subverbose           - int      - verbose/debug flag
 *
 *
 * Right-handed coordinate system, RHCS, utilized in this code is:
 *   +X=north, +Y=east,  +Z=down
 * Differs from left-handed coordinate system, LHCS, of calling program:
 *   +X=east,  +Y=north, +Z=up
 *
 * Input XYZ distances are in km; magnetization_intens is in A/m.
 *
 * Returns total field magnetic anomaly in nT.
 * Return codes: 0=successful completion; <0=ERROR
 *
 * Return codes: 0=successful completion; <0=ERROR
 *
 * Utilizes subroutine: add_tangents()
 *
 * N.B. Demagnetization effects are ignored in this subroutine.
 *
 * FUTURE: Activate the formula for when disturbing field is large
 * FUTURE: with respect to ambient field when this is beneficial.
 *
 */

int rect_prism_mag1(double xstation, double ystation, double zstation,
                    double xmin, double xmax, double ymin, double ymax, double z_deep, double z_shallow,
                    double magnetization_intens, double magnetization_incl, double magnetization_decl,
                    double field_total, double field_incl, double field_decl,
                    double *maganom, int subverbose)

{
    char   subname[NDIMSUBNAME];
    int    i, j, k;
    
    double x[2], y[2], z[2] ;
    double bs[2], hs[2];
    double anom_add=0.;
    
    double g1, g2, g3, t1, t2, t3, w1, w2;
    /*double _g2, __g2, ___g2;*/
    double ft1, ft2;
    double fg1, fg2, fg3;
    double ai, bj, hk, r;
    double pl, pm, pn;
    /*double xa, yb ; */
    double ad, aa, ab, bb, aabb, bd, rr;
    double hx, hy, hz;
    double cos_magnetization_decl, sin_magnetization_decl, cos_magnetization_incl, sin_magnetization_incl;
    double cos_field_decl,  sin_field_decl,  cos_field_incl,  sin_field_incl;
    
    int    isa[2], isb[2];
    int    ija, ijb, isai, isbj;
    int    igna, isgn, igab;
    int    iga[2], igb[2], igh[2];
    int    sgh=1;
    
    iga[0]=igb[0]=igh[0]=-1;
    iga[1]=igb[1]=igh[1]= 1;
    
    
    (void)strcpy(subname, "rect_prism_mag1\0");
    
    if (magnetization_intens == 0.) /* no anomaly */
    {
        if (subverbose > 1) Rprintf("%s: WARNING: magnetization_intens=%f; skipping calculations\n", subname, magnetization_intens);
        return 0;
    }
    if (z_deep == z_shallow) /* no anomaly -- zero thickness */
    {
        if (subverbose > 0) Rprintf("%s: WARNING: z_deep, %f, == z_shallow, %f; skipping calculations\n", subname, z_deep, z_shallow);
        return 0;
    }
    if (xmin == xmax) /* no anomaly -- zero dimension */
    {
        if (subverbose > 0) Rprintf("%s: WARNING: xmin, %f, == xmax, %f; skipping calculations\n", subname, xmin, xmax);
        return 0;
    }
    if (ymin == ymax) /* no anomaly -- zero dimension */
    {
        if (subverbose > 0) Rprintf("%s: WARNING: ymin, %f, == ymax, %f; skipping calculations\n", subname, ymin, ymax);
        return 0;
    }

    /* Convert LHCS to RHCS centered on observation point -- +Z=down */
    y[0] = xmin - xstation;
    x[0] = ymin - ystation;
    y[1] = xmax - xstation;
    x[1] = ymax - ystation;
    z[0] = zstation - z_shallow;
    z[1] = zstation - z_deep;
    if (fabs(x[0]) < eps) x[0] = 0;
    if (fabs(x[1]) < eps) x[1] = 0;
    if (fabs(y[0]) < eps) y[0] = 0;
    if (fabs(y[1]) < eps) y[1] = 0;
    if (fabs(z[0]) < eps) z[0] = 0;
    if (fabs(z[1]) < eps) z[1] = 0;
    
    if (subverbose > 1) {
        Rprintf("%s: DEBUG: x[0], x[1] are: %g %g\n", subname, x[0], x[1]);
        Rprintf("%s: DEBUG: y[0], y[1] are: %g %g\n", subname, y[0], y[1]);
        Rprintf("%s: DEBUG: z[0], z[1] are: %g %g\n", subname, z[0], z[1]);
    }
    
    if (xstation >= xmin && xstation <= xmax &&
        ystation >= ymin && ystation <= ymax && z[0] < 0. && z[1] > 0. ) {
        *maganom = NAN;
        return -1;
    }
    
    if (subverbose > 1) {
       Rprintf("%s: DEBUG: zstation, z_deep, z_shallow are: %g %g %g\n",
                     subname, zstation, z_deep, z_shallow);
       Rprintf("%s: DEBUG: z[0], z[1] are: %g %g\n", subname, z[0], z[1]);
    }
    
    /* components of earth's field */
    sin_field_decl = sin(field_decl*RAD_PER_DEG);
    sin_field_incl = sin(field_incl*RAD_PER_DEG);
    cos_field_decl = cos(field_decl*RAD_PER_DEG);
    cos_field_incl = cos(field_incl*RAD_PER_DEG);
   
    /* components of anomalous magnetization */
    sin_magnetization_decl = sin(magnetization_decl*RAD_PER_DEG);
    sin_magnetization_incl = sin(magnetization_incl*RAD_PER_DEG);
    cos_magnetization_decl = cos(magnetization_decl*RAD_PER_DEG);
    cos_magnetization_incl = cos(magnetization_incl*RAD_PER_DEG);

    pl = magnetization_intens * cos_magnetization_decl * cos_magnetization_incl;
    pm = magnetization_intens * sin_magnetization_decl * cos_magnetization_incl;
    pn = magnetization_intens * sin_magnetization_incl;
    
    for (j=0; j<2; j++) {
        ai = x[j];
        bj = y[j];

        //if (fabs(ai) < eps) ai = 0;
        //if (fabs(bj) < eps) bj = 0;

        if (ai < 0.) {
            isa[j] = -1;
        } else if (ai == 0.) {
            isa[j] =  0;
        } else {
            isa[j] =  1;
        }

        if (bj < 0.) {
            isb[j] = -1;
        } else if (bj == 0.) {
            isb[j] =  0;
        } else {
            isb[j] =  1;
        }
    
        x[j]  = ai;
        y[j]  = bj;
        bs[j] = bj * bj;
        hs[j] = z[j] * z[j];
    } /* end for, j */
    
    ija = (isa[0]+isa[1]) / 2;
    ijb = (isb[0]+isb[1]) / 2;
 
    g1 = g2 = g3 = 1.;
    t1 = t2 = 0.;
    w1 = w2 = 0.;
    
    for (i=0; i<2; i++) {
        ai = x[i];
        isai = isa[i];
        ad = fabs(ai);
        igna = iga[i];
        aa = ai * ai;

        for (j=0; j<2; j++) {
            bj = y[j];
            isbj = isb[j];
            ab = ai * bj;
            bb = bs[j];
            aabb = aa + bb;
            igab = igna * igb[j];
            bd = fabs(bj);

            for (k=0; k<2; k++) {
                isgn = igab * igh[k];
                ft1 = ft2 = 0.;
                hk = z[k];
                rr = aabb + hs[k];
                if (rr == 0.) {
                    *maganom = NAN;
                    return -11;
                }
                r = sqrt(rr);
                fg1 = bj + r;
                fg2 = hk + r;
                fg3 = ai + r;
           
                if (hk == 0.) {
                    if (isbj == 0) {
                        if (ija < 0) {
                            fg3 = 1. / ad;
                        } else if (ija == 0) {
                            *maganom = NAN;
                            return -12;
                        }
                    } /* end if */
           
                    if (isai == 0) {
                        if (ijb < 0) {
                            fg1 = 1. / bd;
                        } else if (ijb == 0) {
                            *maganom = NAN;
                            return -12;
                        }
                    } /* end if */
                } else { /* end if -- hk==0. */
                    if (ab != 0.) {
                        ft1 = bj * hk / (ai * r);
                        ft2 = ai * hk / (bj * r);
                    }
                } /* end else */
           
                if (isgn != 1) { /* negative arctangent and reciprocal log terms */
                    g1 /= fg1;
                    g2 /= fg2;
                    g3 /= fg3;
                    ft1 *= -1.;
                    ft2 *= -1.;
                } else {
                    g1 *= fg1;
                    g2 *= fg2;
                    g3 *= fg3;
                }
           
                if (ab != 0.) {
                    add_tangents(ft1, &w1, &t1);
                    add_tangents(ft2, &w2, &t2);
                }
                //Rcout << "rr=" << rr << " hk=" << hk << " fg2=" << fg2 << " g2=" << g2 << "\n";
            } /* end for, k */
        } /* end for, j */
    } /* end for, i */
    
    /*_g2 = g2;*/
    
    g1 = log(g1);
    g2 = log(g2);
    g3 = log(g3);
    t1 = -1. * (atan(w1) + t1);
    t2 = -1. * (atan(w2) + t2);
    t3 = -1. * (t1 + t2);
    
    hx = sgh * (pl*t1 + pm*g2 + pn*g1);
    hy = sgh * (pl*g2 + pm*t2 + pn*g3);
    hz = sgh * (pl*g1 + pm*g3 + pn*t3);
    
    if (subverbose > 1) {
        Rprintf("%s: DEBUG: pl, pm, pn are: %g %g %g\n", subname, pl, pm, pn);
        Rprintf("%s: DEBUG: g1, g2, g3, t1, t2, t3 are: %g %g %g %g %g %g\n", subname, g1, g2, g3, t1, t2, t3);
        Rprintf("%s: DEBUG: hx, hy, hz are: %g %g %g\n", subname, hx, hy, hz);
    }
    
    /* The 1.e2 multiplier is the result of multiplying a factor Cm=1.e-7 for this SI
     * calculation (see Blakely, 1995, p.67) by a factor of 1.e9 to convert T to nT.
     */
    hx *= 1.e2;
    hy *= 1.e2;
    hz *= 1.e2;
    
    /* Simple formula assumes disturbing field small w.r.t. ambient field. */
    anom_add = hx * cos_field_decl * cos_field_incl +
               hy * sin_field_decl * cos_field_incl +
               hz * sin_field_incl;
    
    /* Formula for disturbing field large w.r.t. ambient field. */
    /* !!!inactive!!!
    Rprintf("1: anom_add = %g\n", anom_add) ;
    anom_add = sqrt( pow((hx + cos_field_decl * cos_field_incl * field_total), 2.) +
                     pow((hy + sin_field_decl * cos_field_incl * field_total), 2.) +
                     pow((hz + sin_field_incl * field_total), 2.) ) - field_total ;
    Rprintf("2: anom_add = %g\n", anom_add) ;
    Rprintf("2: field_total = %g\n", field_total) ;
    */
    
    *maganom = anom_add;
    if(!std::isfinite(*maganom)) {
        *maganom = NAN;
        return -10;
    }
   
    if (subverbose > 1) Rprintf("%s: DEBUG: *maganom is: %f\n", subname, *maganom);
    
    return 0;
} /* end of rect_prism_mag1() */
/************************ rect_prism_mag1() *************************/


/************************ add_tangents() ****************************/
/*
 * Based on the subroutine tandet() from Plouff (1975a).
 *
 * This is a subroutine for adding tangents-of-angles.
 * Negative input_tangent for subtraction.
 * Keeps track of PI & PI/2's in sum_extra so that a simple
 *    atan()of the sum_tangent yields the correct value.
 *
 *  Plouff_Code --> This_Code
 *  -----------     ---------
 *  tandet()    --> add_tangents()
 *  f           --> input_tangent
 *  zze         --> sum_tangent
 *  sume        --> sum_extra
 *
 * Calling variables:
 *      input_tangent   - double   - input tangent value
 *      *sum_tangent    - double * - output sum-of-tangent values
 *      *sum_extra      - double * - output sum of PI or PI/2 factors
 *
 */

void add_tangents(double input_tangent, double *sum_tangent, double *sum_extra)
{
    double q1, q2 ;
    
    q1 = *sum_tangent + input_tangent;
    q2 = 1.0 - *sum_tangent * input_tangent;
    
    if (q2 == 0.0) {
        if (q1 < 0.) {
            *sum_extra -= (M_PI/2.);
        } else {
            *sum_extra += (M_PI/2.);
        }
        *sum_tangent = 0.;
    } else {
        *sum_tangent = q1 / q2;
        if (q2 > 0.) {
            return;
        }
        if (q1 < 0.) {
            *sum_extra -= M_PI;
        } else if (q1 > 0.) {
            *sum_extra += M_PI;
        }
    }
    
    return;
} /* end of add_tangents() */
/************************ add_tangents() ****************************/

/************************ references ********************************/
/*
 * Plouff, D. 1975a. “Derivation of Formulas and FORTRAN Programs to
 * Compute Gravity Anomalies of Prisms.” No. PB-243-526. National Technical
 * Information Service.
 * https://ntrl.ntis.gov/NTRL/dashboard/searchResults/titleDetail/PB243526.xhtml.
 * 
 * Plouff, D. 1975b. “Derivation of Formulas and FORTRAN Programs to Compute
 * Magnetic Anomalies of Prisms.” No. PB-243-525. National Technical
 * Information Service.
 * https://ntrl.ntis.gov/NTRL/dashboard/searchResults/titleDetail/PB243525.xhtml.
 *
 */ 
/************************ references ********************************/
