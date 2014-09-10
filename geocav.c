#include <stdio.h>
#include <math.h>
#include "linus.h"

#include "triangles.h"

#define RADIANS_TO_DEGREES 57.29577951308232e0

/* static double cmatrx[300][300]; */

static long ije[10000];

static double xc1[] = {
 0.35468432, 0.10960348,-0.28694564,-0.28694564, 0.10960348, 0.67765737,
 0.20940764,-0.54823631,-0.54823631, 0.20940764, 0.56409198,-0.43863285,
-0.83518195,-0.07753799, 0.78726083, 0.78726083,-0.07753799,-0.83518195,
-0.43863285, 0.56409198, 0.49635327,-0.65783978,-0.9029206 , 0.09980416,
 0.96460301, 0.96460301, 0.09980416,-0.9029206 ,-0.65783978, 0.49635327,
 0.9029206 ,-0.09980416,-0.96460301,-0.49635327, 0.65783978, 0.65783978,
-0.49635327,-0.96460301,-0.09980416, 0.9029206 , 0.83518195, 0.07753799,
-0.78726083,-0.56409198, 0.43863285, 0.43863285,-0.56409198,-0.78726083,
 0.07753799, 0.83518195, 0.54823631,-0.20940764,-0.67765737,-0.20940764,
 0.54823631, 0.28694564,-0.10960348,-0.35468432,-0.10960348, 0.28694564};

static double yc1[] = {
 0.00000000, 0.33732483, 0.20847821,-0.20847821,-0.33732483, 0.00000000,
 0.64449042, 0.39831701,-0.39831701,-0.64449042, 0.64449042, 0.73564184,
-0.1898388 ,-0.85296863,-0.33732483, 0.33732483, 0.85296863, 0.1898388 ,
-0.73564184,-0.64449042, 0.85296863, 0.73564184,-0.39831701,-0.98181528,
-0.20847821, 0.20847821, 0.98181528, 0.39831701,-0.73564184,-0.85296863,
 0.39831701, 0.98181528, 0.20847821,-0.85296863,-0.73564184, 0.73564184,
 0.85296863,-0.20847821,-0.98181528,-0.39831701, 0.18983880, 0.85296863,
 0.33732483,-0.64449042,-0.73564184, 0.73564184, 0.64449042,-0.33732483,
-0.85296863,-0.1898388 , 0.39831701, 0.64449042, 0.00000000,-0.64449042,
-0.39831701, 0.20847821, 0.33732483, 0.00000000,-0.33732483,-0.20847821};

static double zc1[] = {
 0.93498611, 0.93498611, 0.93498611, 0.93498611, 0.93498611, 0.73537779,
 0.73537779, 0.73537779, 0.73537779, 0.73537779, 0.51617086, 0.51617086,
 0.51617086, 0.51617086, 0.51617086, 0.51617086, 0.51617086, 0.51617086,
 0.51617086, 0.51617086, 0.16148652, 0.16148652, 0.16148652, 0.16148652,
 0.16148652, 0.16148652, 0.16148652, 0.16148652, 0.16148652, 0.16148652,
-0.16148652,-0.16148652,-0.16148652,-0.16148652,-0.16148652,-0.16148652,
-0.16148652,-0.16148652,-0.16148652,-0.16148652,-0.51617086,-0.51617086,
-0.51617086,-0.51617086,-0.51617086,-0.51617086,-0.51617086,-0.51617086,
-0.51617086,-0.51617086,-0.73537779,-0.73537779,-0.73537779,-0.73537779,
-0.73537779,-0.93498611,-0.93498611,-0.93498611,-0.93498611,-0.93498611};


// Returns accessible surface area

//atot = geocav(natoms, atoms, wmax);


double geocav(long ncor, Atom3d *atoms, int winmax)
{
	long i, iij, j, k, l1, n3, n4, nejci, ninf, nsup, ntrian, nts, un3;
	int ires;
	int more;
	double atp, ats, dd, dij2, fndiv, rei, rrej, sre2;
	double xei, xpl, xsl, xsm;
	double yei, ypl, ysl, ysm;
	double zei, zpl, zsl, zsm;
	double dx, dy, dz, dr;
	double atot=0.0e0;
	double pi=3.14159265358979323846e0;
	Atom3d *b, *a;

	ntrian = 1;
	fndiv = 4.0e0 * pi/60.0e0;
	for (i=0; i<ncor; i++) {
		a = &atoms[i];
		a->sasa = 0.0e0;
		ires = a->resnum;
		rei = a->sasarad;
		xei = a->x;
		yei = a->y;
		zei = a->z;
		ats = rei * rei * fndiv;
		iij = nejci = 0;
		for (j=i-1; j>-1; j--) {
			b = &atoms[j];
			if ((ires - b->resnum) > winmax) continue;
			dr = rei + b->sasarad;
			dx = fabs(xei - b->x);
			if (dx > dr) continue;
			dy = fabs(yei - b->y);
			if (dy > dr) continue;
			dz = fabs(zei - b->z);
			if (dz > dr) continue;
			dij2 = dx*dx + dy*dy + dz*dz;
			sre2 = dr * dr;
			if ( dij2 < sre2) {
				ije[iij] = j;
				iij = iij + 1;
				nejci = iij;
			}
		}
		for (j=i+1; j<ncor; j++) {
			b = &atoms[j];
			if ((b->resnum - ires) > winmax) continue;
			dr = rei + b->sasarad;
			dx = fabs(xei - b->x);
			if (dx > dr) continue;
			dy = fabs(yei - b->y);
			if (dy > dr) continue;
			dz = fabs(zei - b->z);
			if (dz > dr) continue;
			sre2 = dr * dr;
			dij2 = dx*dx + dy*dy + dz*dz;
			if ( dij2 < sre2) {
				ije[iij] = j;
				iij = iij + 1;
				nejci = iij;
			}
		}
                nsup = -1;
		un3 = 0;
		for (j=0; j<60; j++) {
			xpl = ypl = zpl = 0.0e0;
			nts = 0;
			ninf = nsup + 1;
			nsup = ninf + ntrian-1;
			for (k=ninf; k<=nsup; k++) {
				xsl = xc1[k]*rei;
				ysl = yc1[k]*rei;
				zsl = zc1[k]*rei;
				xsm = xsl + xei;
				ysm = ysl + yei;
				zsm = zsl + zei;
				l1 = un3;
				more = 1;
				for (n3=l1; n3<nejci; n3++){
					un3 = n3;
					n4 = ije[n3];
					b = &atoms[n4];
					dr = b->sasarad;
					dx = xsm - b->x;
					dy = ysm - b->y;
					dz = zsm - b->z;
					dd = dx*dx + dy*dy + dz*dz;
					rrej = dr * dr;
					if (dd < rrej) {
						more = 0;
						break;
					}
				}
				if (!more)
					continue;
				for (n3=0; n3<l1; n3++) {
					un3 = n3;
					n4 = ije[n3];
					b = &atoms[n4];
					dx = xsm - b->x;
					dy = ysm - b->y;
					dz = zsm - b->z;
					dd = dx*dx + dy*dy + dz*dz;
					dr = b->sasarad;
					rrej = dr * dr;
					if (dd < rrej) {
						more = 0;
						break;
					}
				}
				if (more) {
					xpl += xsl;
					ypl += ysl;
					zpl += zsl;
					nts++;
				}
			}
			if (nts == 0)
				continue;
			else {
				atp = ats * nts;
				atot += atp;
				a->sasa += atp;
			}
		}
	}
	return atot;
}

