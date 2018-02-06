/*
	ENT  --  Entropy calculation and analysis of putative
		 random sequences.

        Designed and implemented by John "Random" Walker in May 1985.

	Multiple analyses of random sequences added in December 1985.

	Bit stream analysis added in September 1997.

	Terse mode output, getopt() command line processing,
	optional stdin input, and HTML documentation added in
	October 1998.
	
	Documentation for the -t (terse output) option added
	in July 2006.
	
	Replaced table look-up for chi square to probability
	conversion with algorithmic computation in January 2008.

	For additional information and the latest version,
	see http://www.fourmilab.ch/random/

*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#else
#include <unistd.h>
#endif

#define isISOspace(x)	((isascii(((unsigned char) (x))) && isspace(((unsigned char) (x)))) || ((x) == 0xA0))
#define isISOalpha(x)	((isoalpha[(((unsigned char) (x))) / 8] & (0x80 >> ((((unsigned char) (x))) % 8))) != 0)
#define isISOupper(x)	((isoupper[(((unsigned char) (x))) / 8] & (0x80 >> ((((unsigned char) (x))) % 8))) != 0)
#define isISOlower(x)	((isolower[(((unsigned char) (x))) / 8] & (0x80 >> ((((unsigned char) (x))) % 8))) != 0)
#define isISOprint(x)   ((((x) >= ' ') && ((x) <= '~')) || ((x) >= 0xA0))
#define toISOupper(x)   (isISOlower(x) ? (isascii(((unsigned char) (x))) ?  \
                            toupper(x) : (((((unsigned char) (x)) != 0xDF) && \
                            (((unsigned char) (x)) != 0xFF)) ? \
			    (((unsigned char) (x)) - 0x20) : (x))) : (x))
#define toISOlower(x)   (isISOupper(x) ? (isascii(((unsigned char) (x))) ?  \
                            tolower(x) : (((unsigned char) (x)) + 0x20)) \
			    : (x))

#define UPDATE  "January 28th, 2008"

#define FALSE 0
#define TRUE  1

#ifdef M_PI
#define PI	 M_PI
#else
#define PI	 3.14159265358979323846
#endif

#define	Z_MAX          6.0            /* maximum meaningful z value */

const unsigned char isoalpha[32] = {
    0,0,0,0,0,0,0,0,127,255,255,224,127,255,255,224,0,0,0,0,0,0,0,0,255,255,
    254,255,255,255,254,255
};

const unsigned char isoupper[32] = {
    0,0,0,0,0,0,0,0,127,255,255,224,0,0,0,0,0,0,0,0,0,0,0,0,255,255,254,254,
    0,0,0,0
};

const unsigned char isolower[32] = {
    0,0,0,0,0,0,0,0,0,0,0,0,127,255,255,224,0,0,0,0,0,0,0,0,0,0,0,1,255,255,
    254,255
};

#define FALSE 0
#define TRUE  1

#define log2of10 3.32192809488736234787

static int binary = FALSE;	   /* Treat input as a bitstream */

static long ccount[256],	   /* Bins to count occurrences of values */
	    totalc = 0; 	   /* Total bytes counted */
static double prob[256];	   /* Probabilities per bin for entropy */

/*  RT_LOG2  --  Calculate log to the base 2  */

static double rt_log2(double x)
{
    return log2of10 * log10(x);
}

#define MONTEN	6		      /* Bytes used as Monte Carlo
					 co-ordinates.	This should be no more
					 bits than the mantissa of your
                                         "double" floating point type. */

static int mp, sccfirst;
static unsigned int monte[MONTEN];
static long inmont, mcount;
static double ccexp, incirc, montex, montey, montepi,
	      scc, sccun, sccu0, scclast, scct1, scct2, scct3,
	      ent, chisq, datasum;

/*  RT_INIT  --  Initialise random test counters.  */

static void rt_init(int binmode)
{
    int i;

    binary = binmode;	       /* Set binary / byte mode */

    /* Initialise for calculations */

    ent = 0.0;		       /* Clear entropy accumulator */
    chisq = 0.0;	       /* Clear Chi-Square */
    datasum = 0.0;	       /* Clear sum of bytes for arithmetic mean */

    mp = 0;		       /* Reset Monte Carlo accumulator pointer */
    mcount = 0; 	       /* Clear Monte Carlo tries */
    inmont = 0; 	       /* Clear Monte Carlo inside count */
    incirc = 65535.0 * 65535.0;/* In-circle distance for Monte Carlo */

    sccfirst = TRUE;	       /* Mark first time for serial correlation */
    scct1 = scct2 = scct3 = 0.0; /* Clear serial correlation terms */

    incirc = pow(pow(256.0, (double) (MONTEN / 2)) - 1, 2.0);

    for (i = 0; i < 256; i++) {
	ccount[i] = 0;
    }
    totalc = 0;
}

/*  RT_ADD  --	Add one or more bytes to accumulation.	*/

static void rt_add(void *buf, int bufl)
{
    unsigned char *bp = buf;
    int oc, c, bean;

    while (bean = 0, (bufl-- > 0)) {
       oc = *bp++;

       do {
	  if (binary) {
	     c = !!(oc & 0x80);
	  } else {
	     c = oc;
	  }
	  ccount[c]++;		  /* Update counter for this bin */
	  totalc++;

	  /* Update inside / outside circle counts for Monte Carlo
	     computation of PI */

	  if (bean == 0) {
	      monte[mp++] = oc;       /* Save character for Monte Carlo */
	      if (mp >= MONTEN) {     /* Calculate every MONTEN character */
		 int mj;

		 mp = 0;
		 mcount++;
		 montex = montey = 0;
		 for (mj = 0; mj < MONTEN / 2; mj++) {
		    montex = (montex * 256.0) + monte[mj];
		    montey = (montey * 256.0) + monte[(MONTEN / 2) + mj];
		 }
		 if ((montex * montex + montey *  montey) <= incirc) {
		    inmont++;
		 }
	      }
	  }

	  /* Update calculation of serial correlation coefficient */

	  sccun = c;
	  if (sccfirst) {
	     sccfirst = FALSE;
	     scclast = 0;
	     sccu0 = sccun;
	  } else {
	     scct1 = scct1 + scclast * sccun;
	  }
	  scct2 = scct2 + sccun;
	  scct3 = scct3 + (sccun * sccun);
	  scclast = sccun;
	  oc <<= 1;
       } while (binary && (++bean < 8));
    }
}

/*  RT_END  --	Complete calculation and return results.  */

static void rt_end(double *r_ent, double *r_chisq, double *r_mean,
	    double *r_montepicalc, double *r_scc)
{
    int i;

    /* Complete calculation of serial correlation coefficient */

    scct1 = scct1 + scclast * sccu0;
    scct2 = scct2 * scct2;
    scc = totalc * scct3 - scct2;
    if (scc == 0.0) {
       scc = -100000;
    } else {
       scc = (totalc * scct1 - scct2) / scc;
    }

    /* Scan bins and calculate probability for each bin and
       Chi-Square distribution.  The probability will be reused
       in the entropy calculation below.  While we're at it,
       we sum of all the data which will be used to compute the
       mean. */
       
    ccexp = totalc / (binary ? 2.0 : 256.0);  /* Expected count per bin */
    for (i = 0; i < (binary ? 2 : 256); i++) {
       double a = ccount[i] - ccexp;;
       
       prob[i] = ((double) ccount[i]) / totalc;       
       chisq += (a * a) / ccexp;
       datasum += ((double) i) * ccount[i];
    }

    /* Calculate entropy */

    for (i = 0; i < (binary ? 2 : 256); i++) {
       if (prob[i] > 0.0) {
	  ent += prob[i] * rt_log2(1 / prob[i]);
       }
    }

    /* Calculate Monte Carlo value for PI from percentage of hits
       within the circle */

    montepi = 4.0 * (((double) inmont) / mcount);

    /* Return results through arguments */

    *r_ent = ent;
    *r_chisq = chisq;
    *r_mean = datasum / totalc;
    *r_montepicalc = montepi;
    *r_scc = scc;
}

/*FUNCTION poz: probability of normal z value */
/*ALGORITHM
	Adapted from a polynomial approximation in:
		Ibbetson D, Algorithm 209
		Collected Algorithms of the CACM 1963 p. 616
	Note:
		This routine has six digit accuracy, so it is only useful for absolute
		z values < 6.  For z values >= to 6.0, poz() returns 0.0.
*/
static double        /*VAR returns cumulative probability from -oo to z */
poz(const double z)  /*VAR normal z value */
{
    double y, x, w;

    if (z == 0.0) {
    	x = 0.0;
    } else {
	y = 0.5 * fabs(z);
	if (y >= (Z_MAX * 0.5)) {
    	    x = 1.0;
	} else if (y < 1.0) {
	   w = y * y;
	   x = ((((((((0.000124818987 * w
		   -0.001075204047) * w +0.005198775019) * w
		   -0.019198292004) * w +0.059054035642) * w
		   -0.151968751364) * w +0.319152932694) * w
		   -0.531923007300) * w +0.797884560593) * y * 2.0;
	} else {
	    y -= 2.0;
	    x = (((((((((((((-0.000045255659 * y
		    +0.000152529290) * y -0.000019538132) * y
		    -0.000676904986) * y +0.001390604284) * y
		    -0.000794620820) * y -0.002034254874) * y
		    +0.006549791214) * y -0.010557625006) * y
		    +0.011630447319) * y -0.009279453341) * y
		    +0.005353579108) * y -0.002141268741) * y
		    +0.000535310849) * y +0.999936657524;
    	}
    }
    return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));
}

/*
	Module:       chisq.c
	Purpose:      compute approximations to chisquare distribution probabilities
	Contents:     pochisq()
	Uses:         poz() in z.c (Algorithm 209)
	Programmer:   Gary Perlman
	Organization: Wang Institute, Tyngsboro, MA 01879
	Copyright:    none
	Tabstops:     4
*/

#define	LOG_SQRT_PI     0.5723649429247000870717135 /* log (sqrt (pi)) */
#define	I_SQRT_PI       0.5641895835477562869480795 /* 1 / sqrt (pi) */
#define	BIGX           20.0         /* max value to represent exp (x) */
#define	ex(x)             (((x) < -BIGX) ? 0.0 : exp(x))

/*FUNCTION pochisq: probability of chi sqaure value */
/*ALGORITHM Compute probability of chi square value.
	Adapted from:
		Hill, I. D. and Pike, M. C.  Algorithm 299
		Collected Algorithms for the CACM 1967 p. 243
	Updated for rounding errors based on remark in
		ACM TOMS June 1985, page 185
*/

static double pochisq(
    	const double ax,    /* obtained chi-square value */
     	const int df	    /* degrees of freedom */
     	)
{
    double x = ax;
    double a, y, s;
    double e, c, z;
    int even;	    	    /* true if df is an even number */

    if (x <= 0.0 || df < 1) {
    	return 1.0;
    }

    a = 0.5 * x;
    even = (2 * (df / 2)) == df;
    if (df > 1) {
    	y = ex(-a);
    }
    s = (even ? y : (2.0 * poz(-sqrt(x))));
    if (df > 2) {
	x = 0.5 * (df - 1.0);
	z = (even ? 1.0 : 0.5);
	if (a > BIGX) {
    	    e = (even ? 0.0 : LOG_SQRT_PI);
    	    c = log(a);
    	    while (z <= x) {
		e = log(z) + e;
		s += ex(c * z - a - e);
		z += 1.0;
    	    }
    	    return (s);
    	} else {
	    e = (even ? 1.0 : (I_SQRT_PI / sqrt(a)));
	    c = 0.0;
	    while (z <= x) {
		    e = e * (a / z);
		    c = c + e;
		    z += 1.0;
    	    }
	    return (c * y + s);
    	}
    } else {
    	return s;
    }
}

/*  GETOPT  --	Dumb version of getopt for brain-dead Windows.  */

#ifdef _WIN32	
static int optind = 1;

static int getopt(int argc, char *argv[], char *opts)
{
    static char *opp = NULL;
    int o;
    
    while (opp == NULL) {
        if ((optind >= argc) || (*argv[optind] != '-')) {
	   return -1;
	}
	opp = argv[optind] + 1;
	optind++;
	if (*opp == 0) {
	    opp = NULL;
	}	
    }
    o = *opp++;
    if (*opp == 0) { 
	opp = NULL;
    }
    return strchr(opts, o) == NULL ? '?' : o;
}
#endif

/*  Main program  */

int entmain(const char* arg)
{
	int i, oc, opt;
	long ccount[256];	      /* Bins to count occurrences of values */
	long totalc = 0;	      /* Total character count */
	char *samp;
	double montepi, chip,
	       scc, ent, mean, chisq;
	FILE *fp = stdin;
	int counts = FALSE,	      /* Print character counts */
	    fold = FALSE,	      /* Fold upper to lower */
	    binary = FALSE,	      /* Treat input as a bitstream */
	    terse = FALSE;	      /* Terse (CSV format) output */

#if 0
        while ((opt = getopt(argc, argv, "bcftu?BCFTU")) != -1) {
	    switch (toISOlower(opt)) {
                 case 'b':
		    binary = TRUE;
		    break;

                 case 'c':
		    counts = TRUE;
		    break;

                 case 'f':
		    fold = TRUE;
		    break;

                 case 't':
		    terse = TRUE;
		    break;

                 case '?':
                 case 'u':
		    return 0;
	    }
	}
	if (optind < argc) {
	   if (optind != (argc - 1)) {
              printf("Duplicate file name.\n");
	      return 2;
	   }
#endif
           if ((fp = fopen(arg, "rb")) == NULL) {
              printf("Cannot open file %s\n", arg);
	      return 2;
	   }
//	}
	
#ifdef _WIN32

	    /** Warning!  On systems which distinguish text mode and
		binary I/O (MS-DOS, Macintosh, etc.) the modes in the open
        	statement for "fp" should have forced the input file into
        	binary mode.  But what if we're reading from standard
		input?  Well, then we need to do a system-specific tweak
        	to make sure it's in binary mode.  While we're at it,
        	let's set the mode to binary regardless of however fopen
		set it.

		The following code, conditional on _WIN32, sets binary
		mode using the method prescribed by Microsoft Visual C 7.0
        	("Monkey C"); this may require modification if you're
		using a different compiler or release of Monkey C.	If
        	you're porting this code to a different system which
        	distinguishes text and binary files, you'll need to add
		the equivalent call for that system. */

	    _setmode(_fileno(fp), _O_BINARY);
#endif

        samp = binary ? "bit" : "byte";
	memset(ccount, 0, sizeof ccount);

	/* Initialise for calculations */

	rt_init(binary);

	/* Scan input file and count character occurrences */

	unsigned char ocb;
	int b;
	unsigned char ob = ocb;
	while ((oc = fgetc(fp)) != EOF) {

	   if (fold && isISOalpha(oc) && isISOupper(oc)) {
	      oc = toISOlower(oc);
	   }
	   ocb = (unsigned char) oc;
	   totalc += binary ? 8 : 1;
	   if (binary) {

	    for (b = 0; b < 8; b++) {
		ccount[ob & 1]++;
		ob >>= 1;
	    }
	   } else {
	       ccount[ocb]++;	      /* Update counter for this bin */
	   }
	   rt_add(&ocb, 1);
	}
	fclose(fp);

	/* Complete calculation and return sequence metrics */

	rt_end(&ent, &chisq, &mean, &montepi, &scc);

	if (terse) {
           printf("0,File-%ss,Entropy,Chi-square,Mean,Monte-Carlo-Pi,Serial-Correlation\n",
              binary ? "bit" : "byte");
           printf("1,%ld,%f,%f,%f,%f,%f\n",
	      totalc, ent, chisq, mean, montepi, scc);
	}

	/* Calculate probability of observed distribution occurring from
	   the results of the Chi-Square test */

    	chip = pochisq(chisq, (binary ? 1 : 255));

	/* Print bin counts if requested */

	if (counts) {
	   if (terse) {
              printf("2,Value,Occurrences,Fraction\n");
	   } else {
              printf("Value Char Occurrences Fraction\n");
	   }
	   for (i = 0; i < (binary ? 2 : 256); i++) {
	      if (terse) {
                 printf("3,%d,%ld,%f\n", i,
		    ccount[i], ((double) ccount[i] / totalc));
	      } else {
		 if (ccount[i] > 0) {
                    printf("%3d   %c   %10ld   %f\n", i,
		       /* The following expression shows ISO 8859-1
			  Latin1 characters and blanks out other codes.
			  The test for ISO space replaces the ISO
			  non-blanking space (0xA0) with a regular
                          ASCII space, guaranteeing it's rendered
                          properly even when the font doesn't contain
			  that character, which is the case with many
			  X fonts. */
                       (!isISOprint(i) || isISOspace(i)) ? ' ' : i,
		       ccount[i], ((double) ccount[i] / totalc));
		 }
	      }
	   }
	   if (!terse) {
              printf("\nTotal:    %10ld   %f\n\n", totalc, 1.0);
	   }
	}

	/* Print calculated results */

	if (!terse) {
           printf("\"%.9f\", ", ent);

#if 0
           printf("Entropy = %f bits per %s.\n", ent, samp);
           printf("\nOptimum compression would reduce the size\n");
           printf("of this %ld %s file by %d percent.\n\n", totalc, samp,
		    (short) ((100 * ((binary ? 1 : 8) - ent) /
			      (binary ? 1.0 : 8.0))));
	   printf(
              "Chi square distribution for %ld samples is %1.2f, and randomly\n",
	      totalc, chisq);
	   if (chip < 0.0001) {
              printf("would exceed this value less than 0.01 percent of the times.\n\n");
	   } else if (chip > 0.9999) {
              printf("would exceed this value more than than 99.99 percent of the times.\n\n");
	   } else {
              printf("would exceed this value %1.2f percent of the times.\n\n",
		 chip * 100);
	   }
	   printf(
              "Arithmetic mean value of data %ss is %1.4f (%.1f = random).\n",
	      samp, mean, binary ? 0.5 : 127.5);
           printf("Monte Carlo value for Pi is %1.9f (error %1.2f percent).\n",
	      montepi, 100.0 * (fabs(PI - montepi) / PI));
           printf("Serial correlation coefficient is ");
	   if (scc >= -99999) {
              printf("%1.6f (totally uncorrelated = 0.0).\n", scc);
	   } else {
              printf("undefined (all values equal!).\n");
	   }
#endif
	}

	return 0;
}
