#ifndef DEFS_H
#define DEFS_H

#define SQR(x) ( (x)*(x) )
#define BD_LAYER 2
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define SIGN(a) ( ((a) < 0.) ? -1. : 1. )
#define TINY_NUMBER __DBL_MIN__
//#define TINY_NUMBER 1.0e-20
#define HUGE_NUMBER 1.0e+20

enum {
nghost = 4
};

#endif
