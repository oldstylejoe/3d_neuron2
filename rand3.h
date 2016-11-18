/***************************************************************************
                          rand3.h  -  description
                             -------------------
    begin                : Tue Apr 18 2000
    copyright            : (C) 2000 by joe snider
    email                : 
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/*random number generator class
 * 
 *I haven't read it, but NR cites
 * Knuth, D.E. 1981, Seminumerical Algorithms, 2nd ed., vol2 of The Art of 
 *     Computer Programming (Reading, MA: Addison-Wessley), chpt 3.2-3.3
 * */

//random number generator yanked from Num Rec ch 7.1
//Uses a subtractive algorithm by Knuth

#ifndef RAND3
#define RAND3

//#include <stdlib.h>

#define MBIG 1000000000
#define MZ 0
#define FAC (1.0/MBIG)
#define MSEED 161803398
#define SEED_BIG 1000000

//According to Knuth, any large MBIG, and any smaller (but still large)
//MSEED can be substituted for the above values.

class CRand3{
 public:
   //nothing special about 7
   CRand3(unsigned long iseed) ;

 public:
   //returns a random deviate between 0 and 1
   double operator() () {
      if (++inext ==56) inext=1;
      if (++inextp == 56) inextp=1;
      long mj = ma[inext]-ma[inextp];
      if (mj < MZ) mj += MBIG;
      ma[inext]=mj;
      return mj*FAC;
   }

   //returns an integer in the range [0,n)
   //appropriate for use in std::random_shuffle
   int operator() (const int& inN) {
      return int(double(inN)*(*this)());
   }

   double GetRand() {return (*this)();}
   
 private:
   int inext,inextp;
   // the value 56 is special and should not be changed
   long ma[56];
};

//construction
CRand3::CRand3(unsigned long iseed = 7) {
   long mj,mk;
   int i,ii,k;
   
   iseed %= SEED_BIG;             //make sure the seed is large enough
   mj=MSEED - long(iseed);
   mj %= MBIG;
   ma[55]=mj;
   mk=1;
   for (i=1;i<=54;i++) {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
   }
   for (k=1;k<=4;k++) {
      for (i=1;i<=55;i++) {
	 ma[i] -= ma[1+(i+30) % 55];
	 if (ma[i] < MZ) ma[i] += MBIG;
      }
   }
   inext=0;
   inextp=31; //The constant 31 is special see Knuth   
}

#endif
