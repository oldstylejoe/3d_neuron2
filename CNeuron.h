//Joe Snider
//7/04
//
//This holds a CBranch object (a tree of neuron branches), and can read them in ...
//Provides interface for reading neurons in from a NeuroLucida file.
//Also able to test if a given box intersects any of the branches.
//
//9/05
//Changed to allow multiple branches per neuron.

#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <stdlib.h>
#include <boost/shared_ptr.hpp>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "point.h"
#include "Branch.h"
#include "rand3.h"
#include "3dvector.h"
#include "histogram.h"
#include "tnt.h"
#include "jama_eig.h"
#include "h_vector.h"

using namespace std;
using namespace TNT;
using namespace JAMA;
using namespace boost;

#ifndef NEURON
#define NEURON

//#define ADD_ERROR_ASC_HACK
//const double ERROR_ASC_HACK = 10.;

//#define TRIM_TO_FIRST_DAUGHTER
//#define TWO_DIMENSIONAL
//#define ROTATE_AROUND_Z
const int KURTOSIS_HISTOGRAM_BINS = 50; //number of bins in kurtosis histogram

const double PI = acos(-1.);

//sorting class
//Type T must be an style style container (vector works).
template <class T>
class CSort {
public:
   bool operator()(const T& inX, const T& inY) const {
      return lexicographical_compare(inX.begin(), inX.end(), inY.begin(), inY.end());
   }
};

//Evaluate the line integral in 3D.
//Used in CNeuron::ProductMoment by the gsl integration routine.
//Takes 9 parameters: power, ax, bx, ay, by, az, bz, first term, second term.
//   first and second term are 00, 01, 02, ... to evaluate xx, xy, xz, ... respectively.
//Used to integrate along a line.
//Be sure to scale the answer by |b-a|, the length of the line.
double LineIntegratorCovariances(double x, void* params) {
   const double n = ((double*) params)[0];
   const int t1 = int( ((double*) params)[7] + 0.5 );
   const int t2 = int( ((double*) params)[8] + 0.5 );
   const double a1 = ((double*) params)[1+2*t1];
   const double b1 = ((double*) params)[2+2*t1];
   const double a2 = ((double*) params)[1+2*t2];
   const double b2 = ((double*) params)[2+2*t2];

   double dRet = pow( ((b1-a1)*x+a1)*((b2-a2)*x+a2), n);
   return dRet;
}

//Evaluate the line integral in 3D.
//Used in CNeuron::ProductMoment by the gsl integration routine.
//Takes 7 parameters: power, ax, bx, ay, by, az, bz.
//Used to integrate along a line.
//Be sure to scale the answer by |b-a|, the length of the line.
double LineFunction3D(double x, void* params) {
   double n = ((double*) params)[0];
   double ax = ((double*) params)[1];
   double bx = ((double*) params)[2];
   double ay = ((double*) params)[3];
   double by = ((double*) params)[4];
   double az = ((double*) params)[5];
   double bz = ((double*) params)[6];
   double dRet = pow( ((bx-ax)*x+ax)*((by-ay)*x+ay)*((bz-az)*x+az), n);
   return dRet;
}

//Evaluate the line integral in 2D.
//Used in CNeuron::ProductMoment by the gsl integration routine.
//Takes 5 parameters: power, ax, bx, ay, by.
//Used to integrate along a line.
//Be sure to scale the answer by |b-a|, the length of the line.
double LineFunction2D(double x, void* params) {
   double n = ((double*) params)[0];
   double ax = ((double*) params)[1];
   double bx = ((double*) params)[2];
   double ay = ((double*) params)[3];
   double by = ((double*) params)[4];
   double dRet = pow( ((bx-ax)*x+ax)*((by-ay)*x+ay), n);
   return dRet;
}

//Evaluate the line integral in 3D in the x direction.
//Used in CNeuron::ProductMoment by the gsl integration routine.
//Takes 7 parameters: power, ax, bx, ay, by, az, bz.
//Used to integrate along a line.
//Be sure to scale the answer by |b-a|, the length of the line.
double LineFunction3D_X(double x, void* params) {
   double n = ((double*) params)[0];
   double ax = ((double*) params)[1];
   double bx = ((double*) params)[2];
   double ay = ((double*) params)[3];
   double by = ((double*) params)[4];
   double az = ((double*) params)[5];
   double bz = ((double*) params)[6];
   double dRet = pow( ((bx-ax)*x+ax), n);
   return dRet;
}

//Evaluate the line integral in 2D in the x direction.
//Used in CNeuron::ProductMoment by the gsl integration routine.
//Takes 5 parameters: power, ax, bx, ay, by.
//Used to integrate along a line.
//Be sure to scale the answer by |b-a|, the length of the line.
double LineFunction2D_X(double x, void* params) {
   double n = ((double*) params)[0];
   double ax = ((double*) params)[1];
   double bx = ((double*) params)[2];
   double ay = ((double*) params)[3];
   double by = ((double*) params)[4];
   double dRet = pow( ((bx-ax)*x+ax), n);
   return dRet;
}

//Evaluate the line integral in 3D in the y direction.
//Used in CNeuron::ProductMoment by the gsl integration routine.
//Takes 7 parameters: power, ax, bx, ay, by, az, bz.
//Used to integrate along a line.
//Be sure to scale the answer by |b-a|, the length of the line.
double LineFunction3D_Y(double x, void* params) {
   double n = ((double*) params)[0];
   double ax = ((double*) params)[1];
   double bx = ((double*) params)[2];
   double ay = ((double*) params)[3];
   double by = ((double*) params)[4];
   double az = ((double*) params)[5];
   double bz = ((double*) params)[6];
   double dRet = pow( ((by-ay)*x+ay), n);
   return dRet;
}

//Evaluate the line integral in 2D in the y direction.
//Used in CNeuron::ProductMoment by the gsl integration routine.
//Takes 5 parameters: power, ax, bx, ay, by.
//Used to integrate along a line.
//Be sure to scale the answer by |b-a|, the length of the line.
double LineFunction2D_Y(double x, void* params) {
   double n = ((double*) params)[0];
   double ax = ((double*) params)[1];
   double bx = ((double*) params)[2];
   double ay = ((double*) params)[3];
   double by = ((double*) params)[4];
   double dRet = pow( ((by-ay)*x+ay), n);
   return dRet;
}

//Evaluate the line integral in 3D in the Z direction.
//Used in CNeuron::SingleMoment by the gsl integration routine.
//Takes 7 parameters: power, ax, bx, ay, by, az, bz.
//Used to integrate along a line.
//Be sure to scale the answer by |b-a|, the length of the line.
double LineFunction3D_Z(double x, void* params) {
   double n = ((double*) params)[0];
   double ax = ((double*) params)[1];
   double bx = ((double*) params)[2];
   double ay = ((double*) params)[3];
   double by = ((double*) params)[4];
   double az = ((double*) params)[5];
   double bz = ((double*) params)[6];
   double dRet = pow( ((bz-az)*x+az), n);
   return dRet;
}

//Evaluate the line integral in 2D.
//Used in CNeuron::SumMoment by the gsl integration routine.
//Takes 5 parameters: power, ax, bx, ay, by.
//Used to integrate along a line.
//Be sure to scale the answer by |b-a|, the length of the line.
double LineSum2D(double x, void* params) {
   double n = 0.5*((double*) params)[0];
   double ax = ((double*) params)[1];
   double bx = ((double*) params)[2];
   double ay = ((double*) params)[3];
   double by = ((double*) params)[4];
   double dRet = pow( ((bx-ax)*x+ax)*((bx-ax)*x+ax) + ((by-ay)*x+ay)*((by-ay)*x+ay), n);
   return dRet;
}

//Evaluate the line integral in 3D.
//Used in CNeuron::SumMoment by the gsl integration routine.
//Takes 7 parameters: power, ax, bx, ay, by, az, bz.
//Used to integrate along a line.
//Be sure to scale the answer by |b-a|, the length of the line.
double LineSum3D(double x, void* params) {
   double n = 0.5*((double*) params)[0];
   double ax = ((double*) params)[1];
   double bx = ((double*) params)[2];
   double ay = ((double*) params)[3];
   double by = ((double*) params)[4];
   double az = ((double*) params)[5];
   double bz = ((double*) params)[6];
   double dRet = pow( ((bx-ax)*x+ax)*((bx-ax)*x+ax) + 
      ((by-ay)*x+ay)*((by-ay)*x+ay) + 
      ((bz-az)*x+az)*((bz-az)*x+az), n);
   return dRet;
}

class CNeuron{
public:
   //typedefs...

public:
   //constructors
   //constructors
   CNeuron() : m_dWidthX(0.), m_dWidthY(0.), m_dWidthZ(0.),
      m_dStandardDeviationX(0.), m_dStandardDeviationY(0.), m_dStandardDeviationZ(0.),
      m_bFirstPass(true)
   {}
   ~CNeuron() {
      //delete the added branches
      Clear();
   }

public:
   //gets and sets
   //not allowed to set the branch externally
   const vector<shared_ptr<CBranch> >& GetBranch() const {return m_vbrNeuron;}

   //not allowed to set the widths and standard deviation.
   double GetWidthX() const {return m_dWidthX;}
   double GetWidthY() const {return m_dWidthY;}
   double GetWidthZ() const {return m_dWidthZ;}
   double GetStandardDeviationX() {return m_dStandardDeviationX;}
   double GetStandardDeviationY() {return m_dStandardDeviationY;}
   double GetStandardDeviationZ() {return m_dStandardDeviationZ;}

public:
   //interface

   //Scale the sub neurons.
   //Does not change the normalized neuron (normalization not affected by scale).
   void Scale(const double& inX, const double& inY, const double& inZ) {
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->Scale(inX, inY, inZ);
      }
   }

   //Return the maximum level over all branches.
   //Defaults to -1 if this is empty.
   int MaxLevel() const {
      int iLevel = -1;
      for(int i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->MaxLevel(iLevel);
      }
      return iLevel;
   }

   //Count the number of intersections of the input neuron with this neuron.
   //Shift and rotate the input neuron inRepeats times, count the intersections, and average.
   //Returns the average and standard error of the number of intersections.
   void CountIntersections( const CNeuron& inNeuron, double& inMean, double& inStdErr, int inRepeats = 100 ) const {
      CRand3 rand3( rand() );

      //get the segments
      vector<vector<double> > vecInSegments;
      vector<shared_ptr<CBranch> >::const_iterator iter = inNeuron.GetBranch().begin();
      vector<shared_ptr<CBranch> >::const_iterator iter_end = inNeuron.GetBranch().end();
      for(; iter != iter_end; ++iter) {
         (*iter)->AllSegments(vecInSegments);
      }
      vector<vector<double> > vecThisSegments;
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->AllSegments(vecThisSegments);
      }

      //Find the centers of mass to use as the origin.
      double dInCOMX, dInCOMY, dInCOMZ;
      vector<vector<double> > vecInPoints;
      TotalMean(vecInSegments, dInCOMX, dInCOMY, dInCOMZ, vecInPoints);
      double dThisCOMX, dThisCOMY, dThisCOMZ;
      vector<vector<double> > vecThisPoints;
      TotalMean(vecThisSegments, dThisCOMX, dThisCOMY, dThisCOMZ, vecThisPoints);

      //Find the min distance to the center of mass to use as shifts.
      double dWidth = 10.e100;
      for(unsigned n = 0; n < vecThisPoints.size(); ++n) {
         dWidth = min(dWidth, (vecThisPoints[n][0]-dThisCOMX)*(vecThisPoints[n][0]-dThisCOMX)+
            (vecThisPoints[n][1]-dThisCOMY)*(vecThisPoints[n][1]-dThisCOMY)+
            (vecThisPoints[n][2]-dThisCOMZ)*(vecThisPoints[n][2]-dThisCOMZ));
      }
      dWidth = sqrt(dWidth);

      //Shift the input segment to have the same center of mass as this.
      RotateShift(0., 0., 0., dThisCOMX-dInCOMX, dThisCOMY-dInCOMY, dThisCOMZ-dInCOMZ, vecInSegments);

      //Do the counting.
      double dX, dY, dZ;
      double dXOld = 0.;
      double dYOld = 0.;
      double dZOld = 0.;
      vector<int> vecIntersections;
      for(int i = 0; i < inRepeats; ++i) {
         //shift (in a box) and rotate
         dX = dWidth*rand3();
         dY = dWidth*rand3();
         dZ = dWidth*rand3();
#ifdef TWO_DIMENSIONAL
         RotateShift(0., 0., 2.*PI*rand3(), dX-dXOld, dY-dYOld, 0., vecInSegments);
#else
         RotateShift(2.*PI*rand3(), PI*rand3(), 2.*PI*rand3(), dX-dXOld, dY-dYOld, dZ-dZOld, vecInSegments);
#endif
         dXOld = dX;
         dYOld = dY;
         dZOld = dZ;
         int iIntersections = 0;
         for(unsigned j = 0; j < vecInSegments.size(); ++j) {
            for(unsigned k = 0; k < m_vbrNeuron.size(); ++k) {
               iIntersections += m_vbrNeuron[k]->CountSegmentIntersections(vecInSegments[j]);
            }
         }
         vecIntersections.push_back(iIntersections);
      }

      //Find the averages.
      inMean = 0.;
      inStdErr = 0.;
      vector<int>::const_iterator iterI = vecIntersections.begin();
      vector<int>::const_iterator iterI_end = vecIntersections.end();
      for(; iterI != iterI_end; ++iterI) {
         inMean += double(*iterI);
      }
      inMean /= double(vecIntersections.size());
      iterI = vecIntersections.begin();
      for(; iterI != iterI_end; ++iterI) {
         inStdErr += (double(*iterI) - inMean)*(double(*iterI) - inMean);
      }
      inStdErr = sqrt(inStdErr)/double(vecIntersections.size());

   }

   //Clear all the branches.
   //This should be the only way that CBranches are ever deleted.
   void Clear() {
      //for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
      //   delete m_vbrNeuron[i];
      //}
      m_vbrNeuron.clear();
      m_vbrNormalizedNeuron.clear();
   }

   //Build a test line.
   //Connects inX1, inY1, inZ1 to inX2, inY2, inZ2 with inSegments segments.
   void BuildLine(const double& inX1, const double& inY1, const double& inZ1,
      const double& inX2, const double& inY2, const double& inZ2, const int& inSegments) {
         //clear anything that's here.
         Clear();
         shared_ptr<CBranch> brInsert( new CBranch(NULL) );
         shared_ptr<CBranch> brInsertN( new CBranch(NULL) );

         //parameterize the line up the line (bad variable names: m is slope, b is offset)
         double m_x = inX2-inX1;
         double m_y = inY2-inY1;
         double m_z = inZ2-inZ1;
         double b_x = inX1;
         double b_y = inY1;
         double b_z = inZ1;
         for(int i = 0; i <= inSegments; ++i) {
            double t = double(i)/double(inSegments);
            brInsert->AddPoint( m_x*t + b_x, m_y*t + b_y, m_z*t + b_z);
            brInsertN->AddPoint( m_x*t + b_x, m_y*t + b_y, m_z*t + b_z);
         }
         m_vbrNeuron.push_back(brInsert);
         m_vbrNormalizedNeuron.push_back(brInsertN);
   }

   //Build a test circle
   //Radius is inR. Offset is inL1 x inL2 x inL3 with inRes points around the circumfrance; rotated by (Euler angles) inPsi, inTheta, and inPhi
   void BuildCircle(const double& inR, const double& inL1, const double& inL2, const double& inL3,
      const int& inRes, const double& inPhi, const double& inTheta, const double& inPsi) {
         //clear anything that's here.
         Clear();
         shared_ptr<CBranch> brInsert( new CBranch(NULL) );

         //set up the rotation matrix
         double a11 = cos(inPsi)*cos(inPhi)-cos(inTheta)*sin(inPhi)*sin(inPsi);
         double a12 = cos(inPsi)*sin(inPhi)+cos(inTheta)*cos(inPhi)*sin(inPsi);
         double a13 = sin(inPsi)*sin(inTheta);
         double a21 = -1.*sin(inPsi)*cos(inPhi)-cos(inTheta)*sin(inPhi)*cos(inPsi);
         double a22 = -1.*sin(inPsi)*sin(inPhi)+cos(inTheta)*cos(inPhi)*cos(inPsi);
         double a23 = cos(inPsi)*sin(inTheta);
         double a31 = sin(inTheta)*sin(inPhi);
         double a32 = -1.*sin(inTheta)*cos(inPhi);
         double a33 = cos(inTheta);
         double x, y, z, t;
         for(int i = 0; i < inRes; ++i) {
            t = double(i)/double(inRes)*2.*PI;
            x = inR*cos(t)+inL1;
            y = inR*sin(t)+inL2;
            z = inL3;
            brInsert->AddPoint( a11*x + a12*y + a13*z,
               a21*x + a22*y + a23*z,
               a31*x + a32*y + a33*z );
         }
         m_vbrNeuron.push_back(brInsert);
   }

   //Build a test cube
   //Size is inL1 x inL2 x inL3 with inRes points on each side; rotated by (Euler angles) inPsi, inTheta, and inPhi
   void BuildCube(const double& inL1, const double& inL2, const double& inL3,
      const int& inRes, const double& inPhi, const double& inTheta, const double& inPsi) {
         //clear anything that's here.
         Clear();

         //set up the rotation matrix
         double a11 = cos(inPsi)*cos(inPhi)-cos(inTheta)*sin(inPhi)*sin(inPsi);
         double a12 = cos(inPsi)*sin(inPhi)+cos(inTheta)*cos(inPhi)*sin(inPsi);
         double a13 = sin(inPsi)*sin(inTheta);
         double a21 = -1.*sin(inPsi)*cos(inPhi)-cos(inTheta)*sin(inPhi)*cos(inPsi);
         double a22 = -1.*sin(inPsi)*sin(inPhi)+cos(inTheta)*cos(inPhi)*cos(inPsi);
         double a23 = cos(inPsi)*sin(inTheta);
         double a31 = sin(inTheta)*sin(inPhi);
         double a32 = -1.*sin(inTheta)*cos(inPhi);
         double a33 = cos(inTheta);
         double x, y, z;
         for(int i = 0; i < inRes; ++i) {
            x = inL1*double(i)/double(inRes);
            for(int j = 0; j < inRes; ++j) {
               y = inL2*double(j)/double(inRes);
               shared_ptr<CBranch> brInsert( new CBranch(NULL) );
               for(int k = 0; k < inRes; ++k) {
                  z = inL3*double(k)/double(inRes);
                  brInsert->AddPoint( a11*x + a12*y + a13*z,
                     a21*x + a22*y + a23*z,
                     a31*x + a32*y + a33*z );
               }
               m_vbrNeuron.push_back(brInsert);
            }
         }
   }

   //Call the external program 'qhull' to find the volume.
   //Assumes that the program 'qhull' is somewhere in the path or the current directory.
   //Also requires a 'cat' command similar to unix to handle piping.
   double QhullVolume() {
      //get all the segments
      vector<vector<double> > vecSegments;
      for(unsigned i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->AllSegments(vecSegments);
      }

      //Save the current set of points to a file (repeats don't bother qhull).
      ofstream ofTempPointsQhull("temp_points_qhull.txt");
      vector<vector<double> >::const_iterator iter = vecSegments.begin();
      vector<vector<double> >::const_iterator iter_end = vecSegments.end();
#ifdef TWO_DIMENSIONAL
      ofTempPointsQhull << "2\n" << 2*vecSegments.size() << "\n";
      for(; iter != iter_end; ++iter) {
         ofTempPointsQhull << iter->at(0) << " " << iter->at(2) << "\n";
         ofTempPointsQhull << iter->at(1) << " " << iter->at(3) << "\n";
      }
#else
      ofTempPointsQhull << "3\n" << 2*vecSegments.size() << "\n";
      for(; iter != iter_end; ++iter) {
         ofTempPointsQhull << iter->at(0) << " " << iter->at(2) << " " << iter->at(4) << "\n";
         ofTempPointsQhull << iter->at(1) << " " << iter->at(3) << " " << iter->at(5) << "\n";
      }
#endif
      ofTempPointsQhull.close();

      //Make the system call to QHull.
      system("cat temp_points_qhull.txt | qhull FA > temp_qhull_volume.txt");

      //Parse the volume out of the file.
      ifstream ifVolume("temp_qhull_volume.txt");
      string strTemp;
      double dVolume = 0.;
      while(ifVolume >> strTemp) {
         if(strTemp.find("volume") != string::npos) {
            if(!(ifVolume >> dVolume)) {
               cout << "Warning in CNeuron::QHullVolume: unable to parse volume...returning 0 and continuing\n" <<flush;
            }
         }
      }
      ifVolume.close();

      return dVolume;
   }

   //Do a Monte Carlo type simulation of the minimum distances to the neuron at inN different radii
   void AverageMinimumDistance(ostream& inOut, const int& inN, const int& inNumTrials) const {
      CRand3 rand3( rand() );

      //get all the segments
      vector<vector<double> > vecSegments;
      for(unsigned i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->AllSegments(vecSegments);
      }

      //the main loop
      double dR, dT, dS, dCosS, dInsertValue;
      double dSlope = 5./double(inN);
      double dOffset = 5./double(inN);
      for(int i = 0; i < inN; ++i) {
         dR = dSlope*double(i) + dOffset;
         vector<double> vecData;
         //get the data
         for(int j = 0; j < inNumTrials; ++j) {
            dT = 2.*PI*rand3();
            dS = rand3();
            dS = 0.5*PI*dS*dS*((rand3()>.5)?1.:-1.);
            dCosS = cos(dS);

#ifdef TWO_DIMENSIONAL
            dInsertValue = m_vbrNormalizedNeuron[0]->MinimumDistance( dR*cos(dT), dR*sin(dT), 0. );
            for(unsigned k = 1; k < m_vbrNormalizedNeuron.size(); ++k) {
               dInsertValue = min( dInsertValue, 
                  m_vbrNormalizedNeuron[k]->MinimumDistance( dR*cos(dT), dR*sin(dT), 0. ) );
            }
#else
            dInsertValue = m_vbrNormalizedNeuron[0]->MinimumDistance( dR*dCosS*sin(dT), 
               dR*dCosS*cos(dT),
               dR*sin(dS) );
            for(int k = 1; k < m_vbrNormalizedNeuron.size(); ++k) {
               dInsertValue = min( dInsertValue, 
                  m_vbrNormalizedNeuron[k]->MinimumDistance( dR*dCosS*sin(dT), 
                  dR*dCosS*cos(dT),
                  dR*sin(dS) ) );
            }
#endif
            vecData.push_back( dInsertValue );
         }

         //dump the mean and std error
         double dMean = 0.;
         vector<double>::const_iterator iterD = vecData.begin();
         vector<double>::const_iterator iterD_end = vecData.end();
         for(; iterD != iterD_end; ++iterD) {
            dMean += *iterD;
         }
         dMean /= double(vecData.size());
         double dStdErr = 0.;
         iterD = vecData.begin();
         iterD_end = vecData.end();
         for(; iterD != iterD_end; ++iterD) {
            dStdErr += (*iterD - dMean)*(*iterD - dMean);
         }
         dStdErr = sqrt(dStdErr);
         dStdErr /= double(vecData.size());

         inOut << dR << " " << dMean << " " << dStdErr << "\n" << flush;
      }
   }

   //Find the box length distribution.
   //Dumps the distribution to the stream inOut.
   void LengthDistribution(ostream& inOut, const int samples = 1000) const {
      //Create an histogram
      //CHistogram<double> hData(0., 2.*max(MaxX(m_vbrNormalizedNeuron),
      //   max(MaxY(m_vbrNormalizedNeuron),
      //   MaxZ(m_vbrNormalizedNeuron))), 1000);
      CHistogram<double> hData(0., 1000., samples);

      for(unsigned i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->LengthDistribution(hData);
      }
      inOut << hData;
   }

   /*//Find the second and fourth moments of a 2-d sample in various directions.
   //Returns geometric mean of kurtosis at minima * kurtosis at maxima of the second moment.
   //Specialized for 2d and mainly just testing for the more general inertial histogram method.
   double RotatedSecondFourthMoments(ostream& inOut, const int& inNumTrials) {
      //get all the segments
      vector<vector<double> > vecSegments;
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->AllSegments(vecSegments);
      }

      //find the center of mass to use as the origin.
      double dCOMX, dCOMY, dCOMZ;
      vector<vector<double> > vecPoints;
      TotalMeanSplit(vecSegments, 100, dCOMX, dCOMY, dCOMZ, vecPoints);

      //Shift so that the center of mass is at the origin.
      ShiftRotatePoints(0., 0., 0., -1.*dCOMX, -1.*dCOMY, -1.*dCOMZ, vecPoints);

      //Do the rotations and find the moments (in the x-direction, WLOG).
      vector<double> vecSecondMoment;
      vector<double> vecFourthMoment;
      vector<double> vecXMax;
      vector<double> vecXMin;
      int iMin = 0;
      double dMin = 1.e10;
      int iMax = 0;
      double dMax = -1.;
      double dX, dY, dM, dTotalMass;
      double dRotateAngle = 2.*acos(-1.)/double(inNumTrials);
      for(int j = 0; j < inNumTrials; ++j) {
         ShiftRotatePoints(dRotateAngle, 0., 0., 0., 0., 0., vecPoints);

         double dXMin = 1.e10;
         double dXMax = -1.e10;
         double dSecondMoment = 0.;
         double dFourthMoment = 0.;
         dTotalMass = 0.;
         vector<vector<double> >::const_iterator iter = vecPoints.begin();
         vector<vector<double> >::const_iterator iter_end = vecPoints.end();
         for(; iter != iter_end; ++iter) {
            dX = iter->at(0);
            dM = iter->at(3);

            dTotalMass += dM;
            dSecondMoment += dX*dX*dM;
            dFourthMoment += dX*dX*dX*dX*dM;

            dXMin = min(dX, dXMin);
            dXMax = max(dX, dXMax);
         }

         //update max and min
         if(dSecondMoment < dMin) {
            iMin = j;
            dMin = dSecondMoment;
         }
         if(dSecondMoment > dMax) {
            iMax = j;
            dMax = dSecondMoment;
         }

         ////rescale to units of the x-range (only for sigma_over_range calcualtion)
         ////Warning: this hack gives the wrong kurtosis
         //dSecondMoment /= (dXMax-dXMin)*(dXMax-dXMin);
         //dFourthMoment /= (dXMax-dXMin)*(dXMax-dXMin)*(dXMax-dXMin)*(dXMax-dXMin);

         //dump and record
         vecSecondMoment.push_back(dSecondMoment/dTotalMass);
         vecFourthMoment.push_back(dFourthMoment/dTotalMass);
         vecXMax.push_back(dXMax);
         vecXMin.push_back(dXMin);
         //inOut << double(j+1)*dRotateAngle << " "
         //   << dSecondMoment/dTotalMass << " "
         //   << dFourthMoment/dTotalMass << "\n";

      }

      //Set the moments and widths with this value.
      m_dWidthX = vecXMax[iMax] - vecXMin[iMax];
      m_dWidthY = vecXMax[iMin] - vecXMin[iMin];
      m_dWidthZ = 0.;
      m_dStandardDeviationX = sqrt(vecSecondMoment[iMax]);
      m_dStandardDeviationY = sqrt(vecSecondMoment[iMin]);
      m_dStandardDeviationZ = 0.;

      //Turn the points into a circle and find the mass at various distances for the circle.
      //First align the max with the x-axis.
      //Note: the +1 part of iMax+1 returns the vecPoints to their start.
      ShiftRotatePoints(2.*double(iMax+1)*acos(-1.)/double(inNumTrials), 0., 0., 0., 0., 0., vecPoints);
      double dXScale = 2./(vecXMax[iMax]-vecXMin[iMax]);
      double dYScale = 2./(vecXMax[iMin]-vecXMin[iMin]);
      //double dXScale = 2./sqrt(vecSecondMoment[iMax]*dTotalMass);
      //double dYScale = 2./sqrt(vecSecondMoment[iMin]*dTotalMass);

      //Scale x by max radius and y by min radius to make a circle of unit radius.
      //and record the distances to a histogram.
      CHistogram<double> histoDistances(0., 2., 20);
      CHistogram<double> histoCount(0., 2., 20);
      histoDistances.zero();
      vector<vector<double> >::const_iterator iter = vecPoints.begin();
      vector<vector<double> >::const_iterator iter_end = vecPoints.end();
      for(; iter != iter_end; ++iter) {
         dX = iter->at(0)*dXScale;
         dY = iter->at(1)*dYScale;
         dM = iter->at(3);
         histoDistances.increment(sqrt(dX*dX+dY*dY), dM);
         histoCount.increment(sqrt(dX*dX+dY*dY), 1.);
         ////dump the rotated image for testing.
         //inOut << dX << " " << dY << "\n";
      }

      //normalize and dump the histogram
      CHistogram<double>::iterator iterH = histoDistances.bin_begin();
      CHistogram<double>::iterator iterHC = histoCount.bin_begin();
      CHistogram<double>::iterator iterH_end = histoDistances.bin_end();
      double dA, dB;
      double dNorm = 0.;
      for(; iterH != iterH_end; ++iterH, ++iterHC) {
         if( *iterHC > 0.5 ) {
            dA = histoDistances.GetBinLeftEdge(iterH);
            dB = histoDistances.GetBinRightEdge(iterH);
            *iterH /= (acos(-1.)*(dB*dB-dA*dA));
            //*iterH /= *iterHC;
            dNorm += *iterH;
         }
      }
      iterH = histoDistances.bin_begin();
      iterH_end = histoDistances.bin_end();
      for(; iterH != iterH_end; ++iterH) {
         inOut << histoDistances.GetBinCenter(iterH) << " "
            << *iterH / dNorm << "\n";
      }

      //return the kurtosis
      //return sqrt( vecFourthMoment[iMin]/vecSecondMoment[iMin]/vecSecondMoment[iMin] *
      //             vecFourthMoment[iMax]/vecSecondMoment[iMax]/vecSecondMoment[iMax] );
      return vecFourthMoment[iMin]/vecSecondMoment[iMin]/vecSecondMoment[iMin];
      //return vecFourthMoment[iMax]/vecSecondMoment[iMax]/vecSecondMoment[iMax];
      //return sqrt(vecSecondMoment[iMin]);
   }
*/
   //Calculate the probability density funciton of the rescaled neuron.
   void InertialHistogram(ostream& inOut) const {
      //get all the segments
      vector<vector<double> > vecSegments;
      for(unsigned i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->AllSegments(vecSegments);
      }

      //find the center of mass to use as the origin.
      double dCOMX, dCOMY, dCOMZ, x, y, z;
      vector<vector<double> > vecPoints;
      TotalMeanSplit(vecSegments, 100, dCOMX, dCOMY, dCOMZ, vecPoints);

      //Create the histogram
      CHistogram<double> histoMass(0., 5., 100);
      histoMass.zero();
      CHistogram<double> histoCount(0., 5., 100);
      histoCount.zero();

      ////Find the second moments in the inertial axis directions, and set up scale factors.
      vector<vector<double> >::const_iterator iter = vecPoints.begin();
      vector<vector<double> >::const_iterator iter_end = vecPoints.end();
      for(; iter != iter_end; ++iter) {
         x = iter->at(0);
         y = iter->at(1);
         z = iter->at(2);
         histoMass.increment(sqrt(x*x+y*y+z*z), iter->at(3));
         histoCount.increment(sqrt(x*x+y*y+z*z), 1.);
      }

      //build the histogram (and normalize).
      double a, b;
      double dNormalization = 0.;
      CHistogram<double>::iterator iterH = histoMass.bin_begin();
      CHistogram<double>::iterator iterHC = histoCount.bin_begin();
      CHistogram<double>::iterator iterH_end = histoMass.bin_end();
      for(; iterH != iterH_end; ++iterH, ++iterHC) {
         if(*iterHC > 0.5) {
            a = histoMass.GetBinLeftEdge(iterH);
            b = histoMass.GetBinRightEdge(iterH);
#ifdef TWO_DIMENSIONAL
            *iterH /= acos(-1.)*(b*b-a*a);
#else
            *iterH /= 4./3.*acos(-1.)*(b*b*b-a*a*a);
#endif
            //*iterH /= double(*iterHC);
            dNormalization += *iterH;
         }
      }

      //dump the histogram
      dNormalization = 1./dNormalization;
      iterH = histoMass.bin_begin();
      iterHC = histoCount.bin_begin();
      iterH_end = histoMass.bin_end();
      for(; iterH != iterH_end; ++iterH, ++iterHC) {
         inOut << histoMass.GetBinCenter(iterH) << " "
            << *iterH * dNormalization << "\n";
      }

   }

   //Find a 2d distribution of the lengths of the segments.
   //Returns the kurtosis (calculated in two dimensions).
   double Distribution2D(ostream& inOut, const int& inBins) const {
      //get all the segments
      vector<vector<double> > vecSegments;
      for(unsigned i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->AllSegments(vecSegments);
      }

      //Set the max size scale, and create histograms for storing the data.
      double dSizeScale = max( MaxX(m_vbrNormalizedNeuron)-MinX(m_vbrNormalizedNeuron), 
         MaxY(m_vbrNormalizedNeuron)-MinY(m_vbrNormalizedNeuron) );
      CHistogram<CHistogram<double> > histoMass(-1.*dSizeScale, dSizeScale, inBins);
      CHistogram<CHistogram<double> > histoCount(-1.*dSizeScale, dSizeScale, inBins);
      CHistogram<CHistogram<double> >::iterator iterHC = histoCount.bin_begin();
      CHistogram<CHistogram<double> >::iterator iterH = histoMass.bin_begin();
      CHistogram<CHistogram<double> >::iterator iterH_end = histoMass.bin_end();
      for(; iterH != iterH_end; ++iterH, ++iterHC) {
         iterH->Set(-1.*dSizeScale, dSizeScale, inBins);
         iterHC->Set(-1.*dSizeScale, dSizeScale, inBins);
      }

      //Break up the segments into steps of (at most) this length.
      double dStepSize = histoMass.GetBinWidth()/10.;  //plenty small

      double dTotalMass = 0.;
      C3DVector vStart, vEnd, vL, vR, vC, vStep;
      vector<vector<double> >::const_iterator iterS = vecSegments.begin();
      vector<vector<double> >::const_iterator iterS_end = vecSegments.end();
      for(; iterS != iterS_end; ++iterS) {
         vStart.SetTo(iterS->at(0), iterS->at(2), 0.);
         vEnd.SetTo(iterS->at(1), iterS->at(3), 0.);
         double dM = (vEnd-vStart).Magnitude();
         int iStepCount = int(dM/dStepSize+1.);
         vStep = vEnd-vStart;
         vStep *= 1./(iStepCount);
         for(int i = 1; i <= iStepCount; ++i) {
            vL = vStart + double(i-1)*vStep;
            vR = vStart + double(i)*vStep;
            vC = 0.5*(vL+vR);
            histoMass[vC.GetX()][vC.GetY()] += (vR-vL).Magnitude();
            histoCount[vC.GetX()][vC.GetY()] += 1.;
         }
         dTotalMass += dM;
      }

      //Normalize and dump to the output stream.
      double dX, dY, dZ;
      double dSecondMoment = 0.;
      double dFourthMoment = 0.;
      iterHC = histoCount.bin_begin();
      iterH = histoMass.bin_begin();
      iterH_end = histoMass.bin_end();
      for(; iterH != iterH_end; ++iterH, ++iterHC) {
         CHistogram<double>::iterator iterHC2 = iterHC->bin_begin();
         CHistogram<double>::iterator iterH2 = iterH->bin_begin();
         CHistogram<double>::iterator iterH2_end = iterH->bin_end();
         for(; iterH2 != iterH2_end; ++iterH2, ++iterHC2) {
            if( *iterHC2 > 0.5 ) {
               dX = histoMass.GetBinCenter(iterH);
               dY = iterH->GetBinCenter(iterH2);
               dZ = *iterH2 / *iterHC2;
               inOut << dX << " " << dY << " " << *iterH2 << "\n";
               dSecondMoment += pow(dX*dY, 2) * dZ;
               dFourthMoment += pow(dX*dY, 4) * dZ;
            }
         }
      }

      return dFourthMoment*dTotalMass/(dSecondMoment*dSecondMoment);

   }

   //Find a 3d distribution of the lengths of the segments.
   void Distribution3D(ostream& inOut, const int& inBins) const {
      //get all the segments
      vector<vector<double> > vecSegments;
      for(unsigned i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->AllSegments(vecSegments);
      }

      //Set the max size scale, and create histograms for storing the data.
      double dSizeScale = max(MaxZ(m_vbrNormalizedNeuron)-MinZ(m_vbrNormalizedNeuron),
         max( MaxX(m_vbrNormalizedNeuron)-MinX(m_vbrNormalizedNeuron), 
         MaxY(m_vbrNormalizedNeuron)-MinY(m_vbrNormalizedNeuron) ) );
      CHistogram<CHistogram<CHistogram<double> > > histoMass(-1.*dSizeScale, dSizeScale, inBins);
      CHistogram<CHistogram<CHistogram<double> > >::iterator iterHX = histoMass.bin_begin();
      CHistogram<CHistogram<CHistogram<double> > >::iterator iterHX_end = histoMass.bin_end();
      for(; iterHX != iterHX_end; ++iterHX) {
         iterHX->Set(-1.*dSizeScale, dSizeScale, inBins);
         CHistogram<CHistogram<double> >::iterator iterHY = iterHX->bin_begin();
         CHistogram<CHistogram<double> >::iterator iterHY_end = iterHX->bin_end();
         for(; iterHY != iterHY_end; ++iterHY) {
            iterHY->Set(-1.*dSizeScale, dSizeScale, inBins);
         }
      }

      //Break up the segments into steps of (at most) this length.
      double dStepSize = histoMass.GetBinWidth()/10.;  //plenty small

      C3DVector vStart, vEnd, vL, vR, vC, vStep;
      vector<vector<double> >::const_iterator iterS = vecSegments.begin();
      vector<vector<double> >::const_iterator iterS_end = vecSegments.end();
      for(; iterS != iterS_end; ++iterS) {
         vStart.SetTo(iterS->at(0), iterS->at(2), iterS->at(4));
         vEnd.SetTo(iterS->at(1), iterS->at(3), iterS->at(5));
         double dM = (vEnd-vStart).Magnitude();
         int iStepCount = int(dM/dStepSize+1.);
         vStep = vEnd-vStart;
         vStep *= 1./(iStepCount);
         for(int i = 1; i <= iStepCount; ++i) {
            vL = vStart + double(i-1)*vStep;
            vR = vStart + double(i)*vStep;
            vC = 0.5*(vL+vR);
            histoMass[vC.GetX()][vC.GetY()][vC.GetZ()] += (vR-vL).Magnitude();
         }
      }

      //Normalize and dump to the output stream.
      iterHX = histoMass.bin_begin();
      iterHX_end = histoMass.bin_end();
      for(; iterHX != iterHX_end; ++iterHX) {
         iterHX->Set(-1.*dSizeScale, dSizeScale, inBins);
         CHistogram<CHistogram<double> >::iterator iterHY = iterHX->bin_begin();
         CHistogram<CHistogram<double> >::iterator iterHY_end = iterHX->bin_end();
         for(; iterHY != iterHY_end; ++iterHY) {
            CHistogram<double>::iterator iterHZ = iterHY->bin_begin();
            CHistogram<double>::iterator iterHZ_end = iterHY->bin_end();
            for(; iterHZ != iterHZ_end; ++iterHZ) {
               inOut << histoMass.GetBinCenter(iterHX) << " "
                  << iterHX->GetBinCenter(iterHY) << " "
                  << iterHY->GetBinCenter(iterHZ) << " "
                  << *iterHZ << "\n";
            }
         }
      }
   }

   //Find the kurtosis of the neuron
   double Kurtosis(ostream& inOut) const {
      //get all the segments
      vector<vector<double> > vecSegments;
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->AllSegments(vecSegments);
      }

      //find the center of mass to use as the origin.
      double dCOMX, dCOMY, dCOMZ;
      vector<vector<double> > vecPoints;
      TotalMean(vecSegments, dCOMX, dCOMY, dCOMZ, vecPoints);

      //Set the max size scale, and create histograms for storing the data.
      double dSizeScale = sqrt( (MaxX(m_vbrNeuron)-MinX(m_vbrNeuron))*(MaxX(m_vbrNeuron)-MinX(m_vbrNeuron)) + 
         (MaxY(m_vbrNeuron)-MinY(m_vbrNeuron))*(MaxY(m_vbrNeuron)-MinY(m_vbrNeuron)) + 
         (MaxZ(m_vbrNeuron)-MinZ(m_vbrNeuron))*(MaxZ(m_vbrNeuron)-MinZ(m_vbrNeuron)) );
      CHistogram<double> histoMass(0., dSizeScale, KURTOSIS_HISTOGRAM_BINS);
      CHistogram<double> histoCount(0., dSizeScale, KURTOSIS_HISTOGRAM_BINS);
      histoMass.zero();
      histoCount.zero();

      //Find the mass m a distance r away from the COM.
      vector<vector<double> >::const_iterator iterS = vecPoints.begin();
      vector<vector<double> >::const_iterator iterS_end = vecPoints.end();
      for(; iterS != iterS_end; ++iterS) {
         double dDist = sqrt( (iterS->at(0)-dCOMX)*(iterS->at(0)-dCOMX) +
            (iterS->at(1)-dCOMY)*(iterS->at(1)-dCOMY) +
            (iterS->at(2)-dCOMZ)*(iterS->at(2)-dCOMZ) );
         histoMass[dDist] += iterS->at(3);//*iterS->at(4);
         histoCount[dDist] += 1.;
      }

      //Normalize, do the sums, and dump to the output stream.
      double dSecondMoment = 0.;
      double dFourthMoment = 0.;
      double dSum = 0.;
      double dTemp;
      CHistogram<double>::iterator iterHC = histoCount.bin_begin();
      CHistogram<double>::iterator iterH = histoMass.bin_begin();
      CHistogram<double>::iterator iterH_end = histoMass.bin_end();
      for(; iterH != iterH_end; ++iterH, ++iterHC) {
         if( *iterHC > 0.5 ) {
            //*iterH /= *iterHC;
            *iterH /= histoMass.GetBinCenter(iterH)*histoMass.GetBinWidth()*2.*acos(-1.);

            dSum += *iterH;
            dTemp = histoMass.GetBinCenter(iterH);
            dTemp *= dTemp;
            dSecondMoment += dTemp* (*iterH);
            dTemp *= dTemp;
            dFourthMoment += dTemp* (*iterH);

            inOut << histoMass.GetBinCenter(iterH) << " " << *iterH << "\n";
            //inOut << histoCount.GetBinCenter(iterHC) << " " << *iterHC << "\n";
         }
      }

      return dSum*dFourthMoment/(dSecondMoment*dSecondMoment);
   }

   //Find the smallest distance to the neuron
   double SmallestDistance(const double& inX, const double& inY, const double& inZ) const {
      double dReturn = m_vbrNeuron[0]->MinimumDistance(inX, inY, inZ);
      for(unsigned i = 1; i < m_vbrNeuron.size(); ++i) {
         dReturn = min(dReturn, m_vbrNeuron[i]->MinimumDistance(inX, inY, inZ));
      }
      return dReturn;
   }

   //Calculate and dump the angle distribution (raw, no averaging)
   void DumpAngles(ostream& inOut) const {
      vector<double> vecAngle;
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->RawAngles(inOut);
      }
   }

   //Calculate and dump the branch lengths (raw, no averaging)
   void DumpLengths(ostream& inOut) const {
      vector<double> vecLength;
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->RawLengths(inOut);
      }
   }

   //save the direction vector
   //Format of output is <r> <x> <y> <z> <n>
   //The error here is proportional to the length of the direction vector
   void DumpCOMAngle(ostream& inOut, const int samples=1000) const {
      //create a histogram
      CHistogram<CPoint> h(0, 1000, samples);
      //initialize the histogram points (under and over flows too)
      CHistogram<CPoint>::iterator iter = h.begin();
      CHistogram<CPoint>::iterator iter_end = h.end();
      for(; iter != iter_end; ++iter) {
         iter->push_back(0);
         iter->push_back(0);
         iter->push_back(0);
         iter->push_back(0);
      }
      for(unsigned i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->GlobalAngle(h);
      }
      
      //dump
      double x, y, z, n;
      iter = h.begin();
      iter_end = h.end();
      for(; iter != iter_end; ++iter) {
         //n = iter->at(3);
         //if(n > 0.5) {
             //only dump if there's some data
             inOut << h.GetBinCenter(iter) << " " << *iter << "\n";
         //}
      }
   }

   //Calculate and dump the number of branches at samples radii
   void DumpBranchCount(ostream& inOut, const int samples = 1000) const {
      //create a histogram
      CHistogram<double> h(0, 1000, samples);
      for(unsigned i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->BranchCount(h);
      }
      inOut << h;
   }

   //calculate and dump the average distance between branch points (line integral along the branch)
   // and the bounding radius.
   //Format: BL, BL_error, Radius
   void DisplayBL_versus_R(ostream& inOut) const {
      vector<double> vecBL;
      unsigned i;
      for(i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->BranchLength(vecBL);
      }

      double dMean = 0.;
      vector<double>::const_iterator iter = vecBL.begin();
      vector<double>::const_iterator iter_end = vecBL.end();
      for(; iter != iter_end; ++iter) {
         dMean += *iter;
      }
      dMean /= double(vecBL.size());
      double dStdErr = 0.;
      iter = vecBL.begin();
      iter_end = vecBL.end();
      for(; iter != iter_end; ++iter) {
         dStdErr += (*iter - dMean)*(*iter - dMean);
      }
      dStdErr = sqrt(dStdErr);
      dStdErr /= double(vecBL.size());

      inOut << dMean << " "
         << dStdErr << " ";

      //get all the segments
      vector<vector<double> > vecSegments;
      for(i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->AllSegments(vecSegments);
      }

      //find the average displacement squared
      double dCMX, dCMY, dCMZ;
      vector<vector<double> > vecPoints;
      TotalMean(vecSegments, dCMX, dCMY, dCMZ, vecPoints);

      vector<double> vecR2;
      for(i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->AverageDisplacementSquared(dCMX, dCMY, dCMZ, vecR2);
      }

      dMean = 0.;
      iter = vecR2.begin();
      iter_end = vecR2.end();
      for(; iter != iter_end; ++iter) {
         dMean += *iter;
      }
      dMean /= double(vecR2.size());
      dStdErr = 0.;
      iter = vecR2.begin();
      iter_end = vecR2.end();
      for(; iter != iter_end; ++iter) {
         dStdErr += (*iter - dMean)*(*iter - dMean);
      }
      dStdErr = sqrt(dStdErr);
      dStdErr /= double(vecR2.size());

      inOut << dMean << " "
         << dStdErr << "\n";

   }

   //Dump the product of the moments in all directions.
   //Dumps powers from inLow to inHigh (inclusive) with step inStep.
   //Assumes all powers are the same (just because I'm lazy).
   //Used for scaling.
   //Evaluates sum_{arbor}(x*y*z)^p*l(x,y,z)
   //   Note, no division by the total length; this has units.
   //Returns the last moment evaluated, rest go to inOut.
   double ProductMoment(ostream& inOut, 
      const int inLow, const int inHigh, const int inStep) const {
      //get all the segments
      vector<vector<double> > vecSegments;
      for(unsigned i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->AllSegments(vecSegments);
      }

      //Create some things needed for the gsl integrator.
      gsl_integration_workspace* w =
         gsl_integration_workspace_alloc(1000);  //1000 is plenty of room

      double dResult, dError;
      double dAlpha[8];

      //Easiest just to have both functions available, this takes no time.
      gsl_function LF_2D;
      LF_2D.function = &LineFunction2D;
      LF_2D.params = &dAlpha;
      gsl_function LF_3D;
      LF_3D.function = &LineFunction3D;
      LF_3D.params = &dAlpha;

      //Do the integration.
      double dSum;
      vector<vector<double> >::const_iterator iter;
      vector<vector<double> >::const_iterator iter_end;
      for(int i = inLow; i <= inHigh; i += inStep) {
         dSum = 0.;
         iter = vecSegments.begin();
         iter_end = vecSegments.end();
         for(; iter != iter_end; ++iter) {
            dAlpha[0] = (double)i;
            for(int j = 0; j < 7; ++j) {
               dAlpha[j+1] = iter->at(j);
            }
#ifdef TWO_DIMENSIONAL
            gsl_integration_qag(&LF_2D, 0, 1, 0, 1.e-7, 1000, GSL_INTEG_GAUSS61, 
               w, &dResult, &dError);
#else
            gsl_integration_qag(&LF_3D, 0, 1, 0, 1.e-7, 1000, GSL_INTEG_GAUSS61, 
               w, &dResult, &dError);
#endif
            //cout << "CNeuron.h gh1 " << i << " " << dSum << "\n" << flush;
            //cout << "CNeuron.h gh1.5 " << *iter << "\n" << flush;
            //cout << "CNeuron.h gh2 " << SegmentLength(*iter) << " " << dResult << " " << dError << "\n" << flush;
            dSum += SegmentLength(*iter)*dResult;
            //cout << "CNeuron.h gh5\n" << flush;
         }
         inOut << dSum << " " << flush;
      }
      inOut << "\n" << flush;

      //clean up the workspace
      gsl_integration_workspace_free(w);

      return dSum;
   }

   //Dump the product of the moments in one direction.
   //Dumps powers from inLow to inHigh (inclusive) with step inStep.
   //Used for scaling.
   //Evaluates sum_{arbor}x^p*l(x,y,z) (or y^p or z^p)
   //   Note, no division by the total length; this has units.
   //Returns the last moment evaluated, rest go to inOut.
   double SingleMoment(ostream& inOut, 
      const int inLow, const int inHigh, const int inStep, const int inDir) const {
      //get all the segments
      vector<vector<double> > vecSegments;
      for(unsigned i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->AllSegments(vecSegments);
      }

      //Create some things needed for the gsl integrator.
      gsl_integration_workspace* w =
         gsl_integration_workspace_alloc(1000);  //1000 is plenty of room

      double dResult, dError;
      double dAlpha[8];

      //Easiest just to have both functions available, this takes no time.
      gsl_function LF_2D;
      LF_2D.params = &dAlpha;
      gsl_function LF_3D;
      LF_3D.params = &dAlpha;
      if(inDir == 0) {
         LF_2D.function = &LineFunction2D_X;
         LF_3D.function = &LineFunction3D_X;
      } else if(inDir == 1) {
         LF_2D.function = &LineFunction2D_Y;
         LF_3D.function = &LineFunction3D_Y;
      } else if(inDir == 2) {
#ifdef TWO_DIMENSIONAL
         cerr << "Error in CNeuron::SingleMoment: Unable to integrate over z in two dimensions...returning\n" << flush;
         return -1.;
#endif
         LF_3D.function = &LineFunction3D_Z;
      } else {
         cerr << "Unknown direction in CNeuron::SingleMoment...returning...\n" << flush;
         return -1.;
      }

      //Do the integration.
      double dSum;
      vector<vector<double> >::const_iterator iter;
      vector<vector<double> >::const_iterator iter_end;
      for(int i = inLow; i <= inHigh; i += inStep) {
         dSum = 0.;
         iter = vecSegments.begin();
         iter_end = vecSegments.end();
         for(; iter != iter_end; ++iter) {
            dAlpha[0] = (double)i;
            for(int j = 0; j < 7; ++j) {
               dAlpha[j+1] = iter->at(j);
            }
#ifdef TWO_DIMENSIONAL
            gsl_integration_qag(&LF_2D, 0, 1, 0, 1.e-7, 1000, GSL_INTEG_GAUSS61, 
               w, &dResult, &dError);
#else
            gsl_integration_qag(&LF_3D, 0, 1, 0, 1.e-7, 1000, GSL_INTEG_GAUSS61, 
               w, &dResult, &dError);
#endif
            //cout << "CNeuron.h gh1 " << i << " " << dSum << "\n" << flush;
            //cout << "CNeuron.h gh1.5 " << *iter << "\n" << flush;
            //cout << "CNeuron.h gh2 " << SegmentLength(*iter) << " " << dResult << " " << dError << "\n" << flush;
            dSum += SegmentLength(*iter)*dResult;
            //cout << "CNeuron.h gh5\n" << flush;
         }
         inOut << dSum << " " << flush;
      }
      inOut << "\n" << flush;

      //clean up the workspace
      gsl_integration_workspace_free(w);

      return dSum;
   }

   //Dump the sum moments in all directions.
   //Dumps powers from inLow to inHigh (inclusive) with step inStep.
   //Assumes all powers are the same (just because I'm lazy).
   //Used for scaling.
   //Evaluates sum_{arbor}(x*x+y*y+z*x)^p/2*l(x,y,z)
   //   Note, no division by the total length; this has units.
   //Returns the last moment evaluated, rest go to inOut.
   double SumMoment(ostream& inOut, 
      const int inLow, const int inHigh, const int inStep) const {
      //get all the segments
      vector<vector<double> > vecSegments;
      for(unsigned i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->AllSegments(vecSegments);
      }

      //Create some things needed for the gsl integrator.
      gsl_integration_workspace* w =
         gsl_integration_workspace_alloc(1000);  //1000 is plenty of room

      double dResult, dError;
      double dAlpha[8];

      //Easiest just to have both functions available, this takes no time.
      gsl_function LF_2D;
      LF_2D.function = &LineSum2D;
      LF_2D.params = &dAlpha;
      gsl_function LF_3D;
      LF_3D.function = &LineSum3D;
      LF_3D.params = &dAlpha;

      //Do the integration.
      double dSum;
      vector<vector<double> >::const_iterator iter;
      vector<vector<double> >::const_iterator iter_end;
      for(int i = inLow; i <= inHigh; i += inStep) {
         dSum = 0.;
         iter = vecSegments.begin();
         iter_end = vecSegments.end();
         for(; iter != iter_end; ++iter) {
            dAlpha[0] = (double)i;
            for(int j = 0; j < 7; ++j) {
               dAlpha[j+1] = iter->at(j);
            }
#ifdef TWO_DIMENSIONAL
            gsl_integration_qag(&LF_2D, 0, 1, 0, 1.e-7, 1000, GSL_INTEG_GAUSS61, 
               w, &dResult, &dError);
#else
            gsl_integration_qag(&LF_3D, 0, 1, 0, 1.e-7, 1000, GSL_INTEG_GAUSS61, 
               w, &dResult, &dError);
#endif
            //cout << "CNeuron.h gh1 " << i << " " << dSum << "\n" << flush;
            //cout << "CNeuron.h gh1.5 " << *iter << "\n" << flush;
            //cout << "CNeuron.h gh2 " << SegmentLength(*iter) << " " << dResult << " " << dError << "\n" << flush;
            dSum += SegmentLength(*iter)*dResult;
            //cout << "CNeuron.h gh5\n" << flush;
         }
         inOut << dSum << " " << flush;
      }
      inOut << "\n" << flush;

      //clean up the workspace
      gsl_integration_workspace_free(w);

      return dSum;
   }

   //Use the sqrt of the determinant of the moment of inertia tensor as the volume.
   double InertialVolume() const {
      if(m_vbrNeuron.size() == 0) {
         return -1.;
      }

      //get all the segments
      vector<vector<double> > vecSegments;
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->AllSegments(vecSegments);
      }

      //create and initialize the moment of inertia tensor (just a vector<vector<double> >)
      vector<vector<double> > vecInertia;
      for(int i = 0; i < 3; ++i) {
         vector<double> vecInsert(3, 0.);
         vecInertia.push_back(vecInsert);
      }

      //find the center of mass to use as the origin (shift theorm may be faster, but I'd be more likely to screw it up)
      double dCOMX, dCOMY, dCOMZ;
      vector<vector<double> > vecPoints;
      TotalMean(vecSegments, dCOMX, dCOMY, dCOMZ, vecPoints);

      //do the sum to find the moment of inertia
      //note: some bad variable names (x, y, z are the position of the point, m its mass (length))
      double m, x, y, z;
      double dM = 0.;
      vector<vector<double> >::const_iterator iter = vecPoints.begin();
      vector<vector<double> >::const_iterator iter_end = vecPoints.end();
      for(; iter != iter_end; ++iter) {
         x = (*iter)[0] - dCOMX;
         y = (*iter)[1] - dCOMY;
         z = (*iter)[2] - dCOMZ;
         m = (*iter)[3];
         dM += m;
         vecInertia[0][0] += m*(y*y+z*z);
         vecInertia[0][1] -= m*x*y;
         vecInertia[0][2] -= m*x*z;
         vecInertia[1][1] += m*(x*x+z*z);
         vecInertia[1][0] -= m*y*x;
         vecInertia[1][2] -= m*y*z;
         vecInertia[2][2] += m*(x*x+y*y);
         vecInertia[2][0] -= m*z*x;
         vecInertia[2][1] -= m*z*y;
      }
#ifdef TWO_DIMENSIONAL
      return sqrt((vecInertia[0][0]*vecInertia[1][1]-vecInertia[0][1]*vecInertia[1][0])/(dM*dM));
#else
      return sqrt(Determinant(vecInertia)/(dM*dM*dM));
#endif
   }

   //Return the total length of the neuron.
   double TotalLength() const {
      double dLength = 0.;
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->TotalLength(dLength);
      }
      return dLength;
   }

   //Return the variance (related to the volume?) of the segments making up the neuron.
   double Variance() const {
      double dReturn = 0.;
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         cerr << "Warning: variance calculation not updated as of 9/05. Don't trust the result\n" << flush;
         m_vbrNeuron[i]->TotalVariance(dReturn);
      }
      return dReturn;
   }

   //Test the self similarity of the neuron by intersecting it with sub-volumes 
   //   inTrials times and measuring the length included in the sub-volume.
   //For now the test volumes are all rectangular prisms
   //  so I can use some of the box counting code.
   //Does a raw dump of: test volume, length in intersection.
   void CheckSelfSimilar(ostream& inOut, const int& inTrials) const {
      CRand3 rand3(rand());

      //A bounding box. Pad a little in the z-direction in 2D.
#ifdef TWO_DIMENSIONAL
      C3DBox boxBound(MinX(m_vbrNeuron), MinY(m_vbrNeuron), MinZ(m_vbrNeuron)-1.,
         MaxX(m_vbrNeuron), MaxY(m_vbrNeuron), MaxZ(m_vbrNeuron)+1.);
#else
      C3DBox boxBound(MinX(m_vbrNeuron), MinY(m_vbrNeuron), MinZ(m_vbrNeuron),
         MaxX(m_vbrNeuron), MaxY(m_vbrNeuron), MaxZ(m_vbrNeuron));
#endif

      //Get all the line segments
      vector<vector<double> > vecSegments;
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->AllSegments(vecSegments);
      }
      vector<vector<double> >::const_iterator iterS;
      vector<vector<double> >::const_iterator iterS_end = vecSegments.end();

      //The main loop
      for(int i = 0; i < inTrials; ++i) {
         C3DBox boxTest = boxBound.RandomSubBox(rand3);
         double dIntsersectionLengthSum = 0.;
         double dDumbIntsersectionLengthSum = 0.;
         //Rotate the neuron (Euler angles)
#ifdef TWO_DIMENSIONAL
         boxTest.SetAB(boxTest.GetA().GetX(), boxTest.GetA().GetY(), -1.,
            boxTest.GetB().GetX(), boxTest.GetB().GetY(), 1.);
         RotateShift(2.*PI*rand3(), 0., 0., 0., 0., 0., vecSegments);
#else
         RotateShift(2.*PI*rand3(), PI*rand3(), 2.*PI*rand3(), 0., 0., 0., vecSegments);
#endif
         for(iterS = vecSegments.begin(); iterS != iterS_end; ++iterS) {
            C3DVector vecL(iterS->at(0), iterS->at(2), iterS->at(4));
            C3DVector vecR(iterS->at(1), iterS->at(3), iterS->at(5));
            double dIntersectionLength = boxTest.IntersectionLength(vecL, vecR);
            if( dIntersectionLength > 0 ) {
               dIntsersectionLengthSum += dIntersectionLength;
            }
            double dDumbIntersectionLength = boxTest.DumbIntersectionLength(vecL, vecR);
            if( dDumbIntersectionLength > 0 ) {
               dDumbIntsersectionLengthSum += dDumbIntersectionLength;
            }
         }
         inOut << boxTest.Volume() << " " 
            << dIntsersectionLengthSum << " "
            << dDumbIntsersectionLengthSum << "\n";
      }
   }

   //Test the self similarity of the neuron by intersecting it with sub-volumes 
   //   inTrials times and measuring the moments.
   //For now the test volumes are all rectangular prisms
   //  so I can use some of the box counting code.
   //Does a raw dump of all the moments
   void CheckSelfSimilarMoments(ostream& inOut, const int& inTrials, const int inMaxMoment = 21) const {
      CRand3 rand3(rand());

      //A bounding box. Pad a little in the z-direction in 2D.
#ifdef TWO_DIMENSIONAL
      C3DBox boxBound(MinX(m_vbrNormalizedNeuron), MinY(m_vbrNormalizedNeuron), MinZ(m_vbrNormalizedNeuron)-1.,
         MaxX(m_vbrNormalizedNeuron), MaxY(m_vbrNormalizedNeuron), MaxZ(m_vbrNormalizedNeuron)+1.);
#else
      C3DBox boxBound(MinX(m_vbrNormalizedNeuron), MinY(m_vbrNormalizedNeuron), MinZ(m_vbrNormalizedNeuron),
         MaxX(m_vbrNormalizedNeuron), MaxY(m_vbrNormalizedNeuron), MaxZ(m_vbrNormalizedNeuron));
#endif

      //Get all the line segments
      vector<vector<double> > vecSegments;
      for(unsigned i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->AllSegments(vecSegments);
      }
      vector<vector<double> > vecPoints;
      double dCOMX, dCOMY, dCOMZ;
      TotalMeanSplit(vecSegments, 100, dCOMX, dCOMY, dCOMZ, vecPoints);
      vector<vector<double> >::const_iterator iterS;
      vector<vector<double> >::const_iterator iterS_end = vecPoints.end();

      //The main loop
      for(int i = 0; i < inTrials; ++i) {
         C3DBox boxTest = boxBound.RandomSubBox(rand3);
         C3DVector vCenter = boxTest.Center();
         vector<double> vecMoments(inMaxMoment,0.);
         for(iterS = vecPoints.begin(); iterS != iterS_end; ++iterS) {
            C3DVector vTemp(iterS->at(0), iterS->at(1), iterS->at(2));
            if( boxTest.IsPointInBox(vTemp) ) {
               vTemp -= vCenter;
               double dXYZ = vTemp.GetX()*vTemp.GetY();
#ifndef TWO_DIMENSIONAL
               dXYZ *= vTemp.GetZ();
#endif
               for(int j = 0; j < inMaxMoment; ++j) {
                  vecMoments[j] += iterS->at(3)*pow(dXYZ, j);
               }
            }
         }
         for(int j = 0; j < inMaxMoment; ++j) {
            inOut << vecMoments[j] << " ";
         }
         inOut << "\n";
      }
   }

   //This version of box count gets all the segments and maps them into boxes.
   //Output is dumped to the output stream.
   void BoxCount2(ostream& inOut, const unsigned& inLow, const unsigned& inStep, const unsigned& inHigh, const int inRepeat = 1) const {
      CRand3 rand3(rand());

      //Get all the line segments
      vector<vector<double> > vecSegments;
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->AllSegments(vecSegments);
      }

      //get the box sizes
      double dGlobalXMin = MinX(m_vbrNeuron) - 0.1*fabs(MinX(m_vbrNeuron));
      double dGlobalXMax = MaxX(m_vbrNeuron) + 0.1*fabs(MaxX(m_vbrNeuron));
      double dGlobalYMin = MinY(m_vbrNeuron) - 0.1*fabs(MinY(m_vbrNeuron));
      double dGlobalYMax = MaxY(m_vbrNeuron) + 0.1*fabs(MaxY(m_vbrNeuron));
      double dGlobalZMin = MinZ(m_vbrNeuron) - 0.1*fabs(MinZ(m_vbrNeuron));
      double dGlobalZMax = MaxZ(m_vbrNeuron) + 0.1*fabs(MaxZ(m_vbrNeuron));
      double dMinStep = max(dGlobalXMax-dGlobalXMin, max(dGlobalYMax-dGlobalYMin, dGlobalZMax-dGlobalZMin))/double(inHigh);

      vector<vector<double> >::iterator iterS = vecSegments.begin();
      vector<vector<double> >::iterator iterS_end = vecSegments.end();
      map<double, vector<int> > mapBoxCount;
      vector<int> vecBlank;

      //the main loop
      for(unsigned i = inLow; i <= inHigh; i += inStep) {
         //map in which to store the occupied boxes
         map<vector<int>, unsigned, CSort<vector<int> > > mapOccupied;

         //set the size of a box
         double dLocalStep = dMinStep*double(i);

         //Do repeated box counts with randomly (not really, but it shouldn't matter) rotated segments.
         //This should help get the area estimate.
         for(int j = 0; j < inRepeat; ++j) {
            //Set up the rotation matrix in terms of Euler angles.
#ifdef TWO_DIMENSIONAL
            double dPhi = 2.*PI*rand3();
            double dTheta = 0.0;
            double dPsi = 0.0;
            double dDX = dMinStep*(2.*rand3()-1.);
            double dDY = dMinStep*(2.*rand3()-1.);
            double dDZ = 0.0;
#else
            double dPhi = 2.*PI*rand3();
            double dTheta = PI*rand3();
            double dPsi = 2.*PI*rand3();
            double dDX = dMinStep*(2.*rand3()-1.);
            double dDY = dMinStep*(2.*rand3()-1.);
            double dDZ = dMinStep*(2.*rand3()-1.);
#endif
            double a11 = cos(dPsi)*cos(dPhi)-cos(dTheta)*sin(dPhi)*sin(dPsi);
            double a12 = cos(dPsi)*sin(dPhi)+cos(dTheta)*cos(dPhi)*sin(dPsi);
            double a13 = sin(dPsi)*sin(dTheta);
            double a21 = -1.*sin(dPsi)*cos(dPhi)-cos(dTheta)*sin(dPhi)*cos(dPsi);
            double a22 = -1.*sin(dPsi)*sin(dPhi)+cos(dTheta)*cos(dPhi)*cos(dPsi);
            double a23 = cos(dPsi)*sin(dTheta);
            double a31 = sin(dTheta)*sin(dPhi);
            double a32 = -1.*sin(dTheta)*cos(dPhi);
            double a33 = cos(dTheta);

            //Do the rotation.
            iterS = vecSegments.begin();
            iterS_end = vecSegments.end();
            for(; iterS != iterS_end; ++iterS) {
               double x1 = iterS->at(0);
               double x2 = iterS->at(1);
               double y1 = iterS->at(2);
               double y2 = iterS->at(3);
               double z1 = iterS->at(4);
               double z2 = iterS->at(5);
               iterS->at(0) = a11*x1 + a12*y1 + a13*z1 + dDX;
               iterS->at(1) = a11*x2 + a12*y2 + a13*z2 + dDX;
               iterS->at(2) = a21*x1 + a22*y1 + a23*z1 + dDY;
               iterS->at(3) = a21*x2 + a22*y2 + a23*z2 + dDY;
               iterS->at(4) = a31*x1 + a32*y1 + a33*z1 + dDZ;
               iterS->at(5) = a31*x2 + a32*y2 + a33*z2 + dDZ;
            }

            //Do the box counting.
            int iBoxCount = DoBoxCount(dLocalStep, vecSegments.begin(), vecSegments.end());
            mapBoxCount.insert(map<double, vector<int> >::value_type(dLocalStep, vecBlank)).first->second.push_back(iBoxCount);
         }
      }

      //Average and dump everything.
      map<double, vector<int> >::const_iterator iterC = mapBoxCount.begin();
      map<double, vector<int> >::const_iterator iterC_end = mapBoxCount.end();
      for(; iterC != iterC_end; ++iterC) {
         double dMean = 0.;
         vector<int>::const_iterator iterB = iterC->second.begin();
         vector<int>::const_iterator iterB_end = iterC->second.end();
         for(; iterB != iterB_end; ++iterB) {
            dMean += double(*iterB);
         }
         dMean /= double(iterC->second.size());

         double dStdErr = 0.;
         iterB = iterC->second.begin();
         for(; iterB != iterB_end; ++iterB) {
            dStdErr += (double(*iterB) - dMean)*(double(*iterB) - dMean);
         }
         dStdErr = sqrt(dStdErr);
         dStdErr /= double(iterC->second.size());

         inOut << iterC->first << " "
            << dMean << " "
            << dStdErr << "\n";
      }
   }

   //This version of box count gets all the segments and maps them into boxes.
   //Data is stored into the inReturn map with format: log<box size> log<occupied count>.
   void BoxCount2(const unsigned& inLow, const unsigned& inStep, const unsigned& inHigh, map<double, double>& inReturn) const {
      //Get all the line segments
      vector<vector<double> > vecSegments;
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->AllSegments(vecSegments);
      }

      //get the box sizes
      double dGlobalXMin = MinX(m_vbrNeuron) - 0.1*fabs(MinX(m_vbrNeuron));
      double dGlobalXMax = MaxX(m_vbrNeuron) + 0.1*fabs(MaxX(m_vbrNeuron));
      double dGlobalYMin = MinY(m_vbrNeuron) - 0.1*fabs(MinY(m_vbrNeuron));
      double dGlobalYMax = MaxY(m_vbrNeuron) + 0.1*fabs(MaxY(m_vbrNeuron));
      double dGlobalZMin = MinZ(m_vbrNeuron) - 0.1*fabs(MinZ(m_vbrNeuron));
      double dGlobalZMax = MaxZ(m_vbrNeuron) + 0.1*fabs(MaxZ(m_vbrNeuron));
      double dMinStep = max(dGlobalXMax-dGlobalXMin, max(dGlobalYMax-dGlobalYMin, dGlobalZMax-dGlobalZMin))/double(inHigh);

      //the main loop
      for(unsigned i = inLow; i <= inHigh; i += inStep) {
         //set the size of a box
         double dLocalStep = dMinStep*double(i);

         //get the box count
         int iBoxCount = DoBoxCount(dLocalStep, vecSegments.begin(), vecSegments.end());

         //store the current level
         inReturn[log(dLocalStep)] = log( double((iBoxCount>0)?iBoxCount:1) );
      }
   }


   //display the neuron (mostly testing)
   void DisplayNeuron(ostream& inOut) const {
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->RecursiveDisplay(inOut);
      }
   }

   //display the neuron (mostly testing)
   void UnformattedDisplayNeuron(ostream& inOut) const {
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->RawRecursiveDisplay(inOut);
      }
   }

   //dump the segments that make up the neuron for later drawing
   void DisplaySegments(ostream& inOut) const {
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         //get all the segments
         vector<vector<double> > vecSegments;
         m_vbrNeuron[i]->AllSegments(vecSegments);

         vector<vector<double> >::const_iterator iter = vecSegments.begin();
         vector<vector<double> >::const_iterator iter_end = vecSegments.end();
         for(; iter != iter_end; ++iter) {
            inOut << iter->at(0) << " "
               << iter->at(2) << " "
               << iter->at(4) << " "
               << iter->at(1) << " "
               << iter->at(3) << " "
               << iter->at(5) << " "
               << iter->at(6) << "\n";
         }
      }
   }

   //dump the segments that make up the neuron for later drawing
   void DisplayNormalizedSegments(ostream& inOut) const {
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         //get all the segments
         vector<vector<double> > vecSegments;
         m_vbrNormalizedNeuron[i]->AllSegments(vecSegments);

         vector<vector<double> >::const_iterator iter = vecSegments.begin();
         vector<vector<double> >::const_iterator iter_end = vecSegments.end();
         for(; iter != iter_end; ++iter) {
            inOut << iter->at(0) << " "
               << iter->at(2) << " "
               << iter->at(4) << " "
               << iter->at(1) << " "
               << iter->at(3) << " "
               << iter->at(5) << " "
               << iter->at(6) << "\n";
         }
      }
   }

   //Find a sub-arbor and write it to a stream for later analysis.
   //Should be called after a merge (only looks at the branch in m_vbrNeuron[0]).
   //In 3D, only writes anything if the sub-arbor avoids the boundaries in the z-direction.
   //  Hopefully, this will account for cutoff.
   //z-cuttoff defaults to 200 units (should be microns from the data).
   void SubArbor(ostream& inOut, string inType = "", double inZCutoff = 200.) const {
      CRand3 rand3( rand() );

      //Find the center (not COM) in the z-direction.
      double dCenterZ = 0.5*(m_vbrNeuron[0]->GetMaxZ()+m_vbrNeuron[0]->GetMinZ());

      //Select a random point from all the points.
      vector<vector<double> > vecvecAllSegments;
      m_vbrNeuron[0]->AllSegments(vecvecAllSegments);
      int iRandomPoint = rand3(vecvecAllSegments.size());
      CPoint pntTrim;
      pntTrim.push_back(vecvecAllSegments[iRandomPoint][0]);
      pntTrim.push_back(vecvecAllSegments[iRandomPoint][2]);
      pntTrim.push_back(vecvecAllSegments[iRandomPoint][4]);

      //Find the new root branch.
      shared_ptr<CBranch> brpIntersect = m_vbrNeuron[0]->DoesPointIntersect(pntTrim, .1);
      if( !brpIntersect->IsNull() ) {
#ifndef TWO_DIMENSIONAL
         if( brpIntersect->GetMaxZ() > dCenterZ + inZCutoff ||
            brpIntersect->GetMinZ() < dCenterZ - inZCutoff ) {
               cout << "CNeuron.h gh1 " << brpIntersect->GetMaxZ() << " "
                  << brpIntersect->GetMinZ() << " " << dCenterZ << "\n" << flush;
               return;
         }
#endif
         //Write the branch
         inOut << "( (Color 1 )\n  (" << inType << ")\n";
         brpIntersect->RecursiveTrimmedDisplay(inOut);
         inOut << ")\n\n";
      } else {
         cout << "Error in CNeuron::SubBranch: Random point does not intersect...returning\n" << flush;
         return;
      }

   }

   //Attempt to merge the seperate arbors.
   //Repeatedly merges arbors onto the first one.
   void MergeArbors() {
      if(m_vbrNeuron.size() > 0) {
         DoMergeArbors(m_vbrNeuron);
         m_vbrNeuron.at(0)->TrimDaughter();
         m_vbrNeuron.at(0)->SetLevel(0);
         m_vbrNeuron.at(0)->SortDaughters();
         m_vbrNeuron.at(0)->UpdateAllMaxMin();

         DoMergeArbors(m_vbrNormalizedNeuron);
         m_vbrNormalizedNeuron.at(0)->TrimDaughter();
         m_vbrNormalizedNeuron.at(0)->SetLevel(0);
         m_vbrNormalizedNeuron.at(0)->SortDaughters();
         m_vbrNormalizedNeuron.at(0)->UpdateAllMaxMin();
      }
   }

   //Trim the first arbor in the list up to the first daughter.
   //The transformation looks like:
   //
   //                    /-----         /-----
   //         ----------/          ->  /
   //                   \              \
   //                    \----          \----
   //
   //This assumes that the root branch is the dangling part (not necessarily true).
   void TrimToFirstDaughter() {
      if(m_vbrNeuron.size() > 0) {
         m_vbrNeuron[0]->TrimToFirstDaughter();
      }
   }

   //Trim the first arbor in the list (call after a merge and it should work best).
   void Trim(const CPoint& inPoint) {
      //Find the new root branch.
      shared_ptr<CBranch> brpIntersect = m_vbrNeuron[0]->DoesPointIntersect(inPoint, .1);
      if( !brpIntersect->IsNull() ) {
         //Set the new branch as a root branch.
         //This will delete any old data.
         brpIntersect->SetRoot();

         //Just set the first member of vbrNeuron to be the new root branch.
         //The fancy new smart pointers (thank you boost! ) should handle memory.
         m_vbrNeuron[0] = brpIntersect;
      } else {
         cout << "The input point does not intersect\n" << flush;
      }
   }

   //Write data in swc format. This is easier and more realiable than .asc.
   void WriteSWCData(ostream& inOut, const int& inType) const {
      time_t timeCurrent;
      time( &timeCurrent );
      inOut << "# Read by 3d_neuron\n";
      inOut << "# Date: " << ctime(&timeCurrent);
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         m_vbrNeuron[i]->RecursiveDisplaySWC(inOut, inType);
      }
   }

   //Write the data to a form that neuroleucidia will recognize (or at least this).
   void WriteData(ostream& inOut, string strType = "") const {
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         inOut << "( (Color " << i << ")\n  (" << strType << ")\n";
         m_vbrNeuron[i]->RecursiveDisplay(inOut);
         inOut << ")\n\n";
      }
   }

   //Write trimmed data to a form that neuroleucidia will recognize (or at least this).
   void WriteTrimmedData(ostream& inOut, string strType = "") const {
      for(unsigned i = 0; i < m_vbrNeuron.size(); ++i) {
         inOut << "( (Color " << i << ")\n  (" << strType << ")\n";
         m_vbrNeuron[i]->RecursiveTrimmedDisplay(inOut);
         inOut << ")\n\n";
      }
   }

   //Read in the data from an input stream.
   //Parses the stream into different branches as marked by an opening '('
   // through a closing ')'
   //The lines containing the '(' and the ')' are stripped.
   //Blank lines and ones beginning with ';' are stripped.
   void ReadData(istream& inFile, string strType = "") {
      string strLine;
      while( getline(inFile, strLine, '\n') ) {
         //skip if the first element is a ';' or the line is empty
         if(strLine[0] != ';' && strLine.find_first_not_of(' ') != string::npos) {
            int iParenthesisCount = ParenthesisCount(strLine);
            if(iParenthesisCount > 0) {
               //A new branch.
               vector<string> vecLines;
               vecLines.push_back(strLine);
               while(iParenthesisCount > 0 && getline(inFile, strLine, '\n')) {
                  if(strLine.find_first_not_of(' ') != string::npos) {
                     iParenthesisCount += ParenthesisCount(strLine);
                     vecLines.push_back(strLine);
                  }
               }
               vecLines.pop_back();

               //Only look at trees with a specific id, i.e. "Axon", ...
               if( vecLines.size() > 2 ) {
                  bool bRightType = false;
                  //Check all the lines.
                  //Could be trouble if the tree has the key word embedded in it.
                  for(int i = 0; i < vecLines.size(); ++i) {
                     bRightType |= (vecLines[i].find(strType) != string::npos);
                  }
                  if( bRightType ) {
                     ParseData(vecLines.begin(), vecLines.end(), m_vbrNeuron);
                     ParseData(vecLines.begin(), vecLines.end(), m_vbrNormalizedNeuron);
                  }
               }
            }
         }
      }

      //Normalize.
      Normalize();
   }

   //Read in the data from an input stream with SWC format.
   //Blank lines and ones beginning with '#' are stripped.
   void ReadSWCData(istream& inFile, int inType) {
      map<int, CPoint> mapPoints;
      CPoint pntCurrentPoint;
      pntCurrentPoint.resize(3);
      set<int> setAddedPoints;

      shared_ptr<CBranch> brpParent(new CBranch(-1));
      shared_ptr<CBranch> brpCurrent(new CBranch(-1));
      shared_ptr<CBranch> brpNParent(new CBranch(-1));
      shared_ptr<CBranch> brpNCurrent(new CBranch(-1));

      //Structure to record a map from ID to branch.
      map<int, shared_ptr<CBranch> > mapID_Branch;
      map<int, shared_ptr<CBranch> > mapID_NBranch;

      int iID, iType, iParent;
      int iOldType = -1;
      int iCount = 0;
      double dX, dY, dZ, dW;
      string strLine;
      while( inFile ) {
         //skip if the first element is a '#'
         if( inFile.peek() == '#' ) {
            if( !getline(inFile, strLine, '\n') ) {
               cerr << "Fatal parse error (1): stopping\n" << flush;
               exit(-1);
            }
         } else {
            if(inFile >> iID >> iType >> dX >> dY >> dZ >> dW >> iParent) {
               //Record the point no matter what (just in case).
               pntCurrentPoint[0] = dX;
               pntCurrentPoint[1] = dY;
#ifdef TWO_DIMENSIONAL
               pntCurrentPoint[2] = 0.;
               dZ = 0.;
#else
               pntCurrentPoint[2] = dZ;
#endif
               mapPoints[iID] = pntCurrentPoint;
               if(iType == inType || iType == 1) {
                  setAddedPoints.insert(iID);
                  ++iCount;
                  if(setAddedPoints.find(iParent) == setAddedPoints.end()) {
                     //The parent was not added so this is a new branch.
                     shared_ptr<CBranch> brpParentTemp(new CBranch(-1));
                     brpParent = brpParentTemp;
                     shared_ptr<CBranch> brpCurrentTemp(new CBranch(0));
                     brpCurrent = brpCurrentTemp;
                     brpCurrent->SetRoot();
                     m_vbrNeuron.push_back(brpCurrent);

                     shared_ptr<CBranch> brpNParentTemp(new CBranch(-1));
                     brpNParent = brpNParentTemp;
                     shared_ptr<CBranch> brpNCurrentTemp(new CBranch(0));
                     brpNCurrent = brpNCurrentTemp;
                     brpNCurrent->SetRoot();
                     m_vbrNormalizedNeuron.push_back(brpNCurrent);

                     brpCurrent->AddPoint(dX, dY, dZ);
                     brpNCurrent->AddPoint(dX, dY, dZ);
                  } else if (iID == iParent + 1) {
                     //new point on the same branch
                     brpCurrent->AddPoint(dX, dY, dZ);
                     brpNCurrent->AddPoint(dX, dY, dZ);
                  } else {
                     //Create a daughter branch.
                     //Find the branch that intersects this point's parent.
                     CPoint pntParent = mapPoints[iParent];
                     brpParent = mapID_Branch[iParent];
                     brpNParent = mapID_NBranch[iParent];
                     if(brpParent->IsNull() || brpNParent->IsNull()) {
                        cerr << "Error: unable to find parent. Stopping read\n" << flush;
                        return;
                     } else {
                        //the raw arbor
                        shared_ptr<CBranch> brpCurrentTemp(new CBranch(brpParent));
                        brpCurrent = brpCurrentTemp;
                        brpParent->AddDaughter(brpCurrent);
                        //add the branchpoint from the parent as the first point of this branch
                        brpCurrent->AddPoint(pntParent[0], pntParent[1], pntParent[2]);
                        brpCurrent->AddPoint(dX, dY, dZ);
                        //cout << dX << " " << dY << " " << dZ << "\n" << flush;

                        //the normalized arbor
                        shared_ptr<CBranch> brpNCurrentTemp(new CBranch(brpNParent));
                        brpNCurrent = brpNCurrentTemp;
                        brpNParent->AddDaughter(brpNCurrent);
                        //add the branchpoint from the parent as the first point of this branch
                        brpNCurrent->AddPoint(pntParent[0], pntParent[1], pntParent[2]);
                        brpNCurrent->AddPoint(dX, dY, dZ);
                     }
                  }
                  //Record the branch of this ID.
                  mapID_Branch[iID] = brpCurrent;
                  mapID_NBranch[iID] = brpNCurrent;
               }
            }
         }
      }

      if(iCount < 2) {
         //Nothing of use was read.
         m_vbrNeuron.clear();
         m_vbrNormalizedNeuron.clear();
      } else {
         //Normalize.
         Normalize();
      }
   }

private:
   //helpers

   //Return the length of the input segment.
   double SegmentLength(const vector<double>& inS) const {
#ifdef TWO_DIMENSIONAL
      return sqrt( (inS[1]-inS[0])*(inS[1]-inS[0]) +
         (inS[3]-inS[2])*(inS[3]-inS[2]) );
#else
      return sqrt( (inS[1]-inS[0])*(inS[1]-inS[0]) +
         (inS[3]-inS[2])*(inS[3]-inS[2]) +
         (inS[5]-inS[4])*(inS[5]-inS[4]) );
#endif
   }

   //Find the minimum and maximum values of various coordinates.
   double MinX(const vector<shared_ptr<CBranch> >& inNeuron) const {
      if(inNeuron.size() == 0) {
         return 0.;
      }
      double dReturn = inNeuron[0]->GetMinX();
      for(unsigned i = 1; i < inNeuron.size(); ++i) {
         dReturn = min(dReturn, inNeuron[i]->GetMinX());
      }
      return dReturn;
   }
   double MinY(const vector<shared_ptr<CBranch> >& inNeuron) const {
      if(inNeuron.size() == 0) {
         return 0.;
      }
      double dReturn = inNeuron[0]->GetMinY();
      for(unsigned i = 1; i < inNeuron.size(); ++i) {
         dReturn = min(dReturn, inNeuron[i]->GetMinY());
      }
      return dReturn;
   }
   double MinZ(const vector<shared_ptr<CBranch> >& inNeuron) const {
      if(inNeuron.size() == 0) {
         return 0.;
      }
      double dReturn = inNeuron[0]->GetMinZ();
      for(unsigned i = 1; i < inNeuron.size(); ++i) {
         dReturn = min(dReturn, inNeuron[i]->GetMinZ());
      }
      return dReturn;
   }
   double MaxX(const vector<shared_ptr<CBranch> >& inNeuron) const {
      if(inNeuron.size() == 0) {
         return 0.;
      }
      double dReturn = inNeuron[0]->GetMaxX();
      for(unsigned i = 1; i < inNeuron.size(); ++i) {
         dReturn = max(dReturn, inNeuron[i]->GetMaxX());
      }
      return dReturn;
   }
   double MaxY(const vector<shared_ptr<CBranch> >& inNeuron) const {
      if(inNeuron.size() == 0) {
         return 0.;
      }
      double dReturn = inNeuron[0]->GetMaxY();
      for(unsigned i = 1; i < inNeuron.size(); ++i) {
         dReturn = max(dReturn, inNeuron[i]->GetMaxY());
      }
      return dReturn;
   }
   double MaxZ(const vector<shared_ptr<CBranch> >& inNeuron) const {
      if(inNeuron.size() == 0) {
         return 0.;
      }
      double dReturn = inNeuron[0]->GetMaxZ();
      for(unsigned i = 1; i < inNeuron.size(); ++i) {
         dReturn = max(dReturn, inNeuron[i]->GetMaxZ());
      }
      return dReturn;
   }

   //Normalize the neurons such that it is alligned with the axis.
   //Modifies m_vbrNormalizedNeuron.
   //Changed 10/06, Joe Snider: Scaling was a bad idea, rotating and shifting are good.
   //Changed 2/08 Joe Snider: this will now recurse up to inLimit times to improve the accuracy.
   //                         Rotate by v transpose now.
   //                         Correctly shifts to COM.
   //                         Always runs at least twice.
   void Normalize(int inLimit = 10) {
      //merge the arbors so they never end at a branch point
      if(m_bFirstPass) {
         //merge so there's only 1 arbor.
         MergeArbors();
#ifdef TRIM_TO_FIRST_DAUGHTER
         TrimToFirstDaughter();
#endif
      }

      gsl_set_error_handler_off();
      unsigned i;
      //get all the segments
      vector<vector<double> > vecSegments;
      for(i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->AllSegments(vecSegments);
      }
      if( vecSegments.size() < 2 ) {return;}

      //find the center of mass to use as the origin.
      double dCOMX, dCOMY, dCOMZ;
      vector<vector<double> > vecPoints;
      TotalMean(vecSegments, dCOMX, dCOMY, dCOMZ, vecPoints);

      //Shift the segments to the COM
      for(int i = 0; i < vecSegments.size(); ++i) {
         vecSegments[0][0] -= dCOMX;
         vecSegments[0][1] -= dCOMX;
         vecSegments[0][2] -= dCOMY;
         vecSegments[0][3] -= dCOMY;
         vecSegments[0][4] -= dCOMZ;
         vecSegments[0][5] -= dCOMZ;
      }

      //cout << "CNeuron.h gh0 " << dCOMX << " " << dCOMY << " " << dCOMZ << "\n\n" << flush;

      //Create some things needed for the gsl integrator.
      gsl_integration_workspace* w =
         gsl_integration_workspace_alloc(1000);  //1000 is plenty of room
      double dResult, dError;
      double dAlpha[9];
      gsl_function F1;
      F1.function = &LineIntegratorCovariances;
      F1.params = &dAlpha;

      //Find the moment of inertia tensor.
      //Uses the TNT library from NIST.
      Array2D<double> a2dInertiaTensor(3,3);
      a2dInertiaTensor = 0.0;
      double m, xx, xy, xz, yy, yz, zz;
      double dM = 0.;
      vector<vector<double> >::const_iterator iter;
      vector<vector<double> >::const_iterator iter_end;
      dAlpha[0] = 1.;
      iter = vecSegments.begin();
      iter_end = vecSegments.end();
      for(; iter != iter_end; ++iter) {
         m = SegmentLength(*iter);
         for(int j = 0; j < 7; ++j) {
            dAlpha[j+1] = iter->at(j);
         }
         //do each integral
         int status = 0;
         dAlpha[7] = 0; dAlpha[8] = 0;
         status += gsl_integration_qag(&F1, 0, 1, 0, 1.e-7, 1000, GSL_INTEG_GAUSS61, 
            w, &xx, &dError);
         dAlpha[7] = 0; dAlpha[8] = 1;
         status += gsl_integration_qag(&F1, 0, 1, 0, 1.e-7, 1000, GSL_INTEG_GAUSS61, 
            w, &xy, &dError);
         dAlpha[7] = 0; dAlpha[8] = 2;
         status += gsl_integration_qag(&F1, 0, 1, 0, 1.e-7, 1000, GSL_INTEG_GAUSS61, 
            w, &xz, &dError);
         dAlpha[7] = 1; dAlpha[8] = 1;
         status += gsl_integration_qag(&F1, 0, 1, 0, 1.e-7, 1000, GSL_INTEG_GAUSS61, 
            w, &yy, &dError);
         dAlpha[7] = 1; dAlpha[8] = 2;
         status += gsl_integration_qag(&F1, 0, 1, 0, 1.e-7, 1000, GSL_INTEG_GAUSS61, 
            w, &yz, &dError);
         dAlpha[7] = 2; dAlpha[8] = 2;
         status += gsl_integration_qag(&F1, 0, 1, 0, 1.e-7, 1000, GSL_INTEG_GAUSS61, 
            w, &zz, &dError);

         if(status != 0) {
            cerr << "Warning: an integral failed on the segment:\n      ("
               << iter->at(0) << ", "
               << iter->at(2) << ", "
               << iter->at(4) << ") - ("
               << iter->at(1) << ", "
               << iter->at(3) << ", "
               << iter->at(5) << ")\n continuing with best estimate\n" << flush;
         }

         a2dInertiaTensor[0][0] += m*(yy+zz);
         a2dInertiaTensor[0][1] -= m*xy;
         a2dInertiaTensor[0][2] -= m*xz;
         a2dInertiaTensor[1][1] += m*(xx+zz);
         a2dInertiaTensor[1][0] -= m*xy;
         a2dInertiaTensor[1][2] -= m*yz;
         a2dInertiaTensor[2][2] += m*(xx+yy);
         a2dInertiaTensor[2][0] -= m*xz;
         a2dInertiaTensor[2][1] -= m*yz;
         dM += m;
      }

      /*cout << "\nCNeuron.h gh1 " 
         << a2dInertiaTensor[0][0] << " "
         << a2dInertiaTensor[0][1] << " "
         << a2dInertiaTensor[0][2] << " "
         << "\n" << flush;
      cout << "CNeuron.h gh2 " 
         << a2dInertiaTensor[1][0] << " "
         << a2dInertiaTensor[1][1] << " "
         << a2dInertiaTensor[1][2] << " "
         << "\n" << flush;
      cout << "CNeuron.h gh3 " 
         << a2dInertiaTensor[2][0] << " "
         << a2dInertiaTensor[2][1] << " "
         << a2dInertiaTensor[2][2] << " "
         << "\n\n" << flush;*/

      //clean up the workspace
      gsl_integration_workspace_free(w);

      //Do the eigenvalue problem (note: moment of inertia is hermitian).
      //Not the best variable names, but there'll be alot of V[1][2]*x+..., so it's best.
      Eigenvalue<double> E1(a2dInertiaTensor);
      Array2D<double> D(3,3);
      Array2D<double> V(3,3);
      E1.getD(D);
      E1.getV(V);

      /*cout << "\nCNeuron.h gh4 " 
         << V[0][0] << " "
         << V[0][1] << " "
         << V[0][2] << " "
         << "\n" << flush;
      cout << "CNeuron.h gh5 " 
         << V[1][0] << " "
         << V[1][1] << " "
         << V[1][2] << " "
         << "\n" << flush;
      cout << "CNeuron.h gh6 " 
         << V[2][0] << " "
         << V[2][1] << " "
         << V[2][2] << " "
         << "\n\n" << flush;*/

      //a2dInertiaTensor = matmult( a2dInertiaTensor, V );
      //a2dInertiaTensor = matmult( transpose(V), a2dInertiaTensor );
      a2dInertiaTensor = matmult( a2dInertiaTensor, transpose(V) );
      a2dInertiaTensor = matmult( V, a2dInertiaTensor );
      /*cout << "\nCNeuron.h gh7 " 
         << a2dInertiaTensor[0][0] << " "
         << a2dInertiaTensor[0][1] << " "
         << a2dInertiaTensor[0][2] << " "
         << "\n" << flush;
      cout << "CNeuron.h gh8 " 
         << a2dInertiaTensor[1][0] << " "
         << a2dInertiaTensor[1][1] << " "
         << a2dInertiaTensor[1][2] << " "
         << "\n" << flush;
      cout << "CNeuron.h gh9 " 
         << a2dInertiaTensor[2][0] << " "
         << a2dInertiaTensor[2][1] << " "
         << a2dInertiaTensor[2][2] << " "
         << "\n\n" << flush;*/

#ifdef ROTATE_AROUND_Z
      for(i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->Shift(-1.*dCOMX, -1.*dCOMY, -1.*dCOMZ);
      }
      double dBFAngle = Test2DRotation();
      V[0][0] = cos(dBFAngle);
      V[1][1] = cos(dBFAngle);
      V[0][1] = -1.*sin(dBFAngle);
      V[1][0] = sin(dBFAngle);
      V[2][2] = 1.;
      V[2][0] = 0.;
      V[2][1] = 0.;
      V[0][2] = 0.;
      V[1][2] = 0.;
      for(i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->Rotate( transpose(V) );
      }
      //stop recursion, will send a warning
      cerr << "Warning: rotating about only the z-axis, will send a false recursion warning\n" << flush;
      inLimit = 0;
#else //ROTATE_AROUND_Z
      //Shift and rotate the arbors (in that order).
      for(i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->Shift(-1.*dCOMX, -1.*dCOMY, -1.*dCOMZ);
         m_vbrNormalizedNeuron[i]->Rotate( transpose(V) );
      }
#endif
      
      //Recurse to get round-off errors
      m_bFirstPass = false;
      if(inLimit > 0) {
         if(inLimit == 10) {
            Normalize(inLimit-1);
         }
         //cout << "CNeuron.h gh11 ........................."
         //   << (fabs(a2dInertiaTensor[0][1]) + fabs(a2dInertiaTensor[0][2]) + fabs(a2dInertiaTensor[1][2])) << " " 
         //   << a2dInertiaTensor[0][0] + a2dInertiaTensor[1][1] + a2dInertiaTensor[2][2] << "\n" << flush;
         if( (fabs(a2dInertiaTensor[0][1]) + fabs(a2dInertiaTensor[0][2]) + fabs(a2dInertiaTensor[1][2]))/
             (a2dInertiaTensor[0][0] + a2dInertiaTensor[1][1] + a2dInertiaTensor[2][2]) > 1.e-6 ) {
            //cout << "CNeuron.h gh12 .........................\n" << flush;
            Normalize(inLimit-1);
         }
      }
      if( inLimit < 1 ) {
         cerr << "Warning in CNeuron::Normalize(): recursion limit reached...continuing with best estimate\n" << flush;
      }

      //////////////////////////////////////////////////////////begin hack/////////////////////////////////////////////
      ////Make a non-linear transformation to destroy scaling.
      ////Warn that it's happening because this is a hack.
      //cerr << "Warning: CNeuron::Normalize uses the non-linear transformation hack\n" << flush;
      //CRand3 rand3(rand());
      //double dScaleHack = 1000.*(1.5*rand3()+5.)/TotalLength();
      //for(i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
      //   //m_vbrNormalizedNeuron[i]->NonlinearHack(TotalLength());
      //   m_vbrNormalizedNeuron[i]->Scale(dScaleHack, dScaleHack, dScaleHack);
      //}
      //////////////////////////////////////////////////////////end hack//////////////////////////////////////////////

      //Store some convenient but costly numbers.
      m_dWidthX = MaxX(m_vbrNormalizedNeuron) - MinX(m_vbrNormalizedNeuron);
      m_dStandardDeviationX = sqrt(.5*(-1.*D[0][0] + D[1][1] + D[2][2])/dM);
      m_dWidthY = MaxY(m_vbrNormalizedNeuron) - MinY(m_vbrNormalizedNeuron);
      m_dStandardDeviationY = sqrt(.5*(D[0][0] - D[1][1] + D[2][2])/dM);
#ifdef TWO_DIMENSIONAL
      m_dWidthZ = 0.;
      m_dStandardDeviationZ = 0.;
#else
      m_dWidthZ = MaxZ(m_vbrNormalizedNeuron) - MinZ(m_vbrNormalizedNeuron);
      m_dStandardDeviationZ = sqrt(.5*(D[0][0] + D[1][1] - D[2][2])/dM);
#endif
   }

   //Attempt to merge the seperate arbors.
   //Repeatedly merges arbors onto the first one.
   void DoMergeArbors(vector<shared_ptr<CBranch> >& inNeuron) {
      double dMergeDistance = 1.;
      while(inNeuron.size() > 1) {
         map<double, int> mapDistID;
         //find an arbor that intersects
         for(unsigned i = 1; i < inNeuron.size(); ++i) {
            if(inNeuron[i]->GetLine().size() > 1) {
               //Check if the arbor intersects the first one.
               CPoint pSubBranchRoot = inNeuron[i]->GetLine().front();
               mapDistID[inNeuron[0]->MinimumDistance(pSubBranchRoot[0], pSubBranchRoot[1], pSubBranchRoot[2])] = i;
            } else {
               cout << "Warning in CNeuron::DoMergeArbors: " << i << " has no points...removing it\n" << flush;
               inNeuron.erase(inNeuron.begin() + i);
               i = 0;
            }
         }

         //empty arbors may have been removed.
         if(inNeuron.size() > 1) {
            //Merge the closest arbor
            //Find the closest branch (this could be optimized, so could MinimumDistance for that matter).
            int iClosestID = mapDistID.begin()->second;
            double dClosestDist = mapDistID.begin()->first;
            CPoint pSubBranchRoot = inNeuron[iClosestID]->GetLine().front();
            shared_ptr<CBranch> brpIntersect = inNeuron[0]->DoesPointIntersect(pSubBranchRoot, dClosestDist+1.e-6);
            if(brpIntersect->IsNull()) {
               cout << "Error in CNeuron::DoMergeArbors: unable to find intersecting branch...halting merge\n" << flush;
               return;
            }

            //Find the intersection point.
            CPoint pClosest = brpIntersect->GetLine().FindClosestPoint(pSubBranchRoot);

            ////Insert the intersection point into the current line.
            //inNeuron[0]->InsertPoint(pClosest);
            inNeuron[iClosestID]->AddPointFront(pClosest[0], pClosest[1], pClosest[2]);

            //Add the old branch as a daughter of the first branch
            brpIntersect->AddDaughter(inNeuron[iClosestID]);

            //Remove the pointer to the old arbor from the inNeuron list.
            inNeuron.erase(inNeuron.begin() + iClosestID);
         }

      }

      //Update the max and mins and the level.
      if(inNeuron.size() > 0) {
         inNeuron[0]->UpdateAllMaxMin();
         inNeuron[0]->SetLevel(0);
      }
   }

   //Do the box counting.
   //Grid is centered at (0,0,0) and has unit vectors (inStep,0,0),...
   //Type T must be castable to an stl style iterator of vector<double> (vector<vector<double> >::const_iterator works).
   //Returns the number of occupied boxes of size inStep.
   template<class T>
   int DoBoxCount(const double& inStep, T inSegmentsBegin, T inSegmentsEnd) const {
      //box for testing.
      C3DBox boxTest;
      //map in which to store the occupied boxes.
      //Lower left corner divided by inStep and rounded to double equality problems.
      map<vector<int>, unsigned, CSort<vector<int> > > mapOccupied;
      vector<int> vecInsert(3);

      //loop over all the segments and find the boxes they intersect
      for(; inSegmentsBegin != inSegmentsEnd; ++inSegmentsBegin) {
         //Generate all the boxes that may intersect this segment and test them.
         //2D case handled slightly differently.
#ifdef TWO_DIMENSIONAL
         C3DVector vecL(inSegmentsBegin->at(0)/inStep, inSegmentsBegin->at(2)/inStep, 0.);
         C3DVector vecR(inSegmentsBegin->at(1)/inStep, inSegmentsBegin->at(3)/inStep, 0.);
         int iXMin = Round(min(vecL.GetX(), vecR.GetX()))-1;
         int iXMax = Round(max(vecL.GetX(), vecR.GetX()));
         int iYMin = Round(min(vecL.GetY(), vecR.GetY()))-1;
         int iYMax = Round(max(vecL.GetY(), vecR.GetY()));
         for(int iX = iXMin; iX <= iXMax; ++iX) {
            for(int iY = iYMin; iY <= iYMax; ++iY) {
               boxTest.SetAB(iX, iY, -.1, iX+1, iY+1, .1);
               if( boxTest.DoesSegmentIntersect(vecL, vecR) ) {
                  vecInsert[0] = iX;
                  vecInsert[1] = iY;
                  vecInsert[2] = 0;
                  mapOccupied.insert(map<vector<int>, unsigned, CSort<vector<int> > >::value_type(vecInsert, 0)).first->second += 1;
               }
            }
         }
#else
         C3DVector vecL(inSegmentsBegin->at(0)/inStep, 
            inSegmentsBegin->at(2)/inStep, 
            inSegmentsBegin->at(4)/inStep);
         C3DVector vecR(inSegmentsBegin->at(1)/inStep, 
            inSegmentsBegin->at(3)/inStep, 
            inSegmentsBegin->at(5)/inStep);
         int iXMin = Round(min(vecL.GetX(), vecR.GetX()))-1;
         int iXMax = Round(max(vecL.GetX(), vecR.GetX()));
         int iYMin = Round(min(vecL.GetY(), vecR.GetY()))-1;
         int iYMax = Round(max(vecL.GetY(), vecR.GetY()));
         int iZMin = Round(min(vecL.GetZ(), vecR.GetZ()))-1;
         int iZMax = Round(max(vecL.GetZ(), vecR.GetZ()));
         for(int iX = iXMin; iX <= iXMax; ++iX) {
            for(int iY = iYMin; iY <= iYMax; ++iY) {
               for(int iZ = iZMin; iZ <= iZMax; ++iZ) {
                  boxTest.SetAB(iX, iY, iZ, iX+1, iY+1, iZ+1);
                  if( boxTest.DoesSegmentIntersect(vecL, vecR) ) {
                     vecInsert[0] = iX;
                     vecInsert[1] = iY;
                     vecInsert[2] = iZ;
                     mapOccupied.insert(map<vector<int>, unsigned, CSort<vector<int> > >::value_type(vecInsert, 0)).first->second += 1;
                  }
               }
            }
         }
#endif
      }

      //////////////////////////////////////////test/////////////////////////////////////////////////
      //if( fabs(0.15-inStep) < 0.01 ) {
      //   cout << "CNeuron.h gh1\n" << flush;
      //   ofstream ofTest("test_box_count.txt");
      //   map<vector<int>, unsigned, CSort<vector<int> > >::const_iterator iter = mapOccupied.begin();
      //   map<vector<int>, unsigned, CSort<vector<int> > >::const_iterator iter_end = mapOccupied.end();
      //  for(; iter != iter_end; ++iter) {
      //      ofTest << iter->first.at(0) << " "
      //        << iter->first.at(1) << " "
      //        << iter->first.at(2) << "\n";
      //   }
      //   ofTest.close();
      //}
      //////////////////////////////////////////end test/////////////////////////////////////////////

      return (int)mapOccupied.size();
      //return int(GetOccupiedCount(mapOccupied)+0.5);
   }

   //Do the actual parsing.
   //T is convertable to stl style iterators through an array of strings.
   template<class T>
   void ParseData(T inBegin, T inEnd, vector<shared_ptr<CBranch> >& inStorage) {
      CRand3 rand3( rand() );   //used by the error hack, doesn't hurt
      shared_ptr<CBranch> brpParent(new CBranch(-1));
      shared_ptr<CBranch> brpCurrent(new CBranch(-1));
      char cPosition;
      string strLine;
      size_t uFirstChar, uSecondChar, uLeft, uRight;
      double dX, dY, dZ;
      for(; inBegin != inEnd; ++inBegin) {
         strLine = *inBegin;
         uFirstChar = strLine.find_first_not_of(' ');
         cPosition = strLine[uFirstChar];
         uSecondChar = strLine.find_first_not_of(' ', uFirstChar+1);
         if(cPosition == '(') {
            if(uFirstChar >= strLine.size()-2 || uFirstChar == 0) {
               //new node
               if( brpCurrent->IsNull() ) {
                  //Create the root branch.
                  //shared_ptr<CBranch> brpParentTemp(new CBranch(NULL));
                  //brpParent = brpParentTemp;
                  shared_ptr<CBranch> brpCurrentTemp(new CBranch(NULL));
                  brpCurrent = brpCurrentTemp;
                  inStorage.push_back(brpCurrent);
               } else {
                  //Create a daughter branch. 
                  brpParent = brpCurrent;
                  shared_ptr<CBranch> brpCurrentTemp(new CBranch(brpParent));
                  brpCurrent = brpCurrentTemp;
                  brpParent->AddDaughter(brpCurrent);
                  //add the branchpoint from the parent as the first point of this branch
                  brpCurrent->AddPoint(brpParent->GetLine().back()[0], 
                     brpParent->GetLine().back()[1], 
                     brpParent->GetLine().back()[2]);
               }
            } else if (strLine[uSecondChar] == '-' ||
               strLine[uSecondChar] == '0' ||
               strLine[uSecondChar] == '1' ||
               strLine[uSecondChar] == '2' ||
               strLine[uSecondChar] == '3' ||
               strLine[uSecondChar] == '4' ||
               strLine[uSecondChar] == '5' ||
               strLine[uSecondChar] == '6' ||
               strLine[uSecondChar] == '7' ||
               strLine[uSecondChar] == '8' ||
               strLine[uSecondChar] == '9') {
                  //point in the node
                  uLeft = uSecondChar;
                  uRight = strLine.find_first_of(' ', uLeft);
                  dX = atof(strLine.substr(uLeft, uRight-uLeft+1).data());
                  //cout << "CNeuron.h gh5 " << strLine.substr(uLeft, uRight-uLeft+1) << " " << dX << "\n" << flush;

                  uLeft = strLine.find_first_not_of(' ', uRight);
                  uRight = strLine.find_first_of(' ', uLeft);
                  dY = atof(strLine.substr(uLeft, uRight-uLeft+1).data());
                  //cout << "CNeuron.h gh6 " << strLine.substr(uLeft, uRight-uLeft+1) << " " << dY << "\n" << flush;

                  uLeft = strLine.find_first_not_of(' ', uRight);
                  uRight = strLine.find_first_of(' ', uLeft);
                  dZ = atof(strLine.substr(uLeft, uRight-uLeft+1).data());
                  //cout << "CNeuron.h gh7 " << strLine.substr(uLeft, uRight-uLeft+1) << " " << dZ << "\n" << flush;
#ifdef ADD_ERROR_ASC_HACK
                  dX += ERROR_ASC_HACK*(2.*rand3()-1.);
                  dY += ERROR_ASC_HACK*(2.*rand3()-1.);
                  dZ += ERROR_ASC_HACK*(2.*rand3()-1.);
#endif
#ifdef TWO_DIMENSIONAL
                  brpCurrent->AddPoint(dX, dY, 0.);
#else
                  brpCurrent->AddPoint(dX, dY, dZ);
#endif
            } else {
               //some irritating label embedded in the data
               //skip until the parenthesis close
               string strSkipLine;
               int iPLevel = ParenthesisCount(strLine);
               if(iPLevel != 0) {
                  do {
                     ++inBegin;
                     if(inBegin == inEnd) {
                        cerr << "Error: unknown parse error.\nExiting with prejudice\n" << flush;
                        exit(-1);
                     }
                     strSkipLine = *inBegin;
                     iPLevel += ParenthesisCount(strSkipLine);
                  } while(iPLevel > 0);
               }
            }
         } else if (cPosition == '|') {
            //new split off the same parent (no more points to the current one after this)
            shared_ptr<CBranch> brpCurrentTemp(new CBranch(brpParent));
            brpCurrent = brpCurrentTemp;
            brpParent->AddDaughter(brpCurrent);
            brpCurrent->AddPoint(brpParent->GetLine().back()[0], 
               brpParent->GetLine().back()[1], 
               brpParent->GetLine().back()[2]);
         } else if (cPosition == ')') {
            //go up hill
            if(brpParent != NULL) {
               brpCurrent = brpParent;
               brpParent = brpCurrent->GetParent();
               //inFile.getline(cstrJunk, 1000, '\n');
            } else {
               cerr << "Warning: mismatched parenthesis (too many ')')...attempting to continue\n" << flush;
            }
         }
      }
   }

   //Calculate the occupied count by weighting surface boxes by 1/2.
   double GetOccupiedCount(const map<vector<int>, unsigned, CSort<vector<int> > >& inBoxes) const {
      //loop over the occupied boxes and test to see if the neighbors are also occupied
      unsigned long ulOccupiedCount = 0;
      vector<int> vecCheck(3);
      bool bNeighborEmpty;
      map<vector<int>, unsigned, CSort<vector<int> > >::const_iterator iterOC = inBoxes.begin();
      map<vector<int>, unsigned, CSort<vector<int> > >::const_iterator iterOC_end = inBoxes.end();
      for(; iterOC != iterOC_end; ++iterOC) {
         //check all neighbors
         //note: its best to just check them all because it minimizes the "if"'s
         bNeighborEmpty = false;
#ifdef TWO_DIMENSIONAL
         for(int i = 0; i < 2; ++i) {
#else
         for(int i = 0; i < 3; ++i) {
#endif
            vecCheck[0] = iterOC->first[0];
            vecCheck[1] = iterOC->first[1];
            vecCheck[2] = iterOC->first[2];
            for(int j = -1; j < 2; j += 2) {
               vecCheck[i] = iterOC->first[i] + j;
               bNeighborEmpty |= (inBoxes.find(vecCheck) == iterOC_end);
            }
         }

         //if all where occupied add 2 else add 1 (will divide by 2 later)
         ulOccupiedCount += (bNeighborEmpty?1:2);
      }

      return 0.5*double(ulOccupiedCount);
   }

   //count the parenthesis in a string: "(" - ")"
   int ParenthesisCount(const string& inStr) const {
      int iPLevel = 0;
      string::const_iterator iter = inStr.begin();
      string::const_iterator iter_end = inStr.end();
      for(; iter != iter_end; ++iter) {
         if(*iter == '(') {iPLevel += 1;}
         if(*iter == ')') {iPLevel -= 1;}
      }
      return iPLevel;
   }

   //return the determinant of the input matrix (specialized for 3x3)
   double Determinant(const vector<vector<double> >& inM) const {
      return inM[0][0]*inM[1][1]*inM[2][2] + 
         inM[0][1]*inM[1][2]*inM[2][0] +
         inM[0][2]*inM[1][0]*inM[2][1] -
         inM[2][0]*inM[1][1]*inM[0][2] -
         inM[2][1]*inM[1][2]*inM[0][0] -
         inM[2][2]*inM[1][0]*inM[0][1];
   }

   //Return the center of mass of all of the segments stored in inSegments.
   //Each segment is approximated as a point mass at the center of the segment.
   //Uses mass of the segment ~ length.
   //TODO: this may be improved by giving the mass some dependence on level.
   //Modifies inX, inY, inZ to be the coordinates of the COM.
   //Modifies inPoints to contain the coordinates of the centers of the lines and their masses.
   //        Format is inPoint[i][0] = x, [1] = y, [2] = z, [3] = m, [4] = depth.
   void TotalMean(const vector<vector<double> >& inSegments, double& inX, double& inY, double& inZ, vector<vector<double> >& inPoints) const {
      inX = 0.;
      inY = 0.;
      inZ = 0.;
      double dTotalMass = 0.;
      double dMass, dX, dY, dZ;

      vector<vector<double> >::const_iterator iterS = inSegments.begin();
      vector<vector<double> >::const_iterator iterS_end = inSegments.end();
      for(; iterS != iterS_end; ++iterS) {
         dX = iterS->at(0)-iterS->at(1);
         dY = iterS->at(2)-iterS->at(3);
         dZ = iterS->at(4)-iterS->at(5);
         dMass = sqrt( dX*dX + dY*dY + dZ*dZ );

         inX += dMass * 0.5*(iterS->at(0)+iterS->at(1));
         inY += dMass * 0.5*(iterS->at(2)+iterS->at(3));
         inZ += dMass * 0.5*(iterS->at(4)+iterS->at(5));

         dTotalMass += dMass;

         vector<double> vecInsert(5);
         vecInsert[0] = 0.5*(iterS->at(0)+iterS->at(1));
         vecInsert[1] = 0.5*(iterS->at(2)+iterS->at(3));
         vecInsert[2] = 0.5*(iterS->at(4)+iterS->at(5));
         vecInsert[3] = dMass;
         vecInsert[4] = iterS->at(6);
         inPoints.push_back(vecInsert);
      }

      inX /= dTotalMass;
      inY /= dTotalMass;
      inZ /= dTotalMass;
   }

   //Return the center of mass of all of the segments stored in inSegments.
   //Each segment is approximated as a point mass at the center of the segment.
   //Uses mass of the segment ~ length.
   //TODO: this may be improved by giving the mass some dependence on level.
   //Modifies inX, inY, inZ to be the coordinates of the COM.
   //Modifies inPoints to contain the coordinates of the centers of the lines and their masses split up into inSplits equal parts.
   //        Format is inPoint[i][0] = x, [1] = y, [2] = z, [3] = m, [4] = depth.
   void TotalMeanSplit(const vector<vector<double> >& inSegments, const int& inSplits, double& inX, double& inY, double& inZ, vector<vector<double> >& inPoints) const {
      inX = 0.;
      inY = 0.;
      inZ = 0.;
      double dTotalMass = 0.;
      double dMass, dX, dY, dZ;

      vector<vector<double> >::const_iterator iterS = inSegments.begin();
      vector<vector<double> >::const_iterator iterS_end = inSegments.end();
      for(; iterS != iterS_end; ++iterS) {
         dX = iterS->at(0)-iterS->at(1);
         dY = iterS->at(2)-iterS->at(3);
         dZ = iterS->at(4)-iterS->at(5);
         dMass = sqrt( dX*dX + dY*dY + dZ*dZ );

         inX += dMass * 0.5*(iterS->at(0)+iterS->at(1));
         inY += dMass * 0.5*(iterS->at(2)+iterS->at(3));
         inZ += dMass * 0.5*(iterS->at(4)+iterS->at(5));

         dTotalMass += dMass;

         //split up the segments and add the centers (bad variable names).
         //Slope is chosen to pick out centers as c = m*(2*i+1) + b
         dMass /= double(inSplits);
         double mx = 1./2./double(inSplits)*(iterS->at(1)-iterS->at(0));
         double bx = iterS->at(0);
         double my = 1./2./double(inSplits)*(iterS->at(3)-iterS->at(2));
         double by = iterS->at(2);
         double mz = 1./2./double(inSplits)*(iterS->at(5)-iterS->at(4));
         double bz = iterS->at(4);
         for(int i = 0; i < inSplits; ++i) {
            vector<double> vecInsert(5);
            vecInsert[0] = mx*double(2*i+1)+bx;
            vecInsert[1] = my*double(2*i+1)+by;
            vecInsert[2] = mz*double(2*i+1)+bz;
            vecInsert[3] = dMass;
            vecInsert[4] = iterS->at(6);
            inPoints.push_back(vecInsert);
         }
      }

      inX /= dTotalMass;
      inY /= dTotalMass;
      inZ /= dTotalMass;
   }

   //Round the input double to the nearest integer.
   int Round(const double& x) const {
      if(x > 0.) {
         return int(x+0.5);
      } else if (x < 0.) {
         return int(x-0.5);
      } else {
         return 0;
      }
      return 0; //should not get here
   }

   //Rotate then shift the input segments (local format is x1, x2, y1, y2, z1, z2).
   //Uses Euler angles.
   //Modifies inSegments.
   void RotateShift(const double& inPhi, const double& inTheta, const double& inPsi,
      const double& inDX, const double& inDY, const double& inDZ,
      vector<vector<double> > & inSegments) const {

         double x1, x2, y1, y2, z1, z2;
         double a11 = cos(inPsi)*cos(inPhi)-cos(inTheta)*sin(inPhi)*sin(inPsi);
         double a12 = cos(inPsi)*sin(inPhi)+cos(inTheta)*cos(inPhi)*sin(inPsi);
         double a13 = sin(inPsi)*sin(inTheta);
         double a21 = -1.*sin(inPsi)*cos(inPhi)-cos(inTheta)*sin(inPhi)*cos(inPsi);
         double a22 = -1.*sin(inPsi)*sin(inPhi)+cos(inTheta)*cos(inPhi)*cos(inPsi);
         double a23 = cos(inPsi)*sin(inTheta);
         double a31 = sin(inTheta)*sin(inPhi);
         double a32 = -1.*sin(inTheta)*cos(inPhi);
         double a33 = cos(inTheta);

         //Do the rotation and shift.
         vector<vector<double> >::iterator iterS = inSegments.begin();
         vector<vector<double> >::iterator iterS_end = inSegments.end();
         for(; iterS != iterS_end; ++iterS) {
            x1 = iterS->at(0);
            x2 = iterS->at(1);
            y1 = iterS->at(2);
            y2 = iterS->at(3);
            z1 = iterS->at(4);
            z2 = iterS->at(5);
            iterS->at(0) = a11*x1 + a12*y1 + a13*z1 + inDX;
            iterS->at(1) = a11*x2 + a12*y2 + a13*z2 + inDX;
            iterS->at(2) = a21*x1 + a22*y1 + a23*z1 + inDY;
            iterS->at(3) = a21*x2 + a22*y2 + a23*z2 + inDY;
            iterS->at(4) = a31*x1 + a32*y1 + a33*z1 + inDZ;
            iterS->at(5) = a31*x2 + a32*y2 + a33*z2 + inDZ;
         }
   }

   //Shift then rotate the input segments (local format is x1, x2, y1, y2, z1, z2).
   //Uses Euler angles.
   //Modifies inSegments.
   void ShiftRotate(const double& inPhi, const double& inTheta, const double& inPsi,
      const double& inDX, const double& inDY, const double& inDZ,
      vector<vector<double> > & inSegments) const {

         double x1, x2, y1, y2, z1, z2;
         double a11 = cos(inPsi)*cos(inPhi)-cos(inTheta)*sin(inPhi)*sin(inPsi);
         double a12 = cos(inPsi)*sin(inPhi)+cos(inTheta)*cos(inPhi)*sin(inPsi);
         double a13 = sin(inPsi)*sin(inTheta);
         double a21 = -1.*sin(inPsi)*cos(inPhi)-cos(inTheta)*sin(inPhi)*cos(inPsi);
         double a22 = -1.*sin(inPsi)*sin(inPhi)+cos(inTheta)*cos(inPhi)*cos(inPsi);
         double a23 = cos(inPsi)*sin(inTheta);
         double a31 = sin(inTheta)*sin(inPhi);
         double a32 = -1.*sin(inTheta)*cos(inPhi);
         double a33 = cos(inTheta);

         //Do the shift.
         vector<vector<double> >::iterator iterS = inSegments.begin();
         vector<vector<double> >::iterator iterS_end = inSegments.end();
         for(; iterS != iterS_end; ++iterS) {
            iterS->at(0) += inDX;
            iterS->at(1) += inDX;
            iterS->at(2) += inDY;
            iterS->at(3) += inDY;
            iterS->at(4) += inDZ;
            iterS->at(5) += inDZ;
         }
         //do the rotation
         iterS = inSegments.begin();
         iterS_end = inSegments.end();
         for(; iterS != iterS_end; ++iterS) {
            x1 = iterS->at(0);
            x2 = iterS->at(1);
            y1 = iterS->at(2);
            y2 = iterS->at(3);
            z1 = iterS->at(4);
            z2 = iterS->at(5);
            iterS->at(0) = a11*x1 + a12*y1 + a13*z1;
            iterS->at(1) = a11*x2 + a12*y2 + a13*z2;
            iterS->at(2) = a21*x1 + a22*y1 + a23*z1;
            iterS->at(3) = a21*x2 + a22*y2 + a23*z2;
            iterS->at(4) = a31*x1 + a32*y1 + a33*z1;
            iterS->at(5) = a31*x2 + a32*y2 + a33*z2;
         }
   }

   //Shift then rotate the input segments (local format is x, y, z, (other stuff)).
   //Uses Euler angles.
   //Modifies inSegments.
   void ShiftRotatePoints(const double& inPhi, const double& inTheta, const double& inPsi,
      const double& inDX, const double& inDY, const double& inDZ,
      vector<vector<double> > & inPoints) const {

         double x, y, z;
         double a11 = cos(inPsi)*cos(inPhi)-cos(inTheta)*sin(inPhi)*sin(inPsi);
         double a12 = cos(inPsi)*sin(inPhi)+cos(inTheta)*cos(inPhi)*sin(inPsi);
         double a13 = sin(inPsi)*sin(inTheta);
         double a21 = -1.*sin(inPsi)*cos(inPhi)-cos(inTheta)*sin(inPhi)*cos(inPsi);
         double a22 = -1.*sin(inPsi)*sin(inPhi)+cos(inTheta)*cos(inPhi)*cos(inPsi);
         double a23 = cos(inPsi)*sin(inTheta);
         double a31 = sin(inTheta)*sin(inPhi);
         double a32 = -1.*sin(inTheta)*cos(inPhi);
         double a33 = cos(inTheta);

         //Do the shift (if it's there).
         vector<vector<double> >::iterator iterP = inPoints.begin();
         vector<vector<double> >::iterator iterP_end = inPoints.end();
         if( fabs(inDX) > 1.e-12 || fabs(inDY) > 1.e-12 || fabs(inDZ) > 1.e-12 ) {
            for(; iterP != iterP_end; ++iterP) {
               iterP->at(0) += inDX;
               iterP->at(1) += inDY;
               iterP->at(2) += inDZ;
            }
         }
         
         //do the rotation
         iterP = inPoints.begin();
         iterP_end = inPoints.end();
         for(; iterP != iterP_end; ++iterP) {
            x = iterP->at(0);
            y = iterP->at(1);
            z = iterP->at(2);
            iterP->at(0) = a11*x + a12*y + a13*z;
            iterP->at(1) = a21*x + a22*y + a23*z;
            iterP->at(2) = a31*x + a32*y + a33*z;
         }
   }

   //Test that the brute force 2d rotation is the same as the inertia one.
   //returns the angle for the two d rotation that maximizes <x^2>.
   double Test2DRotation() const {
      //get all the segments
      vector<vector<double> > vecSegments;
      for(unsigned i = 0; i < m_vbrNormalizedNeuron.size(); ++i) {
         m_vbrNormalizedNeuron[i]->AllSegments(vecSegments);
      }

      //Create some things needed for the gsl integrator.
      gsl_integration_workspace* w =
         gsl_integration_workspace_alloc(1000);  //1000 is plenty of room

      double dResult, dError;
      double dAlpha[8];

      //Easiest just to have both functions available, this takes no time.
      gsl_function LF_2D;
      LF_2D.params = &dAlpha;
      LF_2D.function = &LineFunction2D_X;
      gsl_function LF_3D;
      LF_3D.params = &dAlpha;
      LF_3D.function = &LineFunction3D_X;

      //Do the rotations and find the moments (in the x-direction, WLOG).
      int iMax = 0;
      double dMax = -1.;
      double dX, dM;
      int iNumTrials = 100;
      double dRotateAngle = 2.*acos(-1.)/double(iNumTrials);
      for(int j = 0; j < iNumTrials; ++j) {
         ShiftRotate(dRotateAngle, 0., 0., 0., 0., 0., vecSegments);
         double dSecondMoment = 0.;
         vector<vector<double> >::const_iterator iter = vecSegments.begin();
         vector<vector<double> >::const_iterator iter_end = vecSegments.end();
         for(; iter != iter_end; ++iter) {
            dAlpha[0] = 2.0;
            for(int j = 0; j < 7; ++j) {
               dAlpha[j+1] = iter->at(j);
            }
#ifdef TWO_DIMENSIONAL
            gsl_integration_qag(&LF_2D, 0, 1, 0, 1.e-7, 1000, GSL_INTEG_GAUSS61, 
               w, &dResult, &dError);
#else
            gsl_integration_qag(&LF_3D, 0, 1, 0, 1.e-7, 1000, GSL_INTEG_GAUSS61, 
               w, &dResult, &dError);
#endif
            dSecondMoment += SegmentLength(*iter)*dResult;
         }
         //update max and min
         if(dSecondMoment > dMax) {
            iMax = j;
            dMax = dSecondMoment;
         }
      }

      //clean up the workspace
      gsl_integration_workspace_free(w);

      return dRotateAngle*double(iMax);
   }

private:
   //data members
   vector<shared_ptr<CBranch> > m_vbrNeuron;
   vector<shared_ptr<CBranch> > m_vbrNormalizedNeuron;

   bool m_bFirstPass; //used by Normalize, no public access

   double m_dWidthX;
   double m_dWidthY;
   double m_dWidthZ;
   double m_dStandardDeviationX;
   double m_dStandardDeviationY;
   double m_dStandardDeviationZ;

};

#endif //NEURON
