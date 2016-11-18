//Joe Snider
//7/04
//
//Store the branches of a neuron as a tree with points defining the branch.
//
//Pointer based access to parents and daughters, be sure not to delete them.
//TODO: create smart pointers (overload new and delete) that don't allow deletion while still in the tree.
//DONE: pointers are smart now.

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include "point.h"
#include "line.h"
#include "3dvector.h"
#include "histogram.h"
#include "rand3.h"

using namespace std;
using namespace boost;

#ifndef BRANCH
#define BRANCH

const double LINES_CLOSE = .001; //Should be about twice the size of a spine in local units.
const double POINT_NEAR_DISTANCE = 1.e-6; //should get rid of this

const double MIN_SIGNIFICANT_ANGLE_LENGTH = 2.;  //set the minimum length of a branch to give reliable angles
const int ANGLE_POINTS = 5; //set the number of points away from a vertex to measure the angle
                            //
                            //   --------     <= could use at most 2
                            //       \
                            //        \

class CBranch: public enable_shared_from_this<CBranch> {
public:
   //typedefs...
   typedef shared_ptr<CBranch> BranchPointerType;
   typedef vector<BranchPointerType> DaughtersHolderType;

   typedef CPoint PointType;
   typedef CLine LineType;

public:
   //constructors

   //Warning this does not add the current branch to the parent.
   //  I can't figure out how to make a shared_from_this() before the constructor is finished.
   //  Be sure to set this as a daughter elsewhere.
   CBranch(shared_ptr<CBranch> inParent) {
      //initial values for the max and min
      m_dMinX = 1.e10;
      m_dMaxX = -1.e10;
      m_dMinY = 1.e10;
      m_dMaxY = -1.e10;
      m_dMinZ = 1.e10;
      m_dMaxZ = -1.e10;

      m_Parent = inParent;
      m_iLevel = inParent->GetLevel() + 1;
      //inParent->AddDaughter(shared_from_this());

      UpdateMaxMin(m_dMinX, m_dMaxX, m_dMinY, m_dMaxY, m_dMinZ, m_dMaxZ);
   }

   //Construct a root or null branch.
   //inLevel of 0 indicates a root branch.
   //inLevel < 0 is a null branch.
   CBranch(const int& inLevel) {
      //initial values for the max and min
      m_dMinX = 1.e10;
      m_dMaxX = -1.e10;
      m_dMinY = 1.e10;
      m_dMaxY = -1.e10;
      m_dMinZ = 1.e10;
      m_dMaxZ = -1.e10;

      m_Parent = shared_ptr<CBranch>();
      m_iLevel = inLevel;
   }

   //This is somewhat dangerous.
   //Branches should never be created stand alone; not sure how to enforce this.
   //If this is only a trim, the parent should probably call UpdateMaxMin.
   ~CBranch() {
      ////Delete the daughters.
      //DaughtersHolderType::iterator iter = m_Daughters.begin();
      //DaughtersHolderType::iterator iter_end = m_Daughters.end();
      //for(; iter != iter_end; ++iter) {
      //   delete *iter;
      //}
      //Clear the vector of daughters.
      m_Daughters.clear();
   }

   //Warn about copy construction.
   CBranch(const CBranch& inCopy) {
      cerr << "CBranch.h error: copy constructor called...returning\n" << flush;
   }

public:
   //gets and sets
   //minima and maxima can not be changed except through AddPoint
   //The minima and maxima include all daughters of the current branch.
   double GetMinX() const {return m_dMinX;}
   double GetMaxX() const {return m_dMaxX;}
   double GetMinY() const {return m_dMinY;}
   double GetMaxY() const {return m_dMaxY;}
   double GetMinZ() const {return m_dMinZ;}
   double GetMaxZ() const {return m_dMaxZ;}

   shared_ptr<CBranch> GetParent() const {return m_Parent;}
   void SetParent(shared_ptr<CBranch> inParent) {
      m_Parent = inParent;
      m_iLevel = inParent->GetLevel() + 1;
      //inParent->AddDaughter(shared_from_this());
   }

   const DaughtersHolderType& GetDaughters() const {return m_Daughters;}
   void SetDaughters(const DaughtersHolderType& inDaughters) {
      m_Daughters.clear();
      DaughtersHolderType::const_iterator iter = inDaughters.begin();
      DaughtersHolderType::const_iterator iter_end = inDaughters.end();
      for(; iter != iter_end; ++iter) {
         m_Daughters.push_back(*iter);
      }
   }

   int GetLevel() const {return m_iLevel;}
   void SetLevel(const int inLevel) {
      m_iLevel = inLevel;
      //recursively set the daughter levels.
      DaughtersHolderType::iterator iter = m_Daughters.begin();
      DaughtersHolderType::iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         (*iter)->SetLevel(inLevel+1);
      }
   }

   const LineType& GetLine() const {return m_Line;}
   //LineType& GetLine() {return m_Line;}

public:
   //interface

   //Merge daughters uphill.
   //Takes branches that would have ended at a daughter and turns them into one branch.
   //This is a memory leak for now. I think the daughter could just be deleted, but it's dangerous.
   //This does not preserve divisions (if any).
   //Pointers are smart now, so no memory leak (Joe 3/07).
   //Be sure to update the levels elsewhere.
   void TrimDaughter() {
      DaughtersHolderType::iterator iter = m_Daughters.begin();
      DaughtersHolderType::iterator iter_end = m_Daughters.end();
      if(m_Line.size() > 0) {
         while(iter != iter_end) {
            if( (m_Line.back() - (*iter)->GetLine().front()).LengthSquared() < 1.e-12 ) {
               //Insert the new points
               deque<CPoint>::const_iterator iterL = (*iter)->GetLine().GetTurningPoints().begin();
               deque<CPoint>::const_iterator iterL_end = (*iter)->GetLine().GetTurningPoints().end();
               ++iterL;
               m_Line.AddTurningPoint(iterL, iterL_end);

               //copy the daughter's daughters to this (next step invalidates iterators)
               DaughtersHolderType tempDaughters;
               tempDaughters.assign( (*iter)->m_Daughters.begin(), (*iter)->m_Daughters.end() );

               //remove the daughter before invalidation
               m_Daughters.erase(iter);

               //append the grand daughters
               m_Daughters.insert(m_Daughters.end(), tempDaughters.begin(), tempDaughters.end());

               //Set the parent of the daughters
               UpdateDaughters();

               //restart the search
               iter = m_Daughters.begin();
               iter_end = m_Daughters.end();
            } else {
               //check the next daughter
               ++iter;
            }
         }
      }

      //recurse after everything is merged because of the restarts and
      // in case the line has length zero.
      //  (would be slow to put it in the first loop)
      iter = m_Daughters.begin();
      iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         (*iter)->TrimDaughter();
      }
   }

   //The number of daughters
   int DaughterCount() const {
      return m_Daughters.size();
   }

   //Trim to the first daughter.
   //Keeps a little bit of the tail to prevent problems with later algorithms.
   //  Ensures that the line has at least 2 points.
   //Does nothing if there is no daughter.
   void TrimToFirstDaughter() {
      //Make sure there is a daughter (do nothing otherwise).
      if(m_Daughters.size() > 0) {
         //Find the daughter closest to the start of the branch.
         map<double, int> mapDistance_ID;
         for(int i = 0; i < m_Daughters.size(); ++i) {
            mapDistance_ID[m_Line.DistanceToStart(m_Daughters[i]->GetLine().at(0))] = i;
         }

         //Get the first point after the daughter (should actually be the daughter).
         int iTrimTo = 0;
         for(; iTrimTo < m_Line.size() && m_Line.DistanceToStart(m_Line[iTrimTo]) < mapDistance_ID.begin()->first-1.e-12; ++iTrimTo);

         //Do the trimming (some errors first)
         if(iTrimTo == 0) {
            //nothing to do (not an errror, just a daughter close to the start)
            return;
         } else if(iTrimTo == m_Line.size()) {
            //Error
            cerr << "Error in CBranch::TrimToFirstDaughter: Unable to find trim point .. returning\n" << flush;
            return;
         } else {
            //The point's ok, so decrement it one (make sure there are 2 points) and trim.
            iTrimTo--;
            m_Line.PopFront( iTrimTo );
         }
      }
   }

   //Recursively find the maximium level.
   //Uses a depth first search that could be optimized (hopefully I won't have to bother).
   //inLevel is modified to contain the max level.
   void MaxLevel(int& inLevel) const {
      inLevel = max(inLevel, GetLevel());

      //recurse over the daughters
      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         (*iterD)->MaxLevel(inLevel);
      }
   }

   //Get the first point on the line.
   //Returns an empty point if the line is too short.
   CPoint FirstPoint() const {
      if(m_Line.size() > 0) {
         return m_Line[0];
      }
      CPoint pRet;
      return pRet;
   }

   //Get the last point on the line (makes insertion easier)
   //Returns an empty point if the line is too short.
   CPoint LastPoint() const {
      if(m_Line.size() > 0) {
         return m_Line[m_Line.size()-1];
      }
      CPoint pRet;
      return pRet;
   }

   //Recursively shift by the input ammounts.
   //Just adds the input values to the x,y,z coords.
   void Shift(const double& inX, const double& inY, const double& inZ) {
      //if(m_Line.size() < 1) {return;}
      PointType pShift;
      pShift.push_back(inX);
      pShift.push_back(inY);
      pShift.push_back(inZ);
      DoShift(pShift);

      UpdateAllMaxMin();
   }

   //Recursively scale by the input ammounts in each direction.
   //Just adds the input values to the x,y,z coords.
   void Scale(const double& inX, const double& inY, const double& inZ) {
      PointType pScale;
      pScale.push_back(inX);
      pScale.push_back(inY);
      pScale.push_back(inZ);
      DoScale(pScale);

      UpdateAllMaxMin();
   }

   //Recursively rotate by the input ammounts.
   //Just just apply the input matrix to the points.
   //Warning: bad name; there is no check that inM is a rotation.
   //Type T must have operator [][] access to a double (it's a matrix).
   template<class T>
   void Rotate(const T& inM) {
      DoRotate(inM);
      UpdateAllMaxMin();
   }

   //Count the intersections of the input segment with this branch and all sub branches.
   //Uses the tree structure to make the algorithm O(NlogN) (N is the number of segments in the branch).
   //inSegment stores the endpoints of a segment in format x1, x2, y1, y2, z1, z2.
   //Returns the number of intersections.
   int CountSegmentIntersections(const vector<double>& inSegment) const {
      double x1 = inSegment[0];
      double x2 = inSegment[1];
      double y1 = inSegment[2];
      double y2 = inSegment[3];
      double z1 = inSegment[4];
      double z2 = inSegment[5];
      int iIntersections = 0;
      DoCountSegmentIntersections(x1, x2, y1, y2, z1, z2, iIntersections);
      return iIntersections;
   }

   //Find the minimum distance to the input point
   double MinimumDistance(const double& inX, const double& inY, const double& inZ) const {
      //double dMinimum = 10.*pow(max(GetMaxX()-GetMinX(), max(GetMaxY()-GetMinY(), GetMaxZ()-GetMinZ())), 2);
      double dMinimum = 1.e200;
      DoMinimumDistance(inX, inY, inZ, dMinimum);
      return sqrt(dMinimum);
   }

   //Find the max distance to the input point (a boundary)
   double BoundingRadius(const double& inX, const double& inY, const double& inZ) const {
      double dReturn = -1.;
      DoBoundingRadius(dReturn, inX, inY, inZ);
      return sqrt(dReturn);
   }

   //Recursively calculate the radial distribution of length.
   //Input histogram is modified to store the length distribution.
   //Note that the normalization is such that the sum of the histogram is 
   //   equal to the total length.
   //Breaks up the line of this branch (and all daughters) in inMaxSize pieces
   //   and stores the lengths at the appropriate distance.
   //Note: this is the same regardless of dimension.
   void LengthDistribution(CHistogram<double>& inHisto, const double& inMaxSize = 0.1) {
      //Do nothing if there aren't enough turning points (except recurse).
      if( m_Line.size() > 1 ) {
         //Keep the original divisions to reset
         unsigned uOriginalDivisions = m_Line.size();

         //Divide into the desired sized segments
         unsigned uDivisions = unsigned(m_Line.Length()/inMaxSize+1.);
         m_Line.Divide(uDivisions);

         //Store these lengths
         LineType::const_iterator iterL = m_Line.begin();
         LineType::const_iterator iterR = m_Line.begin();
         ++iterR;
         LineType::const_iterator iter_end = m_Line.end();
         for(; iterR != iter_end; ++iterL, ++iterR) {
            inHisto.increment(0.5*(*iterL + *iterR).Length(), (*iterR - *iterL).Length());
         }

         //Undo the division
         m_Line.Divide(uOriginalDivisions);
      } else {
         cout << "Warning: Line with 0 or 1 points in CBranch...Continuing\n" << flush;
      }

      //Do the recursion
      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         (*iterD)->LengthDistribution(inHisto, inMaxSize);
      }
   }

   //Recursively find the angle formed between this branch and all daughters.
   //Dumps angles to the input stream.
   //Changed to measure the angle locally (3/07 Joe Snider).
   void RawAngles(ostream& inOut) const {
      DaughtersHolderType::const_iterator iter = m_Daughters.begin();
      DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         if( GetLine().size() > ANGLE_POINTS) {
            CPoint pntCurrent;
            for(int i = 0; i < GetLine().size(); ++i) {
               if( (GetLine().at(i) - (*iter)->GetLine().at(0)).LengthSquared() < 1.e-12) {
                  if(i+ANGLE_POINTS < GetLine().size()) {
                     pntCurrent = GetLine().at(i+ANGLE_POINTS) - GetLine().at(i);
                  } else {
                     pntCurrent = GetLine().at(i)-GetLine().at(i);
                  }
               }
            }
            double dCurrentLength = pntCurrent.Length();
            if(dCurrentLength > MIN_SIGNIFICANT_ANGLE_LENGTH && (*iter)->GetLine().size() > ANGLE_POINTS) {
               CPoint pntDaughter = (*iter)->GetLine().at(ANGLE_POINTS) - (*iter)->GetLine().at(0);
               double dDaughterLength = pntDaughter.Length();
               if(dDaughterLength > 1.e-12) {
                  inOut << acos(pntCurrent.Dot(pntDaughter)/(dCurrentLength*dDaughterLength)) << "\n";
               }
            }
         }
         (*iter)->RawAngles(inOut);
      }
   }

   //Recursively find the total direction wrt the radial.
   //Branches are treated as straight lines between digitized points.
   //Added 5/10 (Joe Snider)
   void GlobalAngle(CHistogram<CPoint>& inHisto) const {
      if(GetLine().size() > 1) {
         //find the range
         LineType::const_iterator iterL = m_Line.begin();
         LineType::const_iterator iterR = m_Line.begin();
         ++iterR;
         LineType::const_iterator iter_end = m_Line.end();
         for(; iterR != iter_end; ++iterL, ++iterR) {
            CPoint v = *iterR - *iterL;
            CPoint rHat = *iterL;
            if(v.LengthSquared() > 1.e-6 && rHat.LengthSquared() > 1.e-6) {
               //calculate the direction vector
               //v *= 1./v.Length();
               rHat *= 1./rHat.Length();
               
               if(v.GetDimension() == 2) {
                   //promote to 3d
                   v.push_back(0);
                   rHat.push_back(0);
               }

               //Build the rotation matrix to put rHat in the zHat direction
               //Could simplify the sin(atan(x)) stuff, but not worth it
               //rotate about z
               //TODO: may actually be worth it since atan is blowing up sometimes
               const double phi = atan2(rHat[1], rHat[0]);
               double R1[3][3];
               R1[0][0] = cos(phi); R1[0][1] = sin(phi); R1[0][2] = 0.;
               R1[1][0] = -sin(phi); R1[1][1] = cos(phi); R1[1][2] = 0.;
               R1[2][0] = 0.; R1[2][1] = 0.; R1[2][2] = 1.;

               //rotate about y
               const double theta = acos(rHat[2]);
               double R2[3][3];
               R2[0][0] = cos(theta); R2[0][1] = 0.; R2[0][2] = -sin(theta);
               R2[1][0] = 0.; R2[1][1] = 1.; R2[1][2] = 0.;
               R2[2][0] = sin(theta); R2[2][1] = 0.; R2[1][2] = cos(theta);

               //Apply the rotation to the branch vector
               v.LinearTransform(R1);
               v.LinearTransform(R2);

               //add an extra 1 at the end to count the vectors
               v.push_back(1.);

               //store into the appropriate bins
               CHistogram<CPoint>::iterator low = inHisto.find(iterL->Length());
               CHistogram<CPoint>::iterator high = inHisto.find(iterR->Length());
               *high += v;
               if(high == inHisto.overflow()) {
                  high = inHisto.bin_end();
               }
               for(; low < high; ++low) {
                  *low += v;
               }
            }
         }
      }

      //Do the recursion
      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         (*iterD)->GlobalAngle(inHisto);
      }
   }

   //Recursively count each branch in a range
   //Added 5/10 (Joe Snider)
   void BranchCount(CHistogram<double>& inHisto) const {
      //Do nothing if there aren't enough turning points (except recurse).
      if( m_Line.size() > 1 ) {
         //find the range
         LineType::const_iterator iterL = m_Line.begin();
         LineType::const_iterator iterL_end = m_Line.end();
         double minR2 = iterL->LengthSquared();
         double maxR2 = iterL->LengthSquared();
         ++iterL;
         for(; iterL != iterL_end; ++iterL) {
            const double l2 = iterL->LengthSquared();
            minR2 = (minR2<l2)?minR2:l2;
            maxR2 = (maxR2>l2)?maxR2:l2;
         }
         //increment the appropriate parts of the histogram
         CHistogram<double>::iterator low = inHisto.find(sqrt(minR2));
         CHistogram<double>::iterator high = inHisto.find(sqrt(maxR2));
         *high += 1;
         if(high == inHisto.overflow()) {
            high = inHisto.bin_end();
         }
         for(; low < high; ++low) {
            *low += 1;
         }
      }

      //Do the recursion
      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         (*iterD)->BranchCount(inHisto);
      }
   }

   //Recursively find the lengths of this branch and all daughters.
   //Dump them to the input stream.
   void RawLengths(ostream& inOut) const {
      CPoint pntStart = FirstPoint();

      DaughtersHolderType::const_iterator iter = m_Daughters.begin();
      DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         if(!IsRoot()) {
            inOut << m_Line.DistanceToStart((*iter)->FirstPoint()) - m_Line.DistanceToStart(pntStart) << "\n";
            pntStart = (*iter)->FirstPoint();
         }
         (*iter)->RawLengths(inOut);
      }

      //include the end
      if(DaughterCount() == 0) {
         //this is a leaf just add the length of this branch
         inOut << m_Line.Length() << "\n";
      } else if (!IsRoot()) {
         //add the last daughter to the end
         inOut << m_Line.DistanceToStart(LastPoint()) - m_Line.DistanceToStart(m_Daughters.back()->FirstPoint()) << "\n";
      }
   }

   //Recursively find the lengths of this branch and all daughters.
   //Store them into the input (T must have push_back, std::vector works)
   template<class T>
   void BranchLength(T& inData) const {
      CPoint pntStart = FirstPoint();

      DaughtersHolderType::const_iterator iter = m_Daughters.begin();
      DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         if(!IsRoot()) {
            inData.push_back(m_Line.DistanceToStart((*iter)->FirstPoint()) - m_Line.DistanceToStart(pntStart));
            pntStart = (*iter)->FirstPoint();
         }
         (*iter)->BranchLength(inData);
      }

      //include the end
      if(DaughterCount() == 0) {
         //this is a leaf just add the length of this branch
         inData.push_back(m_Line.Length());
      } else if (!IsRoot()) {
         //add the last daughter to the end
         inData.push_back(m_Line.DistanceToStart(LastPoint()) - m_Line.DistanceToStart(m_Daughters.back()->FirstPoint()));
      }
   }

   //recursively calculate the average displacement squared of the leaves
   //the type T must have push_back defined (std::vector works)
   template<class T>
   void AverageDisplacementSquared(const double& inX, const double& inY, const double& inZ, T& inR2) const {
      //if this is a leaf find the displacement squared
      if(IsLeaf()) {
         if(m_Line.size() > 0) {
            inR2.push_back( DistanceSquared(m_Line.back()[0], m_Line.back()[1], m_Line.back()[2], inX, inY, inZ) );
         }
      } else {
         //call the daughters
         //recursively get the daughter line lengths
         DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
         DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
         for(; iterD != iterD_end; ++iterD) {
            (*iterD)->AverageDisplacementSquared(inX, inY, inZ, inR2);
         }
      }
   }

   //Recursively calculate the mean (center).
   void TotalMean(double& inX, double& inY, double& inZ) const {
      unsigned uN = 0;
      inX = 0.;
      inY = 0.;
      inZ = 0.;
      //call the helper that does the recursion
      DoTotalMean(uN, inX, inY, inZ);
      inX /= double(uN);
      inY /= double(uN);
      inZ /= double(uN);
   }

   //Recursively calculate the total variance.
   void TotalVariance(double& inV) const {
      double dMX = 0.;
      double dMY = 0.;
      double dMZ = 0.;
      TotalMean(dMX, dMY, dMZ);

      //call the helper that does the recursion
      unsigned uN = 0;
      double dVX = 0.;
      double dVY = 0.;
      double dVZ = 0.;
      DoTotalVariance(uN, dVX, dVY, dVZ, dMX, dMY, dMZ);
      inV = sqrt(dVX*dVY*dVZ/double(uN*uN*uN));
   }

   //Return the total length of this and all daughters.
   double TotalLength() const {
      double dTotalLength = 0.;
      TotalLength(dTotalLength);
      return dTotalLength;
   }

   //Recursively calculate the total length of this branch plus all daughters
   void TotalLength(double& inL) const {
      //add up this neuron
      inL += m_Line.Length();

      //recurse over the daughters
      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         (*iterD)->TotalLength(inL);
      }
   }

   //Modify the input vector to contain all the line segments of this branch and all daughters.
   //No checking for duplicates.
   //Format is inSegment[i][0] = x_left, inSegment[i][1] = x_right, ...
   //          inSegment[i][6] = depth
   //Returns the number of segments.
   int AllSegments(vector<vector<double> >& inSegments) const {
      //store the segments for this branch
      if(m_Line.size() > 1) {
         LineType::const_iterator iter_left = m_Line.begin();
         LineType::const_iterator iter_right = m_Line.begin();
         ++iter_right;
         ////only add the first point in the parent (it's duplicated in the daughters)
         //if(m_iLevel == 0) {
         //   vector<double> vecInsert;
         //   vecInsert.push_back(iter_left->operator[](0));
         //   vecInsert.push_back(iter_right->operator[](0));
         //   vecInsert.push_back(iter_left->operator[](1));
         //   vecInsert.push_back(iter_right->operator[](1));
         //   vecInsert.push_back(iter_left->operator[](2));
         //   vecInsert.push_back(iter_right->operator[](2));
         //   inSegments.push_back(vecInsert);
         //}
         //++iter_left;
         //++iter_right;
         LineType::const_iterator iter_end = m_Line.end();
         for(; iter_right != iter_end; ++iter_left, ++iter_right) {
            vector<double> vecInsert;
            vecInsert.push_back(iter_left->at(0));
            vecInsert.push_back(iter_right->at(0));
            vecInsert.push_back(iter_left->at(1));
            vecInsert.push_back(iter_right->at(1));
            vecInsert.push_back(iter_left->at(2));
            vecInsert.push_back(iter_right->at(2));
            vecInsert.push_back(GetLevel());
            inSegments.push_back(vecInsert);
         }
      }
      //recurse over the daughters
      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         (*iterD)->AllSegments(inSegments);
      }

      return (int)inSegments.size();
   }

   //Search this branch and all daughters to see if the input box intersects.
   //Uses the max/min to determine if the branch has to be checked.
   //Returns true (and stops searching) if the branch connects with the input box.
   bool BoxIntersect(const double& inX1, const double& inX2,
                     const double& inY1, const double& inY2,
                     const double& inZ1, const double& inZ2) const {
      //check if the box is inside the current max/min (note: max/min covers all daughters)
      double x, y, z, dLX, dLY, dLZ, dRX, dRY, dRZ, DX, DY, DZ, mx, my, mz;
      if(IsInRange(inX1, inX2, inY1, inY2, inZ1, inZ2)) {
         //box is in range so have to check the branch
         if(m_Line.size() > 1) {
            bool bLeftInBox, bRightInBox;
            LineType::const_iterator iterB_left = m_Line.begin();
            LineType::const_iterator iterB_right = m_Line.begin();
            ++iterB_right;
            LineType::const_iterator iterB_end = m_Line.end();
            for(; iterB_right != iterB_end; ++iterB_left, ++iterB_right) {
               //the endpoints of the current segment ( variable name is L = left, R = right, then X,Y,Z)
               dLX = iterB_left->operator[](0);
               dLY = iterB_left->operator[](1);
               dLZ = iterB_left->operator[](2);
               dRX = iterB_right->operator[](0);
               dRY = iterB_right->operator[](1);
               dRZ = iterB_right->operator[](2);

               //check where the endpoints are in relation to the box (in or out)
               bLeftInBox = PointInBox(dLX, dLY, dLZ, inX1, inX2, inY1, inY2, inZ1, inZ2);
               bRightInBox = PointInBox(dLX, dLY, dLZ, inX1, inX2, inY1, inY2, inZ1, inZ2);

               //the (painful) logic to test intersection
               if( bLeftInBox || bRightInBox ) {
                  //at least one endpoint is in the box, thus the box is occupied so return true and stop
                  return true;
               } else {
                  //Both points out of the box, life is harder.
                  //All I can think of here is to check all the faces. If any of them intersect the line,
                  //  then the line intersects the cube (note: both endpoints known to be outside the cube
                  //  at this point.
                  //This could be a subroutine, but it wouldn't be much prettier.
                  //warning: bad but short variable names to follow (this is ugly anyway):
                  //   (mx = slope in x, bx = offset in x, DX = delta x of the endopints, x = test point, ...)
                  DX = dRX - dLX;
                  DY = dRY - dLY;
                  DZ = dRZ - dLZ;
                  
                  //Check that line is not parallel to the face.
                  //(no intersection by assumption, misses those overlapping the face)
                  if(fabs(DZ) > 1.e-10) {
                     //Face 1 (z=zmin)
                     //find the equation of the current line segment projecting out the z-direction
                     //   slopes handle units ( [x]/[z]*z_face and [y]/[z]*z_face ), offset is the left point
                     mx = DX/DZ;
                     my = DY/DZ;
                     //find the point of intersection with the lower face parallel to z (at z=inZ1)
                     x = mx*inZ1 + dLX;
                     y = my*inZ1 + dLY;
                     //now see if the point intersects with the face
                     if ( (x > inX1) && (x < inX2) && (y > inY1) && (y < inY2) ) {
                        return true;
                     }
                     //the same process 5 more times (slope is the same here)
                     //Face 2 (z=inZ2)
                     x = mx*inZ2 + dLX;
                     y = my*inZ2 + dLY;
                     if ( (x > inX1) && (x < inX2) && (y > inY1) && (y < inY2) ) {
                        return true;
                     }
                  }
                  if(fabs(DY) > 1.e-10) {
                     //Face 3 (y=inY1)
                     mx = DX/DY;
                     mz = DZ/DY;
                     x = mx*inY1 + dLX;
                     z = mz*inY1 + dLZ;
                     if ( (x > inX1) && (x < inX2) && (z > inZ1) && (z < inZ2) ) {
                        return true;
                     }
                     //Face 4 (y=inY2)
                     x = mx*inY2 + dLX;
                     z = mz*inY2 + dLZ;
                     if ( (x > inX1) && (x < inX2) && (z > inZ1) && (z < inZ2) ) {
                        return true;
                     }
                  }
                  if(fabs(DX) > 1.e-10) {
                     //Face 5 (x=inX1)
                     my = DY/DX;
                     mz = DZ/DX;
                     y = my*inX1 + dLY;
                     z = mz*inX1 + dLZ;
                     if ( (y > inY1) && (y < inY2) && (z > inZ1) && (z < inZ2) ) {
                        return true;
                     }
                     //Face 6 (x=inX2)
                     y = my*inX2 + dLY;
                     z = mz*inX2 + dLZ;
                     if ( (y > inY1) && (y < inY2) && (z > inZ1) && (z < inZ2) ) {
                        return true;
                     }
                  }
                  //----end of painful code; if you think of something better, have at it.
               }
            }
         }

         //If we get here, then we have to check all the daughter branches.
         //If any daughter returns true, then we stop and return true.
         DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
         DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
         for(; iterD != iterD_end; ++iterD) {
            if( (*iterD)->BoxIntersect(inX1, inX2, inY1, inY2, inZ1, inZ2) ) {
               return true;
            }
         }

         //No daughters intersected, so we got to a non-intersecting box the hard way.
         return false;
      } else {
         //box is out of range, so don't have to check the branch or any daughters
         return false;
      }
      return false; //should not get here
   }

   //Add a daughter branch to this branch.
   void AddDaughter(shared_ptr<CBranch> inDaughter) {
      m_Daughters.push_back(inDaughter);
      inDaughter->SetParent(shared_from_this());
      UpdateMaxMin(inDaughter->GetMinX(), inDaughter->GetMaxX(),
         inDaughter->GetMinY(), inDaughter->GetMaxY(),
         inDaughter->GetMinZ(), inDaughter->GetMaxZ());
   }

   //Recursively sort the branches so that the daughters come in
   //  order of distance to the start.
   //Call before doing any measurements.
   void SortDaughters() {
      //Do the sorting.
      multimap<double, BranchPointerType> mapSorter;
      DaughtersHolderType::iterator iter = m_Daughters.begin();
      DaughtersHolderType::iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         mapSorter.insert(
            multimap<double, BranchPointerType>::
            value_type(m_Line.DistanceToStart((*iter)->FirstPoint()), *iter) );
      }

      m_Daughters.clear();
      multimap<double, BranchPointerType>::const_iterator iterM = mapSorter.begin();
      multimap<double, BranchPointerType>::const_iterator iterM_end = mapSorter.end();
      for(; iterM != iterM_end; ++iterM) {
         m_Daughters.push_back(iterM->second);
      }
   }

   //Test if this branch is a leaf (no daughters).
   bool IsLeaf() const {
      return (m_Daughters.size() == 0);
   }

   //Test if this branch is the root (m_iLevel == 0).
   bool IsRoot() const {
      return (m_iLevel == 0);
   }

   //Use negetive level values to signal a null branch.
   bool IsNull() const {
      return (m_iLevel < 0);
   }

   //Set this branch as a root branch.
   //Just changes the level and clears the parent pointer.
   void SetRoot() {
      shared_ptr<CBranch> pbrTempParent( new CBranch(-1) );
      m_Parent = pbrTempParent;
      SetLevel(0);
   }

   //Check if the input point intersects this branch.
   //Return a pointer to this if it intersects and NULL if not.
   //Can't be const since this returns a pointer that may change the branch (it's level... at least).
   shared_ptr<CBranch> DoesPointIntersect(const CPoint& inPoint, const double& inNearDistance = POINT_NEAR_DISTANCE) {
      shared_ptr<CBranch> brpReturn(new CBranch(-1));
      DoDoesPointIntersect(inPoint, inNearDistance, brpReturn);
      return brpReturn;
   }

   //Add a point to the front of the line
   void AddPointFront(const double& inX, const double& inY, const double& inZ) {
      m_Line.push_front(inX, inY, inZ);
      UpdateMaxMin(inX, inX, inY, inY, inZ, inZ);
   }

   //Add a point
   void AddPoint(const double& inX, const double& inY, const double& inZ) {
      m_Line.push_back(inX, inY, inZ);
      UpdateMaxMin(inX, inX, inY, inY, inZ, inZ);
   }

   //Save the branch in a Neurolucida like format (swc).
   //Good enough to fool my NL reader, but not tested on NL itself.
   //This is not terribly efficient, but shouldn't matter much.
   //  If it's a bottle neck, then make the search for a parent faster.
   void RecursiveDisplaySWC(ostream& inOut, const int& inType) const {
      vector<CPoint> vecPoints;
      //Starts at the root (by definition)
      LineType::const_iterator iterL = GetLine().begin();
      LineType::const_iterator iterL_end = GetLine().end();
      vecPoints.push_back(*iterL);
      inOut << "1 1 " << *iterL << ((iterL->size() == 2)?" 1 1 -1\n":" 1 -1\n");
      ++iterL;
      for(; iterL != iterL_end; ++iterL) {
         vecPoints.push_back(*iterL);
         inOut << vecPoints.size() << " " << inType << " "
			 << *iterL << ((iterL->size() == 2)?" 1 1 ":" 1 ") 
            << vecPoints.size()-1 << "\n";
      }

      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         //call the daughter version of the display
         (*iterD)->RecursiveDisplaySWC(inOut, inType, vecPoints);
      }
   }

   //Display the number of daughters and points, indented by level (mostly for testing).
   void RecursiveDisplay(ostream& inOut) const {
      //do the indenting
      int i;

      //display the daughter count and the points
      LineType::const_iterator iterL = m_Line.begin();
      LineType::const_iterator iterL_end = m_Line.end();
      if(m_iLevel > 0) {
         ++iterL;
      }
      for(; iterL != iterL_end; ++iterL) {
         for(i = 0; i <= m_iLevel; ++i) {
            inOut << "  ";
         }
         inOut << "(    " << iterL->at(0) << "  " 
            << iterL->at(1) << "  " 
            << iterL->at(2) << "  0.22);\n";
      }

      //call the daughter display
      DaughtersHolderType::const_iterator iter = m_Daughters.begin();
      DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
      if( iter != iter_end ) {
         for(i = 0; i <= m_iLevel; ++i) {
            inOut << "  ";
         }
         inOut << "(\n";
         while(iter != iter_end) {
            (*iter)->RecursiveDisplay(inOut);
            ++iter;
            if( iter != iter_end ) {
               for(i = 0; i <= m_iLevel; ++i) {
                  inOut << "  ";
               }
               inOut << "|\n";
            }
         }
         for(i = 0; i <= m_iLevel; ++i) {
            inOut << "  ";
         }
         inOut << ")\n";
      }
   }

   //Display the number of daughters and points, indented by level (mostly for testing).
   //Skips the initial part of the first branch (trims the dangling line).
   void RecursiveTrimmedDisplay(ostream& inOut) const {
      int iTrimTo = 0;
      if(m_Daughters.size() > 0) {
         //Find the daughter closest to the start of the branch.
         map<double, int> mapDistance_ID;
         for(int i = 0; i < m_Daughters.size(); ++i) {
            double dDummy4 = m_Line.DistanceToStart(m_Daughters[i]->GetLine().at(0));
            mapDistance_ID[m_Line.DistanceToStart(m_Daughters[i]->GetLine().at(0))] = i;
         }

         //Get the first point after the daughter (should actually be the daughter).
         double dDistanceToDaughter = mapDistance_ID.begin()->first-1.e-12;
         for(; iTrimTo < m_Line.size() && m_Line.DistanceToStart(m_Line[iTrimTo]) < dDistanceToDaughter; iTrimTo++) {
            double dDummy = m_Line.DistanceToStart(m_Line[iTrimTo]);
            double dDoubleDummy = 2.;
         }
      }

      //Catch an error.
      if(iTrimTo == m_Line.size()) {
         //Error
         cerr << "Error in CBranch::RecursiveTrimmedDisplay: Unable to find trim point .. skipping trim\n" << flush;
      }

      //do the indenting
      int i;

      //display the daughter count and the points
      LineType::const_iterator iterL = m_Line.begin();
      LineType::const_iterator iterL_end = m_Line.end();
      
      //Skip some
      for(int j = 0; j < iTrimTo-1; j++) {
         iterL++;
      }
      //Display the rest
      for(; iterL != iterL_end; ++iterL) {
         for(i = 0; i <= m_iLevel; ++i) {
            inOut << "  ";
         }
         inOut << "(    " << iterL->at(0) << "  " 
            << iterL->at(1) << "  " 
            << iterL->at(2) << "  0.22);\n";
      }

      //call the daughter display
      DaughtersHolderType::const_iterator iter = m_Daughters.begin();
      DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
      if( iter != iter_end ) {
         for(i = 0; i <= m_iLevel; ++i) {
            inOut << "  ";
         }
         inOut << "(\n";
         while(iter != iter_end) {
            (*iter)->RecursiveDisplay(inOut);
            ++iter;
            if( iter != iter_end ) {
               for(i = 0; i <= m_iLevel; ++i) {
                  inOut << "  ";
               }
               inOut << "|\n";
            }
         }
         for(i = 0; i <= m_iLevel; ++i) {
            inOut << "  ";
         }
         inOut << ")\n";
      }
   }

   //Display the number of daughters and points, indented by level (mostly for testing).
   void RawRecursiveDisplay(ostream& inOut) const {
      //display the daughter count and the points
      LineType::const_iterator iterL = m_Line.begin();
      LineType::const_iterator iterL_end = m_Line.end();
      for(; iterL != iterL_end; ++iterL) {
         inOut << iterL->operator[](0) << " " << iterL->operator[](1) << " " << iterL->operator[](2) << "\n";
      }

      //call the daughter display
      DaughtersHolderType::const_iterator iter = m_Daughters.begin();
      DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         (*iter)->RawRecursiveDisplay(inOut);
      }
   }

   //Display the segments and all recursive daughter segments
   //Format is x1 y1 z1 x2 y2 z2
   void RecursiveDisplaySegments(ostream& inOut, const int& inID) const {
      if(m_Line.size() > 1) {
         //display the daughter count and the points
         LineType::const_iterator iterL = m_Line.begin();
         LineType::const_iterator iterL2 = m_Line.begin();
         iterL2++;
         LineType::const_iterator iterL_end = m_Line.end();
         for(; iterL2 != iterL_end; ++iterL, ++iterL2) {
            inOut << iterL->operator[](0) << " " << iterL->operator[](1) << " " << iterL->operator[](2) << " "
               << iterL2->operator[](0) << " " << iterL2->operator[](1) << " " << iterL2->operator[](2) << " " << inID << "\n";
         }
      }

      //call the daughter display
      DaughtersHolderType::const_iterator iter = m_Daughters.begin();
      DaughtersHolderType::const_iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         (*iter)->RecursiveDisplaySegments(inOut, inID);
      }
   }

   //Update the daughters to set the parent properly.
   void UpdateDaughters() {
      DaughtersHolderType::iterator iter = m_Daughters.begin();
      DaughtersHolderType::iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         (*iter)->SetParent(shared_from_this());
      }
   }

   //Update all the maxima and minima.
   //Call if a daughter branch is deleted or added without the aid of AddPoint.
   //Also called by Shift, Scale, and Rotate.
   void UpdateAllMaxMin() {
      m_dMinX = m_Line.Min(0);
      m_dMinY = m_Line.Min(1);
      m_dMinZ = m_Line.Min(2);
      m_dMaxX = m_Line.Max(0);
      m_dMaxY = m_Line.Max(1);
      m_dMaxZ = m_Line.Max(2);
      UpdateMaxMin(m_dMinX, m_dMaxX, m_dMinY, m_dMaxY, m_dMinZ, m_dMaxZ);

      //daughter recursion
      DaughtersHolderType::iterator iter = m_Daughters.begin();
      DaughtersHolderType::iterator iter_end = m_Daughters.end();
      for(; iter != iter_end; ++iter) {
         (*iter)->UpdateAllMaxMin();
      }
   }


   //Update the maxima and minima recursively over this branch and all parents.
   void UpdateMaxMin(const double& inMinX,
                     const double& inMaxX,
                     const double& inMinY,
                     const double& inMaxY,
                     const double& inMinZ,
                     const double& inMaxZ) {
      //find the local max and min
      m_dMinX = min(m_dMinX, inMinX);
      m_dMinY = min(m_dMinY, inMinY);
      m_dMinZ = min(m_dMinZ, inMinZ);
      m_dMaxX = max(m_dMaxX, inMaxX);
      m_dMaxY = max(m_dMaxY, inMaxY);
      m_dMaxZ = max(m_dMaxZ, inMaxZ);

      //update the parent (unless we're at the root)
      if(!IsRoot()) {
         m_Parent->UpdateMaxMin(m_dMinX, m_dMaxX, m_dMinY, m_dMaxY, m_dMinZ, m_dMaxZ);
      }
   }

private:
   //helper functions

   //Test if the lines (defined by endpoints) intersect.
   //Formula stolen from mathworld.
   bool LineIntersect(const double& a1_x, const double& a1_y, const double& a1_z,
      const double& a2_x, const double& a2_y, const double& a2_z,
      const double& b1_x, const double& b1_y, const double& b1_z,
      const double& b2_x, const double& b2_y, const double& b2_z) const {
      //This is just a horrid formula (bad notation to follow)
      C3DVector a1(a1_x, a1_y, a1_z);
      C3DVector a2(a2_x, a2_y, a2_z);
      C3DVector b1(b1_x, b1_y, b1_z);
      C3DVector b2(b2_x, b2_y, b2_z);

      double s1, s2;
      double dDenom = ((a1-a2).Cross(b1-b2))*((a1-a2).Cross(b1-b2));
      if(dDenom < 1.e-14) {
         //These are parallel
         return false;
      } else {
         s1 = ((a1-a2).Cross(b1-b2))*((b1-b2).Cross(a2-b2)) / dDenom;
         s2 = ((a1-a2).Cross(b1-b2))*((a1-a2).Cross(a2-b2)) / dDenom;
         if( 0. < s1 && s1 < 1. && 0. < s2 && s2 < 1. ) {
            //The nearest point is on the line segment.
            double dDenom2 = ((a2-a1).Cross(b2-b1)).Magnitude();
            if(dDenom2 < 1.e-10) {
               //These are parallel; should have already been caught.
               return false;
            } else {
               double D = fabs((b1-a1)*((a2-a1).Cross(b2-b1)))/dDenom2;
               //Test if the lines get close enough to be considered intersecting.
               return (D < LINES_CLOSE);
            }
         }
      }
      return false; //should not get here
   }


   //Do the intersection recursion.
   //Modifies inIntersections.
   void DoCountSegmentIntersections(const double& x1, 
      const double& x2, 
      const double& y1, 
      const double& y2, 
      const double& z1, 
      const double& z2, 
      int& inIntersections) const {

      //Check if either of the input points are in the current bounding box.
      if( IsInRange(x1, y1, z1) || IsInRange(x2, y2, z2) ) {
         //Count the intersections with the current line.
         if( m_Line.size() > 1 ) {
            LineType::const_iterator iterL = m_Line.begin();
            LineType::const_iterator iterR = m_Line.begin();
            ++iterR;
            LineType::const_iterator iter_end = m_Line.end();
            for(; iterR != iter_end; ++iterL, ++iterR) {
               if( LineIntersect(iterL->at(0), iterL->at(1), iterL->at(2), 
                  iterR->at(0), iterR->at(1), iterR->at(2), 
                  x1, y1, z1, x2, y2, z2) ) {
                     ++inIntersections;
                  }
            }
         }
         
         //Check the daughters.
         DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
         DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
         for(; iterD != iterD_end; ++iterD) {
            (*iterD)->DoCountSegmentIntersections(x1, x2, y1, y2, z1, z2, inIntersections);
         }
      }

   }


   //return the distance squared
   double DistanceSquared(const double& inX1, const double& inY1, const double& inZ1,
      const double& inX2, const double& inY2, const double& inZ2) const {
         return (inX1-inX2)*(inX1-inX2) + (inY1-inY2)*(inY1-inY2) + (inZ1-inZ2)*(inZ1-inZ2);
      }

   //Do the bounding radius recursion
   void DoBoundingRadius(double& inR, const double& inX, const double& inY, const double& inZ) const {
      //check this neuron
      LineType::const_iterator iter = m_Line.begin();
      LineType::const_iterator iter_end = m_Line.end();
      for(; iter != iter_end; ++iter) {
         inR = max(inR, DistanceSquared(inX, inY, inZ, iter->operator[](0), iter->operator[](1), iter->operator[](2)));
      }

      //recurse over the daughters (if necessary)
      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         //no reason to check the daughter if its bounding box is inside the current radius
         if( DistanceSquared(inX, inY, inZ, (*iterD)->GetMinX(), (*iterD)->GetMinY(), (*iterD)->GetMinZ()) > inR ||
             DistanceSquared(inX, inY, inZ, (*iterD)->GetMinX(), (*iterD)->GetMaxY(), (*iterD)->GetMinZ()) > inR ||
             DistanceSquared(inX, inY, inZ, (*iterD)->GetMinX(), (*iterD)->GetMinY(), (*iterD)->GetMaxZ()) > inR ||
             DistanceSquared(inX, inY, inZ, (*iterD)->GetMinX(), (*iterD)->GetMaxY(), (*iterD)->GetMaxZ()) > inR ||
             DistanceSquared(inX, inY, inZ, (*iterD)->GetMaxX(), (*iterD)->GetMinY(), (*iterD)->GetMinZ()) > inR ||
             DistanceSquared(inX, inY, inZ, (*iterD)->GetMaxX(), (*iterD)->GetMaxY(), (*iterD)->GetMinZ()) > inR ||
             DistanceSquared(inX, inY, inZ, (*iterD)->GetMaxX(), (*iterD)->GetMinY(), (*iterD)->GetMaxZ()) > inR ||
             DistanceSquared(inX, inY, inZ, (*iterD)->GetMaxX(), (*iterD)->GetMaxY(), (*iterD)->GetMaxZ()) > inR ) {
                (*iterD)->DoBoundingRadius(inR, inX, inY, inZ);
             }
      }
   }

   //Do the total mean recursion.
   void DoTotalMean(unsigned& inN, double& inX, double& inY, double& inZ) const {
      //add up this neuron
      LineType::const_iterator iter = m_Line.begin();
      LineType::const_iterator iter_end = m_Line.end();
      for(; iter != iter_end; ++iter) {
         ++inN;
         inX += iter->operator[](0);
         inY += iter->operator[](1);
         inZ += iter->operator[](2);
      }

      //recurse over the daughters
      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         (*iterD)->DoTotalMean(inN, inX, inY, inZ);
      }
   }

   //Do total variance recursion.
   void DoTotalVariance(unsigned& inN, double& inVX, double& inVY, double& inVZ,
      const double& inMX, const double& inMY, const double& inMZ) const {
      //add up this neuron
      LineType::const_iterator iter = m_Line.begin();
      LineType::const_iterator iter_end = m_Line.end();
      for(; iter != iter_end; ++iter) {
         ++inN;
         inVX += (iter->operator[](0)-inMX)*(iter->operator[](0)-inMX);
         inVY += (iter->operator[](1)-inMY)*(iter->operator[](1)-inMY);
         inVZ += (iter->operator[](2)-inMZ)*(iter->operator[](2)-inMZ);
      }

      //recurse over the daughters
      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         (*iterD)->DoTotalVariance(inN, inVX, inVY, inVZ, inMX, inMY, inMZ);
      }
   }

   //Do the point intersection recursion.
   void DoDoesPointIntersect(const CPoint& inPoint, const double& inNearDistance, shared_ptr<CBranch>& inRet) {
      //Check this line.
      if( m_Line.IsPointNearLine(inPoint, inNearDistance) ) {
         inRet = shared_from_this();
      }

      if(inRet->IsNull()) {
         //Recurse over the daughters.
         DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
         DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
         for(; iterD != iterD_end; ++iterD) {
            (*iterD)->DoDoesPointIntersect(inPoint, inNearDistance, inRet);
         }
      }
   }

   //Return true if the point (first 3 numbers) is in the box (last 6 numbers)
   bool PointInBox(const double& inPX, const double& inPY, const double& inPZ,
                   const double& inBX1, const double& inBX2,
                   const double& inBY1, const double& inBY2,
                   const double& inBZ1, const double& inBZ2) const {
       return (inPX > inBX1) && (inPX < inBX2) && (inPY > inBY1) && (inPY < inBY2) && (inPZ > inBZ1) && (inPZ < inBZ2);
   }

   //return true if the point is in the bounding box
   bool IsInRange(const double& inX, const double& inY, const double& inZ) const {
      return ( (inX > m_dMinX) && (inX < m_dMaxX) ) ||
             ( (inY > m_dMinY) && (inY < m_dMaxY) ) ||
             ( (inZ > m_dMinZ) && (inZ < m_dMaxZ) );
   }

   //return true if the points are in the bounding box
   bool IsInRange(const double& inX1, const double& inX2, 
                  const double& inY1, const double& inY2,
                  const double& inZ1, const double& inZ2) const {
      return ( (inX1 > m_dMinX) && (inX1 < m_dMaxX) ) ||
             ( (inX2 > m_dMinX) && (inX2 < m_dMaxX) ) ||
             ( (inY1 > m_dMinY) && (inY1 < m_dMaxY) ) ||
             ( (inY2 > m_dMinY) && (inY2 < m_dMaxY) ) ||
             ( (inZ1 > m_dMinZ) && (inZ1 < m_dMaxZ) ) ||
             ( (inZ2 > m_dMinZ) && (inZ2 < m_dMaxZ) );
   }

   //Find the minimum distance to the input point
   void DoMinimumDistance(const double& inX, const double& inY, const double& inZ, double& inMD) const {
      LineType::const_iterator iter = m_Line.begin();
      LineType::const_iterator iter_end = m_Line.end();
      for(; iter != iter_end; ++iter) {
         inMD = min(inMD, DistanceSquared( iter->operator[](0), iter->operator[](1), iter->operator[](2), inX, inY, inZ ));
      }

      //recurse over the daughters
      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         (*iterD)->DoMinimumDistance(inX, inY, inZ, inMD);
      }
   }

   //Do the shift recursion.
   void DoShift(const PointType& inShift) {
      m_Line.Shift(inShift);

      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         (*iterD)->DoShift(inShift);
      }
   }

   //Do the scaling recursion.
   void DoScale(const PointType& inScale) {
      m_Line.Scale(inScale);

      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         (*iterD)->DoScale(inScale);
      }
   }

   //Do the rotation recursion.
   //Class T must have access to a double by operator [][] (it's a matrix).
   template<class T>
   void DoRotate(const T& inM) {
      m_Line.Rotate(inM);

      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         (*iterD)->DoRotate(inM);
      }
   }

   //Do recursion
   void RecursiveDisplaySWC(ostream& inOut, const int& inType, vector<CPoint>& inPoints) const {
      //cout << "branch.h gh1 " << GetLine().size() << "\n" << flush;
      //it is possible to have a blank daughter
      if(GetLine().size() > 1) {
         //find the parent id (-1 if it's not found)
         int iParent = -1;
         vector<CPoint>::const_iterator iterP 
            = find(inPoints.begin(), inPoints.end(), GetLine().at(0));
         if(iterP != inPoints.end()) {
            iParent = (int)distance((vector<CPoint>::const_iterator)inPoints.begin(), iterP) + 1;
         } else {
            cout << "Failed to find a parent in CBranch::RecursiveDisplay. Adding a new root and continuing\n" << flush;
         }

         LineType::const_iterator iterL = GetLine().begin();
         LineType::const_iterator iterL_end = GetLine().end();
         //the first point on the line is a duplicate
         ++iterL;

         //the second point has a special parent
         inPoints.push_back(*iterL);
         inOut << inPoints.size() << " " << inType << " "
			 << *iterL << ((iterL->size() == 2)?" 1 1 ":" 1 " )
            << iParent << "\n";
         ++iterL;

         //the rest of the points go in order
         for(; iterL != iterL_end; ++iterL) {
            inPoints.push_back(*iterL);
            inOut << inPoints.size() << " " << inType << " "
				<< *iterL << ((iterL->size() == 2)?" 1 1 ":" 1 ") 
               << inPoints.size()-1 << "\n";
         }
      }

      DaughtersHolderType::const_iterator iterD = m_Daughters.begin();
      DaughtersHolderType::const_iterator iterD_end = m_Daughters.end();
      for(; iterD != iterD_end; ++iterD) {
         //call the daughter version of the display
         (*iterD)->RecursiveDisplaySWC(inOut, inType, inPoints);
      }
   }

private:
   //data members
   shared_ptr<CBranch> m_Parent;
   DaughtersHolderType m_Daughters;
   int m_iLevel;
   LineType m_Line;
   double m_dMinX;
   double m_dMaxX;
   double m_dMinY;
   double m_dMaxY;
   double m_dMinZ;
   double m_dMaxZ;

};

#endif //BRANCH
