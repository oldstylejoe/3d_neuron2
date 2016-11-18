//Joe Snider
//1/05
//
//A line. A specialization of a deque<CPoints>.
//The line is represented as a deque<CPoints> with the format: P1, P2, P3, ... where each
//  point is connected by a straight line.
//The line between (P1,P2), ()... can be split up by calling Divide().

#include <iostream>
#include <math.h>
#include <deque>
#include <vector>

#include "point.h"

#ifndef LINE__
#define LINE__

using namespace std;

const double ON_DISTANCE = 1.e-6; //Sets how close a test point has to be to be "on" the line.

class CLine: public deque<CPoint> {
public:
   //typedefs

public:
   //constructors
   CLine(): deque<CPoint>() {}

   CLine(const CLine& inCopy): deque<CPoint>() {
       *this = inCopy;
   }

   CLine& operator=(const CLine& inCopy) {
       //Get the turning points.
       m_dqTurningPoints.assign(inCopy.m_dqTurningPoints.begin(), inCopy.m_dqTurningPoints.end());
       assign(inCopy.begin(), inCopy.end());

       return *this;
   }

public:
   //gets and sets

   const deque<CPoint>& GetTurningPoints() const {return m_dqTurningPoints;}

public:
   //interface

   //Shift the line.
   //Maintains the number of divisions (approximately, see Divide() for notes).
   void Shift(const CPoint& inShift) {
      deque<CPoint>::iterator iter = m_dqTurningPoints.begin();
      deque<CPoint>::iterator iter_end = m_dqTurningPoints.end();
      for(; iter != iter_end; ++iter) {
         *iter += inShift;
      }
      Divide((unsigned)size());
   }

   //Scale the line in each direction seperately.
   //Maintains the number of divisions (approximately, see Divide() for notes).
   void Scale(const CPoint& inScale) {
      deque<CPoint>::iterator iter = m_dqTurningPoints.begin();
      deque<CPoint>::iterator iter_end = m_dqTurningPoints.end();
      for(; iter != iter_end; ++iter) {
         iter->Scale(inScale);
      }
      Divide((unsigned)size());
   }

   //Rotate the line.
   //Warning no checking that inM is actually a rotation.
   //Maintains the number of divisions (approximately, see Divide() for notes).
   //Type T must have operator [][] access to doubles (it's a matrix).
   template <class T>
   void Rotate(const T& inM) {
      deque<CPoint>::iterator iter = m_dqTurningPoints.begin();
      deque<CPoint>::iterator iter_end = m_dqTurningPoints.end();
      for(; iter != iter_end; ++iter) {
         iter->LinearTransform(inM);
      }
      Divide(size());
   }

   //Find the min (max follows) in the input direction.
   //No testing that inDir exists.
   double Min(const unsigned& inDir) const {
      if(m_dqTurningPoints.size() == 0) {
         return 0.;
      }
      double dRet = m_dqTurningPoints[0].at(inDir);
      for(unsigned i = 1; i < m_dqTurningPoints.size(); ++i) {
         dRet = min(dRet, m_dqTurningPoints[i].at(inDir));
      }
      return dRet;
   }
   double Max(const unsigned& inDir) const {
      if(m_dqTurningPoints.size() == 0) {
         return 0.;
      }
      double dRet = m_dqTurningPoints[0].at(inDir);
      for(unsigned i = 1; i < m_dqTurningPoints.size(); ++i) {
         dRet = max(dRet, m_dqTurningPoints[i].at(inDir));
      }
      return dRet;
   }

   //Two lines are equal if the have the same turning points.
   bool operator==(const CLine& inLine) const {
      if(inLine.m_dqTurningPoints.size() == m_dqTurningPoints.size()) {
         for(unsigned i = 0; i < m_dqTurningPoints.size(); ++i) {
            if( (inLine.m_dqTurningPoints[i] - m_dqTurningPoints[i]).Length() > 1.e-12 ) {
               return false;
            }
         }
         return true;
      }
      return false;
   }

   //Return a point that represents a vector from the start to the end.
   CPoint StartToEnd() const {
      return (back()-front());
   }

   //Pop inN  points off the front of the line.
   //Maintains the approximate number of divisions.
   void PopFront(const int& inN) {
      if(0 < inN && inN < m_dqTurningPoints.size()) {
         unsigned uSize = size();
         double dOldLength = Length();
         for(int i = 0; i < inN; ++i) {
            m_dqTurningPoints.pop_front();
         }
         UnDivide();
         Divide( (unsigned) (double(uSize)*Length()/dOldLength) );
      }
   }

   //Shrink the line to the input point.
   //Will shrink past turning points and possibly remove turning points.
   //Returns true if anything changed
   bool Shrink(const CPoint& inPoint) {
      if( m_dqTurningPoints.size() < 2 ) {
         //Do nothing if there are no points; should not get here.
         cerr << "Warning in CLine::Shrink: not enough points on the line ("
               << m_dqTurningPoints.size() << " are there and need at least 2...returning\n" << flush;
         return false;
      }

      //Find the sub-segment that intersects the input point
      //This is the most common case, so check it first.
      int iIntersectingSegment = -1;
      for(unsigned i = 0; i < m_dqTurningPoints.size() && iIntersectingSegment < 0; ++i) {
         if( inPoint == m_dqTurningPoints[i] ) {
            iIntersectingSegment = i;
         }
      }
      //Did not find it as a turning point check intermediate parts
      if( iIntersectingSegment < 0 ) {
          for(unsigned i = 1; i < m_dqTurningPoints.size() && iIntersectingSegment < 0; ++i) {
            if( DoDoesSphereHover(m_dqTurningPoints[i-1], m_dqTurningPoints[i], inPoint, ON_DISTANCE) ) {
               iIntersectingSegment = i;
               cerr << "Warning in CLine::Shrink: unlikely line search for point " << inPoint 
                  << " line is ";
               for(int i = 0; i < m_dqTurningPoints.size(); ++i) {
                  cerr << m_dqTurningPoints[i] << " " << flush;
               }
               /*//////////////////////////////////testing only///////////////////////////////////
               for(unsigned j = 0; j < m_dqTurningPoints.size(); ++j) {
                  if( inPoint == m_dqTurningPoints[j] ) {cout << "force\n" << flush;}
               }
               //////////////////////////////////end testing only///////////////////////////////*/
            }
         }
      }

      //Do nothing if the input point is not on the line.
      if( iIntersectingSegment < 0 ) {
         cerr << "Warning in CLine::Shrink: could not find segment " << inPoint 
               << " line is ";
         for(int i = 0; i < m_dqTurningPoints.size(); ++i) {
            cerr << m_dqTurningPoints[i] << " " << flush;
         }
         cerr << "\n" << flush;
         return false;
      }

      //Record the number of subsegments and lengths.
      int iOldSegments = (int)size();
      double dOldLength = Length();

      //Set the turning points.
      m_dqTurningPoints.resize(iIntersectingSegment);

      //Add the input point unless it's already there (usually it's already there)
      //Note: could optimize by checking this first, but it's not a bottleneck, and this is safer
      if( m_dqTurningPoints.size() != 0) {
          if( m_dqTurningPoints.back() != inPoint ) {
              m_dqTurningPoints.push_back(inPoint);
          }
      } else {
         m_dqTurningPoints.push_back(inPoint);
      }

      //Redivide such that the new segments lengths are about the same as the old ones.
      Divide( int(double(iOldSegments)*Length()/dOldLength + 0.5) );

      return true;
   }

   //Return the distance along the line to the input point.
   //Returns the total length if the input point is not on the line.
   //"On" is defined as within ON_DISTANCE of the line.
   double DistanceToStart(const CPoint& inPoint) const {
      double dRet = 0.;
      if(m_dqTurningPoints.size() > 1) {
         if(inPoint == m_dqTurningPoints[0]) {return 0.;}
         CPoint pStart;
         CPoint pEnd = m_dqTurningPoints[0];
         int intersectingID = IsPointOnLine(inPoint);
         if(intersectingID < 0) {return Length();}
         for(unsigned i = 0; i < intersectingID; ++i) {
             pStart = m_dqTurningPoints[i];
             pEnd = m_dqTurningPoints[i+1];
             dRet += (pEnd-pStart).Length();
         }
         dRet += (pEnd-inPoint).Length();
      }
      return dRet;
   }

   //Remove any divisions of the line.
   //Does not change the turning points, just the intermediate segments.
   void UnDivide() {
      assign(m_dqTurningPoints.begin(), m_dqTurningPoints.end());
   }

   //Divide the line into (approximately) inNum segments.
   //This splits up the divisions between turning points by length,
   //  so there may be rounding error that adds at most one extra
   //  division per sub-segment.
   //Does not change the length, shape, etc.. of the line.
   void Divide(const unsigned& inNum) {
      //Do nothing if there aren't enough turning points.
      if(m_dqTurningPoints.size() < 2) {
         return;
      }

      //Undo any previous divisions.
      clear();

      //Find the number of divisions of each segment
      // weighted by the length of the segments.
      double dLength = Length();
      vector<int> vecDivisions(m_dqTurningPoints.size()-1);
      for(unsigned i = 0; i < m_dqTurningPoints.size()-1; ++i) {
         //Make sure there's at least one division.
         double dSegLength = (m_dqTurningPoints[i]-m_dqTurningPoints[i+1]).Length();
         //vecDivisions[i] += max(1, int(double(inNum)*dSegLength/dLength + 0.5));
         vecDivisions[i] = max(1, int(double(inNum)*dSegLength/dLength + 0.5));
      }

      //Do the dividing
      for(unsigned j = 0; j < m_dqTurningPoints.size()-1; ++j) {
         //Find the slope and offset for a parametric representation of the line.
         CPoint pntSlope = m_dqTurningPoints[j+1] - m_dqTurningPoints[j];
         CPoint pntOffset = m_dqTurningPoints[j];

         //Insert the new segments. To get N-segments we need to insert N-1 new points.
         deque<CPoint>::push_back(m_dqTurningPoints[j]);
         for(int k = 1; k < vecDivisions[j]; ++k) {
            CPoint pntInsert(pntSlope);
            pntInsert *= double(k)/double(inNum);
            pntInsert += pntOffset;
            deque<CPoint>::push_back(pntInsert);
         }
      }
      //insert the last point
      deque<CPoint>::push_back(m_dqTurningPoints.back());
   }

   //Allow easy addition of common dimensions
   void push_back(const double& inX, const double& inY) {
      CPoint pInsert;
      pInsert.push_back(inX);
      pInsert.push_back(inY);
      AddTurningPoint(pInsert);
   }

   void push_back(const double& inX, const double& inY, const double& inZ) {
      CPoint pInsert;
      pInsert.push_back(inX);
      pInsert.push_back(inY);
      pInsert.push_back(inZ);
      AddTurningPoint(pInsert);
   }
   template <class T>
   void push_back(const T& inX) {
      AddTurningPoint(inX);
   }

   //Insert a new point.
   //Adds a new point to the m_dqTurningPoint list and to the base class.
   //Type T must be castable to a CPoint.
   template <class T>
   void AddTurningPoint(const T& inX) {
      m_dqTurningPoints.push_back(inX);
      deque<CPoint>::push_back(inX);
   }

   //Insert multiple new points (stl iterator style).
   //Adds a new point to the m_dqTurningPoint list and to the base class.
   //Type T must be forward iterator with value castable to a CPoint.
   template <class T>
   void AddTurningPoint(T inBegin, T inEnd) {
      for(; inBegin != inEnd; ++inBegin) {
         m_dqTurningPoints.push_back(*inBegin);
         deque<CPoint>::push_back(*inBegin);
      }
   }

   //Allow easy addition of common dimensions
   void push_front(const double& inX, const double& inY) {
      CPoint pInsert;
      pInsert.push_back(inX);
      pInsert.push_back(inY);
      AddTurningPointFront(pInsert);
   }

   void push_front(const double& inX, const double& inY, const double& inZ) {
      CPoint pInsert;
      pInsert.push_back(inX);
      pInsert.push_back(inY);
      pInsert.push_back(inZ);
      AddTurningPointFront(pInsert);
   }
   template <class T>
   void push_front(const T& inX) {
      AddTurningPointFront(inX);
   }

   //Insert a new point.
   //Adds a new point to the m_dqTurningPoint list and to the base class.
   //Type T must be castable to a CPoint.
   template <class T>
   void AddTurningPointFront(const T& inX) {
      m_dqTurningPoints.push_front(inX);
      deque<CPoint>::push_front(inX);
   }

   //Return the total length of all the segments.
   double Length() const {
      double dRet = 0.;
      if(m_dqTurningPoints.size() > 1) {
         deque<CPoint>::const_iterator iterL = m_dqTurningPoints.begin();
         deque<CPoint>::const_iterator iterR = m_dqTurningPoints.begin();
         ++iterR;
         deque<CPoint>::const_iterator iter_end = m_dqTurningPoints.end();
         for(; iterR != iter_end; ++iterL, ++iterR) {
            dRet += (*iterL-*iterR).Length();
         }
      }
      return dRet;
   }

   //Check if any of the segments intersect the input sphere.
   //Algorithm first checks for a real intersection point of the extended line,
   // and, if there is one, then it checks if the intersection point lies on the segment.
   //The test point is frequently a turning point, so those are checked explicitly.
   bool DoesSphereIntersect(const CPoint& inCenter, const double& inRadius) const {
      if(size() < 2) {
         return false;
      } else {
         //loop over all start and stop points
         const_iterator iterL = m_dqTurningPoints.begin();
         const_iterator iterR = m_dqTurningPoints.begin();
         ++iterR;
         const_iterator iter_end = m_dqTurningPoints.end();
         if(*iterL == inCenter) {return true;}
         for(; iterR != iter_end; ++iterL, ++iterR) {
            if(*iterR == inCenter) {return true;}
            if( DoDoesSphereIntersect(*iterL, *iterR, inCenter, inRadius) ) {
               return true;
            }
         }
      }

      //If every test fails, then we get here and return false
      return false;
   }

   //Check if any of the segments intersect the input sphere and that the center of the input sphere
   //  is above (in a cylinder around) at least one segment.
   //Algorithm first checks for a real intersection point of the extended line,
   // and if there is one, then it checks that the angles formed by the triangle of the center 
   // of the input sphere and two end points of the segment are all less than pi/2. 
   //The test point is frequently a turning point, so those are checked explicitly.
   bool DoesSphereHover(const CPoint& inCenter, const double& inRadius) const {
      if(size() < 2) {
         return false;
      } else {
         //loop over all start and stop points
         const_iterator iterL = m_dqTurningPoints.begin();
         const_iterator iterR = m_dqTurningPoints.begin();
         ++iterR;
         const_iterator iter_end = m_dqTurningPoints.end();
         if(*iterL == inCenter) {return true;}
         for(; iterR != iter_end; ++iterL, ++iterR) {
            if(*iterR == inCenter) {return true;}
            if( DoDoesSphereHover(*iterL, *iterR, inCenter, inRadius) ) {
               return true;
            }
         }
      }

      //If every test fails, then we get here and return false
      return false;
   }
   
   //Return the segment number (0, 1, ...) the input point intersects,
   //or -1 if no intersection.
   int IsPointOnLine(const CPoint& inPoint) const {
       if(m_dqTurningPoints.size() < 2) {return -1;}
       //first check the turning points (frequently the case, so do it first)
       if(inPoint == m_dqTurningPoints[0]) {return 0;}
       for(unsigned i = 1; i < m_dqTurningPoints.size(); ++i) {
           if(inPoint == m_dqTurningPoints[i]) {return i-1;}
       }
       
       //next check between the segments.
       for(unsigned i = 0; i < m_dqTurningPoints.size()-1; ++i) {
          if(DoDoesSphereHover(m_dqTurningPoints[i], m_dqTurningPoints[i+1], inPoint, ON_DISTANCE)) {
             return i;
          }
       }
       
       //not found
       return -1;
   }

   //Return true if the input point is near the line.
   //Just calls DoesSphereIntersect
   bool IsPointNearLine(const CPoint& inPoint, const double& inDistance) const {
      return DoesSphereIntersect(inPoint, inDistance);
   }

   //Find the closest point on the line to the input point.
   CPoint FindClosestPoint(const CPoint& inPoint) const {
      CPoint pReturn;
      if(m_dqTurningPoints.size() > 1) {
         double dDist = 1.e100;
         for(unsigned i = 0; i < m_dqTurningPoints.size()-1; ++i) {
            //Find the perpendicular distance (restricted to on the line) to the input point.
            //bad variable names
            CPoint x0 = m_dqTurningPoints[i];
            CPoint x1 = m_dqTurningPoints[i+1];
            CPoint x10 = x1-x0;
            //restrict to the segment in parametric form.
            double t = ( (inPoint-x0)*x10 )/( x10*x10 );
            t = (t<0.)?0.:t;
            t = (t>1.)?1.:t;

            //the intersection point
            CPoint a = x10*t+x0;
            double dDistTemp = (inPoint-a).LengthSquared();
            if(dDistTemp < dDist) {
               pReturn = a;
               dDist = dDistTemp;
            }
         }
      }
      return pReturn;
   }

private:
   //helpers

   //Helper for Sphere hovering (see DoesSphereHover for comments).
   //Just does the test on the line segment from inStart to inEnd.
   bool DoDoesSphereHover(const CPoint& inStart, const CPoint& inEnd, const CPoint& inCenter, const double& inRadius) const {
      //warning: bad but short variable names
      const double r2 = inRadius*inRadius;
      CPoint x1 = inStart;
      x1 -= inCenter;
      CPoint x2 = inEnd;
      x2 -= inCenter;
      const double x1_2 = x1.Dot(x1);
      const double x2_2 = x2.Dot(x2);
      //Insure that the segment has non-zero length.
      if( !(x1-x2).IsZero() ) {
         //Check if either end is in the circle (that means intersection by definition).
         //Check the discriminant (see derivation elsewhere (page 89 of notebook)).
         const CPoint x2_minus_x1 = x2-x1;
         const double x1_dot_x2_minus_x1 = x1.Dot(x2_minus_x1);
         const double x2_minus_x1_squared = x2_minus_x1.Dot(x2_minus_x1);
         const double D = x1_dot_x2_minus_x1*x1_dot_x2_minus_x1 - x2_minus_x1_squared*(x1_2-r2);

         //check the discriminant to determine if the sphere intersects the extended line
         if( (x1_2 < r2) || (x2_2 < r2) || D > 0. ) {
            //Check that the angles formed by the triangle (x_1, x_2, p) are less than pi/2.
            //Namely
            //                      Input point
            //                       |
            //                       V
            //                    /-p-\                      Input line
            //                /---     ---\                   |
            //            /---             ----\              V
            //  ------(x_1)---------------------(x_2)---------------
            //            Angles have to be less than pi/2 (90 degrees) on both sides.
            //Find the cosine of the angles from the law of cosines, and check the sign.
            //(I know there's no cosine; do the math).
            if( x1_dot_x2_minus_x1 < 0. && x2.Dot(x2_minus_x1) > 0.) {
               return true;
            }
         }
      }

      //All tests failed, return false.
      return false;
   }

   //Helper for Sphere intersection (see DoesSphereIntersect for comments).
   //Just does the test on the line segment from inStart to inEnd.
   bool DoDoesSphereIntersect(const CPoint& inStart, const CPoint& inEnd, const CPoint& inCenter, const double& inRadius) const {
      //warning: bad but short variable names
      double r2 = inRadius*inRadius;
      CPoint x1 = inStart;
      x1 -= inCenter;
      CPoint x2 = inEnd;
      x2 -= inCenter;
      double x1_2 = x1.Dot(x1);
      double x2_2 = x2.Dot(x2);
      //Insure that the segment has non-zero length.
      if( !(x1-x2).IsZero() ) {
         //Second, check if either end is in the circle (that means intersection by definition) 
         if( (x1_2 < r2) || (x2_2 < r2) ) {
            return true;
         } else {
            //have to work a little harder
            //check the discriminant (see derivation elsewhere (page 89 of notebook))
            CPoint x2_minus_x1 = x2-x1;
            double temp1 = -1.*x1.Dot(x2_minus_x1);
            double temp2 = x2_minus_x1.Dot(x2_minus_x1);
            double D = temp1*temp1 - temp2*(x1_2-r2);

            //check the discriminant to determine if the sphere intersects the extended line
            if( D >= 0 ) {
               //Have to actually find the intersection points now, and check if they're on the segment (not just the line)
               //First step is to find a dimension in which the coordinates are not equal.
               //Note: if there is no such dimension, then this segment has length zero, and we would not get here.
               int iDimTest = 0;
               while( x1[iDimTest] == x2[iDimTest] ) {++iDimTest;}
               //Two solutions to the parametric eqn for the line.
               //If either are on the segment, then intersection occurs.
               //Note: only one coordinate needs to be checked since the new point is on the line.
               D = sqrt(D);
               double t1 = (temp1 + D)/temp2;
               double i1 = x1[iDimTest] + (x2[iDimTest]-x1[iDimTest])*t1;
               //Now, find the appropriate boundaries, in the sense of smaller and larger on the real line.
               double x_left, x_right;
               if(x1[iDimTest] < x2[iDimTest]) {
                  x_left = x1[iDimTest];
                  x_right = x2[iDimTest];
               } else {
                  x_left = x2[iDimTest];
                  x_right = x1[iDimTest];
               }
               if( (x_left <= i1) && (i1 <= x_right) ) {
                  return true;
               }
            }
         }
      }

      //no tests passed so return false.
      return false;
   }

private:
   //data members
   deque<CPoint> m_dqTurningPoints;
};

ostream& operator<<(ostream& inOut, const CLine& inLine) {
   CLine::const_iterator iter = inLine.begin();
   CLine::const_iterator iter_end = inLine.end();
   for(; iter != iter_end; ++iter) {
      inOut << *iter << " ";
   }
   return inOut;
}

#endif //LINE__
