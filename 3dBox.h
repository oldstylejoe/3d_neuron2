//Joe Snider
//11/05
//
//A 3d box that can tell if it intersects with various shapes.
//Box is specified by opposite corners (any two will do).
//The box is oriented with the x-y-z axis.
//Corners are stored such that a = (min(x), min(y), min(z)), b = (max(x), max(y), max(z))
//   This allows some optimizations.

#ifndef __3D_BOX
#define __3D_BOX

#include <iostream>
#include <math.h>
#include <vector>
#include <map>

#include "3dvector.h"

using namespace std;

const int DEPTH = 10;

class C3DBox{
public:
   //typedefs...

public:
   //constructors

   C3DBox(): m_vecA(0.,0.,0.), m_vecB(0.,0.,0.) {}
   C3DBox(const double& inA1, const double& inA2, const double& inA3,
      const double& inB1, const double& inB2, const double& inB3) {
         m_vecA.SetTo( min(inA1, inB1), min(inA2, inB2), min(inA3, inB3) );
         m_vecB.SetTo( max(inA1, inB1), max(inA2, inB2), max(inA3, inB3) );
   }
   C3DBox(const C3DVector& inA, const C3DVector& inB): m_vecA(inA), m_vecB(inB) {}

   C3DBox(const C3DBox& inCopy) {
      *this = inCopy;
   }
   C3DBox& operator=(const C3DBox& inCopy) {
      m_vecA = inCopy.GetA();
      m_vecB = inCopy.GetB();
      return *this;
   }

   ~C3DBox() {}

public:
   //gets and sets
   C3DVector GetA() const {return m_vecA;}
   C3DVector GetB() const {return m_vecB;}
   void SetAB(const double& inA1, const double& inA2, const double& inA3,
      const double& inB1, const double& inB2, const double& inB3) {
         m_vecA.SetTo( min(inA1, inB1), min(inA2, inB2), min(inA3, inB3) );
         m_vecB.SetTo( max(inA1, inB1), max(inA2, inB2), max(inA3, inB3) );
   }

public:
   //interface

   //Return the center of the box.
   C3DVector Center() const {
      return 0.5*(m_vecA+m_vecB);
   }

   double DumbIntersectionLength(const C3DVector& inX, const C3DVector& inY) const {
      //If both end points are inside return the length of the line.
      if( IsPointInBox(inX)&&IsPointInBox(inY) ) {
         return (inX-inY).Magnitude();
      }

      //If one is in and one out, return 1/2 the length.
      if( IsPointInBox(inX)||IsPointInBox(inY) ) {
         return 0.5*(inX-inY).Magnitude();
      }

      //handle both out
      if( DoesSegmentIntersect(inX, inY) ) {
         //Check if the segment is long enough to matter
         double dLongestSide = max(fabs(inX.GetX()-inY.GetX()), 
            max(fabs(inX.GetY()-inY.GetY()), fabs(inX.GetZ()-inY.GetZ())));
         if( 3.*(inX-inY).Magnitude() > dLongestSide) {
            //just return the longest side
            return dLongestSide;
         }
      }

      return 0.;
   }

   //Return the length of the segment within the box (may be 0).
   double IntersectionLength(const C3DVector& inX, const C3DVector& inY) const {
      //Get the corners (don't change the order).
      //In order coords are: 000, 001, 010, 011, 100, 101, 110, 111
      vector<C3DVector> vecCorners;
      vecCorners.push_back( m_vecA );
      vecCorners.push_back( C3DVector(m_vecA.GetX(), m_vecA.GetY(), m_vecB.GetZ()) );
      vecCorners.push_back( C3DVector(m_vecA.GetX(), m_vecB.GetY(), m_vecA.GetZ()) );
      vecCorners.push_back( C3DVector(m_vecA.GetX(), m_vecB.GetY(), m_vecB.GetZ()) );
      vecCorners.push_back( C3DVector(m_vecB.GetX(), m_vecA.GetY(), m_vecA.GetZ()) );
      vecCorners.push_back( C3DVector(m_vecB.GetX(), m_vecA.GetY(), m_vecB.GetZ()) );
      vecCorners.push_back( C3DVector(m_vecB.GetX(), m_vecB.GetY(), m_vecA.GetZ()) );
      vecCorners.push_back( m_vecB );

      //store the faces; 6 groups of 3 indices of vecCorners (3 corners makes a face)
      //Warning: the order of vectors in each face matters (makes face's unit vector positive).
      vector<vector<int> > vecFaces(6);
      vecFaces[0].push_back(0); vecFaces[0].push_back(2); vecFaces[0].push_back(1);
      vecFaces[1].push_back(4); vecFaces[1].push_back(6); vecFaces[1].push_back(5);
      vecFaces[2].push_back(0); vecFaces[2].push_back(1); vecFaces[2].push_back(4);
      vecFaces[3].push_back(2); vecFaces[3].push_back(3); vecFaces[3].push_back(6);
      vecFaces[4].push_back(0); vecFaces[4].push_back(4); vecFaces[4].push_back(2);
      vecFaces[5].push_back(1); vecFaces[5].push_back(5); vecFaces[5].push_back(3);

      //If both end points are inside return the length of the line.
      if( IsPointInBox(inX)&&IsPointInBox(inY) ) {
         return (inX-inY).Magnitude();
      }

      bool bJunk1 = IsPointInBox(inX);
      bool bJunk2 = IsPointInBox(inY);

      //Check if the line intersects at all
      if( DoesSegmentIntersect(inX, inY) ) {
         /////////////////////////testing only////////////////////////////
         //ofstream ofDump("error_dump.txt");
         //for(int j = 0; j < 8; ++j) {
         //   ofDump << vecCorners[j] << "\n";
         //}
         //ofDump << inX << "\n" << inY << "\n" << flush;
         ////////////////////////end testing only.///////////////////////

         //Find which faces (or face) are intersected
         int iF1 = -1; //face id 1 and face id 2
         int iF2 = -1;
         for(int i = 0; i < 6 && iF1 < 0 && iF2 < 0; ++i) {
            if( DoesSegmentIntersectFace( 
               vecCorners[vecFaces[i][0]], vecCorners[vecFaces[i][1]], 
               vecCorners[vecFaces[i][2]], inX, inY) ) {
                  if( iF1 < 0 ) {
                     iF1 = i;
                  } else {
                     iF2 = i;
                  }
            }
         }

         //Find the point that face 1 intersects.
         C3DVector vecIntersection1 = 
            FaceLineIntersection(vecCorners[vecFaces[iF1][0]], vecCorners[vecFaces[iF1][1]], 
            vecCorners[vecFaces[iF1][2]], inX, inY);

         //May intersect only one face.
         if(iF2 < 0) {
            //Return distance from interior endpoint to the face
            if( IsPointInBox(inX) ) {
               return (vecIntersection1-inX).Magnitude();
            } else {
               return (vecIntersection1-inY).Magnitude();
            }
         }

         //Get the other face intersection point
         C3DVector vecIntersection2 = 
            FaceLineIntersection(vecCorners[vecFaces[iF2][0]], vecCorners[vecFaces[iF2][1]], 
            vecCorners[vecFaces[iF2][2]], inX, inY);

         //Return the distance between face intersectin points.
         return (vecIntersection2-vecIntersection1).Magnitude();
      }

      //Default to no intersection
      return 0.;
   }

   //Return the volume of the box.
   double Volume() const {
      return fabs( (m_vecA.GetX()-m_vecB.GetX())*
         (m_vecA.GetY()-m_vecB.GetY())*
         (m_vecA.GetZ()-m_vecB.GetZ()) );
   }

   //Create a random sub-box of the input box.
   //Type T must have operator() return a double in [0,1) (CRand3 works)
   template<class T>
   C3DBox RandomSubBox(T& inRand) const {
      //Just pick two points on the line from a->b in each direction.
      double dRX1 = inRand();
      double dRX2 = inRand();
      double dRY1 = inRand();
      double dRY2 = inRand();
      double dRZ1 = inRand();
      double dRZ2 = inRand();
      C3DBox boxReturn( (m_vecB.GetX()-m_vecA.GetX())*min(dRX1,dRX2) + m_vecA.GetX(),
         (m_vecB.GetY()-m_vecA.GetY())*min(dRY1,dRY2) + m_vecA.GetY(),
         (m_vecB.GetZ()-m_vecA.GetZ())*min(dRZ1,dRZ2) + m_vecA.GetZ(),
         (m_vecB.GetX()-m_vecA.GetX())*max(dRX1,dRX2) + m_vecA.GetX(),
         (m_vecB.GetY()-m_vecA.GetY())*max(dRY1,dRY2) + m_vecA.GetY(),
         (m_vecB.GetZ()-m_vecA.GetZ())*max(dRZ1,dRZ2) + m_vecA.GetZ());
      return boxReturn;
   }

   //Test if a point is in the box.
   bool IsPointInBox(const C3DVector& inV) const {
      return m_vecA.GetX() < inV.GetX() && inV.GetX() < m_vecB.GetX() &&
         m_vecA.GetY() < inV.GetY() && inV.GetY() < m_vecB.GetY() &&
         m_vecA.GetZ() < inV.GetZ() && inV.GetZ() < m_vecB.GetZ();
   }

   //Test if the input segment intersects the box.
   bool DoesSegmentIntersect(const double& inX1, const double& inY1, const double& inZ1,
      const double& inX2, const double& inY2, const double& inZ2) const {
         C3DVector vecX(inX1, inY1, inZ1);
         C3DVector vecY(inX2, inY2, inZ2);
         return DoesSegmentIntersect(vecX, vecY);
   }

   bool DoesSegmentIntersect(const C3DVector& inX, const C3DVector& inY) const {
      //Find the corners.
      //Warning: the order matters
      C3DVector vec000 = m_vecA;
      C3DVector vec001(m_vecA.GetX(), m_vecA.GetY(), m_vecB.GetZ());
      C3DVector vec010(m_vecA.GetX(), m_vecB.GetY(), m_vecA.GetZ());
      C3DVector vec011(m_vecA.GetX(), m_vecB.GetY(), m_vecB.GetZ());
      C3DVector vec100(m_vecB.GetX(), m_vecA.GetY(), m_vecA.GetZ());
      C3DVector vec101(m_vecB.GetX(), m_vecA.GetY(), m_vecB.GetZ());
      C3DVector vec110(m_vecB.GetX(), m_vecB.GetY(), m_vecA.GetZ());
      C3DVector vec111 = m_vecB;

      //Considered to intersect if any endpoint is in the box.
      //This could just be put in the bunch of ||'s below, but I want to be sure it runs first.
      if( IsPointInBox(inX)||IsPointInBox(inY) ) {
         return true;
      }

      //Check the 6 faces (checking a sphere around the box is slooow, so just do this)
      //Warning: the order matters
      return DoesSegmentIntersectFace(vec000, vec010, vec001, inX, inY) ||
         DoesSegmentIntersectFace(vec100, vec110, vec101, inX, inY) ||
         DoesSegmentIntersectFace(vec000, vec001, vec100, inX, inY) ||
         DoesSegmentIntersectFace(vec010, vec011, vec110, inX, inY) ||
         DoesSegmentIntersectFace(vec000, vec100, vec010, inX, inY) ||
         DoesSegmentIntersectFace(vec001, vec101, vec011, inX, inY);
   }

   //Test if the input circle intersects the box.
   bool DoesCircleIntersect(const C3DVector& inCenter, const double& inRadius) const {
      //Pad inRadius to the 3 directions and test if inCenter is contained.
      return( m_vecA.GetX()-inRadius < inCenter.GetX() &&
         inCenter.GetX() < m_vecB.GetX()+inRadius &&
         m_vecA.GetY()-inRadius < inCenter.GetY() &&
         inCenter.GetY() < m_vecB.GetY()+inRadius &&
         m_vecA.GetZ()-inRadius < inCenter.GetZ() &&
         inCenter.GetZ() < m_vecB.GetZ()+inRadius );
   }

private:
   //helpers

   //Return the intersection point of a face (a,b,c) and a segment (x,y).
   //Actually finds the intersection point of the plane defined by the face
   //  with the line defined by the segment.
   //Note that this is specialized for faces parallel to the axis, i.e.
   //   with unit vector || to x, y, or z axis.
   C3DVector FaceLineIntersection(const C3DVector& a, const C3DVector& b, const C3DVector& c,
      const C3DVector& x, const C3DVector& y) const {

         //Find the direction of the unit vector of the plane.
         const C3DVector n = (b-a).Cross(c-a);

         //Since this is such a bottle neck, it's worth it to have just these 3 if statements.
         //The code is just repeated for each direction.
         //Remember: faces are oriented with the x,y,z axis and chosen in a special order
         if (n.GetX() > 1.e-12) {
            const double t = (a.GetX() - x.GetX()) / (y.GetX() - x.GetX());
            return (y-x)*t+x;
         } else if (n.GetY() > 1.e-12) {
            const double t = (a.GetY() - x.GetY()) / (y.GetY() - x.GetY());
            return (y-x)*t+x;
         } else {
            const double t = (a.GetZ() - x.GetZ()) / (y.GetZ() - x.GetZ());
            return (y-x)*t+x;
         }

         return C3DVector(0.,0.,0.);  //should not get here
   }

   //Check if a segment intersects a face.
   //The faces are parallel to one of the axis (by definition).
   // so just find the intersection point of the extended line and plane
   // and see if intersection point is on the (restricted) segment and the face.
   //x,y specify the line.
   //a,b,c specify the face (see notes)
   //This is a significant bottleneck, so optimize as well as possible.
   bool DoesSegmentIntersectFace(const C3DVector& a, const C3DVector& b, const C3DVector& c,
      const C3DVector& x, const C3DVector& y) const {
         //Find the direction of the unit vector of the plane.
         const C3DVector n = (b-a).Cross(c-a);

         //Since this is such a bottle neck, it's worth it to have just these 3 if statements.
         //The code is just repeated for each direction.
         //Remember: faces are oriented with the x,y,z axis and chosen in a special order
         if (n.GetX() > 1.e-12) {
            const double t = (a.GetX() - x.GetX()) / (y.GetX() - x.GetX());
            if( 0. < t && t < 1. ) {
               const C3DVector p = (y-x)*t+x;
               return a.GetY() < p.GetY() && p.GetY() < b.GetY() &&
                  a.GetZ() < p.GetZ() && p.GetZ() < c.GetZ();
            }
         } else if (n.GetY() > 1.e-12) {
            const double t = (a.GetY() - x.GetY()) / (y.GetY() - x.GetY());
            if( 0. < t && t < 1. ) {
               const C3DVector p = (y-x)*t+x;
               return a.GetZ() < p.GetZ() && p.GetZ() < b.GetZ() &&
                  a.GetX() < p.GetX() && p.GetX() < c.GetX();
            }
         } else {
            const double t = (a.GetZ() - x.GetZ()) / (y.GetZ() - x.GetZ());
            if( 0. < t && t < 1. ) {
               const C3DVector p = (y-x)*t+x;
               return a.GetX() < p.GetX() && p.GetX() < b.GetX() &&
                  a.GetY() < p.GetY() && p.GetY() < c.GetY();
            }
         }

         return false; //should not get here
   }

   //Check if the point is within inR of the origin.
   bool IsPointInCircle(const C3DVector& inV, const double& inR) const {
      return inV.MagnitudeSquared() < inR*inR;
   }

   //Check if the line connecting the two input points 
   //is inside the sphere of radius inR about the origin.
   //Looks for solution of parameterized line connecting a to b (b+(a-b)t) and intersecting sphere.
   //Also requires 0<t<1 to restrict to the segment connecting a and b.
   //Bad notation.
   bool IsLineInCircle(const C3DVector& a, const C3DVector& b, const double& r) const {
      C3DVector a_minus_b = a-b;
      double D = pow(b*a_minus_b, 2) - a_minus_b.MagnitudeSquared()*(b*b-r*r);
      if(D > 0) {
         D = sqrt(D);
         double t1 = (b*a_minus_b + D)/a_minus_b.MagnitudeSquared();
         double t2 = (b*a_minus_b - D)/a_minus_b.MagnitudeSquared();
         return (0.<t1 && t1<1.) || (0.<t2 && t2<1.);
      } else {
         return false;
      }
      return false; //should not get here
   }

private:
   //data members
   C3DVector m_vecA;
   C3DVector m_vecB;
};

#endif //__3D_BOX
