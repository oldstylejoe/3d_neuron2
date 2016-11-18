//Joe Snider
//10/05
//
//A 3d vector class. Defines cross product, so must be in 3d.

#ifndef joe_3D_VECTOR
#define joe_3D_VECTOR

#include <iostream>
#include <math.h>

using namespace std;

class C3DVector{
public:
   //typedefs

public:
   //constructors
   C3DVector() {}
   C3DVector(const double& inX, const double& inY, const double& inZ): m_dX(inX), m_dY(inY), m_dZ(inZ) {}
   
   ~C3DVector() {}

   C3DVector(const C3DVector& inCopy) {
      *this = inCopy;
   }

   C3DVector& operator=(const C3DVector& inCopy) {
      m_dX = inCopy.GetX();
      m_dY = inCopy.GetY();
      m_dZ = inCopy.GetZ();

      return *this;
   }

public:
   //gets and sets

   double GetX() const {return m_dX;}
   void SetX(const double& inX) {m_dX = inX;}

   double GetY() const {return m_dY;}
   void SetY(const double& inY) {m_dY = inY;}

   double GetZ() const {return m_dZ;}
   void SetZ(const double& inZ) {m_dZ = inZ;}

public:
   //interface

   void SetTo(const double& inX, const double& inY, const double& inZ) {
      m_dX = inX;
      m_dY = inY;
      m_dZ = inZ;
   }

   double Magnitude() const {
      return sqrt(m_dX*m_dX+m_dY*m_dY+m_dZ*m_dZ);
   }

   double MagnitudeSquared() const {
      return m_dX*m_dX+m_dY*m_dY+m_dZ*m_dZ;
   }

   //Will return the z-comp if i is not 0 or 1
   double operator[](const int& i) const {
      if(i == 0) {
         return m_dX;
      } else if(i == 1) {
         return m_dY;
      } else {
         return m_dZ;
      }

      return 0.; //should not get here
   }

   //Will modify the z-comp if i is not 0 or 1
   double& operator[](const int& i) {
      if(i == 0) {
         return m_dX;
      } else if(i == 1) {
         return m_dY;
      } else {
         return m_dZ;
      }
      return m_dZ; //should not get here
   }

   C3DVector operator+(const C3DVector& inX) const {
      C3DVector vReturn(m_dX+inX.GetX(), m_dY+inX.GetY(), m_dZ+inX.GetZ());
      return vReturn;
   }
   C3DVector operator-(const C3DVector& inX) const {
      C3DVector vReturn(m_dX-inX.GetX(), m_dY-inX.GetY(), m_dZ-inX.GetZ());
      return vReturn;
   }
   C3DVector operator*(const double& inScale) const {
      C3DVector vReturn(m_dX*inScale, m_dY*inScale, m_dZ*inScale);
      return vReturn;
   }
   C3DVector operator/(const double& inScale) const {
      C3DVector vReturn(m_dX/inScale, m_dY/inScale, m_dZ/inScale);
      return vReturn;
   }
   C3DVector& operator+=(const C3DVector& inX) {
      m_dX += inX.GetX();
      m_dY += inX.GetY();
      m_dZ += inX.GetZ();
      return *this;
   }
   C3DVector& operator-=(const C3DVector& inX) {
      m_dX -= inX.GetX();
      m_dY -= inX.GetY();
      m_dZ -= inX.GetZ();
      return *this;
   }
   C3DVector& operator*=(const double& inScale) {
      m_dX *= inScale;
      m_dY *= inScale;
      m_dZ *= inScale;
      return *this;
   }
   C3DVector& operator/=(const double& inScale) {
      m_dX /= inScale;
      m_dY /= inScale;
      m_dZ /= inScale;
      return *this;
   }
   double operator*(const C3DVector& inX) const {
      return m_dX*inX.GetX() + m_dY*inX.GetY() + m_dZ*inX.GetZ();
   }
   C3DVector Cross(const C3DVector& inX) const {
      C3DVector vReturn(m_dY*inX.GetZ() - m_dZ*inX.GetY(),
         inX.GetX()*m_dZ - m_dX*inX.GetZ(),
         m_dX*inX.GetY() - inX.GetX()*m_dY);
      return vReturn;
   }
   bool operator==(const C3DVector& inX) const {
      return ( fabs(m_dX-inX.GetX()) < 1.e-12 &&
         fabs(m_dY-inX.GetY()) < 1.e-12 &&
         fabs(m_dZ-inX.GetZ()) < 1.e-12 );
   }
   bool operator!=(const C3DVector& inX) const {
      return !( fabs(m_dX-inX.GetX()) < 1.e-12 &&
         fabs(m_dY-inX.GetY()) < 1.e-12 &&
         fabs(m_dZ-inX.GetZ()) < 1.e-12 );
   }

private:
   //helpers

private:
   //data members
   double m_dX;
   double m_dY;
   double m_dZ;
};

//Have to do scaling on the left externally (i.e. double*vec).
C3DVector operator*(const double& inScale, const C3DVector& inV) {
   return inV*inScale;
}
C3DVector operator/(const double& inScale, const C3DVector& inV) {
   return inV/inScale;
}

//Unformatted streaming.
ostream& operator<<(ostream& inOut, const C3DVector& inV) {
   inOut << inV.GetX() << " " << inV.GetY() << " " << inV.GetZ();
   return inOut;
}

#endif //3D_VECTOR
