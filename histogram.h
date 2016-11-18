//Joe Snider
//10/03
//
//Histogram class
//basically just a hash added onto a vector (with overflow)
//
//stl style is used in the interface 
//   function names (i.e. lower case and underscore)
//
//m_dLow is the low edge of the first bin
//m_dHigh is the high edge of the last bin
//there are m_iBins accessible + 2 over/underflows = m_iBins+2 total bins
//
//The axis type is always double, T specifies what's held in the bins.
//Note, can have a bin type of another histogram.
//
//todo: add reverse iterators

#ifndef MY_HISTOGRAM
#define MY_HISTOGRAM

#include <iostream>
#include <vector>

using namespace std;

//must be able to instantiate a T from 0 and have operator<<
template<class T>
class CHistogram : public vector<T> {
public:
   //typedef
   typedef vector<T> BaseType;
   typedef typename BaseType::iterator iterator;
   typedef typename BaseType::const_iterator const_iterator;

public:
   //constructors
   CHistogram(const double& inLow,
      const double& inHigh,
      const int& inBins):
   vector<T>(2, 0),
      m_dLow(inLow),
      m_dHigh(inHigh),
      m_iBins(inBins),
      m_dNegetiveLow(1.),
      m_dInvBinWidth(1.)
   {
      SetBins();
   }

   //allow default construction
   CHistogram():
   vector<T>(2, 0),
      m_dLow(0.),
      m_dHigh(1.),
      m_iBins(1),
      m_dNegetiveLow(1.),
      m_dInvBinWidth(1.)
   {
      SetBins();
   }

   //allow construction from int
   CHistogram(const int& inDummy):
   vector<T>(2, 0),
      m_dLow(0.),
      m_dHigh(0.),
      m_iBins(0),
      m_dNegetiveLow(1.),
      m_dInvBinWidth(1.)
   {
      SetBins();
   }
   
   //copy constructor
   CHistogram(const CHistogram& inCopy) {
	   *this = inCopy;
	}
	
	CHistogram& operator=(const CHistogram& inCopy) {
	   m_dLow = inCopy.GetLow();
	   m_dHigh = inCopy.GetHigh();
	   m_iBins = inCopy.GetBins();
	   SetBins();
	   const_iterator iterC = inCopy.begin();
	   const_iterator iterC_end = inCopy.end();
      iterator iter = vector<T>::begin();
	   for(; iterC != iterC_end; ++iterC, ++iter) {
		   *iter = *iterC;
		}
		return *this;
   }

   ~CHistogram() {}

public:
   //gets and sets
   //left edge of the first bin
   double GetLow() const {return m_dLow;}
   void SetLow(const double& inLow) {
      m_dLow = inLow;
      SetBins();
   }

   //right edge of the last bin
   double GetHigh() const {return m_dHigh;}
   void SetHigh(const double& inHigh) {
      m_dHigh = inHigh;
      SetBins();
   }

   int GetBins() const {return m_iBins;}
   void SetBins(const int& inBins) {
      m_iBins = inBins;
      SetBins();
   }

   void Set(const double& inLow, const double& inHigh, const int& inBins) {
      m_dLow = inLow;
      m_dHigh = inHigh;
      m_iBins = inBins;
      SetBins();
   }

   double GetBinWidth() const {return m_dBinWidth;}
   //not allowed to set the bin width

public:
   //interface (stl style functions)

   //return iterators to inX   
   iterator find(const double& inX) {
	   return vector<T>::begin() + Hash(inX);
	}
   const_iterator find(const double& inX) const {
	   return vector<T>::begin() + Hash(inX);
	}

	//subscript operator access	
	T operator[](const double& inX) const {
		return *(find(inX));
	}
	T& operator[](const double& inX) {
		return (*(find(inX)));
	}

   //the underflow bin is the second to last bin (always exists)
   iterator underflow() {
      iterator iterRet = vector<T>::begin();
      for(unsigned i = 0; i < vector<T>::size()-2; ++i, ++iterRet);
      return iterRet;
   }
   const_iterator underflow() const {
      const_iterator iterRet = vector<T>::begin();
      for(unsigned i = 0; i < vector<T>::size()-2; ++i, ++iterRet);
      return iterRet;
   }

   //the overflow bin is the last bin (always exists)
   iterator overflow() {
      iterator iterRet = vector<T>::begin();
      for(int i = 0; i < vector<T>::size()-1; ++i, ++iterRet);
      return iterRet;
   }
   const_iterator overflow() const {
      const_iterator iterRet = vector<T>::begin();
      for(int i = 0; i < vector<T>::size()-1; ++i, ++iterRet);
      return iterRet;
   }

   //allow looping the bins
   iterator bin_begin() {
      return vector<T>::begin();
   }
   const_iterator bin_begin() const {
      return vector<T>::begin();
   }
   iterator bin_end() {
      return underflow();
   }
   const_iterator bin_end() const {
      return underflow();
   }

   //allow access to the hash with a testing only warning
   int TestHash(const double& inX) const {
      cout << "Warning: testing only access to hash. " << flush;
      return Hash(inX);
   }

   ///Increment the appropriate bin for the input value.
   void increment(const double& inX) {
      vector<T>::operator[](Hash(inX)) += 1;
   }

   ///Increment the appropriate bin for the input value, by a value.
   void increment(const double& inX, const T& inY) {
      vector<T>::operator[](Hash(inX)) += inY;
   }

   //Zero the histogram including over/underflow (no rebinning)
   void zero() {
      iterator iter = vector<T>::begin();
      iterator iter_end = vector<T>::end();
      for(; iter != iter_end; ++iter) {
         *iter = 0;
      }
   }

   //Prints the bin center then value.
   void print(ostream& inOut) const {
      const_iterator iter = bin_begin();
      const_iterator iter_end = bin_end();
      for(int i = 0; iter != iter_end; ++iter, ++i) {
         inOut << m_dLow+(double(i)+.5)*m_dBinWidth << " "
            << *iter << "\n";
      }
      inOut << flush;
   }
   
   //map from iterator to bin left edge (slow for now)
   double GetBinLeftEdge(const iterator& inBin) const {
	   const_iterator iter = bin_begin();
	   int i = 0;
	   for(; iter != inBin; ++iter, ++i);
	   return m_dLow + double(i)*m_dBinWidth;
	}
   double GetBinLeftEdge(const const_iterator& inBin) const {
	   const_iterator iter = bin_begin();
	   int i = 0;
	   for(; iter != inBin; ++iter, ++i);
	   return m_dLow + double(i)*m_dBinWidth;
	}

   //map from iterator to bin center (slow for now)
   double GetBinCenter(const iterator& inBin) const {
	   const_iterator iter = bin_begin();
	   int i = 0;
	   for(; iter != inBin; ++iter, ++i);
	   return m_dLow + (double(i)+.5)*m_dBinWidth;
	}
   double GetBinCenter(const const_iterator& inBin) const {
	   const_iterator iter = bin_begin();
	   int i = 0;
	   for(; iter != inBin; ++iter, ++i);
	   return m_dLow + (double(i)+.5)*m_dBinWidth;
	}

   //map from iterator to bin right edge (slow for now)
   double GetBinRightEdge(const iterator& inBin) const {
	   const_iterator iter = bin_begin();
	   int i = 0;
	   for(; iter != inBin; ++iter, ++i);
	   return m_dLow + double(i+1)*m_dBinWidth;
	}
   double GetBinRightEdge(const const_iterator& inBin) const {
	   const_iterator iter = bin_begin();
	   int i = 0;
	   for(; iter != inBin; ++iter, ++i);
	   return m_dLow + double(i+1)*m_dBinWidth;
	}
	
	//return the number of bins
	int GetNumBins() const {return (vector<T>::size()-2);}

private:
   //helpers

   ///Update the bins to match the data.
   void SetBins() {
      if(m_iBins > 0) {
         m_dBinWidth = (m_dHigh - m_dLow) / double(m_iBins);
         vector<T>::clear();
         vector<T>::resize(m_iBins+2, 0);
      } else {
         //cerr << "Error: no bins to set bin width...defaulting to 1\n" << flush;
         m_dBinWidth = 1.;
         vector<T>::clear();
         vector<T>::resize(2, 0);
      }
      m_dNegetiveLow = -1.*m_dLow;
      m_dInvBinWidth = 1./m_dBinWidth;
   }

   ///Hash function.
   ///Returns the correct bin number (or over/underflow).
   int Hash(const double& inX) const {
      if (inX < m_dLow) {
         return m_iBins;
      } else if (inX >= m_dHigh) {
         return m_iBins+1;
      } else {
         return int((inX+m_dNegetiveLow)*m_dInvBinWidth);
      }
      return -1; //should not get here
   }

private:
   //data
   double m_dLow;
   double m_dHigh;
   int m_iBins;
   double m_dBinWidth;

   //used by the hash
   double m_dNegetiveLow;
   double m_dInvBinWidth;
};

//an output operator for the values
template<class T>
ostream& operator<<(ostream& inOut, const CHistogram<T>& inH) {
   typename CHistogram<T>::const_iterator iter = inH.bin_begin();
   typename CHistogram<T>::const_iterator iter_end = inH.bin_end();
   for(; iter != iter_end; ++iter) {
      inOut << inH.GetBinCenter(iter) << " " << *iter << "\n";
   }
   inOut << flush;
   return inOut;
}

#endif //MY_HISTOGRAM
