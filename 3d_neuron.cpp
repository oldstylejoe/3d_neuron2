//Joe Snider
//7/04
//
//Test file for reading in 3d neuron data.

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <string>
#include <algorithm>

#include "3dBox.h"
#include "3dvector.h"
#include "CNeuron.h"
#include "Branch.h"
#include "rand3.h"

const int NUM_SUB_ARBORS = 10;

using namespace std;

//Used to convert strings to lower case.
char ConvertLower(char c) {return tolower(c);}

//Do linear regression to find an offset on std::map type iterators.
//The input types must be const forward iterators with (* ).first and .second (like std::map::const_iterator).
//Modifies inOffset and inOffsetError (to the offset and error of the best fit line).
//Note that this assumes a value (probably -1 or -3) of the slope, and does not allow that to be changed.
template<class T>
void LinearRegression(const T& inIterLow, const T& inIterHigh, const double& inSlope, double& inOffset, double& inOffsetError) {
   //find the means
   double dMeanX = 0.;
   double dMeanY = 0.;
   double dCount = 0.;
   T iter = inIterLow;
   T iter_end = inIterHigh;
   for(; iter != iter_end; ++iter) {
//      cout << "3d_neuron.cpp gh1 " << iter->first << " " << iter->second << "\n" << flush;
      dMeanX += iter->first;
      dMeanY += iter->second;
      dCount += 1.;
   }
   dMeanX /= dCount;
   dMeanY /= dCount;
   inOffset = dMeanY - inSlope*dMeanX;

   //find the error
   double dError2X = 0.;
   double dError2Y = 0.;
   iter = inIterLow;
   for(; iter != iter_end; ++iter) {
      dError2X += (iter->first - dMeanX)*(iter->first - dMeanX);
      dError2Y += (iter->second - dMeanY)*(iter->second - dMeanY);
   }
   inOffsetError = sqrt(inSlope*inSlope*dError2X+dError2Y)/dCount;
}

int main(int argc, char* argv[])
{
   srand( (unsigned) time(NULL) );
   CRand3 rand3( rand() );

   vector<string> vecFiles;
   ifstream ifFiles("files.txt");
   if(ifFiles.good()) {
      string strTemp;
      while(ifFiles >> strTemp) {
         if(strTemp[0] != '#') {
            vecFiles.push_back(strTemp);
         }
      }
   } else {
      cerr << "Error: enable to open file list (files.txt)...exiting\n" << flush;
      return -1;
   }
   ifFiles.close();

   vector<double> vecTotalLength;
   vector<double> vecInertialVolume;

   vector<string> vecTypes;
   //vecTypes.push_back("Axon");
   //vecTypes.push_back("Dendrite");
   //vecTypes.push_back("Apical");
   vecTypes.push_back("");


   ofstream ofSSGlobal("self_similar.txt");
   ofstream ofKurtosis("kurtosis.txt");
   ofstream ofKurtosis2D("kurtosis_2d.txt");
   ofstream ofKurtosisDirectional("kurtosis_directional.txt");
   ofstream ofIntersection("intersections.txt");
   ofstream ofTotalLengths("total_lengths.txt");
   ofstream ofBranchAngles("branch_angles.txt");
   ofstream ofBranchLengths("branch_lengths.txt");
   ofstream ofWidths("widths.txt");
   ofstream ofProductMoment("product_moment.txt");
   ofstream ofSumMoment("sum_moment.txt");
   ofstream ofSingleMoment_x("single_moment_x.txt");
   ofstream ofSingleMoment_y("single_moment_y.txt");
   ofstream ofSingleMoment_z("single_moment_z.txt");
   ofstream ofFilesDistribution("files_distribution.txt");
   ofstream ofFilesLengthDistribution("files_length_distribution.txt");
   ofstream ofFilesBranchCount("files_branch_count.txt");
   ofstream ofFilesGlobalAngles("files_global_angles.txt");
   ofstream ofAngleFiles("files_angles.txt");
   ofstream ofLengthFiles("files_lengths.txt");
   ofstream ofAllGlobalAngles("global_angles.txt");
   //loop over the files
   for(unsigned i = 0; i < vecFiles.size(); ++i) {
      for(int j = 0; j < vecTypes.size(); ++j) {
         CNeuron N1;
         string strType = vecTypes[j];
         int iType = 2+j;

         //get the file name
         cout << "Parsing " << vecFiles[i].c_str() << " ..." << flush;
         ifstream ifIn(vecFiles[i].c_str());
         if(ifIn.good()) {
            string tempFileName = vecFiles[i];
            transform(tempFileName.begin(), tempFileName.end(), tempFileName.begin(), ConvertLower);
            if(tempFileName.compare(tempFileName.length()-3, 3, "swc") == 0) {
               if(iType == 2 || iType == 3 || iType == 4) {
                  N1.ReadSWCData(ifIn, iType);
               } else {
                  cerr << "Unknown type...ignoring...\n" << flush;
               }
            } else if(tempFileName.compare(tempFileName.length()-3, 3, "asc") == 0 ||
               tempFileName.compare(tempFileName.length()-7, 7, "asc.txt") == 0) {
               N1.ReadData(ifIn, strType);
            } else {
               cerr << "Unknown extension...ignoring...\n" << flush;
            }
         } else {
            cerr << "Unable to open the input file\n" << flush;
         }
         ifIn.close();
         cout << "done\n" << flush;

         ofstream ofSegments("segments.txt");
         N1.DisplaySegments(ofSegments);
         ofSegments.close();
         ofstream ofNormalizedSegments("normalized_segments.txt");
         N1.DisplayNormalizedSegments(ofNormalizedSegments);
         ofNormalizedSegments.close();

         //check that the neuron read something
         if(N1.GetBranch().size() > 0 && N1.MaxLevel() > 2) {
            ////////////////////////////////////////////////begin hack/////////////////////////////////
            ////Part of the sub arbor scaling hack
            //CPoint pTrimPoint;
            //pTrimPoint.push_back( vecTrimX[(i/2)%6] );
            //pTrimPoint.push_back( vecTrimY[(i/2)%6] );
            //pTrimPoint.push_back( 0. );
            //N1.Trim(pTrimPoint);
            //N1.TrimToFirstDaughter();
            ///////////////////////////////////////////////end hack/////////////////////////////////////

            /*//save to a file
            string strSaveName;
            strSaveName.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strSaveName.append(strType);
            strSaveName.append("_clean.swc");
            ofstream ofWriteTest(strSaveName.c_str());
            N1.WriteSWCData(ofWriteTest, iType);
            ofWriteTest.close();*/

            //ofstream ofPoints("points.txt");
            //N1.UnformattedDisplayNeuron(ofPoints);
            //ofPoints.close();

            //for(int iBoxRepeat = 0; iBoxRepeat < 10; ++iBoxRepeat) {
            //   //set up the output file name
            //   string strOutput;
            //   strOutput.append(vecFiles[i], 0, vecFiles[i].length()-4);
            //   strOutput.append(strType);
            //   strOutput.append("_box_count");
            //   char cstrBufferNumber[100];
            //   sprintf(cstrBufferNumber, "%d", iBoxRepeat);
            //   strOutput.append(cstrBufferNumber);
            //   strOutput.append(".txt");

            //   //do the box counting
            //   ofstream ofOut(strOutput.c_str());
            //   cout << "Doing the box counting..." << flush;
            //   N1.BoxCount2(ofOut, 1, 1, 1000, 10);
            //   cout << "done\nData saved to " << strOutput << "\n" << flush;
            //   ofOut.close();
            //}

            /*//do the check self-similarity test
            cout << "Doing the self similar test..." << flush;
            N1.CheckSelfSimilarMoments(ofSSGlobal,100);
            ofSSGlobal << flush;
            string strOutputSS;
            strOutputSS.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strOutputSS.append(strType);
            strOutputSS.append("_self_similar.txt");
            ofstream ofSS(strOutputSS.c_str());
            N1.CheckSelfSimilar(ofSS, 10);
            ofSS.close();
            cout << "done.\n" << flush;*/

            /*//set up the output file name
            string strOutputMD;
            strOutputMD.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strOutputMD.append(strType);
            strOutputMD.append("_distance_minimum.txt");
            //do the minimum distance calc
            ofstream ofOutMD(strOutputMD.c_str());
            cout << "Doing the minimum distance calculation..." << flush;
            N1.AverageMinimumDistance(ofOutMD, 100, 100);
            cout << "done\nData saved to " << strOutputMD << "\n" << flush;
            ofOutMD.close();*/

            //do the distribution calculation
            string strOutputLD;
            cout << "Doing the length distribution calculation..." << flush;
            strOutputLD.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strOutputLD.append(strType);
            strOutputLD.append("_length_distribution.txt");
            ofstream ofOutLD(strOutputLD.c_str());
            N1.LengthDistribution(ofOutLD, 100);
            ofOutLD.close();
            ofFilesLengthDistribution << strOutputLD << "\n" << flush;
            cout << "done\nData saved to " << strOutputLD << "\n" << flush;

            /*//set up the kurtosis distribution output file name
            string strOutputKD;
            strOutputKD.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strOutputKD.append(strType);
            strOutputKD.append("_kurtosis_distribution.txt");
            //do the kurtosis calculation
            cout << "Doing the kurtosis calculation..." << flush;
            ofstream ofKurtosisDistribution(strOutputKD.c_str());
            ofKurtosis << vecFiles[i].c_str() << " " << N1.Kurtosis(ofKurtosisDistribution) << "\n";
            ofKurtosisDistribution.close();
            cout << "done.\n" << flush;*/

            double dTotalLength = N1.TotalLength();
            double dInertialVolume = N1.InertialVolume();
            double dQhullVolume = N1.QhullVolume();
            string strTotalLengthID;
            strTotalLengthID.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strTotalLengthID.append(strType);
            ofTotalLengths << strTotalLengthID.c_str() << " " 
               << dTotalLength << " "
               << dInertialVolume << " "
               << dQhullVolume << "\n" << flush;
            vecTotalLength.push_back(dTotalLength);
            vecInertialVolume.push_back(dInertialVolume);

            //dump the count at radii
            cout << "Finding and dumping the branch counts..." << flush;
            string strOutputBC;
            strOutputBC.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strOutputBC.append(strType);
            strOutputBC.append("_branch_count.txt");
            //do the minimum distance calc
            ofstream ofOutBC(strOutputBC.c_str());
            N1.DumpBranchCount(ofOutBC);
            ofFilesBranchCount << strOutputBC << "\n" << flush;
            cout << "done\n" << flush;
            ofOutBC.close();

            //dump the angle at radii
            cout << "Finding and dumping the global angles..." << flush;
            string strOutputAG;
            strOutputAG.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strOutputAG.append(strType);
            strOutputAG.append("_angle_global.txt");
            ofstream ofOutAG(strOutputAG.c_str());
            N1.DumpCOMAngle(ofOutAG, 100);
            ofFilesGlobalAngles << strOutputAG << "\n" << flush;
            cout << "done\n" << flush;
            ofOutAG.close();

            /*//dump the 2d distribution
            cout << "Finding and dumping the 2d distribution..." << flush;
            string strOutputDist;
            strOutputDist.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strOutputDist.append(strType);
            strOutputDist.append("_distribution.txt");
            ofFilesDistribution << strOutputDist << "\n" << flush;
            ofstream ofDistribution(strOutputDist.c_str());
            ofKurtosis2D << vecFiles[i].c_str() << " " << N1.Distribution2D(ofDistribution, 200) << "\n";
            ofDistribution.close();
            cout << "done.\n" << flush;*/

            /*//dump the 3d distribution
            cout << "Finding and dumping the 3d distribution..." << flush;
            string strOutputDist;
            strOutputDist.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strOutputDist.append(strType);
            strOutputDist.append("_distribution.txt");
            ofFilesDistribution << strOutputDist << "\n" << flush;
            ofstream ofDistribution(strOutputDist.c_str());
            N1.Distribution2D(ofDistribution, 40);
            ofDistribution.close();
            cout << "done.\n" << flush;*/

            /*//Do the variance in various directions calculation.
            cout << "Finding and dumping the distribution..." << flush;
            string strOutputDirectional;
            strOutputDirectional.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strOutputDirectional.append(strType);
            strOutputDirectional.append("_directional.txt");
            ofstream ofDirectional(strOutputDirectional.c_str());
            ofDirectional.close();
            cout << "done.\n" << flush;*/

            //Save the widths
            string strWidthID;
            strWidthID.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strWidthID.append(strType);
            ofWidths << strWidthID << " "
               << N1.GetWidthX() << " "
               << N1.GetWidthY() << " "
               << N1.GetWidthZ() << " "
               << N1.GetStandardDeviationX() << " "
               << N1.GetStandardDeviationY() << " "
               << N1.GetStandardDeviationZ() << "\n" << flush;

            /*//Find the product moment
            cout << "Finding the product moment..." << flush;
            string strPMID;
            strPMID.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strPMID.append(strType);
            ofProductMoment << strPMID << " ";
            N1.ProductMoment(ofProductMoment, 0, 20, 1);
            ofProductMoment << flush;
            cout << "done\n" << flush;*/

            /*//Find the sum moment
            cout << "Finding the sum moments..." << flush;
            string strSMID;
            strSMID.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strSMID.append(strType);
            ofSumMoment << strSMID << " ";
            N1.SumMoment(ofSumMoment, 0, 20, 1);
            ofSumMoment << flush;
            cout << "done\n" << flush;*/

            /*//Find the single moment
            cout << "Finding the single moments...x..." << flush;
            string strSingleMID_x;
            strSingleMID_x.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strSingleMID_x.append(strType);
            ofSingleMoment_x << strSingleMID_x << " ";
            N1.SingleMoment(ofSingleMoment_x, 0, 20, 1, 0);
            ofSingleMoment_x << flush;
            cout << "y..." << flush;
            string strSingleMID_y;
            strSingleMID_y.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strSingleMID_y.append(strType);
            ofSingleMoment_y << strSingleMID_y << " ";
            N1.SingleMoment(ofSingleMoment_y, 0, 20, 1, 1);
            ofSingleMoment_y << flush;
            cout << "z..." << flush;
            string strSingleMID_z;
            strSingleMID_z.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strSingleMID_z.append(strType);
            ofSingleMoment_z << strSingleMID_z << " ";
            N1.SingleMoment(ofSingleMoment_z, 0, 20, 1, 2);
            ofSingleMoment_z << flush;
            cout << "done\n" << flush;*/

            //Do the branch angle calculation.
            cout << "Finding and dumping the branch angles..." << flush;
            //string strOutputBA;
            //strOutputBA.append(vecFiles[i], 0, vecFiles[i].length()-4);
            //strOutputBA.append(strType);
            //strOutputBA.append("_angles.txt");
            //ofstream ofBA(strOutputBA.c_str());
            //N1.DumpAngles(ofBA);
            N1.DumpAngles(ofBranchAngles);
            ofBranchAngles << flush;
            //ofBA.close();
            //ofAngleFiles << strOutputBA.c_str() << "\n" << flush;
            cout << "done.\n" << flush;

            /*//Do the branch length calculation.
            cout << "Finding and dumping the branch lengths..." << flush;
            string strOutputBL;
            strOutputBL.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strOutputBL.append(strType);
            strOutputBL.append("_lengths.txt");
            ofstream ofBL(strOutputBL.c_str());
            N1.DumpLengths(ofBL);
            N1.DumpLengths(ofBranchLengths);
            ofBL.close();
            ofLengthFiles << strOutputBL.c_str() << "\n" << flush;
            cout << "done.\n" << flush;*/

            /*//set up the output file name
            string strOutputIH;
            strOutputIH.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strOutputIH.append(strType);
            strOutputIH.append("_inertial_histogram.txt");
            //do the inertial histogram calculation.
            ofstream ofOutIH(strOutputIH.c_str());
            cout << "Doing the inertial histogram calculation..." << flush;
            N1.InertialHistogram(ofOutIH);
            cout << "done\nData saved to " << strOutputIH << "\n" << flush;
            ofOutIH.close();*/

            //for(int j = 0; j <= i; ++j) {
            //   cout << "   Intersecting with " << vecFiles[j].c_str() << "..." << flush;
            //   CNeuron N2;
            //   ifstream ifIn2(vecFiles[j].c_str());
            //   if(ifIn2.good()) {
            //      N2.ReadData(ifIn2);
            //      double dIntersection, dIntersectionError;
            //      N1.CountIntersections(N2, dIntersection, dIntersectionError);
            //      ofIntersection << vecFiles[i].c_str() << " "
            //         << vecFiles[j].c_str() << " "
            //         << dIntersection << " "
            //         << dIntersectionError << " "
            //         << vecTotalLength[j] << " " 
            //         << vecInertialVolume[j] << " "
            //         << vecTotalLength[i] << " " 
            //         << vecInertialVolume[i] << "\n" << flush;
            //   } else {
            //      cerr << "Unable to open file...skipping..." << flush;
            //   }
            //   cout << "done\n" << flush;
            //   ifIn2.close();
            //}

            ////do the branch length calculation
            //N1.DisplayBL_versus_R(ofBranchLengths);

            //display the neuron in a format the can be ploted in a drawing program
            string strSegments;
            strSegments.append(vecFiles[i], 0, vecFiles[i].length()-4);
            strSegments.append(strType);
            strSegments.append("_segments.txt");
            ofstream ofSegments2(strSegments.c_str());
            //N1.DisplaySegments(ofSegments2);
            N1.DisplayNormalizedSegments(ofSegments2);
            ofSegments2.close();

            /*//Merge and dump sub arbors
            cout << "Merging arbors ... " << flush;
            N1.MergeArbors();
            cout << "done.\n" << flush;
            cout << "Saving subarbors ... " << flush;
            for(int iSA = 0; iSA < NUM_SUB_ARBORS; ++iSA) {
               string strSubArbors;
               strSubArbors.append(vecFiles[i], 0, vecFiles[i].length()-4);
               strSubArbors.append(strType);
               strSubArbors.append("_subarbor");
               char cstrBufferSubString[1000]; //plenty big
               sprintf(cstrBufferSubString, "%d", iSA);
               strSubArbors.append(cstrBufferSubString);
               strSubArbors.append(".my.asc");
               ofstream ofSA(strSubArbors.c_str());
               N1.SubArbor(ofSA, strType, 400.);
               ofSA.close();
            }
            cout << "done\n" << flush;*/

         }
      }
   }
   ofLengthFiles.close();
   ofAngleFiles.close();
   ofFilesLengthDistribution.close();
   ofFilesDistribution.close();
   ofSingleMoment_z.close();
   ofSingleMoment_y.close();
   ofSingleMoment_x.close();
   ofSumMoment.close();
   ofProductMoment.close();
   ofWidths.close();
   ofBranchLengths.close();
   ofTotalLengths.close();
   ofIntersection.close();
   ofKurtosis.close();
   ofKurtosis2D.close();
   ofKurtosisDirectional.close();
   ofSSGlobal.close();

	return 0;
}

