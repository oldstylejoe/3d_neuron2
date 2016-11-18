{
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetPadBorderMode(2);
   gStyle->SetPadColor(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetTitleColor(0);
   gStyle->SetStatColor(0);

   #include <fstream.h>
   #include <vector.h>

   THStack * hsData = new THStack("hsData", "");

   int ipColor[] = {2,4,5};
   int ipMarker[] = {21,22,23};

   vector<char*> vecTypes;
   vecTypes.push_back("");
   //vecTypes.push_back("Dendrite");
   //vecTypes.push_back("Apical");

   for(int i = 0; i < vecTypes.size(); ++i) {
      TString strDataName("Data");
      TString strCountName("Count");
      TH1D * h1dData = new TH1D("h1dData", "", 40, -400., 400.);
      //TH1D * h1dData = new TH1D("h1dData", "", 100, 0., 100.);
      h1dData->SetName(strDataName+vecTypes[i]);

      //get the data
      ifstream ifFiles("files_angles.txt");
      //ifstream ifFiles("files_lengths.txt");
      while(ifFiles) {
         TString strTemp;
         strTemp.ReadLine(ifFiles);
         if(strTemp.Contains(vecTypes[i]) && strTemp.Length() > 1) {
            cout << "Reading from " << (char*)strTemp << " ..." << flush;

            int iCount = 0;
            double dX, dY, dYE;
            ifstream ifIn((char*)strTemp);
            while(ifIn >> dX) {
               h1dData->Fill(dX*360./(2.*acos(-1.)), 1.);
            }
            ifIn.close();
            cout << "done\n" << flush;
         }
      }
      ifFiles.close();

      h1dData->Sumw2();
      h1dData->Scale(1./h1dData->Integral());
      h1dData->SetMarkerColor(ipColor[i]);
      h1dData->SetMarkerStyle(ipMarker[i]);
      h1dData->SetLineColor(ipColor[i]);

      //TString strFitName("Fitter");
      //TF1 * fFitter = new TF1(strFitName+vecTypes[i], "[0]*exp(-1.*(x/[1])**6.)", 0., 2.);
      //fFitter->SetParameters(0.2, 0.3);
      //fFitter->SetLineColor(ipColor[i]);
      //h1dData->Fit(fFitter, "r+");

      hsData->Add(h1dData, "");

      TString strAverage("average_angles");
      strAverage += vecTypes[i];
      strAverage += ".txt";
      ofstream ofAverage((char*)strAverage);
      for(int j = 1; j <= h1dData->GetNbinsX(); ++j) {
         ofAverage << h1dData->GetBinCenter(j) << " "
            << h1dData->GetBinContent(j) << " "
            << h1dData->GetBinError(j) << "\n";
      }
      ofAverage.close();
   }

   //draw the data
   TCanvas * pCanvas = new TCanvas("pCanvas", "", 1);

   hsData->Draw("b nostack");
   hsData->GetXaxis()->SetTitle("Normalized radius (unitless)");
   hsData->GetXaxis()->SetTitleColor(1);
   hsData->GetXaxis()->CenterTitle();
   hsData->GetYaxis()->SetTitle("Probability");
   hsData->GetYaxis()->CenterTitle();

   pCanvas->Update();
   pCanvas->Draw();
}
