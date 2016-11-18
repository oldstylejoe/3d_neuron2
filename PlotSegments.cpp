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

   //get the data
   TList * listSegments = new TList();

   int iColor;
   double dX1, dY1, dZ1, dX2, dY2, dZ2;
   double dXMin = 1.e10;
   double dYMin = 1.e10;
   double dZMin = 1.e10;
   double dXMax = -1.e10;
   double dYMax = -1.e10;
   double dZMax = -1.e10;

   ifstream ifSegments("xm980804pyrAxon_segments.txt");
   while(ifSegments >> dX1 >> dY1 >> dZ1 >> dX2 >> dY2 >> dZ2 >> iColor) {
      TPolyLine3D * plTemp = new TPolyLine3D(2);
      plTemp->SetPoint(0, dX1, dY1, dZ1);
      plTemp->SetPoint(1, dX2, dY2, dZ2);
      plTemp->SetLineColor(iColor+2);
      listSegments->Add(plTemp);

      dXMin = TMath::Min(dX1, TMath::Min(dX2, dXMin));
      dYMin = TMath::Min(dY1, TMath::Min(dY2, dYMin));
      dZMin = TMath::Min(dZ1, TMath::Min(dZ2, dZMin));
      dXMax = TMath::Max(dX1, TMath::Max(dX2, dXMax));
      dYMax = TMath::Max(dY1, TMath::Max(dY2, dYMax));
      dZMax = TMath::Max(dZ1, TMath::Max(dZ2, dZMax));
   }

   cout << dXMin << " " << dXMax << "\n";
   cout << dYMin << " " << dYMax << "\n";
   cout << dZMin << " " << dZMax << "\n" << flush;

   //draw the data
   TCanvas * pCanvas = new TCanvas("pCanvas", "", 1);
   TView * viewCanvas = new TView(1);
   viewCanvas->SetRange(dXMin - 0.01*TMath::Abs(dXMin),
      dYMin - 0.01*TMath::Abs(dYMin),
      dZMin - 0.01*TMath::Abs(dZMin),
      dXMax + 0.01*TMath::Abs(dXMax),
      dYMax + 0.01*TMath::Abs(dYMax),
      dZMax + 0.01*TMath::Abs(dZMax));
   viewCanvas->SetOutlineToCube();

   listSegments->Draw();

   pCanvas->Update();
   pCanvas->Draw();
}
