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

   char cstrFileName[] = "normalized_segments.txt";

   int iColor;
   double dX1, dY1, dZ1, dX2, dY2, dZ2;
   double dXMin = 1.e10;
   double dYMin = 1.e10;
   double dXMax = -1.e10;
   double dYMax = -1.e10;

   ifstream ifSegments(cstrFileName);
   while(ifSegments >> dX1 >> dY1 >> dZ1 >> dX2 >> dY2 >> dZ2 >> iColor) {
      dXMin = TMath::Min(dX1, TMath::Min(dX2, dXMin));
      dYMin = TMath::Min(dY1, TMath::Min(dY2, dYMin));
      dXMax = TMath::Max(dX1, TMath::Max(dX2, dXMax));
      dYMax = TMath::Max(dY1, TMath::Max(dY2, dYMax));
   }
   ifSegments.close();

   //draw the data
   dXMin -= 10.;
   dYMin -= 10.;
   dXMax += 10.;
   dYMax += 10.;
   TCanvas * pCanvas = new TCanvas("pCanvas", cstrFileName, 500, 500);

   //TH2D * h2Axis = new TH2D("h2Axis", "", 10, dXMin-10., dXMax+10., 10, dYMin-10., dYMax+10.);
   TH2D * h2Axis = new TH2D("h2Axis", "", 10, -20., 20., 10, -20., 20.);
   h2Axis->Draw("");

   ifstream ifSegments2(cstrFileName);
   TLine lineDrawer;
   lineDrawer.SetLineWidth(3);
   while(ifSegments2 >> dX1 >> dY1 >> dZ1 >> dX2 >> dY2 >> dZ2 >> iColor) {
      lineDrawer.DrawLine(dX1, dY1, dX2, dY2);
      //cout << dX1 << " " << dY2 << " " << dX2 << " " << dY2 << "\n";
   }
   lineDrawer.DrawLine(dXMin+5, dYMax+5, dXMin+105, dYMax+5);
   ifSegments2.close();

   ifstream ifAngles("temp2.txt");
   lineDrawer.SetLineColor(2);
   while(ifAngles >> dX1 >> dY1 >> dX2 >> dY2) {
      lineDrawer.DrawLine(dX1, dY1, dX2, dY2);
      //cout << dX1 << " " << dY2 << " " << dX2 << " " << dY2 << "\n";
   }

   pCanvas->Update();
   pCanvas->Draw();
}
