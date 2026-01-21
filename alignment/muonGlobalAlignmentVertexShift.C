#include <MFTTracking/Constants.h>

TFile* fAnalysisResults;

TH1* histDxVsY;

double scaleME = 1;

constexpr int nPoints = 4;

TH1* GetTH1(TFile* f, TString histname)
{
  return (TH1*)f->Get(histname);

  //TString histname = TString::Format("ST%d/DE%d/Occupancy_B_XY_%d", station, de, de);
  TKey *key = f->GetKey(histname);
  std::cout << "histname: " << histname << "  key: " <<key << std::endl;
  if (!key) return NULL;
  return (TH1*)key->ReadObjectAny(TH1::Class());
}


TH2* GetTH2(TFile* f, TString histname)
{
  return (TH2*)f->Get(histname);

  //TString histname = TString::Format("ST%d/DE%d/Occupancy_B_XY_%d", station, de, de);
  TKey *key = f->GetKey(histname);
  std::cout << "histname: " << histname << "  key: " <<key << std::endl;
  if (!key) return NULL;
  return (TH2*)key->ReadObjectAny(TH2::Class());
}


TH3* GetTH3(TFile* f, TString histname)
{
  return (TH3*)f->Get(histname);

  //TString histname = TString::Format("ST%d/DE%d/Occupancy_B_XY_%d", station, de, de);
  TKey *key = f->GetKey(histname);
  std::cout << "histname: " << histname << "  key: " <<key << std::endl;
  if (!key) return NULL;
  return (TH3*)key->ReadObjectAny(TH3::Class());
}

int GetLR(int chamber, int deIndex)
{
  if (deIndex < 0)
    return -1;

  if (chamber >= 1 && chamber <= 4) {
    if (deIndex >= 1 && deIndex <= 2)
      return 0;
    else if (deIndex <= 3)
      return 1;
  } else if (chamber <= 6) {
    if (deIndex >= 5 && deIndex <= 13)
      return 0;
    else if (deIndex <= 17)
      return 1;
  } else {
    if (deIndex >= 7 && deIndex <= 19)
      return 0;
    else if (deIndex <= 25)
      return 1;
  }

  return -1;
}

Double_t DoubleSidedCB2(double x, double mu, double width, double a1, double p1, double a2, double p2)
{
  double u   = (x-mu)/width;
  double A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
  double A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
  double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
  double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

  double result(1);
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u,-p2);
  return result;
}

double DoubleSidedCB(double* x, double *par)
{
  return(par[0] * DoubleSidedCB2(x[0], par[1],par[2],par[3],par[4],par[5],par[6]));
}

double DoubleSidedCBwithLinBgd(double* x, double *par)
{
  return(par[0] * DoubleSidedCB2(x[0], par[1],par[2],par[3],par[4],par[5],par[6]) + par[7] + x[0] * par[8]);
}

TH2* PlotDCAMFT(std::string histName)
{
  std::string fullHistName = std::string("qa-muon/alignment/") + histName;
  TH2* histogram = GetTH2(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogram << std::endl;
  if (!histogram)
    return nullptr;

  histogram->GetYaxis()->SetRangeUser(-0.1, 0.1);
  histogram->Draw("col");

  return histogram;
}

double PhiFitFunc(double *x, double *par) {
  double phi = x[0] * TMath::Pi() / 180.f;
  double fitval = (phi > 0) ? par[0]*TMath::Sin(phi + par[1]) + par[2] : par[0]*TMath::Sin(phi + par[1]) + par[3];
  return fitval;
}

std::pair<double, double> PlotDCAPhiProjection(TH2* histogram2, float yMin, float yMax, int projRebin, TCanvas& c, double phaseMin, double phaseMax, bool printFits = false)
{
  std::string fullHistName = histogram2->GetName();
  TH1* histogramMean = histogram2->ProjectionX();
  histogramMean->SetName(TString::Format("%s-mean", fullHistName.c_str()));
  histogramMean->SetTitle(TString::Format("%s", histogram2->GetTitle()));
  //histogramMean->SetMinimum(yMin);
  //histogramMean->SetMaximum(yMax);

  int entriesMin = histogram2->GetEntries() / 20;
  int peakMin = 10;
  // skip first bin because it is biased?
  float meanMin = 1000000;
  float meanMax = -1000000;
  for (int bin = 1; bin <= histogram2->GetXaxis()->GetNbins(); bin++) {
    TH1* proj = histogram2->ProjectionY(TString::Format("%s-%d", fullHistName.c_str(), bin), bin, bin);
    proj->SetTitle(TString::Format("%s - bin %d", histogram2->GetTitle(), bin));
    proj->Rebin(projRebin);
    proj->SetLineColor(kRed);

    double mean = -10000;
    double meanErr = 0;
    double sigma = -10000;
    double sigmaErr = 0;

    int valuePeak = proj->GetMaximum();
    int binPeak = proj->GetMaximumBin();
    double xPeak = proj->GetXaxis()->GetBinCenter(binPeak);

    std::cout << std::format("\"{}\" -> peak = {}", proj->GetName(), valuePeak) << std::endl;

    if (valuePeak > peakMin) {

      if (printFits) {
        //std::cout << std::format("[\"{}\"] bin {} max={:0.2f}  binPeak={}  xPeak={:0.2f}\n", fullHistName, bin, valuePeak, binPeak, xPeak);
        //std::cout << std::format("[\"{}\"] bin {}  max={:0.2f}  binPeak={}  xPeak={:0.2f}",
        //    fullHistName, bin, valuePeak, binPeak, xPeak) << std::endl;
      }
      //std::cout << std::format("Bin {}  max={}  binPeak={}  xPeak={}\n", bin, valuePeak, binPeak, xPeak);

      TF1 fcb("fcb", DoubleSidedCB, xPeak - 0.1, xPeak + 0.1, 7);
      fcb.SetNpx(1000);
      fcb.SetLineColor(kBlack);
      //fcb.SetParameter(0, valuePeak);
      double par[7];
      par[0]=35000;
      par[1]=xPeak;
      par[2]=0.01;
      par[3]=1;
      par[4]=1;
      par[5]=1;
      par[6]=1;
      fcb.SetParameters(&par[0]);

      //fcb.FixParameter(1, xPeak);
      fcb.SetParLimits(2, 0.01, 0.1);
      //proj->Fit("fcb", "BRN");
      //fcb.ReleaseParameter(1);

      if (printFits) {
        c.cd(3);
        proj->Draw("E");
        //projME->Draw("E same");
        c.cd(4);
      }

      TH1* projCorr = (TH1*)proj->Clone(TString::Format("%s-%d-copy", fullHistName.c_str(), bin));
      //projCorr->Add(projME, -1);

      if (printFits)
        projCorr->Fit("fcb", "BR");
      else
        projCorr->Fit("fcb", "BRQ");

      if (printFits) {
        projCorr->Draw("E");
        c.SaveAs("residuals_AO2D.pdf");
      }

      mean = fcb.GetParameter(1);
      meanErr = fcb.GetParError(1);
      sigma = fcb.GetParameter(2);
      sigmaErr = fcb.GetParError(2);

      if (mean < meanMin) meanMin = mean;
      if (mean > meanMax) meanMax = mean;
    }

    histogramMean->SetBinContent(bin, mean);
    histogramMean->SetBinError(bin, meanErr);
  }

  TF1 phiFit("phiFit", PhiFitFunc, -180.f, 180.f, 4);
  phiFit.SetLineColor(kRed);
  phiFit.SetLineStyle(kDotted);
  //phiFit.FixParameter(0, 0.02);
  phiFit.FixParameter(1, (phaseMax + phaseMin) / 2.f);
  phiFit.SetParLimits(1, phaseMin, phaseMax);
  phiFit.SetParameter(2, 0);
  phiFit.SetParameter(3, 0);
  phiFit.FixParameter(2, 0);
  phiFit.FixParameter(3, 0);
  //histogramMean->Fit("phiFit", "RN");
  //phiFit.ReleaseParameter(0);
  //phiFit.ReleaseParameter(2);
  //phiFit.ReleaseParameter(3);
  histogramMean->Fit("phiFit", "R");
  histogramMean->SetMinimum(meanMin);
  histogramMean->SetMaximum(meanMax);
  histogramMean->Draw("E");

  TLine* line1 = new TLine(histogramMean->GetXaxis()->GetXmin(), 0, histogramMean->GetXaxis()->GetXmax(), 0);
  line1->SetLineStyle(kDashed);
  line1->Draw();

  TLine* line2 = new TLine(0, yMin, 0, yMax);
  line2->SetLineStyle(kDashed);
  line2->Draw();

  std::cout << std::format("Amplitude: {:0.4f}", phiFit.GetParameter(0)) << std::endl;
  std::cout << std::format("Phase:     {:0.4f}", phiFit.GetParameter(1)) << std::endl;
  std::cout << std::format("Offsets:   top={:0.4f} bottom={:0.4f}", phiFit.GetParameter(2), phiFit.GetParameter(3)) << std::endl;

  return {phiFit.GetParameter(0), phiFit.GetParError(0)};
}

std::pair<double, double> PlotDCAPhiProjection(std::string histName, float yMin, float yMax, int projRebin, TCanvas& c, double phaseMin, double phaseMax, bool printFits = false)
{
  std::string fullHistName = std::string("qa-muon/alignment/") + histName;
  TH2* histogram2 = GetTH2(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogram2 << std::endl;
  if (!histogram2)
    return {};
  histogram2->SetName(fullHistName.c_str());

  return PlotDCAPhiProjection(histogram2, yMin, yMax, projRebin, c, phaseMin, phaseMax, printFits);
}

void PlotDCAPhiProjection3D(std::string histName, float yMin, float yMax, int projRebin, TCanvas& c, double phaseMin, double phaseMax, bool printFits = false)
{
  std::string fullHistName = std::string("muon-global-alignment/") + histName;
  TH3* histogram3 = GetTH3(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogram3 << std::endl;
  if (!histogram3)
    return;

  std::vector<double> xv;
  std::vector<double> yv;
  std::vector<double> exv;
  std::vector<double> eyv;

  for (int bin = 1; bin <= histogram3->GetXaxis()->GetNbins(); bin++) {
    //if (bin == 4) continue;
    histogram3->GetXaxis()->SetRange(bin, bin);
    TH2* proj = (TH2*)histogram3->Project3D("zy");
    proj->SetTitle(TString::Format("%s - z bin %d", histogram3->GetTitle(), bin));
    c.Clear();
    proj->Draw("col");
    c.SaveAs("residuals_AO2D.pdf");

    c.Clear();
    auto amplitude = PlotDCAPhiProjection(proj, yMin, yMax, projRebin, c, phaseMin, phaseMax, printFits);
    c.SaveAs("residuals_AO2D.pdf");

    //gr->AddPointError(histogram3->GetXaxis()->GetBinCenter(bin), amplitude.first, 0, amplitude.second);
    std::cout << std::format("Bin #{} -> {}", bin, histogram3->GetXaxis()->GetBinCenter(bin)) << std::endl;
    xv.push_back(histogram3->GetXaxis()->GetBinCenter(bin));
    exv.push_back(0);
    yv.push_back(amplitude.first);
    eyv.push_back(amplitude.second);

  }
  TGraphErrors* gr = new TGraphErrors(xv.size(), xv.data(), yv.data(), exv.data(), eyv.data());
  gr->Draw("A*");
  gr->SetTitle("Modulation amplitude vs. z shift");
  gr->GetXaxis()->SetTitle("z shift (mm)");
  gr->GetYaxis()->SetTitle("amplitude");

  TLine* line1 = new TLine(gr->GetHistogram()->GetXaxis()->GetXmin(), 0, gr->GetHistogram()->GetXaxis()->GetXmax(), 0);
  line1->SetLineStyle(kDashed);
  line1->Draw();

  TLine* line2 = new TLine(0, gr->GetHistogram()->GetMinimum(), 0, gr->GetHistogram()->GetMaximum());
  line2->SetLineStyle(kDashed);
  line2->Draw();

  TF1 linFit("linFit", "pol1");
  linFit.SetLineColor(kRed);
  linFit.SetLineStyle(kDotted);
  gr->Fit("linFit", "B");

  std::cout << std::format("Optimal z shift: {} mm", -linFit.GetParameter(0) / linFit.GetParameter(1)) << std::endl;
}

void muonGlobalAlignmentVertexShift()
{
  fAnalysisResults = new TFile("AnalysisResults.root");
  //fAnalysisResults = new TFile("AnalysisResults/AnalysisResultsFull.root");

  std::array<std::string, 4> quadrants = {"Q0", "Q1", "Q2", "Q3"};
  //std::array<std::string, 1> quadrants = {"Q0"};


  std::array<std::array<std::array<std::pair<double, double>, 2>, 2>, 2> DCAx;
  std::array<std::array<std::array<std::pair<double, double>, 2>, 2>, 2> DCAy;

  std::array<std::array<std::pair<double, double>, 10>, 4> meanDx;
  std::array<std::array<std::pair<double, double>, 10>, 4> meanDy;

  std::array<std::array<TH1*, 4>, 10> dxVsXhistograms;
  std::array<std::array<TH1*, 4>, 10> dxVsYhistograms;
  std::array<std::array<TH1*, 4>, 10> dyVsXhistograms;
  std::array<std::array<TH1*, 4>, 10> dyVsYhistograms;

  std::array<std::array<std::array<std::array<std::pair<double, double>, 10>, 2>, 2>, 2> meanDx_LR_TB_PN;
  std::array<std::array<std::array<std::array<std::pair<double, double>, 10>, 2>, 2>, 2> meanDy_LR_TB_PN;

  std::array<std::array<std::array<std::array<std::pair<double, double>, 10>, 2>, 2>, 2> mchMeanDx_LR_TB_PN;
  std::array<std::array<std::array<std::array<std::pair<double, double>, 10>, 2>, 2>, 2> mchMeanDy_LR_TB_PN;

  std::array<std::array<std::array<TH1*, 2>, 2>, 10> dxVsDEhistograms;
  std::array<std::array<std::array<TH1*, 2>, 2>, 10> dyVsDEhistograms;

  std::array<std::array<std::array<TH1*, 2>, 2>, 10> dxVsPhiHistograms;
  std::array<std::array<std::array<TH1*, 2>, 2>, 10> dyVsPhiHistograms;

  gStyle->SetOptStat(0);
  //gStyle->SetOptStat(1111);
  //gStyle->SetOptFit(1111);
  
  TCanvas c("c", "c", 1200, 800);
  c.SaveAs("residuals_AO2D.pdf(");

  TH1* h1 = GetTH1(fAnalysisResults, "qa-muon/alignment/DCA/vertex_z");
  if (h1) {
    h1->Draw();
    c.SaveAs("residuals_AO2D.pdf");
  }

  TH2* h2 = PlotDCAMFT("DCA/MFT/DCA_y_vs_x");
  if (h2) {
    h2->GetXaxis()->SetRangeUser(-0.1, 0.1);
    h2->GetYaxis()->SetRangeUser(-0.1, 0.1);
    c.SaveAs("residuals_AO2D.pdf");

    h1 = h2->ProjectionX();
    h1->SetTitle("DCA(x)");
    h1->Draw();
    c.SaveAs("residuals_AO2D.pdf");

    h1 = h2->ProjectionY();
    h1->SetTitle("DCA(y)");
    h1->Draw();
    c.SaveAs("residuals_AO2D.pdf");
  }

  c.Clear();
  PlotDCAPhiProjection3D("DCA/MFT/DCA_x_vs_phi_vs_zshift", -0.1, 0.1, 1, c,
      -TMath::Pi()/4.f + TMath::Pi()/2.f, TMath::Pi()/4.f + TMath::Pi()/2.f, false);
  c.SaveAs("residuals_AO2D.pdf");
  PlotDCAPhiProjection3D("DCA/MFT/DCA_y_vs_phi_vs_zshift", -0.1, 0.1, 1, c,
      -TMath::Pi()/4.f, TMath::Pi()/4.f, false);
  c.SaveAs("residuals_AO2D.pdf");

  c.Clear();
  c.SaveAs("residuals_AO2D.pdf)");
}
