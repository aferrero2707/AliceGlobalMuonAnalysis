#include <MFTTracking/Constants.h>

TFile* fAnalysisResults;
std::string pdfFileName;

double epsilon = 1.e-5;
double mchMomMin = 50;

TH1* histDxVsY;

double scaleME = 1;

constexpr int nPoints = 4;

TH1* GetTH1(TFile* f, TString histname)
{
  TH1* result = (TH1*)f->Get(histname);
  std::cout << std::format("{} -> {}", histname.Data(), (void*)result) << std::endl;
  return result;

  //TString histname = TString::Format("ST%d/DE%d/Occupancy_B_XY_%d", station, de, de);
  TKey *key = f->GetKey(histname);
  std::cout << "histname: " << histname << "  key: " <<key << std::endl;
  if (!key) return NULL;
  return (TH1*)key->ReadObjectAny(TH1::Class());
}


TH2* GetTH2(TFile* f, TString histname)
{
  TH2* result = (TH2*)f->Get(histname);
  std::cout << std::format("{} -> {}", histname.Data(), (void*)result) << std::endl;
  return result;

  //TString histname = TString::Format("ST%d/DE%d/Occupancy_B_XY_%d", station, de, de);
  TKey *key = f->GetKey(histname);
  std::cout << "histname: " << histname << "  key: " <<key << std::endl;
  if (!key) return NULL;
  return (TH2*)key->ReadObjectAny(TH2::Class());
}


TH3* GetTH3(TFile* f, TString histname)
{
  TH3* result = (TH3*)f->Get(histname);
  std::cout << std::format("{} -> {}", histname.Data(), (void*)result) << std::endl;
  return result;

  //TString histname = TString::Format("ST%d/DE%d/Occupancy_B_XY_%d", station, de, de);
  TKey *key = f->GetKey(histname);
  std::cout << "histname: " << histname << "  key: " <<key << std::endl;
  if (!key) return NULL;
  return (TH3*)key->ReadObjectAny(TH3::Class());
}

THnSparse* GetTHnSparse(TFile* f, TString histname)
{
  THnSparse* result = (THnSparse*)f->Get(histname);
  std::cout << std::format("{} -> {}", histname.Data(), (void*)result) << std::endl;
  return result;
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


int getChamberIndex(int deId)
{
  return (deId / 100) - 1;
}

int getNumDEinChamber(int chIndex)
{
  int nDE = 0;
  switch (chIndex) {
    case 0:
    case 1:
    case 2:
    case 3:
      nDE = 4;
      break;
    case 4:
    case 5:
      nDE = 18;
      break;
    case 6:
    case 7:
    case 8:
    case 9:
      nDE = 26;
      break;
    default:
      break;
  }
  return nDE;
}

int getDEindexInChamber(int deId)
{
  return (deId - 100) % 100;
}

int getChamberOffset(int chIndex)
{
  int offset = 0;
  for (int c = 0; c < chIndex; c++) {
    offset += getNumDEinChamber(c);
  }
  return offset;
}

int getDEindex(int deId)
{
  auto idx = getDEindexInChamber(deId);
  int offset = getChamberOffset(getChamberIndex(deId));

  return idx + offset;
}

int getDEFromIndex(int index)
{
  int deId = 0;
  for (int chamber = 9; chamber >= 0; chamber--) {
    int offset = getChamberOffset(chamber);
    if (offset > index) {
      continue;
    }

    int indexInChamber = index - offset;
    deId = (chamber + 1) * 100 + indexInChamber;
    break;
  }

  return deId;
}


Double_t VariableWidthGaussian(double* x, double *par)
{
  // distance from mean
  Double_t dx = x[0] - par[1];
  // quadratically variable sigma
  Double_t sigma = par[2] * (par[3] + par[4] * dx + par[5] * dx * dx + par[6] * dx * dx * dx);

  // gaussian + quadratic background
  Double_t result = par[0] * TMath::Exp(-0.5 * (dx * dx) / (sigma * sigma)) + par[7] + par[8] * x[0] + par[9] * x[0] * x[0];
  return result;
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

void PlotDCAProjection(TH2* histogram2, float yMin, float yMax, int projRebin, TCanvas& c, bool printFits = false)
{
  std::string fullHistName = histogram2->GetName();
  TH1* histogramMean = histogram2->ProjectionX();
  histogramMean->SetName(TString::Format("%s-mean", fullHistName.c_str()));
  histogramMean->SetTitle(TString::Format("%s", histogram2->GetTitle()));
  histogramMean->SetMinimum(yMin);
  histogramMean->SetMaximum(yMax);

  int entriesMin = histogram2->GetEntries() / 20;
  int peakMin = 50;
  // skip first bin because it is biased?
  for (int bin = 1; bin <= histogram2->GetXaxis()->GetNbins(); bin++) {
    TH1* proj = histogram2->ProjectionY(TString::Format("%s-%d", fullHistName.c_str(), bin), bin, bin);
    proj->SetTitle(TString::Format("%s - bin %d", histogram2->GetTitle(), bin));
    proj->Rebin(projRebin);
    proj->SetLineColor(kRed);

    /*TH1* projME = histogram2ME->ProjectionY(TString::Format("%s-%d-ME", fullHistNameME.c_str(), bin), bin, bin);
    projME->Rebin(projRebin);

    double integral1 = proj->Integral(1, proj->GetXaxis()->FindBin(-0.2));
    double integral2 = proj->Integral(proj->GetXaxis()->FindBin(0.2), proj->GetXaxis()->GetNbins());
    double integralME1 = projME->Integral(1, projME->GetXaxis()->FindBin(-0.2));
    double integralME2 = projME->Integral(projME->GetXaxis()->FindBin(0.2), projME->GetXaxis()->GetNbins());
    double scaleME = (integral1 + integral2) / (integralME1 + integralME2);
    projME->Scale(scaleME);*/

    double mean = -10000;
    double meanErr = 0;
    double sigma = -10000;
    double sigmaErr = 0;

    int valuePeak = proj->GetMaximum();
    int binPeak = proj->GetMaximumBin();
    double xPeak = proj->GetXaxis()->GetBinCenter(binPeak);

    if (valuePeak > peakMin) {

      if (printFits) {
        std::cout << std::format("[\"{}\"] bin {} max={}  binPeak={}  xPeak={}\n", fullHistName, bin, valuePeak, binPeak, xPeak);
        //std::cout << std::format("[\"{}\"] bin {}  max={:0.2f}  binPeak={}  xPeak={:0.2f}",
        //    fullHistName, bin, valuePeak, binPeak, xPeak) << std::endl;
      }
      //std::cout << std::format("Bin {}  max={}  binPeak={}  xPeak={}\n", bin, valuePeak, binPeak, xPeak);

      TF1 fcb("fcb", DoubleSidedCB, -0.1, 0.1, 7);
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
        c.SaveAs(pdfFileName.c_str());
      }

      mean = fcb.GetParameter(1);
      meanErr = fcb.GetParError(1);
      sigma = fcb.GetParameter(2);
      sigmaErr = fcb.GetParError(2);

      //mean = fgaus.GetParameter(1);
      //meanErr = fgaus.GetParError(1);
      //sigma = fgaus.GetParameter(2);
      //sigmaErr = fgaus.GetParError(2);
    }

    histogramMean->SetBinContent(bin, mean);
    histogramMean->SetBinError(bin, meanErr);
  }

  TF1 linFit("linFit", "pol1");
  linFit.SetLineColor(kRed);
  linFit.SetLineStyle(kDotted);
  histogramMean->Fit("linFit", "Q");
  histogramMean->Draw("E");

  TLine* line1 = new TLine(histogramMean->GetXaxis()->GetXmin(), 0, histogramMean->GetXaxis()->GetXmax(), 0);
  line1->SetLineStyle(kDashed);
  line1->Draw();

  TLine* line2 = new TLine(0, yMin, 0, yMax);
  line2->SetLineStyle(kDashed);
  line2->Draw();

  std::cout << std::format("Slope: {:0.4f} mm / 10 m", linFit.GetParameter(1) * 1000 * 10) << std::endl;
  //c.SaveAs(pdfFileName.c_str());
  //histogramSigma->Draw("E");
  //c.SaveAs(pdfFileName.c_str());

  histDxVsY = histogramMean;
}

void PlotDCAProjection(std::string histName, float yMin, float yMax, int projRebin, TCanvas& c, bool printFits = false)
{
  std::string fullHistName = std::string("qa-muon/alignment/") + histName;
  TH2* histogram2 = GetTH2(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogram2 << std::endl;
  if (!histogram2)
    return;

  return PlotDCAProjection(histogram2, yMin, yMax, projRebin, c, printFits);
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
  int peakMin = 50;
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
        c.SaveAs(pdfFileName.c_str());
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
  std::string fullHistName = std::string("qa-muon/alignment/") + histName;
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
    c.SaveAs(pdfFileName.c_str());

    c.Clear();
    auto amplitude = PlotDCAPhiProjection(proj, yMin, yMax, projRebin, c, phaseMin, phaseMax, printFits);
    c.SaveAs(pdfFileName.c_str());

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

std::pair<double, double> PlotDCAMCH(TH1* histogram, int rebin = 1)
{
  //std::string fullHistName = std::string("qa-muon/alignment/") + histName;
  //TH1* histogram = GetTH1(fAnalysisResults, fullHistName);
  //std::cout << fullHistName << " -> " << histogram << std::endl;
  if (!histogram)
    return {nan(""), nan("")};
  histogram->Rebin(rebin);

  histogram->Draw("E");

  int valuePeak = histogram->GetMaximum();
  int binPeak = histogram->GetMaximumBin();
  double xPeak = histogram->GetXaxis()->GetBinCenter(binPeak);

  if (valuePeak < 10)
    return {nan(""), nan("")};

  bool cbFit = false;

  if (cbFit) {
    TF1 fcb("fcb", DoubleSidedCB, -6, 6, 7);
    fcb.SetNpx(1000);
    fcb.SetLineColor(kBlack);
    //fcb.SetParameter(0, valuePeak);
    double par[7];
    par[0]=35000;
    par[1]=xPeak;
    par[2]=0.5;
    par[3]=1;
    par[4]=1;
    par[5]=1;
    par[6]=1;
    fcb.SetParameters(&par[0]);
    fcb.FixParameter(1, xPeak);
    fcb.SetParLimits(2, 0.5, 10);
    histogram->Fit("fcb", "BRN");
    fcb.ReleaseParameter(1);
    histogram->Fit("fcb", "BRN");
    histogram->Fit("fcb", "BR");

    return {fcb.GetParameter(1), fcb.GetParError(1)};
  } else {
    TF1 fgaus2("fgaus2", "gaus(0)+pol2(3)");
    fgaus2.SetNpx(1000);
    fgaus2.SetLineColor(kRed);
    fgaus2.SetParameter(0, valuePeak / 10);
    //fgaus2.SetParLimits(0, 0, valuePeak * 10);
    fgaus2.SetParameter(1, xPeak);
    fgaus2.FixParameter(1, xPeak);
    fgaus2.SetParameter(2, 1);
    fgaus2.SetParLimits(2, 0, 10);
    histogram->Fit("fgaus2", "BQN");
    fgaus2.ReleaseParameter(1);
    fgaus2.SetParLimits(1, -10, 10);
    histogram->Fit("fgaus2", "BQN");

    TF1 fgaus("fgaus", "gaus(0)", fgaus2.GetParameter(1) - 1.5 * fgaus2.GetParameter(2), fgaus2.GetParameter(1) + 1.5 * fgaus2.GetParameter(2));
    fgaus.SetNpx(1000);
    fgaus.SetLineColor(kBlack);
    //fgaus.SetParLimits(0, 0, 100000000000.0);
    fgaus.SetParameter(0, fgaus2.GetParameter(0));
    fgaus.SetParameter(1, fgaus2.GetParameter(1));
    fgaus.SetParameter(2, fgaus2.GetParameter(2));
    fgaus.SetParLimits(2, 0, 100);
    //fgaus.SetParameter(3, fgaus2.GetParameter(0)/10);
    //fgaus.FixParameter(3, 0);
    //fgaus.SetParameter(4, fgaus2.GetParameter(1));
    //fgaus.SetParameter(5, fgaus2.GetParameter(2)*5);
    //fgaus.SetParLimits(5, fgaus2.GetParameter(2)*2, 100);
    histogram->Fit("fgaus", "RQB+");

    return {fgaus.GetParameter(1), fgaus.GetParError(1)};
  }
}

void PlotDCAMCHvsMomentum(THnSparse* hnx, THnSparse* hny, TCanvas& c)
{
  // axes assignments:
  // 0: momentum
  // 1: quadrant
  // 2: sign
  // 3: DCA
  int pAxisId = 0;
  int quadrantAxisId = 1;
  int signAxisId = 2;
  int dcaAxisId = 3;

  int momRebin = 1;

  if (!hnx || !hny) return;

  c.Clear();
  c.cd();

  std::array<THnSparse*, 2> hn{ hnx, hny };

  std::array<std::string, 2> coordinates{ "x", "y" };

  std::array<std::string, 4> quadrants = {"Q0", "Q1", "Q2", "Q3"};

  int colors[4] = {kBlue, kRed, kOrange, kCyan};
  //int markers[4] = {kStar, kCircle, kMultiply, kFullDotLarge};
  int markers[2] = {kStar, kCircle};
  int lineStyles[2] = {kSolid, kDashed};


  std::vector<double> momenta;
  // dca[coordinate][sign][quadrant][momentum]
  std::array<std::array<std::array<std::vector<double>, 4>, 2>, 2> momentum;
  std::array<std::array<std::array<std::vector<double>, 4>, 2>, 2> momentumErr;
  std::array<std::array<std::array<std::vector<double>, 4>, 2>, 2> dca;
  std::array<std::array<std::array<std::vector<double>, 4>, 2>, 2> dcaErr;

  // loop on coordinates
  for (int ci = 0; ci < 2; ci++) {
    auto coordinate = coordinates[ci];

    // (sign == 0 ? "positive" : "negative")
    TMultiGraph* mg = new TMultiGraph();
    mg->SetTitle(Form("DCA(%s) vs. momentum;p (GeV/c); DCA (cm)", coordinate.c_str()));
    TGraphErrors* gr = new TGraphErrors();
    gr->AddPoint(0, 1000);
    gr->AddPoint(hn[ci]->GetAxis(pAxisId)->GetXmax(), 1000);
    mg->Add(gr,"lp");

    TLegend* legend = new TLegend(0.6, 0.8, 0.9, 0.9);
    legend->SetNColumns(4);

    // loop on charge sign
    for (int sign = 0; sign < 2; sign++) {
      hn[ci]->GetAxis(signAxisId)->SetRange(sign + 1, sign + 1);
      // loop on quadrant
      for (int quadrant = 0; quadrant < 4; quadrant++) {
        hn[ci]->GetAxis(quadrantAxisId)->SetRange(quadrant + 1, quadrant + 1);
        // loop on momentum
        for (int mom = 0; mom < hn[ci]->GetAxis(pAxisId)->GetNbins(); mom += momRebin) {
          hn[ci]->GetAxis(pAxisId)->SetRange(mom + 1, mom + momRebin);
          //hn[ci]->GetAxis(pAxisId)->SetRange(5, -1);

          auto* proj = hn[ci]->Projection(dcaAxisId);
          proj->SetTitle(std::format("MCH DCA({}), Q{} {}", coordinate, quadrant, (sign == 0 ? "+" : "-")).c_str());
          auto result = PlotDCAMCH(proj);
          delete proj;
          if (std::isnan(result.first) || std::isnan(result.second)) {
            continue;
          }

          double avgMom = (hn[ci]->GetAxis(pAxisId)->GetBinLowEdge(mom + 1) +
                           hn[ci]->GetAxis(pAxisId)->GetBinUpEdge(mom + momRebin)) / 2.f;
          momentum[ci][sign][quadrant].push_back(avgMom);
          momentumErr[ci][sign][quadrant].push_back(0);
          dca[ci][sign][quadrant].push_back(result.first);
          dcaErr[ci][sign][quadrant].push_back(result.second);
        }
      }

      for (int quadrant = 0; quadrant < 4; quadrant++) {
        if (dca[ci][sign][quadrant].empty()) {
          continue;
        }
        std::cout << std::format("Q{}  dca[0]={}", quadrant, dca[ci][sign][quadrant][0]) << std::endl;
        TGraphErrors* gr = new TGraphErrors(momentum[ci][sign][quadrant].size(),
                                            momentum[ci][sign][quadrant].data(),
                                            dca[ci][sign][quadrant].data(),
                                            momentumErr[ci][sign][quadrant].data(),
                                            dcaErr[ci][sign][quadrant].data());
        gr->SetLineColor(colors[quadrant]);
        gr->SetMarkerColor(colors[quadrant]);
        gr->SetMarkerStyle(markers[sign]);
        gr->SetMarkerSize(2);
        gr->SetLineStyle(lineStyles[sign]);
        mg->Add(gr,"pl");
        auto* entry = legend->AddEntry(gr, (quadrants[quadrant] + (sign == 0 ? " (+)" : " (-)")).c_str(), "P");
        entry->SetTextColor(colors[quadrant]);
      }
    }
    mg->Draw("a");
    mg->SetMinimum(-1);
    mg->SetMaximum(1);
    legend->Draw();

    c.SaveAs(pdfFileName.c_str());
  }
}

#define __SIMPLE_FIT

std::pair<double, double> PlotDXY(TH1* proj, TCanvas& c, std::string pdfName = "", bool printFits = false)
{
  c.Clear();
  proj->SetLineColor(kRed);
  proj->Draw("E");

  TF1 fgaus("fgaus", "gausn(0)+pol2(3)");
  int valuePeak = proj->GetMaximum();
  if (valuePeak < 10)
    return {nan(""), nan("")};

  int binPeak = proj->GetMaximumBin();
  double xPeak = proj->GetXaxis()->GetBinCenter(binPeak);
  fgaus.SetNpx(1000);
  fgaus.SetLineColor(kBlack);
  fgaus.SetParameter(0, valuePeak / 10);
  fgaus.SetParLimits(0, 0, valuePeak * 10);
  fgaus.SetParameter(1, xPeak);
  fgaus.SetParameter(2, 1);
  fgaus.SetParLimits(2, 0, 10);
  fgaus.SetParameter(3, 0);
  fgaus.SetParameter(4, 0);
  fgaus.SetParameter(5, 0);
  proj->Fit("fgaus", "BQN");

#ifdef SIMPLE_FIT
  xPeak = fgaus.GetParameter(1);
  TF1 fgaus2("fgaus2", "gausn(0)", xPeak - fgaus.GetParameter(2), xPeak + fgaus.GetParameter(2));
  fgaus2.SetNpx(1000);
  fgaus2.SetLineColor(kBlue);
  fgaus2.SetParameter(0, fgaus.GetParameter(0));
  fgaus2.SetParameter(1, xPeak);
  fgaus2.SetParameter(2, fgaus.GetParameter(2));
  proj->Fit("fgaus2", "BRQ");
#else
  TF1 fgaus2("fgaus2", VariableWidthGaussian, proj->GetXaxis()->GetXmin(), proj->GetXaxis()->GetXmax(), 10);
  fgaus2.SetNpx(1000);
  fgaus2.SetLineColor(kBlue);
  fgaus2.SetParameter(0, fgaus.GetParameter(0));
  fgaus2.SetParameter(1, xPeak);
  fgaus2.SetParameter(2, fgaus.GetParameter(2));
  fgaus2.SetParameter(3, 1);
  fgaus2.FixParameter(4, 0);
  fgaus2.FixParameter(5, 0);
  fgaus2.FixParameter(6, 0);
  fgaus2.SetParameter(7, 0);
  fgaus2.SetParameter(8, 0);
  fgaus2.SetParameter(9, 0);
  proj->Fit("fgaus2", "BRQN");
  //fgaus2.ReleaseParameter(4);
  //fgaus2.ReleaseParameter(5);
  //fgaus2.ReleaseParameter(6);
  proj->Fit("fgaus2", "BRQ");
#endif

/*
  TF1 fcb("fcb","[0]*ROOT::Math::crystalball_function(x, [1], [2], [3], [4]) + [5]");
  fcb.SetParameters(100, 0.6, -2.13903e+06, 1, xPeak);
  fcb.SetParLimits(3, 0, 100);
  fcb.SetNpx(1000);
  fcb.SetLineColor(kBlack);
  //proj2->Fit("fcb", "BNQ");
  //proj2->Fit("fcb", "BQ");

  proj2->Draw("E");
 */
  if (!pdfName.empty()) {
    proj->GetXaxis()->SetRangeUser(-10.0, 10.0);
    c.SaveAs(pdfName.c_str());
  }

  //PlotDXYProjection(fullHistName.c_str(), fullHistNameME.c_str(), histogram2, histogram2ME, scaleME, -5.0, 5.0, 1, c, true, printFits);

  return {fgaus2.GetParameter(1), fgaus2.GetParError(1)};
  //return {fcb.GetParameter(4), fgaus.GetParError(4)};
}

void PlotChamberCoverage(THnSparse* hnx, int chamber, TCanvas& c)
{
  if (hnx->GetNdimensions() < 5) {
    return;
  }

  int chIndex = chamber - 1;

  int residualAxisId = 0;
  int deAxisId = 1;
  int quadrantAxisId = 2;
  int signAxisId = 3;
  int pAxisId = 4;

  int momRebin = 2;

  if (!hnx) return;

  c.Clear();
  c.cd();

  int colors[4] = {kBlue, kRed, kOrange, kCyan};
  //int markers[4] = {kStar, kCircle, kMultiply, kFullDotLarge};
  int markers[2] = {kStar, kCircle};
  int lineStyles[2] = {kSolid, kDashed};


  hnx->GetAxis(residualAxisId)->SetRange(0, -1);
  hnx->GetAxis(signAxisId)->SetRange(0, -1);

  // select DEs in chamber
  int deIndexMin = getChamberOffset(chIndex);
  int deIndexMax = deIndexMin + getNumDEinChamber(chIndex) - 1;
  hnx->GetAxis(deAxisId)->SetRange(deIndexMin + 1, deIndexMax + 1);

  // select momentum
  if (pAxisId >= 0) {
    // p >= xx GeV/c
    int momBin = hnx->GetAxis(pAxisId)->FindBin(mchMomMin + epsilon);
    int momBinMax = hnx->GetAxis(pAxisId)->GetNbins();
    hnx->GetAxis(pAxisId)->SetRange(momBin, momBinMax);
  }

  // loop on quadrants
  for (int quadrant = 0; quadrant < 4; quadrant++) {

    // select quadrant
    hnx->GetAxis(quadrantAxisId)->SetRange(quadrant + 1, quadrant + 1);

    TH1D* proj[2];
    // loop on charge sign
    for (int sign = 0; sign < 2; sign++) {
      hnx->GetAxis(signAxisId)->SetRange(sign + 1, sign + 1);

      proj[sign] = hnx->Projection(deAxisId);

      // set DE names in X axis
      for (int debin = 1; debin <= proj[sign]->GetXaxis()->GetNbins(); debin++) {
        int deId = getDEFromIndex(deIndexMin + debin - 1);
        proj[sign]->GetXaxis()->SetBinLabel(debin, TString::Format("DE%d", deId));
      }
      //proj->SetName(std::format("d{}_vs_de_{}_{}", (ci == 0 ? "x" : "y"), deIndex+1, (sign == 0 ? "positive" : "negative")).c_str());
      //proj->SetTitle(std::format("#Delta{}, DE{} {}", (ci == 0 ? "x" : "y"), deId, (sign == 0 ? "positive" : "negative")).c_str());
    }

    proj[0]->SetLineColor(colors[quadrant]);
    proj[0]->SetLineStyle(lineStyles[0]);
    proj[1]->SetLineColor(colors[quadrant]);
    proj[1]->SetLineStyle(lineStyles[1]);

    proj[0]->Draw("hist");
    proj[1]->Draw("hist same");

    c.SaveAs(pdfFileName.c_str());
  }
}


void muonGlobalAlignmentChamberCoverage(int chamber = 7, const char* _rootFileName = "AnalysisResults.root")
{
  fAnalysisResults = new TFile(_rootFileName);
  pdfFileName = std::format("CH{}-coverage.pdf", chamber);

  std::string taskPath = "muon-global-alignment"; // = (wagonId > 0) ? std::format("muon-global-alignment_id{}", wagonId) : "muon-global-alignment";
  //taskPath = "muon-global-alignment_id47994";
  //taskPath = "muon-global-alignment_p-15-pt-2_id47994";
  //taskPath = "muon-global-alignment_p-20-pt-2_id47994";
  //taskPath = "muon-global-alignment_p-50-pt-4_id47994";

  //taskPath = "muon-global-alignment_id48150";
  //taskPath = "muon-global-alignment_id48242";
  //taskPath = "muon-global-alignment_id48243";

  std::array<std::string, 4> quadrants = {"Q0", "Q1", "Q2", "Q3"};

  gStyle->SetOptStat(0);
  //gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  
  TCanvas c("c", "c", 1200, 800);
  c.SaveAs((pdfFileName + "(").c_str());

  PlotChamberCoverage(
      GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/dx_vs_de_corr").c_str()),
      chamber, c
  );

  c.Clear();
  c.SaveAs((pdfFileName + ")").c_str());
}
