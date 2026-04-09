#include <MFTTracking/Constants.h>

TFile* fAnalysisResults;
std::string pdfFileName;

double epsilon = 1.e-5;
double mchMomMin = 50;

double chamberZ[10] = {526.16, 545.24, 676.4, 695.4, 967.5,
                              998.5, 1276.5, 1307.5, 1406.6, 1437.6};

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
  //if (valuePeak < 10) {
  //  return {nan(""), nan("")};
  //}
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

void PlotChamberResidualvsMomentum(THnSparse* hnx, THnSparse* hny, TCanvas& c)
{
  if (hnx->GetNdimensions() < 5) {
    return;
  }

  int residualAxisId = 0;
  int deAxisId = 1;
  int quadrantAxisId = 2;
  int signAxisId = 3;
  int pAxisId = 4;

  int momRebin = 2;

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


  // loop on coordinates
  for (int ci = 0; ci < 2; ci++) {
    auto coordinate = coordinates[ci];
    // loop on chamber ID
    for (int chamber = 0; chamber < 10; chamber++) {
      if (chamber > 0 && chamber < 9) continue;
      int deIdMin = getChamberOffset(chamber);
      int deIdMax = deIdMin + getNumDEinChamber(chamber) - 1;
      hn[ci]->GetAxis(deAxisId)->SetRange(deIdMin + 1, deIdMax + 1);

      // dca[coordinate][sign][quadrant][momentum]
      std::array<std::array<std::vector<double>, 4>, 2> momentum;
      std::array<std::array<std::vector<double>, 4>, 2> momentumErr;
      std::array<std::array<std::vector<double>, 4>, 2> dca;
      std::array<std::array<std::vector<double>, 4>, 2> dcaErr;

      TMultiGraph* mg = new TMultiGraph();
      mg->SetTitle(Form("CH%d #Delta(%s) vs. momentum;p (GeV/c); #Delta(%s) (cm)", chamber+1, coordinate.c_str(), coordinate.c_str()));
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

            auto* proj = hn[ci]->Projection(residualAxisId);
            proj->SetTitle(std::format("MCH DCA({}), Q{} {}", coordinate, quadrant, (sign == 0 ? "+" : "-")).c_str());
            auto result = PlotDXY(proj, c);
            if (std::isnan(result.first) || std::isnan(result.second)) {
              continue;
            }
            if (false && ci == 0 && chamber == 9 && sign == 0 && quadrant == 0) {
              proj->Draw();
              c.SaveAs(pdfFileName.c_str());
            }

            double avgMom = (hn[ci]->GetAxis(pAxisId)->GetBinLowEdge(mom + 1) +
                hn[ci]->GetAxis(pAxisId)->GetBinUpEdge(mom + momRebin)) / 2.f;
            momentum[sign][quadrant].push_back(avgMom);
            momentumErr[sign][quadrant].push_back(0);
            dca[sign][quadrant].push_back(result.first);
            dcaErr[sign][quadrant].push_back(result.second);

            delete proj;
          }
        }

        for (int quadrant = 0; quadrant < 4; quadrant++) {
          if (dca[sign][quadrant].empty()) {
            std::cout << std::format("CH{} Q{} ({}) d{}[0]=NA",
                chamber+1, quadrant, (sign == 0 ? "+" : "-"), coordinate) << std::endl;
            continue;
          }
          std::cout << std::format("CH{} Q{} ({}) d{}[0]={}",
              chamber+1, quadrant, (sign == 0 ? "+" : "-"), coordinate, dca[sign][quadrant][0]) << std::endl;
          TGraphErrors* gr = new TGraphErrors(momentum[sign][quadrant].size(),
              momentum[sign][quadrant].data(),
              dca[sign][quadrant].data(),
              momentumErr[sign][quadrant].data(),
              dcaErr[sign][quadrant].data());
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
}

void PlotDXYProjection(const char* fullHistName, const char* fullHistNameME, TH2* histogram2, TH2* histogram2ME, double scaleME, float yMin, float yMax, int projRebin, TCanvas& c, bool subtractBackground, bool printFits = false)
{
  c.Clear();

  TH1* histogramMean = histogram2->ProjectionX();
  histogramMean->SetName(TString::Format("%s-mean", fullHistName));
  histogramMean->SetTitle(TString::Format("%s", histogram2->GetTitle()));
  histogramMean->SetMinimum(yMin);
  histogramMean->SetMaximum(yMax);
  TH1* histogramSigma = histogram2->ProjectionX();
  histogramSigma->SetName(TString::Format("%s-sigma", fullHistName));
  histogramSigma->SetTitle(TString::Format("%s (sigma)", histogram2->GetTitle()));
  histogramSigma->SetMinimum(0);
  histogramSigma->SetMaximum(yMax);

  int entriesMin = histogram2->GetEntries() / 20;
  // skip first bin because it is biased?
  for (int bin = 1; bin <= histogram2->GetXaxis()->GetNbins(); bin++) {
    TH1* proj = histogram2->ProjectionY(TString::Format("%s-%d", fullHistName, bin), bin, bin);
    TH1* projME = histogram2ME->ProjectionY(TString::Format("%s-%d", fullHistNameME, bin), bin, bin);

    proj->Rebin(projRebin);
    projME->Rebin(projRebin);

    projME->Scale(scaleME);
    proj->SetTitle(TString::Format("%s - bin %d", histogram2->GetTitle(), bin));
    double mean = -10000;
    double meanErr = 0;
    double sigma = -10000;
    double sigmaErr = 0;
    if (proj->GetEntries() > entriesMin) {

      if (printFits && subtractBackground) {
        proj->SetLineColor(kRed);
        proj->Draw("E");
        projME->Draw("L same");
        c.SaveAs(pdfFileName.c_str());
      }

      if (subtractBackground)
        proj->Add(projME, -1);

      //proj->Rebin(projRebin);
      int valuePeak = proj->GetMaximum();
      int binPeak = proj->GetMaximumBin();
      double xPeak = proj->GetXaxis()->GetBinCenter(binPeak);
      if (printFits) {
        std::cout << std::format("[\"{}\"] bin {} max={}  binPeak={}  xPeak={}\n", fullHistName, bin, valuePeak, binPeak, xPeak);
        //std::cout << std::format("[\"{}\"] bin {}  max={:0.2f}  binPeak={}  xPeak={:0.2f}",
        //    fullHistName, bin, valuePeak, binPeak, xPeak) << std::endl;
      }
      //std::cout << std::format("Bin {}  max={}  binPeak={}  xPeak={}\n", bin, valuePeak, binPeak, xPeak);

      TF1 fgaus2("fgaus2", "gaus(0)+pol0(3)");
      fgaus2.SetNpx(1000);
      fgaus2.SetLineColor(kBlack);
      fgaus2.SetParameter(0, valuePeak / 10);
      //fgaus2.SetParLimits(0, 0, valuePeak * 10);
      fgaus2.SetParameter(1, xPeak);
      fgaus2.FixParameter(1, xPeak);
      fgaus2.SetParameter(2, 1);
      fgaus2.SetParLimits(2, 0, 10);
      fgaus2.SetParameter(3, 0);
      fgaus2.SetParameter(4, 0);
      proj->Fit("fgaus2", "BQN");
      fgaus2.ReleaseParameter(1);
      fgaus2.SetParLimits(1, -10, 10);
      proj->Fit("fgaus2", "BQN");

      TF1 fgaus("fgaus", "gaus", fgaus2.GetParameter(1) - 2.0 * fgaus2.GetParameter(2), fgaus2.GetParameter(1) + 2.0 * fgaus2.GetParameter(2));
      fgaus.SetNpx(1000);
      fgaus.SetLineColor(kBlack);
      //fgaus.SetParLimits(0, 0, 100000000000.0);
      fgaus.SetParameter(0, fgaus2.GetParameter(0));
      fgaus.SetParameter(1, fgaus2.GetParameter(1));
      fgaus.SetParameter(2, fgaus2.GetParameter(2));
      fgaus.SetParLimits(2, 0, 100);
      proj->Fit("fgaus", "RQBN");

      TF1 fcb("fcb", DoubleSidedCBwithLinBgd, -30, 30, 9);
      fcb.SetNpx(1000);
      fcb.SetLineColor(kBlack);
      //fcb.SetParameter(0, valuePeak);
      double par[9];
      par[0]=fgaus2.GetParameter(0);
      par[1]=fgaus2.GetParameter(1);
      par[2]=fgaus2.GetParameter(2);
      par[3]=1;
      par[4]=1;
      par[5]=1;
      par[6]=1;
      par[7]=0;
      par[8]=0;
      fcb.SetParameters(&par[0]);

      fcb.SetParLimits(2, 0.5, 10);
      fcb.FixParameter(3, 1);
      fcb.FixParameter(4, 1);
      fcb.FixParameter(5, 1);
      fcb.FixParameter(6, 1);
      //histogram->Fit("fcb", "BRN");
      //fcb.ReleaseParameter(1);
      proj->Fit("fcb", "BRQN");
      fcb.ReleaseParameter(3);
      fcb.ReleaseParameter(4);
      fcb.ReleaseParameter(5);
      fcb.ReleaseParameter(6);
      if (printFits)
        proj->Fit("fcb", "BR");
      else
        proj->Fit("fcb", "BRQ");

      if (printFits) {
        proj->Draw("E");
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

    histogramSigma->SetBinContent(bin, sigma);
    histogramSigma->SetBinError(bin, sigmaErr);
  }
  //histogramMean->Draw("E");
  //c.SaveAs(pdfFileName.c_str());
  //histogramSigma->Draw("E");
  //c.SaveAs(pdfFileName.c_str());

  histDxVsY = histogramMean;
}

std::array<std::pair<double, double>, 2> PlotDXYvsDE(std::string histName, int chamber, TCanvas& c, bool printFits = false)
{
  std::array<std::pair<double, double>, 2> result;

  c.Clear();
  c.Divide(2, 2);

  std::string fullHistName = std::string("qa-muon/alignmentsame-event/Residuals/") + histName;
  TH2* histogram2 = GetTH2(fAnalysisResults, fullHistName);
  //std::cout << fullHistName << " -> " << histogram2 << std::endl;
  if (!histogram2) return {};

  c.cd(1);
  histogram2->Draw("col");

  std::string fullHistNameME = std::string("qa-muon/alignmentmixed-event/Residuals/") + histName;
  TH2* histogram2ME = GetTH2(fAnalysisResults, fullHistNameME);

  c.cd(2);
  histogram2ME->Draw("col");

  TH1* proj[2];
  TH1* projME[2];

  proj[0] = histogram2->ProjectionY((fullHistName + "_pyL").c_str());
  projME[0] = histogram2ME->ProjectionY((fullHistNameME + "_pyL").c_str());
  proj[1] = histogram2->ProjectionY((fullHistName + "_pyR").c_str());
  projME[1] = histogram2ME->ProjectionY((fullHistNameME + "_pyR").c_str());

  for (int lr = 0; lr < 2; lr++) {
    proj[lr]->Reset();
    proj[lr]->SetLineColor(kRed);
    projME[lr]->Reset();
  }

  for (int bin = 1; bin <= histogram2->GetXaxis()->GetNbins(); bin++) {
    TH1* p = histogram2->ProjectionY(TString::Format("%s-%d", fullHistName.c_str(), bin), bin, bin);
    TH1* pME = histogram2ME->ProjectionY(TString::Format("%s-%d", fullHistNameME.c_str(), bin), bin, bin);

    int lr = GetLR(chamber, bin - 1);
    //std::cout << std::format("CH{} DE{} LR={}", chamber, bin - 1, lr) << std::endl;
    if (lr < 0)
      continue;

    proj[lr]->Add(p);
    projME[lr]->Add(pME);
  }

  for (int lr = 0; lr < 2; lr++) {
    proj[lr]->Rebin(4);
    projME[lr]->Rebin(4);

    double integral1 = proj[lr]->Integral(1, proj[lr]->GetXaxis()->FindBin(-15));
    double integral2 = proj[lr]->Integral(proj[lr]->GetXaxis()->FindBin(15), proj[lr]->GetXaxis()->GetNbins());
    double integralME1 = projME[lr]->Integral(1, projME[lr]->GetXaxis()->FindBin(-15));
    double integralME2 = projME[lr]->Integral(projME[lr]->GetXaxis()->FindBin(15), projME[lr]->GetXaxis()->GetNbins());
    double scaleME = (integral1 + integral2) / (integralME1 + integralME2);

    projME[lr]->Scale(scaleME);

    c.cd(3 + lr);
    TH1* proj2 = (TH1*)proj[lr]->Clone((fullHistName + "_py_" + std::to_string(lr)).c_str());
    proj2->Add(projME[lr], -1);
    proj2->SetLineColor(kRed);
    int valuePeak = proj2->GetMaximum();
    int binPeak = proj2->GetMaximumBin();
    double xPeak = proj2->GetXaxis()->GetBinCenter(binPeak);
    TF1 fgaus("fgaus", "gausn(0)+pol0(3)");
    fgaus.SetNpx(1000);
    fgaus.SetLineColor(kBlack);
    fgaus.SetParameter(0, valuePeak / 10);
    fgaus.SetParLimits(0, 0, valuePeak * 10);
    fgaus.SetParameter(1, xPeak);
    fgaus.SetParameter(2, 1);
    fgaus.SetParLimits(2, 0, 10);
    fgaus.SetParameter(3, 0);
    fgaus.SetParameter(4, 0);
    proj2->Fit("fgaus", "BQN");

    TF1 fcb("fcb", DoubleSidedCBwithLinBgd, proj2->GetXaxis()->GetXmin(), proj2->GetXaxis()->GetXmax(), 9);
    fcb.SetNpx(1000);
    fcb.SetLineColor(kBlack);
    //fcb.SetParameter(0, valuePeak);
    double par[9];
    par[0]=fgaus.GetParameter(0);
    par[1]=fgaus.GetParameter(1);
    par[2]=fgaus.GetParameter(2);
    par[3]=1;
    par[4]=1;
    par[5]=1;
    par[6]=1;
    par[7]=0;
    par[8]=0;
    fcb.SetParameters(&par[0]);

    fcb.FixParameter(1, fgaus.GetParameter(1));
    fcb.SetParLimits(2, fgaus.GetParameter(2) * 0.5, fgaus.GetParameter(2) * 2);
    proj2->Fit("fcb", "BRNQ");
    fcb.ReleaseParameter(1);
    proj2->Fit("fcb", "BRQ");

    proj2->Draw("E");
    //proj[lr]->Draw("E");

    //proj[lr]->Draw("E");
    //projME[lr]->Draw("E same");

    //result[lr].first = fgaus.GetParameter(1);
    //result[lr].second = fgaus.GetParError(1);
    result[lr].first = fcb.GetParameter(1);
    result[lr].second = fcb.GetParError(1);
  }

  c.SaveAs(pdfFileName.c_str());

  PlotDXYProjection(fullHistName.c_str(), fullHistNameME.c_str(), histogram2, histogram2ME, scaleME, -5.0, 5.0, 8, c, false, printFits);

  return result;
  //return {fcb.GetParameter(4), fgaus.GetParError(4)};
}

void PlotZTrend(int n, double* xv, std::array<std::array<std::pair<double, double>, 10>, 4>& values, std::array<std::pair<double, double>, 4>& dca, const char* title, double ymin, double ymax, TCanvas& c)
{
  double exv[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  double yv[10];
  double eyv[10];

  std::array<std::string, 4> quadrants = {"Q0", "Q1", "Q2", "Q3"};

  int colors[4] = {kBlue, kRed, kOrange, kCyan};
  int markers[4] = {kStar, kCircle, kMultiply, kFullDotLarge};

  c.Clear();
  c.cd();

  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(title);
  TGraphErrors* gr = new TGraphErrors();
  gr->AddPoint(0, 1000);
  gr->AddPoint(1500, 1000);
  mg->Add(gr,"lp");

  for (int j = 0; j < quadrants.size(); j++) {
    TGraphErrors* gr = new TGraphErrors();
    gr->SetLineColor(colors[j]);
    gr->SetMarkerColor(colors[j]);
    gr->SetMarkerStyle(kDot);
    gr->SetMarkerSize(0);
    gr->SetLineStyle(kDashed);
    gr->AddPoint(0, dca[j].first);
    gr->AddPoint(xv[0], values[j][0].first);
    mg->Add(gr,"l");
  }

  for (int j = 0; j < quadrants.size(); j++) {
    TGraphErrors* gr = new TGraphErrors();
    gr->SetLineColor(colors[j]);
    gr->SetMarkerColor(colors[j]);
    gr->SetMarkerStyle(markers[2]);
    gr->SetMarkerSize(2);
    gr->AddPoint(0, dca[j].first);
    mg->Add(gr,"p");
  }

  TLegend* legend = new TLegend(0.6, 0.8, 0.9, 0.9);
  legend->SetNColumns(4);
  for (int j = 0; j < quadrants.size(); j++) {
    for (int i = 0; i < 10; i++) {
      yv[i] = values[j][i].first;
      eyv[i] = values[j][i].second;
    }
    TGraphErrors* gr = new TGraphErrors(n, xv, yv, exv, eyv);
    gr->SetLineColor(colors[j]);
    gr->SetMarkerColor(colors[j]);
    gr->SetMarkerStyle(markers[1]);
    gr->SetMarkerSize(2);
    mg->Add(gr,"lp");
    auto* entry = legend->AddEntry(gr, (quadrants[j] + " ").c_str(), "P");
    entry->SetTextColor(colors[j]);
  }
  mg->Draw("a");
  mg->SetMinimum(ymin);
  mg->SetMaximum(ymax);
  legend->Draw();



  c.SaveAs(pdfFileName.c_str());
}

void PlotDECorrelations(TH2* hpos, TH2* hneg, TCanvas& c)
{
  if (!hpos || !hneg) return;

  c.cd();

  // set DE names in axis
  for (int xbin = 1; xbin <= hpos->GetXaxis()->GetNbins(); xbin++) {
    int deId = getDEFromIndex(xbin - 1);
    hpos->GetXaxis()->SetBinLabel(xbin, TString::Format("DE%d", deId));
    hneg->GetXaxis()->SetBinLabel(xbin, TString::Format("DE%d", deId));
  }
  for (int ybin = 1; ybin <= hpos->GetYaxis()->GetNbins(); ybin++) {
    int deId = getDEFromIndex(ybin - 1);
    hpos->GetYaxis()->SetBinLabel(ybin, TString::Format("DE%d", deId));
    hneg->GetYaxis()->SetBinLabel(ybin, TString::Format("DE%d", deId));
  }

  TH2* hvec[2]{hpos, hneg};
  for (int chIndex = 0; chIndex < 10; chIndex++) {
    for (int i = 0; i < 2; i++) {
      std::string histTitle = hvec[i]->GetTitle();
      int deMin = getChamberOffset(chIndex);
      int deMax = deMin + getNumDEinChamber(chIndex);
      hvec[i]->GetXaxis()->SetRangeUser(deMin, deMax);
      hvec[i]->SetTitle(std::format("CH{} {}", chIndex+1, histTitle).c_str());
      hvec[i]->Draw("col");
      c.SaveAs(pdfFileName.c_str());
    }
  }
}

void PlotTrackSlopes(TH2* hxpos, TH2* hxneg, TH2* hypos, TH2* hyneg, TCanvas& c)
{
  if (!hxpos || !hxneg || !hypos || !hyneg) return;

  c.cd();
  c.Clear();
  c.Divide(2,2);

  c.cd(1);
  TH2* h2 = hxpos;
  h2->RebinX(10);
  h2->GetXaxis()->SetRangeUser(-0.2, 0.2);
  h2->GetYaxis()->SetRangeUser(-0.02, 0.02);
  h2->Draw("col");

  c.cd(2);
  h2 = hxneg;
  h2->RebinX(10);
  h2->GetXaxis()->SetRangeUser(-0.2, 0.2);
  h2->GetYaxis()->SetRangeUser(-0.02, 0.02);
  h2->Draw("col");

  c.cd(3);
  h2 = hypos;
  h2->RebinX(10);
  h2->GetXaxis()->SetRangeUser(-0.2, 0.2);
  h2->GetYaxis()->SetRangeUser(-0.02, 0.02);
  h2->Draw("col");

  c.cd(4);
  h2 = hyneg;
  h2->RebinX(10);
  h2->GetXaxis()->SetRangeUser(-0.2, 0.2);
  h2->GetYaxis()->SetRangeUser(-0.02, 0.02);
  h2->Draw("col");

  c.SaveAs(pdfFileName.c_str());
  c.Clear();
}

void PlotChamberResidualsAndDCA(THnSparse* hDCAx, THnSparse* hDCAy, THnSparse* hDEx, THnSparse* hDEy, TCanvas& c)
{
  std::array<std::array<std::pair<double, double>, 4>, 2> DCAx;
  std::array<std::array<std::pair<double, double>, 4>, 2> DCAy;

  std::array<std::array<std::array<std::pair<double, double>, 10>, 4>, 2> meanDx;
  std::array<std::array<std::array<std::pair<double, double>, 10>, 4>, 2> meanDy;

  std::array<std::string, 4> quadrants = {"Q0", "Q1", "Q2", "Q3"};

  c.Clear();
  c.Divide(2, 2);

  std::cout << "Plotting MCH DCA(x)" << std::endl;
  auto* hn = hDCAx;
  if (hn) {
    // axes assignments:
    // 0: momentum
    // 1: quadrant
    // 2: sign
    // 3: DCA
    int pAxisId = 0;
    int quadrantAxisId = 1;
    int signAxisId = 2;
    int dcaAxisId = 3;

    // p >= xx GeV/c
    int momBin = hn->GetAxis(pAxisId)->FindBin(mchMomMin + epsilon);
    int momBinMax = hn->GetAxis(pAxisId)->GetNbins();
    hn->GetAxis(pAxisId)->SetRange(momBin, momBinMax);
    for (int k = 0; k < 2; k++) {
      hn->GetAxis(signAxisId)->SetRange(k + 1, k + 1);
      for (int q = 0; q < quadrants.size(); q++) {
        hn->GetAxis(quadrantAxisId)->SetRange(q + 1, q + 1);
        if (q == 0) c.cd(2);
        if (q == 1) c.cd(1);
        if (q == 2) c.cd(3);
        if (q == 3) c.cd(4);

        auto* proj = hn->Projection(dcaAxisId);
        proj->SetTitle(std::format("MCH DCA(x), Q{} {}", q, (k == 0 ? "+" : "-")).c_str());
        DCAx[k][q] = PlotDCAMCH(proj);
      }
      c.SaveAs(pdfFileName.c_str());
    }
  }

  std::cout << "Plotting MCH DCA(y)" << std::endl;
  hn = hDCAy;
  if (hn) {
    // axes assignments:
    // 0: momentum
    // 1: quadrant
    // 2: sign
    // 3: DCA
    int pAxisId = 0;
    int quadrantAxisId = 1;
    int signAxisId = 2;
    int dcaAxisId = 3;

    // p >= xx GeV/c
    int momBin = hn->GetAxis(pAxisId)->FindBin(mchMomMin + epsilon);
    int momBinMax = hn->GetAxis(pAxisId)->GetNbins();
    hn->GetAxis(pAxisId)->SetRange(momBin, momBinMax);
    for (int k = 0; k < 2; k++) {
      hn->GetAxis(signAxisId)->SetRange(k + 1, k + 1);
      for (int q = 0; q < quadrants.size(); q++) {
        hn->GetAxis(quadrantAxisId)->SetRange(q + 1, q + 1);
        if (q == 0) c.cd(2);
        if (q == 1) c.cd(1);
        if (q == 2) c.cd(3);
        if (q == 3) c.cd(4);

        auto* proj = hn->Projection(dcaAxisId);
        proj->SetTitle(std::format("MCH DCA(y), Q{} {}", q, (k == 0 ? "+" : "-")).c_str());
        DCAy[k][q] = PlotDCAMCH(proj);
        delete proj;
      }
      c.SaveAs(pdfFileName.c_str());
    }
  }

  std::cout << "Plotting MCH DCA(x/y) vs. momentum" << std::endl;
  PlotDCAMCHvsMomentum(hDCAx, hDCAy, c);


  TCanvas c2("c2", "c2", 1200, 800);
  c2.SaveAs("residuals_CH.pdf(");

  // axes assignments (old version)
  // 0: DE index
  // 1: quadrant
  // 2: sign
  // 3: residual

  // axes assignments (new version)
  // 0: residual
  // 1: DE index
  // 2: quadrant
  // 3: sign
  // 4: momentum
  hn = hDEx;
  if (hn) {
    int residualAxisId = 0;
    int deAxisId = 1;
    int quadrantAxisId = 2;
    int signAxisId = 3;
    int pAxisId = 4;

    if (hn->GetNdimensions() == 4) {
      residualAxisId = 3;
      deAxisId = 0;
      quadrantAxisId = 1;
      signAxisId = 2;
      pAxisId = -1;
    }

    hn->GetAxis(0)->SetRange(0, -1);
    hn->GetAxis(1)->SetRange(0, -1);
    hn->GetAxis(2)->SetRange(0, -1);
    hn->GetAxis(3)->SetRange(0, -1);
    if (pAxisId >= 0) {
      // p >= xx GeV/c
      int momBin = hn->GetAxis(pAxisId)->FindBin(mchMomMin + epsilon);
      int momBinMax = hn->GetAxis(pAxisId)->GetNbins();
      hn->GetAxis(pAxisId)->SetRange(momBin, momBinMax);
    }

    // loop on chamber ID
    for (int i = 0; i < 10; i++) {
      int deIdMin = getChamberOffset(i);
      int deIdMax = deIdMin + getNumDEinChamber(i) - 1;
      hn->GetAxis(deAxisId)->SetRange(deIdMin + 1, deIdMax + 1);
      // loop on quadrant
      for (int j = 0; j < hn->GetAxis(quadrantAxisId)->GetNbins(); j++) {
        hn->GetAxis(quadrantAxisId)->SetRange(j + 1, j + 1);
        // loop on charge
        for (int k = 0; k < hn->GetAxis(signAxisId)->GetNbins(); k++) {
          hn->GetAxis(signAxisId)->SetRange(k + 1, k + 1);

          auto* proj = hn->Projection(residualAxisId);
          proj->SetName(std::format("dx_vs_chamber_{}_{}_{}", j, i+1, (k == 0 ? "positive" : "negative")).c_str());
          proj->SetTitle(std::format("#Deltax, Q{} CH{} {}", j, i+1, (k == 0 ? "positive" : "negative")).c_str());
          meanDx[k][j][i] = PlotDXY(proj, c2, "residuals_CH.pdf");
        }
      }
    }
  }

  hn = hDEy;
  if (hn) {
    int residualAxisId = 0;
    int deAxisId = 1;
    int quadrantAxisId = 2;
    int signAxisId = 3;
    int pAxisId = 4;

    if (hn->GetNdimensions() == 4) {
      residualAxisId = 3;
      deAxisId = 0;
      quadrantAxisId = 1;
      signAxisId = 2;
      pAxisId = -1;
    }

    hn->GetAxis(0)->SetRange(0, -1);
    hn->GetAxis(1)->SetRange(0, -1);
    hn->GetAxis(2)->SetRange(0, -1);
    hn->GetAxis(3)->SetRange(0, -1);
    if (pAxisId >= 0) {
      // p >= xx GeV/c
      int momBin = hn->GetAxis(pAxisId)->FindBin(mchMomMin + epsilon);
      int momBinMax = hn->GetAxis(pAxisId)->GetNbins();
      hn->GetAxis(pAxisId)->SetRange(momBin, momBinMax);
    }

    // loop on chamber ID
    for (int i = 0; i < 10; i++) {
      int deIdMin = getChamberOffset(i);
      int deIdMax = deIdMin + getNumDEinChamber(i) - 1;
      hn->GetAxis(deAxisId)->SetRange(deIdMin + 1, deIdMax + 1);
      // loop on quadrant
      for (int j = 0; j < hn->GetAxis(quadrantAxisId)->GetNbins(); j++) {
        hn->GetAxis(quadrantAxisId)->SetRange(j + 1, j + 1);
        // loop on charge
        for (int k = 0; k < hn->GetAxis(signAxisId)->GetNbins(); k++) {
          hn->GetAxis(signAxisId)->SetRange(k + 1, k + 1);

          auto* proj = hn->Projection(residualAxisId);
          proj->SetName(std::format("dy_vs_chamber_{}_{}_{}", j, i+1, (k == 0 ? "positive" : "negative")).c_str());
          proj->SetTitle(std::format("#Deltay, Q{} CH{} {}", j, i+1, (k == 0 ? "positive" : "negative")).c_str());
          meanDy[k][j][i] = PlotDXY(proj, c2, "residuals_CH.pdf");
        }
      }
    }
  }

  c2.Clear();
  c2.SaveAs("residuals_CH.pdf)");

/**/

  double xv[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

  PlotZTrend(10, chamberZ, meanDx[0], DCAx[0], "#Delta(x) vs. chamber z (positive);chamber z (cm); #Delta(x) (cm)", -5.0, 5.0, c);
  PlotZTrend(10, chamberZ, meanDx[0], DCAx[0], "#Delta(x) vs. chamber z (positive);chamber z (cm); #Delta(x) (cm)", -1.0, 1.0, c);
  PlotZTrend(10, chamberZ, meanDx[0], DCAx[0], "#Delta(x) vs. chamber z (positive);chamber z (cm); #Delta(x) (cm)", -0.5, 0.5, c);
  PlotZTrend(10, chamberZ, meanDx[1], DCAx[1], "#Delta(x) vs. chamber z (negative);chamber z (cm); #Delta(x) (cm)", -5.0, 5.0, c);
  PlotZTrend(10, chamberZ, meanDx[1], DCAx[1], "#Delta(x) vs. chamber z (negative);chamber z (cm); #Delta(x) (cm)", -1.0, 1.0, c);
  PlotZTrend(10, chamberZ, meanDx[1], DCAx[1], "#Delta(x) vs. chamber z (negative);chamber z (cm); #Delta(x) (cm)", -0.5, 0.5, c);

  PlotZTrend(10, chamberZ, meanDy[0], DCAy[0], "#Delta(y) vs. chamber z (positive);chamber z (cm); #Delta(y) (cm)", -5.0, 5.0, c);
  PlotZTrend(10, chamberZ, meanDy[0], DCAy[0], "#Delta(y) vs. chamber z (positive);chamber z (cm); #Delta(y) (cm)", -1.0, 1.0, c);
  PlotZTrend(10, chamberZ, meanDy[0], DCAy[0], "#Delta(y) vs. chamber z (positive);chamber z (cm); #Delta(y) (cm)", -0.5, 0.5, c);
  PlotZTrend(10, chamberZ, meanDy[1], DCAy[1], "#Delta(y) vs. chamber z (negative);chamber z (cm); #Delta(y) (cm)", -5.0, 5.0, c);
  PlotZTrend(10, chamberZ, meanDy[1], DCAy[1], "#Delta(y) vs. chamber z (negative);chamber z (cm); #Delta(y) (cm)", -1.0, 1.0, c);
  PlotZTrend(10, chamberZ, meanDy[1], DCAy[1], "#Delta(y) vs. chamber z (negative);chamber z (cm); #Delta(y) (cm)", -0.5, 0.5, c);

  PlotChamberResidualvsMomentum(hDEx, hDEy, c);
}

#define DE_GROUPING_V2

void GetDEShifts(THnSparse* hDEx, THnSparse* hDEy, std::map<int, std::array<std::array<double, 2>, 2>>& shifts, TCanvas& c, std::string residualsDEFileName, std::string residualsGroupsFileName)
{
  //-------------------
  // DE residuals
  //-------------------

#ifdef DE_GROUPING_V1
  std::vector<std::tuple<std::string, int, std::string, int, std::vector<int>>> deGroups {
    {"DE100", 1, "Q0", 0, {100}},
    {"DE101", 1, "Q1", 1, {101}},
    {"DE102", 1, "Q2", 2, {102}},
    {"DE103", 1, "Q3", 3, {103}},
    {"DE200", 2, "Q0", 0, {200}},
    {"DE201", 2, "Q1", 1, {201}},
    {"DE202", 2, "Q2", 2, {202}},
    {"DE203", 2, "Q3", 3, {203}},
    {"DE300", 3, "Q0", 0, {300}},
    {"DE301", 3, "Q1", 1, {301}},
    {"DE302", 3, "Q2", 2, {302}},
    {"DE303", 3, "Q3", 3, {303}},
    {"DE400", 4, "Q0", 0, {400}},
    {"DE401", 4, "Q1", 1, {401}},
    {"DE402", 4, "Q2", 2, {402}},
    {"DE403", 4, "Q3", 3, {403}},
    {"CH5L", 5, "L", 4, {505, 506, 507, 508, 509, 510, 511, 512, 513}},
    {"CH5R", 5, "R", 5, {500, 501, 502, 503, 504, 514, 515, 516, 517}},
    {"CH6L", 6, "L", 4, {605, 606, 607, 608, 609, 610, 611, 612, 613}},
    {"CH6R", 6, "R", 5, {600, 601, 602, 603, 604, 614, 615, 616, 617}},
    {"CH7L", 7, "L", 4, {707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719}},
    {"CH7R", 7, "R", 5, {700, 701, 702, 703, 704, 706, 707, 720, 721, 722, 723, 724, 725}},
    {"CH8L", 8, "L", 4, {807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819}},
    {"CH8R", 8, "R", 5, {800, 801, 802, 803, 804, 806, 807, 820, 821, 822, 823, 824, 825}},
    {"CH9L", 9, "L", 4, {907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919}},
    {"CH9R", 9, "R", 5, {900, 901, 902, 903, 904, 906, 907, 920, 921, 922, 923, 924, 925}},
    {"CH10L", 10, "L", 4, {1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019}},
    {"CH10R", 10, "R", 5, {1000, 1001, 1002, 1003, 1004, 1006, 1007, 1020, 1021, 1022, 1023, 1024, 1025}}
  };

  std:vector<std::string> deGroupNames{
    "Q0", "Q1", "Q2", "Q3", "L", "R"
  };
#endif


#ifdef DE_GROUPING_V2
  std::vector<std::tuple<std::string, int, std::string, int, std::vector<int>>> deGroups {
    {"DE100", 1, "Q0", 0, {100}},
    {"DE101", 1, "Q1", 1, {101}},
    {"DE102", 1, "Q2", 2, {102}},
    {"DE103", 1, "Q3", 3, {103}},
    {"DE200", 2, "Q0", 0, {200}},
    {"DE201", 2, "Q1", 1, {201}},
    {"DE202", 2, "Q2", 2, {202}},
    {"DE203", 2, "Q3", 3, {203}},
    {"DE300", 3, "Q0", 0, {300}},
    {"DE301", 3, "Q1", 1, {301}},
    {"DE302", 3, "Q2", 2, {302}},
    {"DE303", 3, "Q3", 3, {303}},
    {"DE400", 4, "Q0", 0, {400}},
    {"DE401", 4, "Q1", 1, {401}},
    {"DE402", 4, "Q2", 2, {402}},
    {"DE403", 4, "Q3", 3, {403}},

    {"CH5LT", 5, "LT", 4, {505, 506, 507, 508}},
    {"CH5LC", 5, "LC", 5, {509}},
    {"CH5LB", 5, "LB", 6, {510, 511, 512, 513}},

    {"CH5RT", 5, "RT", 7, {501, 502, 503, 504}},
    {"CH5RC", 5, "RC", 8, {500}},
    {"CH5RB", 5, "RB", 9, {514, 515, 516, 517}},

    {"CH6LT", 6, "LT", 4, {605, 606, 607, 608}},
    {"CH6LC", 6, "LC", 5, {609}},
    {"CH6LB", 6, "LB", 6, {610, 611, 612, 613}},

    {"CH6RT", 6, "RT", 7, {601, 602, 603, 604}},
    {"CH6RC", 6, "RC", 8, {600}},
    {"CH6RB", 6, "RB", 9, {614, 615, 616, 617}},

    {"CH7LT", 7, "LT", 4, {707, 708, 709, 710, 711, 712}},
    {"CH7LC", 7, "LC", 5, {713}},
    {"CH7LB", 7, "LB", 6, {714, 715, 716, 717, 718, 719}},

    {"CH7RT", 7, "RT", 7, {701, 702, 703, 704, 705, 706, 707}},
    {"CH7RC", 7, "RC", 8, {700}},
    {"CH7RB", 7, "RB", 9, {720, 721, 722, 723, 724, 725}},

    {"CH8LT", 8, "LT", 4, {807, 808, 809, 810, 811, 812}},
    {"CH8LC", 8, "LC", 5, {813}},
    {"CH8LB", 8, "LB", 6, {814, 815, 816, 817, 818, 819}},

    {"CH8RT", 8, "RT", 7, {801, 802, 803, 804, 805, 806, 807}},
    {"CH8RC", 8, "RC", 8, {800}},
    {"CH8RB", 8, "RB", 9, {820, 821, 822, 823, 824, 825}},

    {"CH9LT", 9, "LT", 4, {907, 908, 909, 910, 911, 912}},
    {"CH9LC", 9, "LC", 5, {913}},
    {"CH9LB", 9, "LB", 6, {914, 915, 916, 917, 918, 919}},

    {"CH9RT", 9, "RT", 7, {901, 902, 903, 904, 905, 906, 907}},
    {"CH9RC", 9, "RC", 8, {900}},
    {"CH9RB", 9, "RB", 9, {920, 921, 922, 923, 924, 925}},

    {"CH10LT", 10, "LT", 4, {1007, 1008, 1009, 1010, 1011, 1012}},
    {"CH10LC", 10, "LC", 5, {1013}},
    {"CH10LB", 10, "LB", 6, {1014, 1015, 1016, 1017, 1018, 1019}},

    {"CH10RT", 10, "RT", 7, {1001, 1002, 1003, 1004, 1005, 1006, 1007}},
    {"CH10RC", 10, "RC", 8, {1000}},
    {"CH10RB", 10, "RB", 9, {1020, 1021, 1022, 1023, 1024, 1025}}
  };

  std:vector<std::string> deGroupNames{
    "Q0", "Q1", "Q2", "Q3", "LT", "LC", "LB", "RT", "RC", "RB"
  };
#endif

  std::array<std::vector<std::array<TH1*, 2>>, 2> hDEResiduals;

  // one map per coordinate
  // map index is the group name (3rd column)
  // each element of the map is a matrix of 2x4 vectors, with z coordinates + errors and residuals + errors for each charge sign
  std::array<std::unordered_map<std::string, std::array<std::array<std::vector<double>, 4>, 2>>, 2> groupResiduals;
  std::array<std::unordered_map<std::string, std::array<std::vector<double>, 4>>, 2> groupResidualsAverage;
  std::array<std::unordered_map<std::string, std::array<std::vector<double>, 4>>, 2> groupResidualsDelta;

  std::string coordinates[2] = {"x", "y"};

  TCanvas c3("c3", "c3", 1200, 800);
  c3.SaveAs((residualsDEFileName + "(").c_str());

  TCanvas c5("c5", "c5", 1200, 800);
  c5.SaveAs((residualsGroupsFileName + "(").c_str());


  int colors[10] = {kBlue, kRed+1, kOrange, kCyan+1, kBlack, kGreen+2, kMagenta+1, kYellow+1, kRed-6, kBlue-6};
  int markers[2] = {kStar, kCircle};
  int lineStyles[2] = {kSolid, kDashed};



  // axis assignments (old version with 4 axes):
  // 0: DE index
  // 1: quadrant
  // 2: charge sign
  // 3: residual

  // axis assignments (new version with 5 axes):
  // 0: residual
  // 1: DE index
  // 2: quadrant
  // 3: charge sign
  // 4: momentum

  // -----------------
  // X direction

  THnSparse* hDE[2] = { hDEx, hDEy };
  if (hDE[0] && hDE[1]) {

    int residualAxisId = 0;
    int deAxisId = 1;
    int quadrantAxisId = 2;
    int signAxisId = 3;
    int pAxisId = 4;

    if (hDE[0]->GetNdimensions() == 4) {
      residualAxisId = 3;
      deAxisId = 0;
      quadrantAxisId = 1;
      signAxisId = 2;
      pAxisId = -1;
    }

    for (int xy = 0; xy < 2; xy++) {
      hDE[xy]->GetAxis(0)->SetRange(0, -1);
      hDE[xy]->GetAxis(1)->SetRange(0, -1);
      hDE[xy]->GetAxis(2)->SetRange(0, -1);
      hDE[xy]->GetAxis(3)->SetRange(0, -1);
      int momBin = 0;
      int momBinMax = -1;
      if (pAxisId >= 0) {
        // p >= xx GeV/c
        momBin = hDE[xy]->GetAxis(pAxisId)->FindBin(mchMomMin + epsilon);
        momBinMax = hDE[xy]->GetAxis(pAxisId)->GetNbins();
        hDE[xy]->GetAxis(pAxisId)->SetRange(momBin, momBinMax);
      }

      TH1* histDxyVsDE[2] = {nullptr, nullptr};
      for (int k = 0; k < hDE[xy]->GetAxis(signAxisId)->GetNbins(); k++) {
        //hDE[xy]->GetAxis(2)->SetRange(k + 1, k + 1);
        //hDE[xy]->GetAxis(0)->SetRange(1, hDE[xy]->GetAxis(0)->GetNbins());
        histDxyVsDE[k] = hDE[xy]->Projection(deAxisId);
        histDxyVsDE[k]->Reset();
        histDxyVsDE[k]->SetName(std::format("d{}_vs_de_{}", (xy == 0 ? "x" : "y"), (k == 0 ? "positive" : "negative")).c_str());
        histDxyVsDE[k]->SetTitle(std::format("#Delta{} vs. DE", (xy == 0 ? "x" : "y")).c_str());
      }
      for (int i = 0; i < hDE[xy]->GetAxis(deAxisId)->GetNbins(); i++) {
        hDE[xy]->GetAxis(deAxisId)->SetRange(i + 1, i + 1);
        for (int k = 0; k < hDE[xy]->GetAxis(signAxisId)->GetNbins(); k++) {
          hDE[xy]->GetAxis(signAxisId)->SetRange(k + 1, k + 1);

          auto* proj = hDE[xy]->Projection(residualAxisId);
          proj->SetName(std::format("d{}_vs_de_{}_{}", (xy == 0 ? "x" : "y"), i+1, (k == 0 ? "positive" : "negative")).c_str());
          proj->SetTitle(std::format("#Delta{}, DE{} {}", (xy == 0 ? "x" : "y"), getDEFromIndex(i), (k == 0 ? "positive" : "negative")).c_str());
          std::cout << std::format("[TOTO] DE{}({}) momRange={},{} signBin={} deBin={}  entries={}", getDEFromIndex(i), (xy == 0 ? "x" : "y"), momBin, momBinMax, k+1, i+1, proj->GetEntries()) << std::endl;
          auto mean = PlotDXY(proj, c3, residualsDEFileName);
          histDxyVsDE[k]->SetBinContent(i + 1, mean.first);
          histDxyVsDE[k]->SetBinError(i + 1, mean.second);
        }
      }
      c.Clear();
      c.cd();
      histDxyVsDE[0]->SetLineColor(kBlue);
      histDxyVsDE[0]->SetMinimum(-1.0);
      histDxyVsDE[0]->SetMaximum(1.0);
      histDxyVsDE[0]->Draw("E");
      histDxyVsDE[1]->SetLineColor(kRed);
      histDxyVsDE[1]->Draw("E same");
      c.SaveAs(pdfFileName.c_str());

      // set DE names in X axis
      for (int xbin = 1; xbin <= histDxyVsDE[0]->GetXaxis()->GetNbins(); xbin++) {
        int deId = getDEFromIndex(xbin - 1);
        histDxyVsDE[0]->GetXaxis()->SetBinLabel(xbin, TString::Format("DE%d", deId));
      }

      std::string histTitle = histDxyVsDE[0]->GetTitle();
      for (int chIndex = 0; chIndex < 10; chIndex++) {
        int deMin = getChamberOffset(chIndex);
        int deMax = deMin + getNumDEinChamber(chIndex);
        histDxyVsDE[0]->GetXaxis()->SetRangeUser(deMin, deMax);
        histDxyVsDE[0]->SetTitle((histTitle + "(CH" + std::to_string(chIndex+1) + ")").c_str());
        histDxyVsDE[0]->SetMinimum(-5.0);
        histDxyVsDE[0]->SetMaximum(5.0);
        c.SaveAs(pdfFileName.c_str());

        histDxyVsDE[0]->SetMinimum(-1.0);
        histDxyVsDE[0]->SetMaximum(1.0);
        c.SaveAs(pdfFileName.c_str());

        histDxyVsDE[0]->SetMinimum(-0.5);
        histDxyVsDE[0]->SetMaximum(0.5);
        c.SaveAs(pdfFileName.c_str());

        if (chIndex < 4) {
          histDxyVsDE[0]->SetMinimum(-0.2);
          histDxyVsDE[0]->SetMaximum(0.2);
          c.SaveAs(pdfFileName.c_str());
        }
      }
    }

    //   std::array<std::vector<std::array<TH1*, 2>>, 2> hDEResiduals;
    //hn->GetAxis(1)->SetRange(1, hn->GetAxis(1)->GetNbins());
    for (const auto& group : deGroups) {
      auto groupId = std::get<2>(group);
      int chamberId = std::get<1>(group) - 1;

      // loop over coordinate
      for (int xy = 0; xy < 2; xy++) {
        hDEResiduals[xy].emplace_back(std::array<TH1*, 2>{nullptr, nullptr});
        // loop over charge sign
        for (int charge = 0; charge < 2; charge++) {
          hDE[xy]->GetAxis(signAxisId)->SetRange(charge + 1, charge + 1);
          // loop over DE ids
          for (auto deId : std::get<4>(group)) {
            auto deIndex = getDEindex(deId);
            hDE[xy]->GetAxis(deAxisId)->SetRange(deIndex + 1, deIndex + 1);
            TH1* proj = hDE[xy]->Projection(residualAxisId);
            if (hDEResiduals[xy].back()[charge]) {
              hDEResiduals[xy].back()[charge]->Add(proj);
              delete proj;
            } else {
              proj->SetTitle(std::format("{} - #Delta{}, {}", std::get<0>(group), (xy == 0 ? "x" : "y"), (charge == 0 ? "positive" : "negative")).c_str());
              hDEResiduals[xy].back()[charge] = proj;
            }
          }

          std::cout << std::format("hDEResiduals[xy].size(): {}", hDEResiduals[xy].size()) << std::endl;

          // compute average shift for this group, charge and coordinate
          auto mean = PlotDXY(hDEResiduals[xy].back()[charge], c5, residualsGroupsFileName);
          // assign the shift to all the Detection Eelements in the group
          if (!std::isnan(mean.first) && !std::isnan(mean.second)) {
            for (auto deId : std::get<4>(group)) {
              shifts[deId][xy][charge] = mean.first;
            }

            groupResiduals[xy][groupId][charge][0].push_back(chamberZ[chamberId]);
            groupResiduals[xy][groupId][charge][1].push_back(0);
            groupResiduals[xy][groupId][charge][2].push_back(mean.first);
            groupResiduals[xy][groupId][charge][3].push_back(mean.second);

            std::cout << std::format("Shifts for group {} ({} {}) -> {:+0.3f}",
                std::get<0>(group), (xy == 0 ? "x" : "y"), (charge == 0 ? "+" : "-"), mean.first, mean.second) << std::endl;
          }
        }

        if (groupResiduals[xy][groupId][0][2].empty()) {
          continue;
        }

        std::cout << "Computing average and delta" << std::endl;
        std::cout << std::format("groupResiduals[xy][groupId][0][2].size(): {}", groupResiduals[xy][groupId][0][2].size()) << std::endl;
        double residualPos = groupResiduals[xy][groupId][0][2].back();
        double residualNeg = groupResiduals[xy][groupId][1][2].back();
        double residualPosErr = groupResiduals[xy][groupId][0][3].back();
        double residualNegErr = groupResiduals[xy][groupId][1][3].back();
        double residualAverage = (residualPos + residualNeg) / 2.0;
        double residualDelta = residualPos - residualNeg;
        double residualErr = std::sqrt(residualPosErr*residualPosErr + residualNegErr*residualNegErr);

        groupResidualsAverage[xy][groupId][0].push_back(chamberZ[chamberId]);
        groupResidualsAverage[xy][groupId][1].push_back(0);
        groupResidualsAverage[xy][groupId][2].push_back(residualAverage);
        groupResidualsAverage[xy][groupId][3].push_back(residualErr);

        groupResidualsDelta[xy][groupId][0].push_back(chamberZ[chamberId]);
        groupResidualsDelta[xy][groupId][1].push_back(0);
        groupResidualsDelta[xy][groupId][2].push_back(residualDelta);
        groupResidualsDelta[xy][groupId][3].push_back(residualErr);
        std::cout << "Computing average and delta done." << std::endl;
      }
    }

    // loop over coordinate
    for (int xy = 0; xy < 2; xy++) {
      // plot the trends for this coordinate
      TMultiGraph* mg = new TMultiGraph();
      mg->SetTitle(Form("#Delta(%s) for groups;z (cm); #Delta(%s) (cm)", coordinates[xy].c_str(), coordinates[xy].c_str()));
      TGraphErrors* gr = new TGraphErrors();
      gr->AddPoint(500, 1000);
      gr->AddPoint(2000, 1000);
      mg->Add(gr,"lp");

      TLegend* legend = new TLegend(0.7, 0.1, 0.98, 0.9);
      legend->SetNColumns(2);

      // loop over groups
      int groupId = 0;
      //for (const auto& [groupName, residuals] : groupResiduals[xy]) {
      for (const auto& groupName : deGroupNames) {
        const auto& residuals = groupResiduals[xy][groupName];
        if (residuals[0].empty()) {
          continue;
        }
        // loop over charge sign
        for (int charge = 0; charge < 2; charge++) {
          std::cout << std::format("Plotting {}{}", groupName, (charge == 0 ? " (+)" : " (-)")) << std::endl;
          TGraphErrors* gr = new TGraphErrors(residuals[charge][0].size(),
              residuals[charge][0].data(),
              residuals[charge][2].data(),
              residuals[charge][1].data(),
              residuals[charge][3].data());
          gr->SetLineColor(colors[groupId]);
          gr->SetMarkerColor(colors[groupId]);
          gr->SetMarkerStyle(markers[charge]);
          gr->SetMarkerSize(2);
          gr->SetLineStyle(lineStyles[charge]);
          mg->Add(gr,"pl");
          auto* entry = legend->AddEntry(gr, (groupName + (charge == 0 ? " (+)" : " (-)")).c_str(), "P");
          entry->SetTextColor(colors[groupId]);

        }
        groupId += 1;
      }

      c5.cd();
      //mg->Draw("a");
      //mg->SetMinimum(-5);
      //mg->SetMaximum(5);
      //legend->Draw();
      //c5.SaveAs(residualsGroupsFileName.c_str());

      mg->Draw("a");
      mg->SetMinimum(-1);
      mg->SetMaximum(1);
      legend->Draw();
      c5.SaveAs(residualsGroupsFileName.c_str());
    }

    // loop over coordinate
    for (int xy = 0; xy < 2; xy++) {
      // plot the trends for this coordinate
      TMultiGraph* mg = new TMultiGraph();
      mg->SetTitle(Form("#Delta(%s) charge average for groups;z (cm); #Delta(%s) (cm)", coordinates[xy].c_str(), coordinates[xy].c_str()));
      TGraphErrors* gr = new TGraphErrors();
      gr->AddPoint(500, 1000);
      gr->AddPoint(2000, 1000);
      mg->Add(gr,"lp");

      TLegend* legend = new TLegend(0.7, 0.1, 0.98, 0.9);
      legend->SetNColumns(1);

      // loop over groups
      int groupId = 0;
      for (const auto& groupName : deGroupNames) {
        const auto& residuals = groupResidualsAverage[xy][groupName];
        if (residuals[0].empty()) {
          continue;
        }
        TGraphErrors* gr = new TGraphErrors(residuals[0].size(),
            residuals[0].data(),
            residuals[2].data(),
            residuals[1].data(),
            residuals[3].data());
        gr->SetLineColor(colors[groupId]);
        gr->SetMarkerColor(colors[groupId]);
        gr->SetMarkerStyle(markers[0]);
        gr->SetMarkerSize(2);
        gr->SetLineStyle(lineStyles[0]);
        mg->Add(gr,"pl");
        auto* entry = legend->AddEntry(gr, groupName.c_str(), "P");
        entry->SetTextColor(colors[groupId]);

        groupId += 1;
      }

      c5.cd();
      mg->Draw("a");
      mg->SetMinimum(-1);
      mg->SetMaximum(1);
      legend->Draw();
      c5.SaveAs(residualsGroupsFileName.c_str());
    }

    // loop over coordinate
    for (int xy = 0; xy < 2; xy++) {
      // plot the trends for this coordinate
      TMultiGraph* mg = new TMultiGraph();
      mg->SetTitle(Form("#Delta(%s) charge difference for groups;z (cm); #Delta(%s) (cm)", coordinates[xy].c_str(), coordinates[xy].c_str()));
      TGraphErrors* gr = new TGraphErrors();
      gr->AddPoint(500, 1000);
      gr->AddPoint(2000, 1000);
      mg->Add(gr,"lp");

      TLegend* legend = new TLegend(0.7, 0.1, 0.98, 0.9);
      legend->SetNColumns(1);

      // loop over groups
      int groupId = 0;
      for (const auto& groupName : deGroupNames) {
        const auto& residuals = groupResidualsDelta[xy][groupName];
        if (residuals[0].empty()) {
          continue;
        }
        TGraphErrors* gr = new TGraphErrors(residuals[0].size(),
            residuals[0].data(),
            residuals[2].data(),
            residuals[1].data(),
            residuals[3].data());
        gr->SetLineColor(colors[groupId]);
        gr->SetMarkerColor(colors[groupId]);
        gr->SetMarkerStyle(markers[0]);
        gr->SetMarkerSize(2);
        gr->SetLineStyle(lineStyles[0]);
        mg->Add(gr,"pl");
        auto* entry = legend->AddEntry(gr, groupName.c_str(), "P");
        entry->SetTextColor(colors[groupId]);

        groupId += 1;
      }

      c5.cd();
      mg->Draw("a");
      mg->SetMinimum(-1);
      mg->SetMaximum(1);
      legend->Draw();
      c5.SaveAs(residualsGroupsFileName.c_str());
    }
  }

  c3.Clear();
  c3.SaveAs((residualsDEFileName + ")").c_str());

  c5.Clear();
  c5.SaveAs((residualsGroupsFileName + ")").c_str());
}


void muonGlobalAlignmentMftMchResiduals(const char* _rootFileName = "AnalysisResults.root", const char* _pdfFileName = "mftMchResiduals.pdf")
{
  //fAnalysisResults = new TFile("AnalysisResults.root");
  //fAnalysisResults = new TFile("AnalysisResults/AnalysisResultsFull.root");
  fAnalysisResults = new TFile(_rootFileName);
  pdfFileName = _pdfFileName;

  //int wagonId = -1;

  // LHC24am without MFT corrections
  //int wagonId = 45637;
  // LHC24an with MFT corrections
  //int wagonId = 47795;

  // LHC24an without MFT corrections
  //int wagonId = 45637;
  // LHC24an with MFT corrections, p > 30 GeV/c
  //int wagonId = 47795;
  // LHC24an with MFT corrections, p > 15 GeV/c
  //int wagonId = 47994;

  // ppRef without MFT corrections
  //int wagonId = 47792;
  // ppRef with MFT corrections
  //int wagonId = 47794;

  // ppRef without MFT corrections - low momentum version
  //int wagonId = 47890;
  // ppRef with MFT corrections
  //int wagonId = 47888;

  std::string taskPath = "muon-global-alignment"; // = (wagonId > 0) ? std::format("muon-global-alignment_id{}", wagonId) : "muon-global-alignment";
  //taskPath = "muon-global-alignment_id47994";
  //taskPath = "muon-global-alignment_p-15-pt-2_id47994";
  //taskPath = "muon-global-alignment_p-20-pt-2_id47994";
  //taskPath = "muon-global-alignment_p-50-pt-4_id47994";

  //taskPath = "muon-global-alignment_id48150";
  //taskPath = "muon-global-alignment_id48242";
  //taskPath = "muon-global-alignment_id48243";

  std::array<std::string, 4> quadrants = {"Q0", "Q1", "Q2", "Q3"};
  //std::array<std::string, 1> quadrants = {"Q0"};


  std::array<std::array<std::pair<double, double>, 4>, 2> DCAx;
  std::array<std::array<std::pair<double, double>, 4>, 2> DCAy;

  std::array<std::array<std::array<std::pair<double, double>, 10>, 4>, 2> meanDx;
  std::array<std::array<std::array<std::pair<double, double>, 10>, 4>, 2> meanDy;

  std::array<std::array<std::pair<double, double>, 154>, 2> meanDxVsDE;
  std::array<std::array<std::pair<double, double>, 154>, 2> meanDyVsDE;

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
  gStyle->SetOptFit(1111);
  
  TCanvas c("c", "c", 1200, 800);
  c.SaveAs((pdfFileName + "(").c_str());

  TCanvas c4("c4", "c4", 1200, 800);
  c4.SaveAs("residuals_tracks.pdf(");

  TCanvas c3("c3", "c3", 1200, 800);
  c3.SaveAs("residuals_DE.pdf(");

  TCanvas c5("c5", "c5", 1200, 800);
  c5.SaveAs("residuals_groups.pdf(");

  c.cd();
  TH1* h1 = GetTH1(fAnalysisResults, (taskPath + "/DCA/vertex_z").c_str());
  if (h1) {
    h1->Draw();
    c.SaveAs(pdfFileName.c_str());
  }

  TH2* h2 = GetTH2(fAnalysisResults, (taskPath + "/DCA/MCH/DCA_y_vs_x").c_str());
  if (h2) {
    h2->Draw("colz");
    c.SaveAs(pdfFileName.c_str());
  }
  //c.SaveAs((pdfFileName + ")").c_str());
  //return;


  std::cout << "Plotting MCH track slopes" << std::endl;
  TH2* hxpos = GetTH2(fAnalysisResults, (taskPath + "/residuals/track_dslopex_vs_slopex_pos").c_str());
  TH2* hxneg = GetTH2(fAnalysisResults, (taskPath + "/residuals/track_dslopex_vs_slopex_neg").c_str());
  TH2* hypos = GetTH2(fAnalysisResults, (taskPath + "/residuals/track_dslopey_vs_slopey_pos").c_str());
  TH2* hyneg = GetTH2(fAnalysisResults, (taskPath + "/residuals/track_dslopey_vs_slopey_neg").c_str());
  PlotTrackSlopes(hxpos, hxneg, hypos, hyneg, c);

  std::cout << "Plotting DE correlations" << std::endl;
  TH2* hpos = GetTH2(fAnalysisResults, (taskPath + "/residuals/DE_correlation_pos").c_str());
  TH2* hneg = GetTH2(fAnalysisResults, (taskPath + "/residuals/DE_correlation_neg").c_str());
  PlotDECorrelations(hpos, hneg, c);

  auto* hn = GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/track_dx").c_str());
  if (hn) {
    TH2* hResidualXp = new TH2F("hResidualXp", "Track #Deltax, positive",
        hn->GetAxis(0)->GetNbins(), hn->GetAxis(0)->GetXmin(), hn->GetAxis(0)->GetXmax(),
        hn->GetAxis(1)->GetNbins(), hn->GetAxis(1)->GetXmin(), hn->GetAxis(1)->GetXmax());
    TH2* hResidualXn = new TH2F("hResidualXn", "Track #Deltax, negative",
        hn->GetAxis(0)->GetNbins(), hn->GetAxis(0)->GetXmin(), hn->GetAxis(0)->GetXmax(),
        hn->GetAxis(1)->GetNbins(), hn->GetAxis(1)->GetXmin(), hn->GetAxis(1)->GetXmax());
    std::array<TH2*, 2> hResidualX{hResidualXp, hResidualXn};
    //hResidualX->Reset();
    for (int i = 0; i < hn->GetAxis(0)->GetNbins(); i++) {
      hn->GetAxis(0)->SetRange(i + 1, i + 1);
      for (int j = 0; j < hn->GetAxis(1)->GetNbins(); j++) {
        hn->GetAxis(1)->SetRange(j + 1, j + 1);
        for (int k = 0; k < hn->GetAxis(2)->GetNbins(); k++) {
          hn->GetAxis(2)->SetRange(k + 1, k + 1);

          hResidualX[k]->SetBinContent(i + 1, j + 1, -100);
          hResidualX[k]->SetBinError(i + 1, j + 1, 0);
          auto* proj = hn->Projection(3);
          if (proj->GetEntries() < 10) continue;
          proj->SetName(std::format("track_dx_{}_{}_{}", j, i+1, (k == 0 ? "positive" : "negative")).c_str());
          proj->SetTitle(std::format("Track #Deltax, X={} Y={} Q={}", i, j, (k == 0 ? "positive" : "negative")).c_str());
          auto mean = PlotDXY(proj, c4, "residuals_tracks.pdf");
          hResidualX[k]->SetBinContent(i + 1, j + 1, mean.first);
          hResidualX[k]->SetBinError(i + 1, j + 1, mean.second);
        }
      }
    }
    c.cd();
    hResidualX[0]->SetMinimum(-1.1);
    hResidualX[0]->SetMaximum(1.1);
    hResidualX[0]->Draw("colz");
    c.SaveAs(pdfFileName.c_str());
    hResidualX[1]->SetMinimum(-1.1);
    hResidualX[1]->SetMaximum(1.1);
    hResidualX[1]->Draw("colz");
    c.SaveAs(pdfFileName.c_str());
  }

  hn = GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/track_dy").c_str());
  if (hn) {
    TH2* hResidualYp = new TH2F("hResidualXp", "Track #Deltay, positive",
        hn->GetAxis(0)->GetNbins(), hn->GetAxis(0)->GetXmin(), hn->GetAxis(0)->GetXmax(),
        hn->GetAxis(1)->GetNbins(), hn->GetAxis(1)->GetXmin(), hn->GetAxis(1)->GetXmax());
    TH2* hResidualYn = new TH2F("hResidualXn", "Track #Deltay, negative",
        hn->GetAxis(0)->GetNbins(), hn->GetAxis(0)->GetXmin(), hn->GetAxis(0)->GetXmax(),
        hn->GetAxis(1)->GetNbins(), hn->GetAxis(1)->GetXmin(), hn->GetAxis(1)->GetXmax());
    std::array<TH2*, 2> hResidualY{hResidualYp, hResidualYn};
    for (int i = 0; i < hn->GetAxis(0)->GetNbins(); i++) {
      hn->GetAxis(0)->SetRange(i + 1, i + 1);
      for (int j = 0; j < hn->GetAxis(1)->GetNbins(); j++) {
        hn->GetAxis(1)->SetRange(j + 1, j + 1);
        for (int k = 0; k < hn->GetAxis(2)->GetNbins(); k++) {
          hn->GetAxis(2)->SetRange(k + 1, k + 1);

          hResidualY[k]->SetBinContent(i + 1, j + 1, -100);
          hResidualY[k]->SetBinError(i + 1, j + 1, 0);
          auto* proj = hn->Projection(3);
          if (proj->GetEntries() < 10) continue;
          proj->SetName(std::format("track_dy_{}_{}_{}", j, i+1, (k == 0 ? "positive" : "negative")).c_str());
          proj->SetTitle(std::format("Track #Deltay, X={} Y={} Q={}", i, j, (k == 0 ? "positive" : "negative")).c_str());
          auto mean = PlotDXY(proj, c4, "residuals_tracks.pdf");
          hResidualY[k]->SetBinContent(i + 1, j + 1, mean.first);
          hResidualY[k]->SetBinError(i + 1, j + 1, mean.second);
        }
      }
    }
    c.cd();
    hResidualY[0]->SetMinimum(-1.1);
    hResidualY[0]->SetMaximum(1.1);
    hResidualY[0]->Draw("colz");
    c.SaveAs(pdfFileName.c_str());
    hResidualY[1]->SetMinimum(-1.1);
    hResidualY[1]->SetMaximum(1.1);
    hResidualY[1]->Draw("colz");
    c.SaveAs(pdfFileName.c_str());
  }

  c4.Clear();
  c4.SaveAs("residuals_tracks.pdf)");



  c.Clear();
  c.SaveAs(pdfFileName.c_str());
  PlotChamberResidualsAndDCA(
      GetTHnSparse(fAnalysisResults, (taskPath + "/DCA/MCH/DCA_x_vs_sign_vs_quadrant_vs_mom").c_str()),
      GetTHnSparse(fAnalysisResults, (taskPath + "/DCA/MCH/DCA_y_vs_sign_vs_quadrant_vs_mom").c_str()),
      GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/dx_vs_de").c_str()),
      GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/dy_vs_de").c_str()),
      c
  );

  std::map<int, std::array<std::array<double, 2>, 2>> deShifts;
  GetDEShifts(
      GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/dx_vs_de").c_str()),
      GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/dy_vs_de").c_str()),
      deShifts, c, "residuals_de.pdf", "residuals_groups.pdf"
  );

  std::map<int, std::array<std::array<double, 2>, 2>> deShiftsCorr;
  c.Clear();
  c.SaveAs(pdfFileName.c_str());
  PlotChamberResidualsAndDCA(
      GetTHnSparse(fAnalysisResults, (taskPath + "/DCA/MCH/DCA_x_vs_sign_vs_quadrant_vs_mom_corr").c_str()),
      GetTHnSparse(fAnalysisResults, (taskPath + "/DCA/MCH/DCA_y_vs_sign_vs_quadrant_vs_mom_corr").c_str()),
      GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/dx_vs_de_corr").c_str()),
      GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/dy_vs_de_corr").c_str()),
      c
  );

  GetDEShifts(
      GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/dx_vs_de_corr").c_str()),
      GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/dy_vs_de_corr").c_str()),
      deShiftsCorr, c, "residuals_de_corr.pdf", "residuals_groups_corr.pdf"
  );

  // update alignment corrections
  std::ofstream correctionsJson("corrections.json");
  std::ofstream halfCorrectionsJson("half-corrections.json");
  correctionsJson << "{" << std::endl;
  halfCorrectionsJson << "{" << std::endl;
  // loop over DE ids
  for (auto [deId, shift] : deShifts) {
    std::cout << std::format("Shifts for DE{}:\n  X[+, -, avg]: {:+0.3f}, {:+0.3f}, {:+0.3f}\n  Y[+, -, avg]: {:+0.3f}, {:+0.3f}, {:+0.3f}",
        deId, shift[0][0], shift[0][1], (shift[0][0] + shift[0][1]) / 2.0,
        shift[1][0], shift[1][1], (shift[1][0] + shift[1][1]) / 2.0) << std::endl;
    auto shiftHalf = shift;
    // add residual shift after alignment corrections
    const auto& shiftCorr = deShiftsCorr.find(deId);
    if (shiftCorr != deShiftsCorr.end()) {
      std::cout << std::format("  Shifts corrections for DE{}:\n    X[+, -]: {:+0.3f}, {:+0.3f}\n    Y[+, -]: {:+0.3f}, {:+0.3f}",
          deId, shiftCorr->second[0][0], shiftCorr->second[0][1],
          shiftCorr->second[1][0], shiftCorr->second[1][1]) << std::endl;
      // loop over coordinate
      for (int xy = 0; xy < 2; xy++) {
        // loop over charge sign
        for (int charge = 0; charge < 2; charge++) {
          shift[xy][charge] += shiftCorr->second[xy][charge];
          shiftHalf[xy][charge] += shiftCorr->second[xy][charge] / 2.0;
        }
      }
      std::cout << std::format("  Final shifts for DE{}:\n    X[+, -, avg]: {:+0.3f}, {:+0.3f}, {:+0.3f}\n    Y[+, -, avg]: {:+0.3f}, {:+0.3f}, {:+0.3f}",
          deId, shift[0][0], shift[0][1], (shift[0][0] + shift[0][1]) / 2.0,
          shift[1][0], shift[1][1], (shift[1][0] + shift[1][1]) / 2.0) << std::endl;
    }

    double averageCorrectionX = -1.0 * (shift[0][0] + shift[0][1]) / 2.0;
    double averageCorrectionY = -1.0 * (shift[1][0] + shift[1][1]) / 2.0;
    correctionsJson << "  \"" << deId << "\": {" << std::endl;
    correctionsJson << "    \"x\": " << averageCorrectionX << ","
        << " \"y\": " << averageCorrectionY << ","
        << " \"z\": 0, \"yaw\": 0, \"pitch\": 0, \"roll\": 0" << std::endl;
    correctionsJson << "  }," << std::endl;

    double averageHalfCorrectionX = -1.0 * (shiftHalf[0][0] + shiftHalf[0][1]) / 2.0;
    double averageHalfCorrectionY = -1.0 * (shiftHalf[1][0] + shiftHalf[1][1]) / 2.0;
    halfCorrectionsJson << "  \"" << deId << "\": {" << std::endl;
    halfCorrectionsJson << "    \"x\": " << averageHalfCorrectionX << ","
        << " \"y\": " << averageHalfCorrectionY << ","
        << " \"z\": 0, \"yaw\": 0, \"pitch\": 0, \"roll\": 0" << std::endl;
    halfCorrectionsJson << "  }," << std::endl;
  }
  correctionsJson << "  \"0\": {}" << std::endl;
  correctionsJson << "}" << std::endl;
  halfCorrectionsJson << "  \"0\": {}" << std::endl;
  halfCorrectionsJson << "}" << std::endl;


/*
  c.Clear();
  c.Divide(2, 2);

  std::cout << "Plotting MCH DCA(x)" << std::endl;
  auto* hn = GetTHnSparse(fAnalysisResults, (taskPath + "/DCA/MCH/DCA_x_vs_sign_vs_quadrant_vs_mom_corr").c_str());
  if (hn) {
    // axes assignments:
    // 0: momentum
    // 1: quadrant
    // 2: sign
    // 3: DCA
    int pAxisId = 0;
    int quadrantAxisId = 1;
    int signAxisId = 2;
    int dcaAxisId = 3;

    // p >= xx GeV/c
    int momBin = hn->GetAxis(pAxisId)->FindBin(mchMomMin + epsilon);
    int momBinMax = hn->GetAxis(pAxisId)->GetNbins();
    hn->GetAxis(pAxisId)->SetRange(momBin, momBinMax);
    for (int k = 0; k < 2; k++) {
      hn->GetAxis(signAxisId)->SetRange(k + 1, k + 1);
      for (int q = 0; q < quadrants.size(); q++) {
        hn->GetAxis(quadrantAxisId)->SetRange(q + 1, q + 1);
        if (q == 0) c.cd(2);
        if (q == 1) c.cd(1);
        if (q == 2) c.cd(3);
        if (q == 3) c.cd(4);

        auto* proj = hn->Projection(dcaAxisId);
        proj->SetTitle(std::format("MCH DCA(x), Q{} {}", q, (k == 0 ? "+" : "-")).c_str());
        DCAx[k][q] = PlotDCAMCH(proj);
      }
      c.SaveAs(pdfFileName.c_str());
    }

    // DCA as function of momentum
  } else {
    hn = GetTHnSparse(fAnalysisResults, (taskPath + "/DCA/MCH/DCA_x_vs_sign_vs_quadrant_vs_vz").c_str());
    for (int k = 0; k < 2; k++) {
      hn->GetAxis(2)->SetRange(k + 1, k + 1);
      for (int q = 0; q < quadrants.size(); q++) {
        hn->GetAxis(1)->SetRange(q + 1, q + 1);
        if (q == 0) c.cd(2);
        if (q == 1) c.cd(1);
        if (q == 2) c.cd(3);
        if (q == 3) c.cd(4);

        auto* proj = hn->Projection(3);
        proj->SetTitle(std::format("MCH DCA(x), Q{} {}", q, (k == 0 ? "+" : "-")).c_str());
        DCAx[k][q] = PlotDCAMCH(proj);
        delete proj;
      }
      c.SaveAs(pdfFileName.c_str());
    }
  }

  std::cout << "Plotting MCH DCA(y)" << std::endl;
  hn = GetTHnSparse(fAnalysisResults, (taskPath + "/DCA/MCH/DCA_y_vs_sign_vs_quadrant_vs_mom_corr").c_str());
  if (hn) {
    // axes assignments:
    // 0: momentum
    // 1: quadrant
    // 2: sign
    // 3: DCA
    int pAxisId = 0;
    int quadrantAxisId = 1;
    int signAxisId = 2;
    int dcaAxisId = 3;

    // p >= xx GeV/c
    int momBin = hn->GetAxis(pAxisId)->FindBin(mchMomMin + epsilon);
    int momBinMax = hn->GetAxis(pAxisId)->GetNbins();
    hn->GetAxis(pAxisId)->SetRange(momBin, momBinMax);
    for (int k = 0; k < 2; k++) {
      hn->GetAxis(signAxisId)->SetRange(k + 1, k + 1);
      for (int q = 0; q < quadrants.size(); q++) {
        hn->GetAxis(quadrantAxisId)->SetRange(q + 1, q + 1);
        if (q == 0) c.cd(2);
        if (q == 1) c.cd(1);
        if (q == 2) c.cd(3);
        if (q == 3) c.cd(4);

        auto* proj = hn->Projection(dcaAxisId);
        proj->SetTitle(std::format("MCH DCA(y), Q{} {}", q, (k == 0 ? "+" : "-")).c_str());
        DCAy[k][q] = PlotDCAMCH(proj);
        delete proj;
      }
      c.SaveAs(pdfFileName.c_str());
    }
  } else {
    hn = GetTHnSparse(fAnalysisResults, (taskPath + "/DCA/MCH/DCA_y_vs_sign_vs_quadrant_vs_vz").c_str());
    for (int k = 0; k < 2; k++) {
      hn->GetAxis(2)->SetRange(k + 1, k + 1);
      for (int q = 0; q < quadrants.size(); q++) {
        hn->GetAxis(1)->SetRange(q + 1, q + 1);
        if (q == 0) c.cd(2);
        if (q == 1) c.cd(1);
        if (q == 2) c.cd(3);
        if (q == 3) c.cd(4);

        auto* proj = hn->Projection(3);
        proj->SetTitle(std::format("MCH DCA(y), Q{} {}", q, (k == 0 ? "+" : "-")).c_str());
        DCAy[k][q] = PlotDCAMCH(proj);
      }
      c.SaveAs(pdfFileName.c_str());
    }
  }

  std::cout << "Plotting MCH DCA(x/y) vs. momentum" << std::endl;
  PlotDCAMCHvsMomentum(GetTHnSparse(fAnalysisResults, (taskPath + "/DCA/MCH/DCA_x_vs_sign_vs_quadrant_vs_mom_corr").c_str()),
                       GetTHnSparse(fAnalysisResults, (taskPath + "/DCA/MCH/DCA_y_vs_sign_vs_quadrant_vs_mom_corr").c_str()),
                       c);
*/
/*
  //c.Clear();
  //c.SaveAs((pdfFileName + ")").c_str());
  //return;

  // axes assignments (old version)
  // 0: DE index
  // 1: quadrant
  // 2: sign
  // 3: residual

  // axes assignments (new version)
  // 0: residual
  // 1: DE index
  // 2: quadrant
  // 3: sign
  // 4: momentum
  hn = GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/dx_vs_de_corr").c_str());
  if (hn) {
    int residualAxisId = 0;
    int deAxisId = 1;
    int quadrantAxisId = 2;
    int signAxisId = 3;
    int pAxisId = 4;

    if (hn->GetNdimensions() == 4) {
      residualAxisId = 3;
      deAxisId = 0;
      quadrantAxisId = 1;
      signAxisId = 2;
      pAxisId = -1;
    }

    hn->GetAxis(0)->SetRange(0, -1);
    hn->GetAxis(1)->SetRange(0, -1);
    hn->GetAxis(2)->SetRange(0, -1);
    hn->GetAxis(3)->SetRange(0, -1);
    if (pAxisId >= 0) {
      // p >= xx GeV/c
      int momBin = hn->GetAxis(pAxisId)->FindBin(mchMomMin + epsilon);
      int momBinMax = hn->GetAxis(pAxisId)->GetNbins();
      hn->GetAxis(pAxisId)->SetRange(momBin, momBinMax);
    }

    // loop on chamber ID
    for (int i = 0; i < 10; i++) {
      int deIdMin = getChamberOffset(i);
      int deIdMax = deIdMin + getNumDEinChamber(i) - 1;
      hn->GetAxis(deAxisId)->SetRange(deIdMin + 1, deIdMax + 1);
      // loop on quadrant
      for (int j = 0; j < hn->GetAxis(quadrantAxisId)->GetNbins(); j++) {
        hn->GetAxis(quadrantAxisId)->SetRange(j + 1, j + 1);
        // loop on charge
        for (int k = 0; k < hn->GetAxis(signAxisId)->GetNbins(); k++) {
          hn->GetAxis(signAxisId)->SetRange(k + 1, k + 1);

          auto* proj = hn->Projection(residualAxisId);
          proj->SetName(std::format("dx_vs_chamber_{}_{}_{}", j, i+1, (k == 0 ? "positive" : "negative")).c_str());
          proj->SetTitle(std::format("#Deltax, Q{} CH{} {}", j, i+1, (k == 0 ? "positive" : "negative")).c_str());
          meanDx[k][j][i] = PlotDXY(proj, c2, "residuals_CH.pdf");
        }
      }
    }
  }

  hn = GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/dy_vs_de_corr").c_str());
  if (hn) {
    int residualAxisId = 0;
    int deAxisId = 1;
    int quadrantAxisId = 2;
    int signAxisId = 3;
    int pAxisId = 4;

    if (hn->GetNdimensions() == 4) {
      residualAxisId = 3;
      deAxisId = 0;
      quadrantAxisId = 1;
      signAxisId = 2;
      pAxisId = -1;
    }

    hn->GetAxis(0)->SetRange(0, -1);
    hn->GetAxis(1)->SetRange(0, -1);
    hn->GetAxis(2)->SetRange(0, -1);
    hn->GetAxis(3)->SetRange(0, -1);
    if (pAxisId >= 0) {
      // p >= xx GeV/c
      int momBin = hn->GetAxis(pAxisId)->FindBin(mchMomMin + epsilon);
      int momBinMax = hn->GetAxis(pAxisId)->GetNbins();
      hn->GetAxis(pAxisId)->SetRange(momBin, momBinMax);
    }

    // loop on chamber ID
    for (int i = 0; i < 10; i++) {
      int deIdMin = getChamberOffset(i);
      int deIdMax = deIdMin + getNumDEinChamber(i) - 1;
      hn->GetAxis(deAxisId)->SetRange(deIdMin + 1, deIdMax + 1);
      // loop on quadrant
      for (int j = 0; j < hn->GetAxis(quadrantAxisId)->GetNbins(); j++) {
        hn->GetAxis(quadrantAxisId)->SetRange(j + 1, j + 1);
        // loop on charge
        for (int k = 0; k < hn->GetAxis(signAxisId)->GetNbins(); k++) {
          hn->GetAxis(signAxisId)->SetRange(k + 1, k + 1);

          auto* proj = hn->Projection(residualAxisId);
          proj->SetName(std::format("dy_vs_chamber_{}_{}_{}", j, i+1, (k == 0 ? "positive" : "negative")).c_str());
          proj->SetTitle(std::format("#Deltay, Q{} CH{} {}", j, i+1, (k == 0 ? "positive" : "negative")).c_str());
          meanDy[k][j][i] = PlotDXY(proj, c2, "residuals_CH.pdf");
        }
      }
    }
  }

  c2.Clear();
  c2.SaveAs("residuals_CH.pdf)");



  double xv[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  double chamberZ[10] = {526.16, 545.24, 676.4, 695.4, 967.5,
                                998.5, 1276.5, 1307.5, 1406.6, 1437.6};

  PlotZTrend(10, chamberZ, meanDx[0], DCAx[0], "#Delta(x) vs. chamber z (positive);chamber z (cm); #Delta(x) (cm)", -5.0, 5.0, c);
  PlotZTrend(10, chamberZ, meanDx[0], DCAx[0], "#Delta(x) vs. chamber z (positive);chamber z (cm); #Delta(x) (cm)", -1.0, 1.0, c);
  PlotZTrend(10, chamberZ, meanDx[0], DCAx[0], "#Delta(x) vs. chamber z (positive);chamber z (cm); #Delta(x) (cm)", -0.5, 0.5, c);
  PlotZTrend(10, chamberZ, meanDx[1], DCAx[1], "#Delta(x) vs. chamber z (negative);chamber z (cm); #Delta(x) (cm)", -5.0, 5.0, c);
  PlotZTrend(10, chamberZ, meanDx[1], DCAx[1], "#Delta(x) vs. chamber z (negative);chamber z (cm); #Delta(x) (cm)", -1.0, 1.0, c);
  PlotZTrend(10, chamberZ, meanDx[1], DCAx[1], "#Delta(x) vs. chamber z (negative);chamber z (cm); #Delta(x) (cm)", -0.5, 0.5, c);

  PlotZTrend(10, chamberZ, meanDy[0], DCAy[0], "#Delta(y) vs. chamber z (positive);chamber z (cm); #Delta(y) (cm)", -5.0, 5.0, c);
  PlotZTrend(10, chamberZ, meanDy[0], DCAy[0], "#Delta(y) vs. chamber z (positive);chamber z (cm); #Delta(y) (cm)", -1.0, 1.0, c);
  PlotZTrend(10, chamberZ, meanDy[0], DCAy[0], "#Delta(y) vs. chamber z (positive);chamber z (cm); #Delta(y) (cm)", -0.5, 0.5, c);
  PlotZTrend(10, chamberZ, meanDy[1], DCAy[1], "#Delta(y) vs. chamber z (negative);chamber z (cm); #Delta(y) (cm)", -5.0, 5.0, c);
  PlotZTrend(10, chamberZ, meanDy[1], DCAy[1], "#Delta(y) vs. chamber z (negative);chamber z (cm); #Delta(y) (cm)", -1.0, 1.0, c);
  PlotZTrend(10, chamberZ, meanDy[1], DCAy[1], "#Delta(y) vs. chamber z (negative);chamber z (cm); #Delta(y) (cm)", -0.5, 0.5, c);

  PlotChamberResidualvsMomentum(GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/dx_vs_de_corr").c_str()),
                                GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/dy_vs_de_corr").c_str()),
                                c);
*/

/*
  //-------------------
  // DE residuals
  //-------------------

  std::vector<std::tuple<std::string, std::string, std::string, std::vector<int>>> deGroups {
    {"DE100", "CH1", "Q0", {100}},
    {"DE101", "CH1", "Q1", {101}},
    {"DE102", "CH1", "Q2", {102}},
    {"DE103", "CH1", "Q3", {103}},
    {"DE200", "CH2", "Q0", {200}},
    {"DE201", "CH2", "Q1", {201}},
    {"DE202", "CH2", "Q2", {202}},
    {"DE203", "CH2", "Q3", {203}},
    {"DE300", "CH3", "Q0", {300}},
    {"DE301", "CH3", "Q1", {301}},
    {"DE302", "CH3", "Q2", {302}},
    {"DE303", "CH3", "Q3", {303}},
    {"DE400", "CH4", "Q0", {400}},
    {"DE401", "CH4", "Q1", {401}},
    {"DE402", "CH4", "Q2", {402}},
    {"DE403", "CH4", "Q3", {403}},
    {"CH5L", "CH5", "L", {505, 506, 507, 508, 509, 510, 511, 512, 513}},
    {"CH5R", "CH5", "R", {500, 501, 502, 503, 504, 514, 515, 516, 517}},
    {"CH6L", "CH6", "L", {605, 606, 607, 608, 609, 610, 611, 612, 613}},
    {"CH6R", "CH6", "R", {600, 601, 602, 603, 604, 614, 615, 616, 617}},
    {"CH7L", "CH7", "L", {707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719}},
    {"CH7R", "CH7", "R", {700, 701, 702, 703, 704, 706, 707, 720, 721, 722, 723, 724, 725}},
    {"CH8L", "CH8", "L", {807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819}},
    {"CH8R", "CH8", "R", {800, 801, 802, 803, 804, 806, 807, 820, 821, 822, 823, 824, 825}},
    {"CH9L", "CH9", "L", {907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919}},
    {"CH9R", "CH9", "R", {900, 901, 902, 903, 904, 906, 907, 920, 921, 922, 923, 924, 925}},
    {"CH10L", "CH10", "L", {1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019}},
    {"CH10R", "CH10", "R", {1000, 1001, 1002, 1003, 1004, 1006, 1007, 1020, 1021, 1022, 1023, 1024, 1025}}
  };

  std::array<std::vector<std::array<TH1*, 2>>, 2> hDEResiduals;

  std::ofstream corrections_x("shifts-x.txt");
  std::ofstream corrections_y("shifts-y.txt");

  // axis assignments (old version with 4 axes):
  // 0: DE index
  // 1: quadrant
  // 2: charge sign
  // 3: residual

  // axis assignments (new version with 5 axes):
  // 0: residual
  // 1: DE index
  // 2: quadrant
  // 3: charge sign
  // 4: momentum

  // -----------------
  // X direction

  THnSparse* hDE[2] = {
      GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/dx_vs_de_corr").c_str()),
      GetTHnSparse(fAnalysisResults, (taskPath + "/residuals/dy_vs_de_corr").c_str())
  };
  if (hDE[0] && hDE[1]) {

    std::ofstream correctionsJson("corrections.json");
    correctionsJson << "{" << std::endl;

    int residualAxisId = 0;
    int deAxisId = 1;
    int quadrantAxisId = 2;
    int signAxisId = 3;
    int pAxisId = 4;

    if (hDE[0]->GetNdimensions() == 4) {
      residualAxisId = 3;
      deAxisId = 0;
      quadrantAxisId = 1;
      signAxisId = 2;
      pAxisId = -1;
    }

    for (int xy = 0; xy < 2; xy++) {
      hDE[xy]->GetAxis(0)->SetRange(0, -1);
      hDE[xy]->GetAxis(1)->SetRange(0, -1);
      hDE[xy]->GetAxis(2)->SetRange(0, -1);
      hDE[xy]->GetAxis(3)->SetRange(0, -1);
      if (pAxisId >= 0) {
        // p >= xx GeV/c
        int momBin = hDE[xy]->GetAxis(pAxisId)->FindBin(mchMomMin + epsilon);
        int momBinMax = hDE[xy]->GetAxis(pAxisId)->GetNbins();
        hDE[xy]->GetAxis(pAxisId)->SetRange(momBin, momBinMax);
      }

      TH1* histDxyVsDE[2] = {nullptr, nullptr};
      for (int k = 0; k < hDE[xy]->GetAxis(signAxisId)->GetNbins(); k++) {
        //hDE[xy]->GetAxis(2)->SetRange(k + 1, k + 1);
        //hDE[xy]->GetAxis(0)->SetRange(1, hDE[xy]->GetAxis(0)->GetNbins());
        histDxyVsDE[k] = hDE[xy]->Projection(deAxisId);
        histDxyVsDE[k]->Reset();
        histDxyVsDE[k]->SetName(std::format("d{}_vs_de_{}", (xy == 0 ? "x" : "y"), (k == 0 ? "positive" : "negative")).c_str());
        histDxyVsDE[k]->SetTitle(std::format("#Delta{} vs. DE", (xy == 0 ? "x" : "y")).c_str());
      }
      for (int i = 0; i < hDE[xy]->GetAxis(deAxisId)->GetNbins(); i++) {
        hDE[xy]->GetAxis(deAxisId)->SetRange(i + 1, i + 1);
        for (int k = 0; k < hDE[xy]->GetAxis(signAxisId)->GetNbins(); k++) {
          hDE[xy]->GetAxis(signAxisId)->SetRange(k + 1, k + 1);

          auto* proj = hDE[xy]->Projection(residualAxisId);
          proj->SetName(std::format("d{}_vs_de_{}_{}", (xy == 0 ? "x" : "y"), i+1, (k == 0 ? "positive" : "negative")).c_str());
          proj->SetTitle(std::format("#Delta{}, DE{} {}", (xy == 0 ? "x" : "y"), getDEFromIndex(i), (k == 0 ? "positive" : "negative")).c_str());
          auto mean = PlotDXY(proj, c3, "residuals_DE.pdf");
          histDxyVsDE[k]->SetBinContent(i + 1, mean.first);
          histDxyVsDE[k]->SetBinError(i + 1, mean.second);
        }
      }
      c.Clear();
      c.cd();
      histDxyVsDE[0]->SetLineColor(kBlue);
      histDxyVsDE[0]->SetMinimum(-1.0);
      histDxyVsDE[0]->SetMaximum(1.0);
      histDxyVsDE[0]->Draw("E");
      histDxyVsDE[1]->SetLineColor(kRed);
      histDxyVsDE[1]->Draw("E same");
      c.SaveAs(pdfFileName.c_str());

      // set DE names in X axis
      for (int xbin = 1; xbin <= histDxyVsDE[0]->GetXaxis()->GetNbins(); xbin++) {
        int deId = getDEFromIndex(xbin - 1);
        histDxyVsDE[0]->GetXaxis()->SetBinLabel(xbin, TString::Format("DE%d", deId));
      }

      std::string histTitle = histDxyVsDE[0]->GetTitle();
      for (int chIndex = 0; chIndex < 10; chIndex++) {
        int deMin = getChamberOffset(chIndex);
        int deMax = deMin + getNumDEinChamber(chIndex);
        histDxyVsDE[0]->GetXaxis()->SetRangeUser(deMin, deMax);
        histDxyVsDE[0]->SetTitle((histTitle + "(CH" + std::to_string(chIndex+1) + ")").c_str());
        histDxyVsDE[0]->SetMinimum(-5.0);
        histDxyVsDE[0]->SetMaximum(5.0);
        c.SaveAs(pdfFileName.c_str());

        histDxyVsDE[0]->SetMinimum(-1.0);
        histDxyVsDE[0]->SetMaximum(1.0);
        c.SaveAs(pdfFileName.c_str());

        histDxyVsDE[0]->SetMinimum(-0.5);
        histDxyVsDE[0]->SetMaximum(0.5);
        c.SaveAs(pdfFileName.c_str());

        if (chIndex < 4) {
          histDxyVsDE[0]->SetMinimum(-0.2);
          histDxyVsDE[0]->SetMaximum(0.2);
          c.SaveAs(pdfFileName.c_str());
        }
      }
    }

    //   std::array<std::vector<std::array<TH1*, 2>>, 2> hDEResiduals;
    //hn->GetAxis(1)->SetRange(1, hn->GetAxis(1)->GetNbins());
    for (const auto& group : deGroups) {
      for (int xy = 0; xy < 2; xy++) {
        hDEResiduals[xy].emplace_back(std::array<TH1*, 2>{nullptr, nullptr});
        // loop over charge sign
        for (int charge = 0; charge < 2; charge++) {
          hDE[xy]->GetAxis(signAxisId)->SetRange(charge + 1, charge + 1);
          // loop over DE ids
          for (auto deId : std::get<3>(group)) {
            auto deIndex = getDEindex(deId);
            hDE[xy]->GetAxis(deAxisId)->SetRange(deIndex + 1, deIndex + 1);
            TH1* proj = hDE[xy]->Projection(residualAxisId);
            if (hDEResiduals[xy].back()[charge]) {
              hDEResiduals[xy].back()[charge]->Add(proj);
              delete proj;
            } else {
              proj->SetTitle(std::format("{} - #Deltax, {}", std::get<0>(group), (charge == 0 ? "positive" : "negative")).c_str());
              hDEResiduals[xy].back()[charge] = proj;
            }
          }
        }
      }

      // loop over charge sign
      double averageCorrectionX = 0;
      for (int charge = 0; charge < 2; charge++) {
        auto mean = PlotDXY(hDEResiduals[0].back()[charge], c5, "residuals_groups.pdf");
        corrections_x << std::format("{} ({}): dx = {:0.3f} +/- {:0.3f}",
            std::get<0>(group), (charge == 0 ? "positive" : "negative"), mean.first, mean.second) << std::endl;

        averageCorrectionX += -mean.first;
      }
      averageCorrectionX /= 2;

      // loop over charge sign
      double averageCorrectionY = 0;
      for (int charge = 0; charge < 2; charge++) {
        auto mean = PlotDXY(hDEResiduals[1].back()[charge], c5, "residuals_groups.pdf");
        corrections_y << std::format("{} ({}): dx = {:0.3f} +/- {:0.3f}",
            std::get<0>(group), (charge == 0 ? "positive" : "negative"), mean.first, mean.second) << std::endl;

        averageCorrectionY += -mean.first;
      }
      averageCorrectionY /= 2;

      // loop over DE ids
      for (auto deId : std::get<3>(group)) {
        correctionsJson << "  \"" << deId << "\": {" << std::endl;
        correctionsJson << "    \"x\": " << averageCorrectionX << ","
            << " \"y\": " << averageCorrectionY << ","
            << " \"z\": 0, \"yaw\": 0, \"pitch\": 0, \"roll\": 0" << std::endl;
        correctionsJson << "  }," << std::endl;
      }
    }
    correctionsJson << "  \"last\": {}" << std::endl;
    correctionsJson << "}" << std::endl;
  }

  c3.Clear();
  c3.SaveAs("residuals_DE.pdf)");

  c5.Clear();
  c5.SaveAs("residuals_groups.pdf)");
*/

  c.Clear();
  c.SaveAs((pdfFileName + ")").c_str());
}
