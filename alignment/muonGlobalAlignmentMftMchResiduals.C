#include <MFTTracking/Constants.h>

TFile* fAnalysisResults;
std::string pdfFileName;

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

std::pair<double, double> PlotDCAMCH(TH1* histogram)
{
  //std::string fullHistName = std::string("qa-muon/alignment/") + histName;
  //TH1* histogram = GetTH1(fAnalysisResults, fullHistName);
  //std::cout << fullHistName << " -> " << histogram << std::endl;
  if (!histogram)
    return {};
  histogram->Rebin(4);

  histogram->Draw("E");

  int valuePeak = histogram->GetMaximum();
  int binPeak = histogram->GetMaximumBin();
  double xPeak = histogram->GetXaxis()->GetBinCenter(binPeak);

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
/*
  fcb.FixParameter(1, xPeak);
  fcb.SetParLimits(2, 0.5, 10);
  histogram->Fit("fcb", "BRN");
  fcb.ReleaseParameter(1);
  histogram->Fit("fcb", "BRN");
  histogram->Fit("fcb", "BR");
*/
/**/
  TF1 fgaus2("fgaus2", "gaus(0)");
  fgaus2.SetNpx(1000);
  fgaus2.SetLineColor(kBlack);
  fgaus2.SetParameter(0, valuePeak / 10);
  //fgaus2.SetParLimits(0, 0, valuePeak * 10);
  fgaus2.SetParameter(1, xPeak);
  fgaus2.FixParameter(1, xPeak);
  fgaus2.SetParameter(2, 1);
  fgaus2.SetParLimits(2, 0, 10);
  histogram->Fit("fgaus2", "BQN");
  fgaus2.ReleaseParameter(1);
  fgaus2.SetParLimits(1, -10, 10);
  histogram->Fit("fgaus2", "LBQN");

  TF1 fgaus("fgaus", "gaus(0)+gaus(3)", fgaus2.GetParameter(1) - 1.0 * fgaus2.GetParameter(2), fgaus2.GetParameter(1) + 1.0 * fgaus2.GetParameter(2));
  fgaus.SetNpx(1000);
  fgaus.SetLineColor(kBlack);
  //fgaus.SetParLimits(0, 0, 100000000000.0);
  fgaus.SetParameter(0, fgaus2.GetParameter(0));
  fgaus.SetParameter(1, fgaus2.GetParameter(1));
  fgaus.SetParameter(2, fgaus2.GetParameter(2));
  fgaus.SetParLimits(2, 0, 100);
  fgaus.SetParameter(3, fgaus2.GetParameter(0)/10);
  fgaus.FixParameter(3, 0);
  fgaus.SetParameter(4, fgaus2.GetParameter(1));
  fgaus.SetParameter(5, fgaus2.GetParameter(2)*5);
  fgaus.SetParLimits(5, fgaus2.GetParameter(2)*2, 100);
  histogram->Fit("fgaus", "RQB");
/**/
  histogram->Draw("E");

  //return {fcb.GetParameter(1), fcb.GetParError(1)};
  return {fgaus.GetParameter(1), fgaus.GetParError(1)};
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

std::pair<double, double> PlotDXY(TH1* proj, TCanvas& c, std::string pdfName = "", bool printFits = false)
{
  c.Clear();
  proj->SetLineColor(kRed);
  proj->Draw("E");

  TF1 fgaus("fgaus", "gausn(0)+pol2(3)");
  int valuePeak = proj->GetMaximum();
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
  proj->Fit("fgaus", "BQ");
/*
  TF1 fcb("fcb","[0]*ROOT::Math::crystalball_function(x, [1], [2], [3], [4]) + [5]");
  fcb.SetParameters(100, 0.6, -2.13903e+06, 1, xPeak);
  fcb.SetParLimits(3, 0, 100);
  fcb.SetNpx(1000);
  fcb.SetLineColor(kBlack);
  //proj2->Fit("fcb", "BNQ");
  //proj2->Fit("fcb", "BQ");

  proj2->GetXaxis()->SetRangeUser(-15.0, 15.0);
  proj2->Draw("E");
 */
  if (!pdfName.empty()) {
    c.SaveAs(pdfName.c_str());
  }

  //PlotDXYProjection(fullHistName.c_str(), fullHistNameME.c_str(), histogram2, histogram2ME, scaleME, -5.0, 5.0, 1, c, true, printFits);

  return {fgaus.GetParameter(1), fgaus.GetParError(1)};
  //return {fcb.GetParameter(4), fgaus.GetParError(4)};
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

void muonGlobalAlignmentMftMchResiduals(const char* _rootFileName = "AnalysisResults.root", const char* _pdfFileName = "mftMchResiduals.pdf")
{
  //fAnalysisResults = new TFile("AnalysisResults.root");
  //fAnalysisResults = new TFile("AnalysisResults/AnalysisResultsFull.root");
  fAnalysisResults = new TFile(_rootFileName);
  pdfFileName = _pdfFileName;

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

  TCanvas c2("c2", "c2", 1200, 800);
  c2.SaveAs("residuals_CH.pdf(");

  TCanvas c3("c3", "c3", 1200, 800);
  c3.SaveAs("residuals_DE.pdf(");

  TCanvas c4("c4", "c4", 1200, 800);
  c4.SaveAs("residuals_tracks.pdf(");

  TCanvas c5("c5", "c5", 1200, 800);
  c4.SaveAs("residuals_groups.pdf(");

  TH1* h1 = GetTH1(fAnalysisResults, "muon-global-alignment/DCA/vertex_z");
  if (h1) {
    h1->Draw();
    c.SaveAs(pdfFileName.c_str());
  }

  c.Clear();
  c.Divide(2, 2);

  auto* hn = GetTHnSparse(fAnalysisResults, "muon-global-alignment/DCA/MCH/DCA_x_vs_sign_vs_quadrant_vs_vz");
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
    }
    c.SaveAs(pdfFileName.c_str());
  }

  hn = GetTHnSparse(fAnalysisResults, "muon-global-alignment/DCA/MCH/DCA_y_vs_sign_vs_quadrant_vs_vz");
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

  //c.Clear();
  //c.SaveAs((pdfFileName + ")").c_str());
  //return;

  hn = GetTHnSparse(fAnalysisResults, "muon-global-alignment/residuals/track_dx");
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

  hn = GetTHnSparse(fAnalysisResults, "muon-global-alignment/residuals/track_dy");
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

  hn = GetTHnSparse(fAnalysisResults, "muon-global-alignment/residuals/dx_vs_chamber");
  for (int i = 0; i < hn->GetAxis(0)->GetNbins(); i++) {
    hn->GetAxis(0)->SetRange(i + 1, i + 1);
    for (int j = 0; j < hn->GetAxis(1)->GetNbins(); j++) {
      hn->GetAxis(1)->SetRange(j + 1, j + 1);
      for (int k = 0; k < hn->GetAxis(2)->GetNbins(); k++) {
        hn->GetAxis(2)->SetRange(k + 1, k + 1);

        auto* proj = hn->Projection(3);
        proj->SetName(std::format("dx_vs_chamber_{}_{}_{}", j, i+1, (k == 0 ? "positive" : "negative")).c_str());
        proj->SetTitle(std::format("#Deltax, Q{} CH{} {}", j, i+1, (k == 0 ? "positive" : "negative")).c_str());
        meanDx[k][j][i] = PlotDXY(proj, c2, "residuals_CH.pdf");
      }
    }
  }

  hn = GetTHnSparse(fAnalysisResults, "muon-global-alignment/residuals/dy_vs_chamber");
  for (int i = 0; i < hn->GetAxis(0)->GetNbins(); i++) {
    hn->GetAxis(0)->SetRange(i + 1, i + 1);
    for (int j = 0; j < hn->GetAxis(1)->GetNbins(); j++) {
      hn->GetAxis(1)->SetRange(j + 1, j + 1);
      for (int k = 0; k < hn->GetAxis(2)->GetNbins(); k++) {
        hn->GetAxis(2)->SetRange(k + 1, k + 1);

        auto* proj = hn->Projection(3);
        proj->SetName(std::format("dy_vs_chamber_{}_{}_{}", j, i+1, (k == 0 ? "positive" : "negative")).c_str());
        proj->SetTitle(std::format("#Deltay, Q{} CH{} {}", j, i+1, (k == 0 ? "positive" : "negative")).c_str());
        meanDy[k][j][i] = PlotDXY(proj, c2, "residuals_CH.pdf");
      }
    }
  }
/**/

  double xv[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  double defaultChamberZ[10] = {526.16, 545.24, 676.4, 695.4, 967.5,
                                998.5, 1276.5, 1307.5, 1406.6, 1437.6};

  PlotZTrend(10, defaultChamberZ, meanDx[0], DCAx[0], "#Delta(x) vs. chamber z (positive);chamber z (cm); #Delta(x) (cm)", -5.0, 5.0, c);
  PlotZTrend(10, defaultChamberZ, meanDx[0], DCAx[0], "#Delta(x) vs. chamber z (positive);chamber z (cm); #Delta(x) (cm)", -1.0, 1.0, c);
  PlotZTrend(10, defaultChamberZ, meanDx[0], DCAx[0], "#Delta(x) vs. chamber z (positive);chamber z (cm); #Delta(x) (cm)", -0.5, 0.5, c);
  PlotZTrend(10, defaultChamberZ, meanDx[1], DCAx[1], "#Delta(x) vs. chamber z (negative);chamber z (cm); #Delta(x) (cm)", -5.0, 5.0, c);
  PlotZTrend(10, defaultChamberZ, meanDx[1], DCAx[1], "#Delta(x) vs. chamber z (negative);chamber z (cm); #Delta(x) (cm)", -1.0, 1.0, c);
  PlotZTrend(10, defaultChamberZ, meanDx[1], DCAx[1], "#Delta(x) vs. chamber z (negative);chamber z (cm); #Delta(x) (cm)", -0.5, 0.5, c);

  PlotZTrend(10, defaultChamberZ, meanDy[0], DCAy[0], "#Delta(y) vs. chamber z (positive);chamber z (cm); #Delta(y) (cm)", -5.0, 5.0, c);
  PlotZTrend(10, defaultChamberZ, meanDy[0], DCAy[0], "#Delta(y) vs. chamber z (positive);chamber z (cm); #Delta(y) (cm)", -1.0, 1.0, c);
  PlotZTrend(10, defaultChamberZ, meanDy[0], DCAy[0], "#Delta(y) vs. chamber z (positive);chamber z (cm); #Delta(y) (cm)", -0.5, 0.5, c);
  PlotZTrend(10, defaultChamberZ, meanDy[1], DCAy[1], "#Delta(y) vs. chamber z (negative);chamber z (cm); #Delta(y) (cm)", -5.0, 5.0, c);
  PlotZTrend(10, defaultChamberZ, meanDy[1], DCAy[1], "#Delta(y) vs. chamber z (negative);chamber z (cm); #Delta(y) (cm)", -1.0, 1.0, c);
  PlotZTrend(10, defaultChamberZ, meanDy[1], DCAy[1], "#Delta(y) vs. chamber z (negative);chamber z (cm); #Delta(y) (cm)", -0.5, 0.5, c);


  //-------------------
  // DE residuals
  //-------------------

  std::vector<std::pair<std::string, std::vector<int>>> deGroups {
    {"DE100", {100}},
    {"DE101", {101}},
    {"DE102", {102}},
    {"DE103", {103}},
    {"DE200", {200}},
    {"DE201", {201}},
    {"DE202", {202}},
    {"DE203", {203}},
    {"DE300", {300}},
    {"DE301", {301}},
    {"DE302", {302}},
    {"DE303", {303}},
    {"DE400", {400}},
    {"DE401", {401}},
    {"DE402", {402}},
    {"DE403", {403}},
    {"CH5L", {505, 506, 507, 508, 509, 510, 511, 512, 513}},
    {"CH5R", {500, 501, 502, 503, 504, 514, 515, 516, 517}},
    {"CH6L", {605, 606, 607, 608, 609, 610, 611, 612, 613}},
    {"CH6R", {600, 601, 602, 603, 604, 614, 615, 616, 617}},
    {"CH7L", {707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719}},
    {"CH7R", {700, 701, 702, 703, 704, 706, 707, 720, 721, 722, 723, 724, 725}},
    {"CH8L", {807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819}},
    {"CH8R", {800, 801, 802, 803, 804, 806, 807, 820, 821, 822, 823, 824, 825}},
    {"CH9L", {907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919}},
    {"CH9R", {900, 901, 902, 903, 904, 906, 907, 920, 921, 922, 923, 924, 925}},
    {"CH10L", {1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019}},
    {"CH10R", {1000, 1001, 1002, 1003, 1004, 1006, 1007, 1020, 1021, 1022, 1023, 1024, 1025}}
  };

  std::array<std::vector<std::array<TH1*, 2>>, 2> hDEResiduals;

  std::ofstream corrections_x("corrections-x.txt");

  // axis assignments:
  // 0: DE index
  // 1: quadrant
  // 2: charge sign
  // 3: residual

  // -----------------
  // X direction

  hn = GetTHnSparse(fAnalysisResults, "muon-global-alignment/residuals/dx_vs_de");
  TH1* histDxVsDE[2] = {nullptr, nullptr};
  for (int k = 0; k < hn->GetAxis(2)->GetNbins(); k++) {
    //hn->GetAxis(2)->SetRange(k + 1, k + 1);
    //hn->GetAxis(0)->SetRange(1, hn->GetAxis(0)->GetNbins());
    histDxVsDE[k] = hn->Projection(0);
    histDxVsDE[k]->Reset();
    histDxVsDE[k]->SetName(std::format("dx_vs_de_{}", (k == 0 ? "positive" : "negative")).c_str());
    histDxVsDE[k]->SetTitle("#Deltax vs. DE");
  }
  for (int i = 0; i < hn->GetAxis(0)->GetNbins(); i++) {
    hn->GetAxis(0)->SetRange(i + 1, i + 1);
    for (int k = 0; k < hn->GetAxis(2)->GetNbins(); k++) {
      hn->GetAxis(2)->SetRange(k + 1, k + 1);

      auto* proj = hn->Projection(3);
      proj->SetName(std::format("dx_vs_de_{}_{}", i+1, (k == 0 ? "positive" : "negative")).c_str());
      proj->SetTitle(std::format("#Deltax, DE{} {}", getDEFromIndex(i), (k == 0 ? "positive" : "negative")).c_str());
      auto mean = PlotDXY(proj, c3, "residuals_DE.pdf");
      histDxVsDE[k]->SetBinContent(i + 1, mean.first);
      histDxVsDE[k]->SetBinError(i + 1, mean.second);
    }
  }
  c.Clear();
  c.cd();
  histDxVsDE[0]->SetLineColor(kBlue);
  histDxVsDE[0]->SetMinimum(-1.0);
  histDxVsDE[0]->SetMaximum(1.0);
  histDxVsDE[0]->Draw("E");
  histDxVsDE[1]->SetLineColor(kRed);
  histDxVsDE[1]->Draw("E same");
  c.SaveAs(pdfFileName.c_str());

  // set DE names in X axis
  for (int xbin = 1; xbin <= histDxVsDE[0]->GetXaxis()->GetNbins(); xbin++) {
    int deId = getDEFromIndex(xbin - 1);
    histDxVsDE[0]->GetXaxis()->SetBinLabel(xbin, TString::Format("DE%d", deId));
  }

  std::string histTitle = histDxVsDE[0]->GetTitle();
  for (int chIndex = 0; chIndex < 10; chIndex++) {
    int deMin = getChamberOffset(chIndex);
    int deMax = deMin + getNumDEinChamber(chIndex);
    histDxVsDE[0]->GetXaxis()->SetRangeUser(deMin, deMax);
    histDxVsDE[0]->SetTitle((histTitle + "(CH" + std::to_string(chIndex+1) + ")").c_str());
    histDxVsDE[0]->SetMinimum(-1.0);
    histDxVsDE[0]->SetMaximum(1.0);
    c.SaveAs(pdfFileName.c_str());

    histDxVsDE[0]->SetMinimum(-0.5);
    histDxVsDE[0]->SetMaximum(0.5);
    c.SaveAs(pdfFileName.c_str());

    if (chIndex < 4) {
      histDxVsDE[0]->SetMinimum(-0.1);
      histDxVsDE[0]->SetMaximum(0.1);
      c.SaveAs(pdfFileName.c_str());
    }
  }

  //   std::array<std::vector<std::array<TH1*, 2>>, 2> hDEResiduals;
  hn->GetAxis(1)->SetRange(1, hn->GetAxis(1)->GetNbins());
  for (const auto& group : deGroups) {
    hDEResiduals[0].emplace_back(std::array<TH1*, 2>{nullptr, nullptr});
    // loop over charge sign
    for (int charge = 0; charge < 2; charge++) {
      hn->GetAxis(2)->SetRange(charge + 1, charge + 1);
      // loop over DE ids
      for (auto deId : group.second) {
        auto deIndex = getDEindex(deId);
        hn->GetAxis(0)->SetRange(deIndex + 1, deIndex + 1);
        TH1* proj = hn->Projection(3);
        if (hDEResiduals[0].back()[charge]) {
          hDEResiduals[0].back()[charge]->Add(proj);
          delete proj;
        } else {
          proj->SetTitle(std::format("{} - #Deltax, {}", group.first, (charge == 0 ? "positive" : "negative")).c_str());
          hDEResiduals[0].back()[charge] = proj;
        }
      }
      auto mean = PlotDXY(hDEResiduals[0].back()[charge], c5, "residuals_groups.pdf");
      corrections_x << std::format("{} ({}): dx = {:0.3f} +/- {:0.3f}",
          group.first, (charge == 0 ? "positive" : "negative"), mean.first, mean.second) << std::endl;
    }
  }

  // -----------------
  // Y direction

  std::ofstream corrections_y("corrections-y.txt");

  hn = GetTHnSparse(fAnalysisResults, "muon-global-alignment/residuals/dy_vs_de");
  TH1* histDyVsDE[2] = {nullptr, nullptr};
  for (int k = 0; k < hn->GetAxis(2)->GetNbins(); k++) {
    hn->GetAxis(2)->SetRange(k + 1, k + 1);
    hn->GetAxis(0)->SetRange(1, hn->GetAxis(0)->GetNbins());
    histDyVsDE[k] = hn->Projection(0);
    histDyVsDE[k]->Reset();
    histDyVsDE[k]->SetName(std::format("dy_vs_de_{}", (k == 0 ? "positive" : "negative")).c_str());
    histDyVsDE[k]->SetTitle("#Deltay vs. DE");
    for (int i = 0; i < hn->GetAxis(0)->GetNbins(); i++) {
      hn->GetAxis(0)->SetRange(i + 1, i + 1);

      auto* proj = hn->Projection(3);
      proj->SetName(std::format("dy_vs_de_{}_{}", i+1, (k == 0 ? "positive" : "negative")).c_str());
      proj->SetTitle(std::format("#Deltay, DE{} {}", i, (k == 0 ? "positive" : "negative")).c_str());
      auto mean = PlotDXY(proj, c3, "residuals_DE.pdf");
      histDyVsDE[k]->SetBinContent(i + 1, mean.first);
      histDyVsDE[k]->SetBinError(i + 1, mean.second);
    }
  }
  c.Clear();
  c.cd();
  histDyVsDE[0]->SetLineColor(kBlue);
  histDyVsDE[0]->SetMinimum(-1.0);
  histDyVsDE[0]->SetMaximum(1.0);
  histDyVsDE[0]->Draw("E");
  histDyVsDE[1]->SetLineColor(kRed);
  histDyVsDE[1]->Draw("E same");
  c.SaveAs(pdfFileName.c_str());

  // set DE names in X axis
  for (int xbin = 1; xbin <= histDyVsDE[0]->GetXaxis()->GetNbins(); xbin++) {
    int deId = getDEFromIndex(xbin - 1);
    histDyVsDE[0]->GetXaxis()->SetBinLabel(xbin, TString::Format("DE%d", deId));
  }

  histTitle = histDyVsDE[0]->GetTitle();
  for (int chIndex = 0; chIndex < 10; chIndex++) {
    int deMin = getChamberOffset(chIndex);
    int deMax = deMin + getNumDEinChamber(chIndex);
    histDyVsDE[0]->GetXaxis()->SetRangeUser(deMin, deMax);
    histDyVsDE[0]->SetTitle((histTitle + "(CH" + std::to_string(chIndex+1) + ")").c_str());
    histDyVsDE[0]->SetMinimum(-1.0);
    histDyVsDE[0]->SetMaximum(1.0);
    c.SaveAs(pdfFileName.c_str());

    histDyVsDE[0]->SetMinimum(-0.5);
    histDyVsDE[0]->SetMaximum(0.5);
    c.SaveAs(pdfFileName.c_str());

    if (chIndex < 4) {
      histDyVsDE[0]->SetMinimum(-0.1);
      histDyVsDE[0]->SetMaximum(0.1);
      c.SaveAs(pdfFileName.c_str());
    }
  }

  //   std::array<std::vector<std::array<TH1*, 2>>, 2> hDEResiduals;
  hn->GetAxis(1)->SetRange(1, hn->GetAxis(1)->GetNbins());
  for (const auto& group : deGroups) {
    hDEResiduals[0].emplace_back(std::array<TH1*, 2>{nullptr, nullptr});
    // loop over charge sign
    for (int charge = 0; charge < 2; charge++) {
      hn->GetAxis(2)->SetRange(charge + 1, charge + 1);
      // loop over DE ids
      for (auto deId : group.second) {
        auto deIndex = getDEindex(deId);
        hn->GetAxis(0)->SetRange(deIndex + 1, deIndex + 1);
        TH1* proj = hn->Projection(3);
        if (hDEResiduals[0].back()[charge]) {
          hDEResiduals[1].back()[charge]->Add(proj);
          delete proj;
        } else {
          proj->SetTitle(std::format("{} - #Deltay, {}", group.first, (charge == 0 ? "positive" : "negative")).c_str());
          hDEResiduals[1].back()[charge] = proj;
        }
      }
      auto mean = PlotDXY(hDEResiduals[0].back()[charge], c5, "residuals_groups.pdf");
      corrections_y << std::format("{} ({}): dy = {:0.3f} +/- {:0.3f}",
          group.first, (charge == 0 ? "positive" : "negative"), mean.first, mean.second) << std::endl;
    }
  }


  c.Clear();
  c.SaveAs((pdfFileName + ")").c_str());

  c3.Clear();
  c2.SaveAs("residuals_CH.pdf)");

  c3.Clear();
  c3.SaveAs("residuals_DE.pdf)");

  c5.Clear();
  c5.SaveAs("residuals_groups.pdf)");
}
