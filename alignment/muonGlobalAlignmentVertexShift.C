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
  return (THnSparse*)f->Get(histname);
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

    //std::cout << std::format("\"{}\" -> peak = {}", proj->GetName(), valuePeak) << std::endl;

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
      par[2]=0.02;
      par[3]=1;
      par[4]=1;
      par[5]=1;
      par[6]=1;
      fcb.SetParameters(&par[0]);

      fcb.FixParameter(1, xPeak);
      fcb.FixParameter(2, 0.02);
      fcb.SetParLimits(3, 0.01f, 2.f);
      fcb.SetParLimits(5, 0.01f, 2.f);
      fcb.SetParLimits(4, 0.5f, 10.f);
      fcb.SetParLimits(6, 0.5f, 10.f);
      proj->Fit("fcb", "BRNQ");
      fcb.SetParLimits(2, 0.01, 0.1);
      fcb.ReleaseParameter(1);
      const char* fitOpt = printFits ? "BRS" : "BRQS";
      TFitResultPtr fitResult = proj->Fit("fcb", fitOpt);
      if (fitResult->Status() != 0) {
        fitResult->Print();
      }

      if (printFits) {
        proj->Draw("E");
        c.SaveAs(pdfFileName.c_str());
      }

      if (fitResult->Status() == 0) {
        mean = fcb.GetParameter(1);
        meanErr = fcb.GetParError(1);
        sigma = fcb.GetParameter(2);
        sigmaErr = fcb.GetParError(2);
      } else {
        mean = 0;
        meanErr = 1.e6;
        sigma = 1.e6;
        sigmaErr = 1.e6;
      }

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
  TFitResultPtr fitResult = histogramMean->Fit("phiFit", "BQRS");
  if (fitResult->Status() != 0) {
    fitResult->Print();
  }
  //histogramMean->SetMinimum(meanMin);
  //histogramMean->SetMaximum(meanMax);
  histogramMean->SetMinimum(-0.1);
  histogramMean->SetMaximum(0.1);
  histogramMean->Draw("E");

  TLine* line1 = new TLine(histogramMean->GetXaxis()->GetXmin(), 0, histogramMean->GetXaxis()->GetXmax(), 0);
  line1->SetLineStyle(kDashed);
  line1->Draw();

  TLine* line2 = new TLine(0, yMin, 0, yMax);
  line2->SetLineStyle(kDashed);
  line2->Draw();

  //std::cout << std::format("Amplitude: {:0.4f}", phiFit.GetParameter(0)) << std::endl;
  //std::cout << std::format("Phase:     {:0.4f}", phiFit.GetParameter(1)) << std::endl;
  //std::cout << std::format("Offsets:   top={:0.4f} bottom={:0.4f}", phiFit.GetParameter(2), phiFit.GetParameter(3)) << std::endl;

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
    //if (bin != 16) continue;
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
    //std::cout << std::format("Bin #{} -> {}", bin, histogram3->GetXaxis()->GetBinCenter(bin)) << std::endl;
    xv.push_back(histogram3->GetXaxis()->GetBinCenter(bin));
    exv.push_back(0);
    yv.push_back(amplitude.first);
    eyv.push_back(amplitude.second);

  }
  TGraphErrors* gr = new TGraphErrors(xv.size(), xv.data(), yv.data(), exv.data(), eyv.data());
  gr->Draw("A*");
  gr->SetTitle("Modulation amplitude vs. z shift");
  gr->GetXaxis()->SetTitle("z shift (mm)");
  gr->GetYaxis()->SetTitle("modulation");

  TLine* line1 = new TLine(gr->GetHistogram()->GetXaxis()->GetXmin(), 0, gr->GetHistogram()->GetXaxis()->GetXmax(), 0);
  line1->SetLineStyle(kDashed);
  line1->Draw();

  TLine* line2 = new TLine(0, gr->GetHistogram()->GetMinimum(), 0, gr->GetHistogram()->GetMaximum());
  line2->SetLineStyle(kDashed);
  line2->Draw();

  TF1 linFit("linFit", "pol1");
  linFit.SetLineColor(kRed);
  linFit.SetLineStyle(kDotted);
  TFitResultPtr fitResult = gr->Fit("linFit", "BQS");
  if (fitResult->Status() != 0) {
    fitResult->Print();
  }

  std::cout << std::format("Optimal z shift: {} cm", -linFit.GetParameter(0) / linFit.GetParameter(1) / 10.f) << std::endl;
}

std::pair<double, double> PlotDCASlopeProjection(TH2* histogram2, float yMin, float yMax, int projRebin, TCanvas& c, bool printFits = false)
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

    //std::cout << std::format("\"{}\" -> peak = {}", proj->GetName(), valuePeak) << std::endl;

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
      par[2]=0.02;
      par[3]=1;
      par[4]=1;
      par[5]=1;
      par[6]=1;
      fcb.SetParameters(&par[0]);

      fcb.FixParameter(1, xPeak);
      fcb.FixParameter(2, 0.02);
      fcb.SetParLimits(3, 0.01f, 2.f);
      fcb.SetParLimits(5, 0.01f, 2.f);
      fcb.SetParLimits(4, 0.5f, 10.f);
      fcb.SetParLimits(6, 0.5f, 10.f);
      proj->Fit("fcb", "BRNQ");
      fcb.SetParLimits(2, 0.01, 0.1);
      fcb.ReleaseParameter(1);
      const char* fitOpt = printFits ? "BRS" : "BRQS";
      TFitResultPtr fitResult = proj->Fit("fcb", fitOpt);
      if (fitResult->Status() != 0) {
        fitResult->Print();
      }

      if (printFits) {
        proj->Draw("E");
        c.SaveAs(pdfFileName.c_str());
      }

      if (fitResult->Status() == 0) {
        mean = fcb.GetParameter(1);
        meanErr = fcb.GetParError(1);
        sigma = fcb.GetParameter(2);
        sigmaErr = fcb.GetParError(2);
      } else {
        mean = 0;
        meanErr = 1.e6;
        sigma = 1.e6;
        sigmaErr = 1.e6;
      }

      if (mean < meanMin) meanMin = mean;
      if (mean > meanMax) meanMax = mean;
    }

    histogramMean->SetBinContent(bin, mean);
    histogramMean->SetBinError(bin, meanErr);
  }

  TF1 linFit("linFit", "pol1");
  linFit.SetLineColor(kRed);
  linFit.SetLineStyle(kDotted);
  TFitResultPtr fitResult = histogramMean->Fit("linFit", "BQS");
  if (fitResult->Status() != 0) {
    fitResult->Print();
  }

  //histogramMean->SetMinimum(meanMin);
  //histogramMean->SetMaximum(meanMax);
  histogramMean->SetMinimum(-0.1);
  histogramMean->SetMaximum(0.1);
  histogramMean->Draw("E");

  TLine* line1 = new TLine(histogramMean->GetXaxis()->GetXmin(), 0, histogramMean->GetXaxis()->GetXmax(), 0);
  line1->SetLineStyle(kDashed);
  line1->Draw();

  TLine* line2 = new TLine(0, yMin, 0, yMax);
  line2->SetLineStyle(kDashed);
  line2->Draw();

  //std::cout << std::format("Amplitude: {:0.4f}", phiFit.GetParameter(0)) << std::endl;
  //std::cout << std::format("Phase:     {:0.4f}", phiFit.GetParameter(1)) << std::endl;
  //std::cout << std::format("Offsets:   top={:0.4f} bottom={:0.4f}", phiFit.GetParameter(2), phiFit.GetParameter(3)) << std::endl;

  return {linFit.GetParameter(1), linFit.GetParError(1)};
}

std::pair<double, double> PlotDCASlopeProjection(std::string histName, float yMin, float yMax, int projRebin, TCanvas& c, bool printFits = false)
{
  std::string fullHistName = std::string("qa-muon/alignment/") + histName;
  TH2* histogram2 = GetTH2(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogram2 << std::endl;
  if (!histogram2)
    return {};
  histogram2->SetName(fullHistName.c_str());

  return PlotDCASlopeProjection(histogram2, yMin, yMax, projRebin, c, printFits);
}

void PlotDCASlopeProjection3D(std::string histName, float yMin, float yMax, int projRebin, TCanvas& c, bool printFits = false)
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
    //if (bin != 16) continue;
    histogram3->GetXaxis()->SetRange(bin, bin);
    TH2* proj = (TH2*)histogram3->Project3D("zy");
    proj->SetTitle(TString::Format("%s - z bin %d", histogram3->GetTitle(), bin));
    c.Clear();
    proj->Draw("col");
    c.SaveAs(pdfFileName.c_str());

    c.Clear();
    auto amplitude = PlotDCASlopeProjection(proj, yMin, yMax, projRebin, c, printFits);
    c.SaveAs(pdfFileName.c_str());

    //gr->AddPointError(histogram3->GetXaxis()->GetBinCenter(bin), amplitude.first, 0, amplitude.second);
    //std::cout << std::format("Bin #{} -> {}", bin, histogram3->GetXaxis()->GetBinCenter(bin)) << std::endl;
    xv.push_back(histogram3->GetXaxis()->GetBinCenter(bin));
    exv.push_back(0);
    yv.push_back(amplitude.first);
    eyv.push_back(amplitude.second);

  }
  TGraphErrors* gr = new TGraphErrors(xv.size(), xv.data(), yv.data(), exv.data(), eyv.data());
  gr->Draw("A*");
  gr->SetTitle("Correlation amplitude vs. z shift");
  gr->GetXaxis()->SetTitle("z shift (mm)");
  gr->GetYaxis()->SetTitle("correlation");

  TLine* line1 = new TLine(gr->GetHistogram()->GetXaxis()->GetXmin(), 0, gr->GetHistogram()->GetXaxis()->GetXmax(), 0);
  line1->SetLineStyle(kDashed);
  line1->Draw();

  TLine* line2 = new TLine(0, gr->GetHistogram()->GetMinimum(), 0, gr->GetHistogram()->GetMaximum());
  line2->SetLineStyle(kDashed);
  line2->Draw();

  TF1 linFit("linFit", "pol1");
  linFit.SetLineColor(kRed);
  linFit.SetLineStyle(kDotted);
  TFitResultPtr fitResult = gr->Fit("linFit", "BQS");
  if (fitResult->Status() != 0) {
    fitResult->Print();
  }

  std::cout << std::format("Optimal z shift: {} cm", -linFit.GetParameter(0) / linFit.GetParameter(1) / 10.f) << std::endl;
}


std::array<double, 4> FitDCA(TH1* proj, TCanvas& c)
{
  double mean{ nan("") };
  double meanErr{ nan("") };
  double sigma{ nan("") };
  double sigmaErr{ nan("") };

  int peakMin = 10;
  if (proj->GetEntries() < peakMin) {
    return {mean, meanErr, sigma, sigmaErr};
  }

  int valuePeak = proj->GetMaximum();
  int binPeak = proj->GetMaximumBin();
  double xPeak = proj->GetXaxis()->GetBinCenter(binPeak);

  //std::cout << std::format("{}  peak={}", proj->GetName(), valuePeak) << std::endl;
  if (valuePeak < peakMin) {
    return {mean, meanErr, sigma, sigmaErr};
  }

  //return {xPeak, 0, 0, 0};

  if (true) {
    TF1 fgaus2("fgaus2", "gausn", xPeak - 0.02, xPeak + 0.02);
    fgaus2.SetParameter(1, xPeak);
    fgaus2.SetParameter(2, 0.02);
    fgaus2.SetParLimits(2, 0.01, 0.1);
    const char* fitOpt2 = "BQRSN";
    TFitResultPtr fitResult2 = proj->Fit("fgaus2", fitOpt2);
    if (fitResult2.Get() && fitResult2->Status() == 0) {
      xPeak = fgaus2.GetParameter(1);
    }

    TF1 fcb("fgaus", "gausn", xPeak - 0.02, xPeak + 0.02);
    fcb.SetNpx(1000);
    fcb.SetLineColor(kBlack);
    fcb.SetParameter(1, xPeak);
    fcb.SetParameter(2, 0.02);
    fcb.SetParLimits(2, 0.01, 0.1);
    const char* fitOpt = "BQRS";
    TFitResultPtr fitResult = proj->Fit("fgaus", fitOpt);

    if (fitResult.Get()) {
      if (fitResult->Status() == 0) {
        mean = fcb.GetParameter(1);
        meanErr = fcb.GetParError(1);
        sigma = fcb.GetParameter(2);
        sigmaErr = fcb.GetParError(2);
      } else {
        //std::cout << "Entries: " << proj->GetEntries() << std::endl;
        //fitResult->Print();
      }
    }
  } else {
    TF1 fcb("fcb", DoubleSidedCB, -0.1, 0.1, 7);
    fcb.SetNpx(1000);
    fcb.SetLineColor(kBlack);
    //fcb.SetParameter(0, valuePeak);
    double par[7];
    par[0]=35000;
    par[1]=xPeak;
    par[2]=0.02;
    par[3]=1;
    par[4]=1;
    par[5]=1;
    par[6]=1;
    fcb.SetParameters(&par[0]);

    fcb.FixParameter(1, xPeak);
    fcb.FixParameter(2, 0.02);
    fcb.SetParLimits(3, 0.01f, 2.f);
    fcb.SetParLimits(5, 0.01f, 2.f);
    fcb.SetParLimits(4, 0.5f, 10.f);
    fcb.SetParLimits(6, 0.5f, 10.f);
    //proj->Fit("fcb", "BRNQ");
    fcb.ReleaseParameter(2);
    fcb.SetParameter(2, 0.02);
    fcb.SetParLimits(2, 0.01, 0.1);
    fcb.ReleaseParameter(1);
    fcb.SetParameter(1, xPeak);
    const char* fitOpt = "BRQS";
    TFitResultPtr fitResult = proj->Fit("fcb", fitOpt);

    if (fitResult.Get() && fitResult->Status() == 0) {
      mean = fcb.GetParameter(1);
      meanErr = fcb.GetParError(1);
      sigma = fcb.GetParameter(2);
      sigmaErr = fcb.GetParError(2);
    }
  }

  //std::cout << std::format("  {}  peak = {}    mean = {:0.3f} +/- {:0.3f}", proj->GetName(), valuePeak, mean, meanErr) << std::endl;

  //proj->Draw();
  //c.SaveAs(pdfFileName.c_str());

  return {mean, meanErr, sigma, sigmaErr};
}

void PlotDCAXY(std::string coordinate, TCanvas& c)
{
  const int xyRebin = 1;

  std::string fullHistName = std::string("muon-global-alignment/DCA/MFT/DCA_") + coordinate;
  THnSparse* histogramFull = GetTHnSparse(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogramFull << std::endl;
  if (!histogramFull) {
    return;
  }

  // DCA axis assignments:
  // 0 -> DCA
  // 1 -> vertex z shift
  // 2 -> track x
  // 3 -> track y

  c.cd();
  TH2* histDCA = histogramFull->Projection(3, 2);
  if (!histDCA) return;

  for (int i = 0; i < 4; i++) {
    std::cout << std::format("Axis #{}: {}", i, histogramFull->GetAxis(i)->GetTitle()) << std::endl;
  }

  for (int zbin = 1; zbin <= histogramFull->GetAxis(1)->GetNbins(); zbin++) {
    //if (zbin != 5) continue;
    histogramFull->GetAxis(1)->SetRange(zbin, zbin);

    histDCA->Reset();
    histDCA->SetTitle(TString::Format("%s - vz = %0.2f", histogramFull->GetTitle(), histogramFull->GetAxis(1)->GetBinCenter(zbin)));

    TH3* h3 = histogramFull->Projection(0, 2, 3);
    std::cout << "3D projection:" << std::endl
        << "  X axis: " << h3->GetXaxis()->GetTitle() << std::endl
        << "  Y axis: " << h3->GetYaxis()->GetTitle() << std::endl
        << "  Z axis: " << h3->GetZaxis()->GetTitle() << std::endl;

    TH1* proj = (TH1*)h3->ProjectionX();
    proj->Draw();
    c.SaveAs(pdfFileName.c_str());
    //continue;

    TH1* projB = (TH1*)h3->ProjectionX(TString::Format("%s_%d_bottom", h3->GetName(), zbin),
        0, -1, 1, 48);
    //0, -1, 0, h3->GetZaxis()->GetNbins()/2);
    TH1* projT = (TH1*)h3->ProjectionX(TString::Format("%s_%d_top", h3->GetName(), zbin),
        0, -1, 49, 48 * 2);
        //0, -1, h3->GetZaxis()->GetNbins()/2 + 1, -1);

    projB->SetLineColor(kGreen);
    projB->Draw("hist same");
    projT->SetLineColor(kRed);
    projT->Draw("hist same");
    proj->GetXaxis()->SetRangeUser(-0.1, 0.1);
    c.SaveAs(pdfFileName.c_str());

    //for (int ybin = 1; ybin <= histogramFull->GetAxis(3)->GetNbins(); ybin++) {
    for (int ybin = 1; ybin <= h3->GetZaxis()->GetNbins(); ybin += xyRebin) {
      //if (ybin != histogramFull->GetAxis(3)->GetNbins()/2) continue;
      //if (ybin < h3->GetZaxis()->GetNbins()/2) continue;
      //if (ybin < 9 || ybin > 12) continue;
      //if (ybin < 50 || ybin > 90) continue;
      //histogramFull->GetAxis(3)->SetRange(ybin, ybin);
      //std::cout << std::format("Starting row {}", ybin) << std::endl;

      //for (int xbin = 1; xbin <= histogramFull->GetAxis(2)->GetNbins(); xbin++) {
      for (int xbin = 1; xbin <= h3->GetYaxis()->GetNbins(); xbin += xyRebin) {
        //if (xbin != (histogramFull->GetAxis(2)->GetNbins()/2 + 10)) continue;
        //if (xbin < h3->GetYaxis()->GetNbins()/2) continue;
        //if (xbin != 21) continue;
        //if (xbin < 25 || xbin > 125) continue;
        //if (xbin < 14 || xbin > 17) continue;
        //histogramFull->GetAxis(2)->SetRange(xbin, xbin);
        //std::cout << std::format("  column {}", xbin) << std::endl;

        //TH1* proj = histogramFull->Projection(0);
        TH1* proj = (TH1*)h3->ProjectionX(TString::Format("%s_%d_%d_%d", h3->GetName(), zbin, xbin, ybin),
            xbin, xbin + xyRebin - 1, ybin, ybin + xyRebin - 1);
        if (!proj) continue;
        //proj->SetName(TString::Format("%s_%d_%d_%d", proj->GetName(), nclusbin, xbin, ybin));
        //proj->GetXaxis()->SetRangeUser(-0.1f, 0.1f);

        auto fitResult = FitDCA(proj, c);

        for (int ybin2 = ybin; ybin2 < ybin + xyRebin; ybin2++) {
          for (int xbin2 = xbin; xbin2 < xbin + xyRebin; xbin2++) {
            if (std::isnan(fitResult[0]) || std::isnan(fitResult[1])) {
              histDCA->SetBinContent(xbin2, ybin2, -100);
              histDCA->SetBinError(xbin2, ybin2, 0);
            } else {
              histDCA->SetBinContent(xbin2, ybin2, fitResult[0]);
              histDCA->SetBinError(xbin2, ybin2, fitResult[1]);
            }
          }
        }
      }
      std::cout << std::format("Row {} done.", ybin) << std::endl;
    }

    histDCA->SetMinimum(-0.025);
    histDCA->SetMaximum(0.025);
    histDCA->Draw("colz");
    c.SaveAs(pdfFileName.c_str());
  }
}


void muonGlobalAlignmentVertexShift(const char* _rootFileName = "AnalysisResults.root", const char* _pdfFileName = "vertexShift.pdf")
{
  //fAnalysisResults = new TFile("AnalysisResults.root");
  //fAnalysisResults = new TFile("AnalysisResults/AnalysisResultsFull.root");
  fAnalysisResults = new TFile(_rootFileName);
  pdfFileName = _pdfFileName;

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
  c.SaveAs((pdfFileName + "(").c_str());

  TH1* h1 = GetTH1(fAnalysisResults, "qa-muon/alignment/DCA/vertex_z");
  if (h1) {
    h1->Draw();
    c.SaveAs(pdfFileName.c_str());
  }

  TH2* h2 = PlotDCAMFT("DCA/MFT/DCA_y_vs_x");
  if (h2) {
    h2->GetXaxis()->SetRangeUser(-0.1, 0.1);
    h2->GetYaxis()->SetRangeUser(-0.1, 0.1);
    c.SaveAs(pdfFileName.c_str());

    h1 = h2->ProjectionX();
    h1->SetTitle("DCA(x)");
    h1->Draw();
    c.SaveAs(pdfFileName.c_str());

    h1 = h2->ProjectionY();
    h1->SetTitle("DCA(y)");
    h1->Draw();
    c.SaveAs(pdfFileName.c_str());
  }

  c.Clear();
  //PlotDCAXY("x", c);
  //PlotDCAXY("y", c);

  c.Clear();
  //PlotDCAPhiProjection3D("DCA/MFT/DCA_x_vs_phi_vs_zshift", -0.1, 0.1, 1, c,
  //    -TMath::Pi()/4.f + TMath::Pi()/2.f, TMath::Pi()/4.f + TMath::Pi()/2.f, false);
  //c.SaveAs(pdfFileName.c_str());
  PlotDCASlopeProjection3D("DCA/MFT/DCA_x_vs_slopex_vs_zshift", -0.1, 0.1, 1, c, false);
  c.SaveAs(pdfFileName.c_str());

  //PlotDCAPhiProjection3D("DCA/MFT/DCA_y_vs_phi_vs_zshift", -0.1, 0.1, 1, c,
  //    -TMath::Pi()/4.f, TMath::Pi()/4.f, false);
  //c.SaveAs(pdfFileName.c_str());
  PlotDCASlopeProjection3D("DCA/MFT/DCA_y_vs_slopey_vs_zshift", -0.1, 0.1, 1, c, false);
  c.SaveAs(pdfFileName.c_str());

  c.Clear();
  c.SaveAs((pdfFileName + ")").c_str());
}
