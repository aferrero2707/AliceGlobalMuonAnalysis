#include <MFTTracking/Constants.h>

#include <deque>

TFile* fAnalysisResults;
std::string pdfFileName;
std::string taskName = "muon-global-alignment";
//std::string taskName = "muon-global-alignment_id48147";


TH1* GetTH1(TFile* f, TString histname)
{
  TH1* result = (TH1*)f->Get(histname);
  std::cout << "histname: " << histname << " -> " << result << std::endl;
  return result;
}

TH2* GetTH2(TFile* f, TString histname)
{
  TH2* result = (TH2*)f->Get(histname);
  std::cout << "histname: " << histname << " -> " << result << std::endl;
  return result;
}

TH3* GetTH3(TFile* f, TString histname)
{
  TH3* result = (TH3*)f->Get(histname);
  std::cout << "histname: " << histname << " -> " << result << std::endl;
  return result;
}


THnSparse* GetTHnSparse(TFile* f, TString histname)
{
  THnSparse* result = (THnSparse*)f->Get(histname);
  std::cout << "histname: " << histname << " -> " << result << std::endl;
  return result;
}




std::array<bool, 10> GetFiredLayers(int pattern)
{
  std::array<bool, 10> result;
  std::bitset<10> bits{static_cast<unsigned long long>(pattern)};
  for (int bit = 0; bit < 10; bit++) {
    result[bit] = bits.test(bit);
  }
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

double DoubleSidedCBwithQuadBgd(double* x, double *par)
{
  return(par[0] * DoubleSidedCB2(x[0], par[1],par[2],par[3],par[4],par[5],par[6]) + par[7] + x[0] * par[8] + x[0] * x[0] * par[9]);
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
    fgaus2.SetLineColor(kRed);
    const char* fitOpt2 = "BQRSN";
    TFitResultPtr fitResult2 = proj->Fit("fgaus2", fitOpt2);
    if (fitResult2.Get() && fitResult2->Status() == 0) {
      xPeak = fgaus2.GetParameter(1);
    }

    TF1 fcb("fgaus", "gausn", xPeak - 0.025, xPeak + 0.025);
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

std::tuple<float, float, float, float> PlotDCAProjection(TH2* histogram2, float yMin, float yMax, int projRebin, TCanvas& c, bool printFits = false)
{
  std::tuple<float, float, float, float> result{ nan(""), nan(""), nan(""), nan("") };
  std::string fullHistName = histogram2->GetName();
  histogram2->GetYaxis()->UnZoom();

  TH1* histogramMean = histogram2->ProjectionX();
  histogramMean->SetName(TString::Format("%s-mean", fullHistName.c_str()));
  histogramMean->SetTitle(TString::Format("%s", histogram2->GetTitle()));
  histogramMean->SetMinimum(yMin);
  histogramMean->SetMaximum(yMax);

  int nValidBins = 0;
  int entriesMin = histogram2->GetEntries() / 100;
  int peakMin = 50;
  // skip first bin because it is biased?
  for (int bin = 1; bin <= histogram2->GetXaxis()->GetNbins(); bin++) {
    TH1* proj = histogram2->ProjectionY(TString::Format("%s-%d", fullHistName.c_str(), bin), bin, bin);

    //proj->Rebin(projRebin);
    //std::cout << std::format("range: {:0.2f} -> {:0.2f}", histogram2->GetYaxis()->GetXmin(), histogram2->GetYaxis()->GetXmax()) << std::endl;
    //std::cout << std::format("range: {:0.2f} -> {:0.2f}", proj->GetXaxis()->GetXmin(), proj->GetXaxis()->GetXmax()) << std::endl;
    proj->GetXaxis()->SetRangeUser(-0.1f, 0.1f);

    proj->SetTitle(TString::Format("%s - bin %d", histogram2->GetTitle(), bin));
    double mean{ nan("") };
    double meanErr{ nan("") };
    double sigma{ nan("") };
    double sigmaErr{ nan("") };

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
      /*if (printFits)
        proj->Fit("fcb", "BR");
      else
        proj->Fit("fcb", "BRQ");*/
      //std::cout << std::format("Fit status: {}", fitResult->Status()) << std::endl;
      //fitResult->Print();

      if (false && printFits) {
        proj->Draw("E");
        c.SaveAs(pdfFileName.c_str());
      }

      if (fitResult.Get() && fitResult->Status() == 0) {
        mean = fcb.GetParameter(1);
        meanErr = fcb.GetParError(1);
        sigma = fcb.GetParameter(2);
        sigmaErr = fcb.GetParError(2);
      }

      //mean = fgaus.GetParameter(1);
      //meanErr = fgaus.GetParError(1);
      //sigma = fgaus.GetParameter(2);
      //sigmaErr = fgaus.GetParError(2);
    }

    if (std::isnan(mean) || std::isnan(meanErr)) continue;

    histogramMean->SetBinContent(bin, mean);
    histogramMean->SetBinError(bin, meanErr);

    nValidBins += 1;
  }

  TF1 linFit("linFit", "pol1");
  linFit.SetLineColor(kRed);
  linFit.SetLineStyle(kDotted);
  if (nValidBins > 0) {
    std::cout << std::format("Fitting \"{}\" with {} valid bins", fullHistName, nValidBins) << std::endl;
    if (printFits) {
      histogramMean->Fit("linFit", "");
    } else {
      histogramMean->Fit("linFit", "Q");
    }
    std::cout << std::format("Fitting done") << std::endl;

    if (std::isnan(linFit.GetParameter(0)) || std::isnan(linFit.GetParameter(1))) {
      for (int bin = 1; bin < histogramMean->GetXaxis()->GetNbins(); bin++) {
        std::cout << std::format("  bin#{}: {:0.3f} +/- {:0.3f}",
            bin, histogramMean->GetBinContent(bin), histogramMean->GetBinError(bin)) << std::endl;
      }
    }

    result = {linFit.GetParameter(0), linFit.GetParError(0), linFit.GetParameter(1), linFit.GetParError(1)};
  }
  histogramMean->Draw("E");

  TLine* line1 = new TLine(histogramMean->GetXaxis()->GetXmin(), 0, histogramMean->GetXaxis()->GetXmax(), 0);
  line1->SetLineStyle(kDashed);
  line1->Draw();

  TLine* line2 = new TLine(0, yMin, 0, yMax);
  line2->SetLineStyle(kDashed);
  line2->Draw();

  //std::cout << std::format("Slope: {:0.4f} mm / 10 m", linFit.GetParameter(1) * 1000 * 10) << std::endl;
  //c.SaveAs(pdfFileName.c_str());
  //histogramSigma->Draw("E");
  //c.SaveAs(pdfFileName.c_str());

  return result; //{linFit.GetParameter(0), linFit.GetParError(0), linFit.GetParameter(1), linFit.GetParError(1)};
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

void PlotTracks(std::string histName, int numberOfClusters, TCanvas& c)
{
  std::string fullHistName = taskName + "/" + histName;
  THnSparse* histogramFull = GetTHnSparse(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogramFull << std::endl;
  if (!histogramFull) return;

  // Layers pattern axis assignments:
  // 0 -> layers pattern
  // 1 -> track x
  // 2 -> track y
  // 3 -> nClusters
  // 4 -> MFT track type

  int binClus = numberOfClusters - 4;


  c.cd();
  //TH2* trackxy = histogramFull->Projection(3, 2);
  //if (!trackxy) return false;
  //trackxy->Draw("col");
  //c.SaveAs(pdfFileName.c_str());

  histogramFull->GetAxis(3)->SetRange(binClus, binClus);
  histogramFull->GetAxis(4)->SetRange(1, 1);

  TH2* trackxy = histogramFull->Projection(2, 1);
  if (!trackxy) return;
  trackxy->SetTitle(std::format("Track (x,y) with nClus={}", numberOfClusters).c_str());
  trackxy->Draw("col");
  c.SaveAs(pdfFileName.c_str());

  std::array<TH2*, 10> hxy;
  for (int layer = 0; layer < 10; layer++) {
    hxy[layer] = (TH2*)trackxy->Clone(TString::Format("%s_%d", trackxy->GetName(), layer));
    hxy[layer]->SetTitle(std::format("Track (x,y) with nClus={} and layer={} fired", numberOfClusters, layer).c_str());
  }

  for (int ybin = 1; ybin <= histogramFull->GetAxis(2)->GetNbins(); ybin++) {
    histogramFull->GetAxis(2)->SetRange(ybin, ybin);
    //std::cout << std::format("ybin={}", ybin) << std::endl;
    for (int xbin = 1; xbin <= histogramFull->GetAxis(1)->GetNbins(); xbin++) {
      histogramFull->GetAxis(1)->SetRange(xbin, xbin);
      //std::cout << std::format("  xbin={}", ybin) << std::endl;
      TH1* layers = histogramFull->Projection(0);
      if (!layers) continue;
      layers->SetName(TString::Format("%s_%d_%d_%d", layers->GetName(), numberOfClusters, xbin, ybin));

      TH1F h("h", "Fired layers", 10, 0, 10);
      double integral = 0;
      for (int i = 1; i <= layers->GetXaxis()->GetNbins(); i++) {
        auto entries = layers->GetBinContent(i);
        integral += entries;
        auto firedLayers = GetFiredLayers(i - 1);
        for (int bit = 0; bit < 10; bit++) {
          if (firedLayers[bit]) {
            h.Fill(bit, entries);
          }
        }
      }
      if (integral > 0)
        h.Scale(1.f / integral);
      for (int layer = 0; layer < 10; layer++) {
        //std::cout << std::format("    layer {} = {}", layer, h.GetBinContent(layer + 1)) << std::endl;
        hxy[layer]->SetBinContent(xbin, ybin, h.GetBinContent(layer + 1));
      }
    }
  }

  for (int layer = 0; layer < 10; layer++) {
    hxy[layer]->Draw("colz");
    c.SaveAs(pdfFileName.c_str());
  }
}


void PlotMftLayerEfficiencies(TCanvas& c)
{
  if (!GetTH2(fAnalysisResults, taskName + "/DCA/MFT/mftTrackEffNum_0")) return;

  TPaveText* title = new TPaveText(0.1, 0.4, 0.9, 0.6, "NDC");
  c.Clear();
  title->AddText("MFT layer efficiencies");
  title->Draw();
  c.SaveAs(pdfFileName.c_str());

  c.Clear();
  for (int i = 0; i < 10; i++) {
    std::string fullHistName = taskName + "/DCA/MFT/mftTrackEffNum_" + std::to_string(i);
    TH2* histogramNum = GetTH2(fAnalysisResults, fullHistName);
    std::cout << fullHistName << " -> " << histogramNum << std::endl;
    if (!histogramNum) continue;

    fullHistName = taskName + "/DCA/MFT/mftTrackEffDen_" + std::to_string(i);
    TH2* histogramDen = GetTH2(fAnalysisResults, fullHistName);
    std::cout << fullHistName << " -> " << histogramDen << std::endl;
    if (!histogramDen) continue;

    histogramNum->Divide(histogramDen);
    histogramNum->SetMinimum(0);
    histogramNum->SetMaximum(1);

    if (i < 4) {
      float mftLadderWidth = 1.7;
      histogramNum->GetXaxis()->SetRangeUser(-mftLadderWidth * 15.f / 2.f, mftLadderWidth * 15.f / 2.f);
      histogramDen->GetXaxis()->SetRangeUser(-mftLadderWidth * 15.f / 2.f, mftLadderWidth * 15.f / 2.f);
      histogramNum->GetYaxis()->SetRangeUser(-10.f, 10.f);
      histogramDen->GetYaxis()->SetRangeUser(-10.f, 10.f);
    }

    histogramDen->Draw("col");
    c.SaveAs(pdfFileName.c_str());

    histogramNum->Draw("col");
    c.SaveAs(pdfFileName.c_str());
  }
}


bool ProcessDCAvsZ(std::string coordinate, int numberOfClusters, TCanvas& c, TCanvas& c2, bool printFits = false)
{
  std::string fullHistName = taskName + "/DCA/MFT/DCA_" + coordinate;
  THnSparse* histogramFull = GetTHnSparse(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogramFull << std::endl;
  if (!histogramFull) return false;

  fullHistName = taskName + "/DCA/MFT/layers";
  THnSparse* histogramLayers = GetTHnSparse(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogramLayers << std::endl;

  fullHistName = taskName + "/DCA/MFT/trackMomentum";
  THnSparse* histogramMomentum = GetTHnSparse(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogramMomentum << std::endl;

  fullHistName = taskName + "/DCA/MFT/trackChi2";
  THnSparse* histogramChi2 = GetTHnSparse(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogramChi2 << std::endl;

  fullHistName = taskName + "/DCA/MFT/slope_" + coordinate;
  THnSparse* histogramSlope = GetTHnSparse(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogramSlope << std::endl;

  // DCA axis assignments:
  // 0 -> DCA
  // 1 -> vertex z
  // 2 -> track x
  // 3 -> track y
  // 4 -> nClusters
  // 5 -> layers

  // Layers pattern axis assignments:
  // 0 -> layers pattern
  // 1 -> track x
  // 2 -> track y
  // 3 -> nClusters

  // Track slope axis assignments:
  // 0 -> slope
  // 1 -> track x
  // 2 -> track y
  // 3 -> nClusters
  // 4 -> layers

  //int numberOfClusters = 7;
  int binClus = numberOfClusters - 4;

  int layerBinMax = 0; //9;

  c.cd();
  TH2* trackxy = histogramFull->Projection(3, 2);
  if (!trackxy) return false;

  /*
  trackxy->Draw("col");
  c.SaveAs(pdfFileName.c_str());

  for (int nclBin = 1; nclBin <= histogramFull->GetAxis(4)->GetNbins(); nclBin++) {
    histogramFull->GetAxis(4)->SetRange(nclBin, nclBin);
    for (int ttypeBin = 1; ttypeBin <= histogramFull->GetAxis(5)->GetNbins(); ttypeBin++) {
      histogramFull->GetAxis(5)->SetRange(ttypeBin, ttypeBin);
      trackxy = histogramFull->Projection(3, 2);
      if (!trackxy) continue;
      trackxy->SetTitle(std::format("Track (x,y) with nClus={} and trackType={}", nclBin + 4, (ttypeBin == 1 ? "LF" : "CA")).c_str());
      trackxy->Draw("col");
      c.SaveAs(pdfFileName.c_str());
    }
  }
*/
  //histogramFull->GetAxis(4)->SetRange(1, histogramFull->GetAxis(4)->GetNbins());
  histogramFull->GetAxis(4)->SetRange(binClus, binClus);
  //histogramFull->GetAxis(5)->SetRange(1, histogramFull->GetAxis(5)->GetNbins());

  if (histogramMomentum) {
    histogramMomentum->GetAxis(3)->SetRange(binClus, binClus);
  }

  if (histogramLayers) {
    histogramLayers->GetAxis(3)->SetRange(binClus, binClus);
  }

  trackxy->Reset();

  std::array<TH2*, 11> histDCA;
  std::array<TH2*, 11> histSlope;
  std::array<TH2*, 11> histSigma;
  for (int i = 0; i < 11; i++) {
    histDCA[i] = (TH2*)trackxy->Clone();
    histSlope[i] = (TH2*)trackxy->Clone();
    histSigma[i] = (TH2*)trackxy->Clone();
    if (i == 0) {
      histDCA[i]->SetTitle(std::format("DCA({}) offset vs. track XY ({} clusters)", coordinate, numberOfClusters).c_str());
      histSlope[i]->SetTitle(std::format("DCA({}) slope vs. track XY ({} clusters)", coordinate, numberOfClusters).c_str());
      histSigma[i]->SetTitle(std::format("#sigma_{{DCA({})}} vs. track XY ({} clusters)", coordinate, numberOfClusters).c_str());
    } else {
      histDCA[i]->SetTitle(std::format("DCA({}) offset vs. track XY ({} clusters, layer {})", coordinate, numberOfClusters, i - 1).c_str());
      histSlope[i]->SetTitle(std::format("DCA({}) slope vs. track XY ({} clusters, layer {})", coordinate, numberOfClusters, i - 1).c_str());
      histSigma[i]->SetTitle(std::format("#sigma_{{DCA({})}} vs. track XY ({} clusters, layer {})", coordinate, numberOfClusters, i - 1).c_str());
    }
  }

  c2.cd();
  for (int ybin = 1; ybin <= histogramFull->GetAxis(3)->GetNbins(); ybin++) {
    //if (ybin != 11) continue;
    //if (ybin < 9 || ybin > 12) continue;
    if (ybin < 15 || ybin > 20) continue;
    histogramFull->GetAxis(3)->SetRange(ybin, ybin);
    for (int xbin = 1; xbin <= histogramFull->GetAxis(2)->GetNbins(); xbin++) {
      //if (xbin != 21) continue;
      //if (xbin < 20 || xbin > 25) continue;
      if (xbin < 14 || xbin > 17) continue;
      histogramFull->GetAxis(2)->SetRange(xbin, xbin);
      for (int lbin = 0; lbin <= histogramFull->GetAxis(5)->GetNbins(); lbin++) {
        if (lbin > layerBinMax) break;

        if (lbin == 0) {
          histogramFull->GetAxis(5)->SetRange(1, histogramFull->GetAxis(5)->GetNbins());
        } else {
          histogramFull->GetAxis(5)->SetRange(lbin, lbin);
        }

        TH2* dcaVsZ = histogramFull->Projection(0, 1);
        std::cout << std::format("ybin={} xbin={} lbin={} dcaVsZ={}", ybin, xbin, lbin, (void*)dcaVsZ) << std::endl;
        if (!dcaVsZ) continue;
        dcaVsZ->SetName(TString::Format("%s_%d_%d_%d_%d", dcaVsZ->GetName(), numberOfClusters, xbin, ybin, lbin));
        //if (dcaVsZ->GetEntries() < 1000) continue;

        std::string dcaFileName = std::format("mft_dca{}_vs_z_{}clus_{}_{}_l{}.pdf", coordinate, numberOfClusters, xbin, ybin, lbin);

        dcaVsZ->GetYaxis()->SetRangeUser(-0.1, 0.1);
        dcaVsZ->Draw("col");
        c2.SaveAs((dcaFileName+"(").c_str());

        auto linfit = PlotDCAProjection(dcaVsZ, -0.1, 0.1, 1, c2, printFits);
        c2.SaveAs(dcaFileName.c_str());

        dcaVsZ->GetYaxis()->SetRangeUser(-0.1, 0.1);
        dcaVsZ->Draw("col");
        if (!std::isnan(std::get<0>(linfit)) && !std::isnan(std::get<2>(linfit))) {
          TLine* line = new TLine(dcaVsZ->GetXaxis()->GetXmin(),
              dcaVsZ->GetXaxis()->GetXmin() * std::get<2>(linfit) + std::get<0>(linfit),
              dcaVsZ->GetXaxis()->GetXmax(),
              dcaVsZ->GetXaxis()->GetXmax() * std::get<2>(linfit) + std::get<0>(linfit));
          line->SetLineColor(kRed);
          line->SetLineStyle(kDashed);
          line->SetLineWidth(2);
          line->Draw();
          //std::cout << std::format("TOTO bin {}-{}: {:0.2f} +/- {:0.2f}", ybin, xbin, std::get<0>(linfit), std::get<1>(linfit)) << std::endl;
          histDCA[lbin]->SetBinContent(xbin, ybin, std::get<0>(linfit));
          histDCA[lbin]->SetBinError(xbin, ybin, 0.1);

          histSlope[lbin]->SetBinContent(xbin, ybin, std::get<2>(linfit));
          histSlope[lbin]->SetBinError(xbin, ybin, 0.1);
        }
        c2.SaveAs(dcaFileName.c_str());

        TH1* proj = dcaVsZ->ProjectionY(std::format("{}_py", dcaVsZ->GetName()).c_str(),
            dcaVsZ->GetXaxis()->FindBin(-0.5),
            dcaVsZ->GetXaxis()->FindBin(0.5));
        proj->SetLineColor(kRed);
        proj->Draw("E");

        int valuePeak = proj->GetMaximum();
        int binPeak = proj->GetMaximumBin();
        double xPeak = proj->GetXaxis()->GetBinCenter(binPeak);
        int peakMin = 50;

        if (valuePeak > peakMin) {
          TF1 fcb("fcb", DoubleSidedCBwithQuadBgd, -0.5, 0.5, 10);
          fcb.SetNpx(1000);
          fcb.SetLineColor(kBlack);
          //fcb.SetParameter(0, valuePeak);
          double par[10];
          par[0]=35000;
          par[1]=xPeak;
          par[2]=0.01;
          par[3]=1;
          par[4]=1;
          par[5]=1;
          par[6]=1;
          par[7]=0;
          par[8]=0;
          par[9]=0;
          fcb.SetParameters(&par[0]);

          //fcb.FixParameter(1, xPeak);
          fcb.SetParLimits(2, 0.005, 0.5);
          //proj->Fit("fcb", "BRN");
          //fcb.ReleaseParameter(1);
          if (printFits)
            proj->Fit("fcb", "BR");
          else
            proj->Fit("fcb", "BR");

          if (!std::isnan(fcb.GetParameter(2)) && !std::isnan(fcb.GetParError(2))) {
            histSigma[lbin]->SetBinContent(xbin, ybin, fcb.GetParameter(2));
            histSigma[lbin]->SetBinError(xbin, ybin, 0.1);
            std::cout << std::format("peak: {:0.3f} +/- {:0.3f}", fcb.GetParameter(1), fcb.GetParameter(2)) << std::endl;
          }
        }
        c2.SaveAs(dcaFileName.c_str());

        if (histogramMomentum) {
          histogramMomentum->GetAxis(1)->SetRange(xbin, xbin);
          histogramMomentum->GetAxis(2)->SetRange(ybin, ybin);
          histogramMomentum->GetAxis(4)->SetRange(lbin, lbin);
          TH1* h = histogramMomentum->Projection(0);
          if (h) {
            h->GetXaxis()->SetRangeUser(0, 30);
            h->Draw();
            c2.SaveAs(dcaFileName.c_str());
          }
        }

        if (histogramChi2) {
          histogramChi2->GetAxis(1)->SetRange(xbin, xbin);
          histogramChi2->GetAxis(2)->SetRange(ybin, ybin);
          histogramChi2->GetAxis(4)->SetRange(lbin, lbin);
          TH1* h = histogramChi2->Projection(0);
          if (h) {
            h->GetXaxis()->SetRangeUser(0, 30);
            h->Draw();
            c2.SaveAs(dcaFileName.c_str());
          }
        }

        if (false && histogramLayers) {
          histogramLayers->GetAxis(1)->SetRange(xbin, xbin);
          histogramLayers->GetAxis(2)->SetRange(ybin, ybin);
          TH1* layers = histogramLayers->Projection(0);
          if (layers) {
            layers->Draw();
            c2.SaveAs(dcaFileName.c_str());

            TH1F h("h", "Fired layers", 10, 0, 10);
            TH2F h2("h2", "Fired layers correlation", 10, 0, 10, 10, 0, 10);
            double integral = 0;
            for (int i = 1; i <= layers->GetXaxis()->GetNbins(); i++) {
              auto entries = layers->GetBinContent(i);
              integral += entries;
              auto firedLayers = GetFiredLayers(i - 1);
              for (int bit = 0; bit < 10; bit++) {
                if (firedLayers[bit]) {
                  h.Fill(bit, entries);
                  for (int bit2 = 0; bit2 < 10; bit2++) {
                    if (firedLayers[bit2]) {
                      h2.Fill(bit, bit2, entries);
                    }
                  }
                }
              }
            }
            if (integral > 0) {
              h.Scale(1.f / integral);
              h2.Scale(1.f / integral);
            }

            h.SetMinimum(0);
            h.SetMaximum(1);
            h.Draw("HIST");
            c2.SaveAs(dcaFileName.c_str());

            h2.SetMinimum(0);
            h2.SetMaximum(1);
            h2.Draw("colz");
            c2.SaveAs(dcaFileName.c_str());
          }
        }

        if (histogramSlope) {
          histogramSlope->GetAxis(1)->SetRange(xbin, xbin);
          histogramSlope->GetAxis(2)->SetRange(ybin, ybin);
          TH1* slope = histogramSlope->Projection(0);
          if (slope) {
            slope->Draw();
            c2.SaveAs(dcaFileName.c_str());
          }
        }

        c2.Clear();
        c2.SaveAs((dcaFileName+")").c_str());
      }
    }
  }

  c.cd();
  c.Clear();
  for (int layer = 0; layer <= layerBinMax; layer++) {
    if (layer ==1) {
      c.Clear();
      c.Divide(4, 2);
    }

    if (layer > 0) {
      c.cd(layer);
    }
    float min = -0.02;
    float max = 0.02;
    for (int ybin = 1; ybin <= histDCA[layer]->GetYaxis()->GetNbins(); ybin++) {
      for (int xbin = 1; xbin <= histDCA[layer]->GetXaxis()->GetNbins(); xbin++) {
        if (histDCA[layer]->GetBinContent(xbin, ybin) <= min)
          histDCA[layer]->SetBinContent(xbin, ybin, min + 1.e-6);
      }
    }
    histDCA[layer]->SetMinimum(min);
    histDCA[layer]->SetMaximum(max);
    histDCA[layer]->Draw("colz1");
    if (layer == 0 || layer == 8) {
      c.SaveAs(pdfFileName.c_str());
    }
  }

  c.cd();
  c.Clear();
  for (int layer = 0; layer <= layerBinMax; layer++) {
    if (layer ==1) {
      c.Clear();
      c.Divide(4, 2);
    }

    if (layer > 0) {
      c.cd(layer);
    }
    float min = -0.002;
    float max = 0.002;
    for (int ybin = 1; ybin <= histSlope[layer]->GetYaxis()->GetNbins(); ybin++) {
      for (int xbin = 1; xbin <= histSlope[layer]->GetXaxis()->GetNbins(); xbin++) {
        if (histSlope[layer]->GetBinContent(xbin, ybin) <= min)
          histSlope[layer]->SetBinContent(xbin, ybin, min + 1.e-6);
      }
    }
    histSlope[layer]->SetMinimum(min);
    histSlope[layer]->SetMaximum(max);
    histSlope[layer]->Draw("colz1");
    if (layer == 0 || layer == 8) {
      c.SaveAs(pdfFileName.c_str());
    }
  }


  c.cd();
  c.Clear();
  for (int layer = 0; layer <= layerBinMax; layer++) {
    if (layer ==1) {
      c.Clear();
      c.Divide(4, 2);
    }

    if (layer > 0) {
      c.cd(layer);
    }
    float min = 0;
    float max = 0.05;
    histSigma[layer]->SetMinimum(min);
    histSigma[layer]->SetMaximum(max);
    histSigma[layer]->Draw("colz1");
    if (layer == 0 || layer == 8) {
      c.SaveAs(pdfFileName.c_str());
    }
  }

  return true;
}


void ProcessDCAlayers(std::string coordinate, int layer, bool fired, TCanvas& c)
{
  std::string fullHistName =
      taskName + "/DCA/MFT/DCA_" + coordinate + "_fine" +
      (fired ? "_fl" : "_l") + std::to_string(layer);
  THnSparse* histogramFull = GetTHnSparse(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogramFull << std::endl;
  if (!histogramFull) return;

  // DCA axis assignments:
  // 0 -> DCA
  // 1 -> track x
  // 2 -> track y

  c.cd();
  TH2* histDCA = histogramFull->Projection(2, 1);
  if (!histDCA) return;
  histDCA->Reset();
  histDCA->SetTitle(histogramFull->GetTitle());

  for (int ybin = 1; ybin <= histogramFull->GetAxis(2)->GetNbins(); ybin++) {
    for (int xbin = 1; xbin <= histogramFull->GetAxis(1)->GetNbins(); xbin++) {
      histDCA->SetBinContent(xbin, ybin, -100);
      histDCA->SetBinError(xbin, ybin, 0);
    }
  }

  std::mutex rootMutex;
  std::deque<int> rowQueue;

  auto processRow = [&]()
  {
    while(true) {
      rootMutex.lock();
      if (rowQueue.empty()) {
        rootMutex.unlock();
        break;
      }
      int row = rowQueue.front();
      rowQueue.pop_front();
      rootMutex.unlock();

      rootMutex.lock();
      histogramFull->GetAxis(2)->SetRange(row, row);
      TH2* h2 = histogramFull->Projection(0, 1);
      h2->SetName(TString::Format("%s_%d", h2->GetName(), row));
      //std::cout << std::format("Starting row {}", row) << std::endl;
      rootMutex.unlock();

      std::vector<std::pair<double, double>> fitResults;
      for (int xbin = 1; xbin <= histogramFull->GetAxis(1)->GetNbins(); xbin++) {
        //if (xbin != 21) continue;
        if (xbin < 150 || xbin > 230) continue;
        //if (xbin < 14 || xbin > 17) continue;

        TH1* proj = (TH1*)h2->ProjectionY(TString::Format("%s_%d_%d", histogramFull->GetName(), xbin, row), xbin, xbin);
        //proj->SetName(TString::Format("%s_%d_%d", proj->GetName(), xbin, row));
        proj->GetXaxis()->SetRangeUser(-0.1f, 0.1f);

        auto fitResult = FitDCA(proj, c);

        fitResults.emplace_back(std::pair<double, double>{fitResult[0], fitResult[1]});

        /*rootMutex.lock();
        if (std::isnan(fitResult[0]) || std::isnan(fitResult[1])) {
          histDCA->SetBinContent(xbin, row, -100);
          histDCA->SetBinError(xbin, row, 0);
        } else {
          histDCA->SetBinContent(xbin, row, fitResult[0]);
          histDCA->SetBinError(xbin, row, fitResult[1]);
        }
        rootMutex.unlock();*/
      }

      rootMutex.lock();
      for (int xbin = 1; xbin <= histogramFull->GetAxis(1)->GetNbins(); xbin++) {
        float mean = fitResults[xbin - 1].first;
        float error = fitResults[xbin - 1].second;
        if (std::isnan(mean) || std::isnan(error)) {
          histDCA->SetBinContent(xbin, row, -100);
          histDCA->SetBinError(xbin, row, 0);
        } else {
          histDCA->SetBinContent(xbin, row, mean);
          histDCA->SetBinError(xbin, row, error);
        }
      }
      std::cout << std::format("Row {} done.", row) << std::endl;
      rootMutex.unlock();
    }
  };

  /*
  for (int ybin = 1; ybin <= histogramFull->GetAxis(2)->GetNbins(); ybin++) {
    //if (ybin != 11) continue;
    //if (ybin < 9 || ybin > 12) continue;
    if (ybin < 100 || ybin > 180) continue;
    rowQueue.push_back(ybin);
  }

  std::array<std::thread, 10> threads {
    std::thread{processRow}, std::thread{processRow}, std::thread{processRow}, std::thread{processRow}, std::thread{processRow},
    std::thread{processRow}, std::thread{processRow}, std::thread{processRow}, std::thread{processRow}, std::thread{processRow}
  };

  for (auto& t : threads) {
    t.join();
  }
  */

  /**/
  for (int ybin = 1; ybin <= histogramFull->GetAxis(2)->GetNbins(); ybin++) {
    //if (ybin < histogramFull->GetAxis(2)->GetNbins()/2) continue;
    //if (ybin < 9 || ybin > 12) continue;
    //if (ybin < 50 || ybin > 90) continue;
    histogramFull->GetAxis(2)->SetRange(ybin, ybin);
    for (int xbin = 1; xbin <= histogramFull->GetAxis(1)->GetNbins(); xbin++) {
      //if (xbin < histogramFull->GetAxis(1)->GetNbins()/2) continue;
      //if (xbin != 21) continue;
      //if (xbin < 25 || xbin > 125) continue;
      //if (xbin < 14 || xbin > 17) continue;
      histogramFull->GetAxis(1)->SetRange(xbin, xbin);

      TH1* proj = histogramFull->Projection(0);
      if (!proj) continue;
      proj->SetName(TString::Format("%s_%d_%d", proj->GetName(), xbin, ybin));
      proj->GetXaxis()->SetRangeUser(-0.1f, 0.1f);

      auto fitResult = FitDCA(proj, c);

      if (std::isnan(fitResult[0]) || std::isnan(fitResult[1])) {
        histDCA->SetBinContent(xbin, ybin, -100);
        histDCA->SetBinError(xbin, ybin, 0);
      } else {
        histDCA->SetBinContent(xbin, ybin, fitResult[0]);
        histDCA->SetBinError(xbin, ybin, fitResult[1]);
      }
    }
    std::cout << std::format("Row {} done.", ybin) << std::endl;
  }
  /**/

  histDCA->SetMinimum(-0.025);
  histDCA->SetMaximum(0.025);
  histDCA->Draw("colz");
  c.SaveAs(pdfFileName.c_str());
}


void ProcessDCA(std::string coordinate, TCanvas& c, TCanvas& c2)
{
  const int xyRebin = 2;

  std::string fullHistName = taskName + "/DCA/MFT/DCA_" + coordinate;
  THnSparse* histogramFull = GetTHnSparse(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogramFull << std::endl;
  if (!histogramFull) return;

  THnSparse* histogramMom = GetTHnSparse(fAnalysisResults, taskName + "/DCA/MFT/trackMomentum");
  THnSparse* histogramChi2 = GetTHnSparse(fAnalysisResults, taskName + "/DCA/MFT/trackChi2");
  THnSparse* histogramLayers = GetTHnSparse(fAnalysisResults, taskName + "/DCA/MFT/layers");

  // DCA axis assignments:
  // 0 -> DCA
  // 1 -> vertex z
  // 2 -> track x
  // 3 -> track y
  // 4 -> # of clusters

  // Track momentum axis assignments:
  // 0 -> momentum
  // 1 -> track x
  // 2 -> track y
  // 3 -> nClusters

  // Track chi2 axis assignments:
  // 0 -> chi2
  // 1 -> track x
  // 2 -> track y
  // 3 -> nClusters

  // Layers axis assignments:
  // 0 -> layers
  // 1 -> track x
  // 2 -> track y
  // 3 -> nClusters

  c.cd();
  TH2* histDCA = histogramFull->Projection(3, 2);
  if (!histDCA) return;

  for (int i = 0; i < 5; i++) {
    std::cout << std::format("Axis #{}: {}", i, histogramFull->GetAxis(i)->GetTitle()) << std::endl;
  }

  // only vertices around zero
  int vzBinMin = histogramFull->GetAxis(1)->FindBin(-1.f + 1.e-5);
  int vzBinMax = histogramFull->GetAxis(1)->FindBin(+1.f - 1.e-5);
  histogramFull->GetAxis(1)->SetRange(vzBinMin, vzBinMax);

  std::vector<std::pair<int, int>> _nclusRanges {
    {6, 10}, {5, 5}, {6, 6}, {7, 7}, {8, 8}, {9, 9}, {10, 10}
  };
  std::vector<std::pair<int, int>> nclusRanges {
    {6, 10}
  };

  //for (int nclusbin = 1; nclusbin <= histogramFull->GetAxis(4)->GetNbins(); nclusbin++) {
  for (auto [nclusMin, nclusMax] : nclusRanges) {
    histogramFull->GetAxis(4)->SetRange(nclusMin - 4, nclusMax - 4);
    if (histogramMom) histogramMom->GetAxis(3)->SetRange(nclusMin - 4, nclusMax - 4);
    if (histogramChi2) histogramChi2->GetAxis(3)->SetRange(nclusMin - 4, nclusMax - 4);
    if (histogramLayers) histogramLayers->GetAxis(3)->SetRange(nclusMin - 4, nclusMax - 4);

    histDCA->Reset();
    if (nclusMin == nclusMax) {
      histDCA->SetTitle(TString::Format("%s - %d clusters", histogramFull->GetTitle(), nclusMin));
    } else {
      histDCA->SetTitle(TString::Format("%s - %d-%d clusters", histogramFull->GetTitle(), nclusMin, nclusMax));
    }

    //for (int ybin = 1; ybin <= histogramFull->GetAxis(3)->GetNbins(); ybin++) {
    //  for (int xbin = 1; xbin <= histogramFull->GetAxis(2)->GetNbins(); xbin++) {
    //    histDCA->SetBinContent(xbin, ybin, -100);
    //    histDCA->SetBinError(xbin, ybin, 0);
    //  }
    //}

    TH2* hxy = histogramFull->Projection(3, 2);
    hxy->Draw("colz");
    c.SaveAs(pdfFileName.c_str());

    TH3* h3 = histogramFull->Projection(0, 2, 3);
    TH1* proj = (TH1*)h3->ProjectionX();
    std::cout << "3D projection:" << std::endl
        << "  X axis: " << h3->GetXaxis()->GetTitle() << std::endl
        << "  Y axis: " << h3->GetYaxis()->GetTitle() << std::endl
        << "  Z axis: " << h3->GetZaxis()->GetTitle() << std::endl;

    proj->Draw();
    c.SaveAs(pdfFileName.c_str());
    //continue;

    TH3* h3Mom = histogramMom ? histogramMom->Projection(0, 1, 2) : nullptr;
    TH3* h3Chi2 = histogramChi2 ? histogramChi2->Projection(0, 1, 2) : nullptr;
    TH3* h3Layers = histogramLayers ? histogramLayers->Projection(0, 1, 2) : nullptr;

    //for (int ybin = 1; ybin <= histogramFull->GetAxis(3)->GetNbins(); ybin++) {
    for (int ybin = 1; ybin <= h3->GetZaxis()->GetNbins(); ybin += xyRebin) {
      //if (ybin != histogramFull->GetAxis(3)->GetNbins()/2) continue;
      //if (ybin < h3->GetZaxis()->GetNbins()/2) continue;
      //if (ybin < 9 || ybin > 12) continue;
      //if (ybin < 50 || ybin > 90) continue;
      //histogramFull->GetAxis(3)->SetRange(ybin, ybin);
      //std::cout << std::format("Starting row {}", ybin) << std::endl;

      if (histogramMom) histogramMom->GetAxis(2)->SetRange(ybin, ybin);
      if (histogramChi2) histogramChi2->GetAxis(2)->SetRange(ybin, ybin);
      if (histogramLayers) histogramLayers->GetAxis(2)->SetRange(ybin, ybin);

      //for (int xbin = 1; xbin <= histogramFull->GetAxis(2)->GetNbins(); xbin++) {
      for (int xbin = 1; xbin <= h3->GetYaxis()->GetNbins(); xbin += xyRebin) {
        //if (xbin != (histogramFull->GetAxis(2)->GetNbins()/2 + 10)) continue;
        //if (xbin < h3->GetYaxis()->GetNbins()/2) continue;
        //if (xbin != 21) continue;
        //if (xbin < 25 || xbin > 125) continue;
        //if (xbin < 14 || xbin > 17) continue;
        //histogramFull->GetAxis(2)->SetRange(xbin, xbin);
        //std::cout << std::format("  column {}", xbin) << std::endl;

        if (histogramMom) histogramMom->GetAxis(1)->SetRange(xbin, xbin);
        if (histogramChi2) histogramChi2->GetAxis(1)->SetRange(xbin, xbin);
        if (histogramLayers) histogramLayers->GetAxis(1)->SetRange(xbin, xbin);

        //TH1* proj = histogramFull->Projection(0);
        TH1* proj = (TH1*)h3->ProjectionX(TString::Format("%s_%d_%d_%d_%d", h3->GetName(), nclusMin, nclusMax, xbin, ybin),
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
        if (std::isnan(fitResult[0]) || std::isnan(fitResult[1])) continue;
        if (!histogramMom && !histogramChi2 && !histogramLayers) continue;

        std::string dcaFileName = std::format("mft_dca{}-clus_{}_{}-{}_{}.pdf", coordinate, nclusMin, nclusMax, xbin, ybin);
        c2.cd();
        proj->Draw();
        c2.SaveAs((dcaFileName + "(").c_str());
        proj->GetXaxis()->SetRangeUser(-0.1, 0.1);
        proj->Draw();
        c2.SaveAs(dcaFileName.c_str());

        if (h3Mom) {
          TH1* proj2 = (TH1*)h3Mom->ProjectionX(TString::Format("%s_%d_%d_%d_%d", h3Mom->GetName(), nclusMin, nclusMax, xbin, ybin),
              xbin, xbin + xyRebin - 1, ybin, ybin + xyRebin - 1);
          //TH1* proj2 = histogramMom->Projection(0);
          c2.SetLogy(kTRUE);
          proj2->GetXaxis()->SetRangeUser(0, 50);
          proj2->Draw();
          c2.SaveAs(dcaFileName.c_str());
          c2.SetLogy(kFALSE);
          delete proj2;
        }

        if (h3Chi2) {
          TH1* proj2 = (TH1*)h3Chi2->ProjectionX(TString::Format("%s_%d_%d_%d_%d", h3Mom->GetName(), nclusMin, nclusMax, xbin, ybin),
              xbin, xbin + xyRebin - 1, ybin, ybin + xyRebin - 1);
          //TH1* proj2 = histogramChi2->Projection(0);
          c2.SetLogy(kTRUE);
          proj2->GetXaxis()->SetRangeUser(0, 100);
          proj2->Draw();
          c2.SaveAs(dcaFileName.c_str());
          c2.SetLogy(kFALSE);
          delete proj2;
        }

        if (h3Layers) {
          TH1* proj2 = (TH1*)h3Layers->ProjectionX(TString::Format("%s_%d_%d_%d_%d", h3Mom->GetName(), nclusMin, nclusMax, xbin, ybin),
              xbin, xbin + xyRebin - 1, ybin, ybin + xyRebin - 1);
          //TH1* proj2 = histogramLayers->Projection(0);
          proj2->Draw();
          c2.SaveAs(dcaFileName.c_str());
          delete proj2;
        }

        c2.Clear();
        c2.SaveAs((dcaFileName + ")").c_str());
      }
      //std::cout << std::format("Row {} done.", ybin) << std::endl;
    }

    c.cd();
    histDCA->SetMinimum(-0.025);
    histDCA->SetMaximum(0.025);
    histDCA->Draw("colz");
    c.SaveAs(pdfFileName.c_str());

    c2.cd();
    c2.Clear();
    histDCA->Draw("colz");
    c2.SaveAs(std::format("DCA{}.root", coordinate).c_str());

    TGraph grVertical;
    // vertical strip
    for (int ybin = 1; ybin <= h3->GetZaxis()->GetNbins(); ybin++) {
      double total = 0;
      int counts = 0;
      int xbinMin = h3->GetYaxis()->FindBin(-1);
      int xbinMax = h3->GetYaxis()->FindBin(+1);
      //std::cout << std::format("y={}", ybin) << std::endl;
      for (int xbin = xbinMin; xbin <= xbinMax; xbin++) {
        double value = histDCA->GetBinContent(xbin, ybin);
        if (value == -100) continue;
        total += value;
        counts += 1;
        //std::cout << std::format("  x={} value={:0.3f}", xbin, value) << std::endl;
      }
      if (counts < 1) continue;
      //std::cout << std::format("  average={:0.3f}", total / counts) << std::endl;

      grVertical.AddPoint(h3->GetZaxis()->GetBinCenter(ybin), total / counts);
    }
    c.cd();
    grVertical.Draw("AL*");
    grVertical.SetMinimum(-0.025);
    grVertical.SetMaximum(0.025);
    c.SaveAs(pdfFileName.c_str());

    // top/bottom averages
    double total[2]{ 0, 0 };
    int counts[2]{ 0, 0 };
    for (int ybin = 1; ybin <= h3->GetZaxis()->GetNbins(); ybin++) {
      for (int xbin = 1; xbin <= h3->GetYaxis()->GetNbins(); xbin++) {
        double value = histDCA->GetBinContent(xbin, ybin);
        if (value == -100) continue;

        int tb = 0;
        if (ybin < (h3->GetZaxis()->GetNbins() / 2)) {
          tb = 0;
        } else {
          tb = 1;
        }
        total[tb] += value;
        counts[tb] += 1;
      }
    }
    TH1F* htb = new TH1F("htb", TString::Format("DCA (%s) top/bottom", coordinate.c_str()), 2, 0, 2);
    htb->GetXaxis()->SetBinLabel(1, "B");
    htb->GetXaxis()->SetBinLabel(2, "T");
    htb->SetBinContent(1, (counts[0] > 0 ? total[0] / counts[0] : 0));
    htb->SetBinContent(2, (counts[1] > 0 ? total[1] / counts[1] : 0));

    std::cout << std::format("DCA ({}) top={:0.3f} bottom={:0.3f}", coordinate.c_str(),
        (counts[1] > 0 ? total[1] / counts[1] : 0), (counts[0] > 0 ? total[0] / counts[0] : 0)) << std::endl;

    c.cd();
    htb->Draw("hist");
    c.SaveAs(pdfFileName.c_str());
  }
}


void ProcessDCAGlobalFwd(std::string coordinate, TCanvas& c, TCanvas& c2)
{
  const int xyRebin = 2;

  std::string fullHistName = taskName + "/DCA/GlobalFwd/DCA_" + coordinate;
  THnSparse* histogramFull = GetTHnSparse(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogramFull << std::endl;
  if (!histogramFull) return;

  fullHistName = taskName + "/DCA/GlobalFwd/DCAMFT_" + coordinate;
  THnSparse* histogramMFT = GetTHnSparse(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogramMFT << std::endl;

  THnSparse* histogramDeltaP = GetTHnSparse(fAnalysisResults, taskName + "/DCA/GlobalFwd/DeltaP");
  THnSparse* histogramDeltaQ = GetTHnSparse(fAnalysisResults, taskName + "/DCA/GlobalFwd/DeltaQ");

  double topOffset=0.0, bottomOffset=0.0;
  //double topOffset=0.003, bottomOffset=-0.007;

  // DCA axis assignments:
  // 0 -> DCA
  // 1 -> vertex z
  // 2 -> track x
  // 3 -> track y
  // 4 -> # of clusters

  // Delta p/Q axis assignments:
  // 0 -> delta p/Q
  // 1 -> momentum
  // 2 -> track x
  // 3 -> track y
  // 4 -> # of clusters

  c.cd();
  TH2* histDCA = histogramFull->Projection(3, 2);
  if (!histDCA) return;

  TH2* histDCA2 = histogramMFT ? histogramMFT->Projection(3, 2) : nullptr;

  for (int i = 0; i < 5; i++) {
    std::cout << std::format("Axis #{}: {}", i, histogramFull->GetAxis(i)->GetTitle()) << std::endl;
  }

  // only vertices around zero
  int vzBinMin = histogramFull->GetAxis(1)->FindBin(-2.f + 1.e-5);
  int vzBinMax = histogramFull->GetAxis(1)->FindBin(+2.f - 1.e-5);
  histogramFull->GetAxis(1)->SetRange(vzBinMin, vzBinMax);
  if (histogramMFT) histogramMFT->GetAxis(1)->SetRange(vzBinMin, vzBinMax);

  std::vector<std::pair<int, int>> _nclusRanges {
    {6, 10}, {5, 5}, {6, 6}, {7, 7}, {8, 8}, {9, 9}, {10, 10}
  };
  std::vector<std::pair<int, int>> nclusRanges {
    {6, 10}
  };

  //for (int nclusbin = 1; nclusbin <= histogramFull->GetAxis(4)->GetNbins(); nclusbin++) {
  for (auto [nclusMin, nclusMax] : nclusRanges) {
    histogramFull->GetAxis(4)->SetRange(nclusMin - 4, nclusMax - 4);
    if (histogramMFT) histogramMFT->GetAxis(4)->SetRange(nclusMin - 4, nclusMax - 4);
    if (histogramDeltaP) histogramDeltaP->GetAxis(4)->SetRange(nclusMin - 4, nclusMax - 4);
    if (histogramDeltaQ) histogramDeltaQ->GetAxis(4)->SetRange(nclusMin - 4, nclusMax - 4);

    histDCA->Reset();
    if (nclusMin == nclusMax) {
      histDCA->SetTitle(TString::Format("%s - %d clusters", histogramFull->GetTitle(), nclusMin));
    } else {
      histDCA->SetTitle(TString::Format("%s - %d-%d clusters", histogramFull->GetTitle(), nclusMin, nclusMax));
    }
    if (histDCA2) {
      histDCA2->Reset();
      if (nclusMin == nclusMax) {
        histDCA2->SetTitle(TString::Format("%s - %d clusters", histogramMFT->GetTitle(), nclusMin));
      } else {
        histDCA2->SetTitle(TString::Format("%s - %d-%d clusters", histogramMFT->GetTitle(), nclusMin, nclusMax));
      }
    }

    //for (int ybin = 1; ybin <= histogramFull->GetAxis(3)->GetNbins(); ybin++) {
    //  for (int xbin = 1; xbin <= histogramFull->GetAxis(2)->GetNbins(); xbin++) {
    //    histDCA->SetBinContent(xbin, ybin, -100);
    //    histDCA->SetBinError(xbin, ybin, 0);
    //  }
    //}

    TH2* hxy = histogramFull->Projection(3, 2);
    hxy->Draw("colz");
    c.SaveAs(pdfFileName.c_str());

    TH3* h3 = histogramFull->Projection(0, 2, 3);
    TH1* proj = (TH1*)h3->ProjectionX();
    std::cout << "3D projection:" << std::endl
        << "  X axis: " << h3->GetXaxis()->GetTitle() << std::endl
        << "  Y axis: " << h3->GetYaxis()->GetTitle() << std::endl
        << "  Z axis: " << h3->GetZaxis()->GetTitle() << std::endl;

    proj->Draw();
    c.SaveAs(pdfFileName.c_str());
    //continue;

    TH3* h3MFT = nullptr;
    if (histogramMFT) {
      h3MFT = histogramMFT->Projection(0, 2, 3);
      proj = (TH1*)h3MFT->ProjectionX();
      proj->Draw();
      c.SaveAs(pdfFileName.c_str());
    }

    //for (int ybin = 1; ybin <= histogramFull->GetAxis(3)->GetNbins(); ybin++) {
    for (int ybin = 1; ybin <= h3->GetZaxis()->GetNbins(); ybin += xyRebin) {
      //if (ybin != histogramFull->GetAxis(3)->GetNbins()/2) continue;
      //if (ybin < h3->GetZaxis()->GetNbins()/2) continue;
      //if (ybin < 9 || ybin > 12) continue;
      //if (ybin < 50 || ybin > 90) continue;
      //histogramFull->GetAxis(3)->SetRange(ybin, ybin);
      std::cout << std::format("Starting row {}", ybin) << std::endl;

      if (histogramDeltaP) histogramDeltaP->GetAxis(3)->SetRange(ybin, ybin + xyRebin - 1);
      if (histogramDeltaQ) histogramDeltaQ->GetAxis(3)->SetRange(ybin, ybin + xyRebin - 1);

      //for (int xbin = 1; xbin <= histogramFull->GetAxis(2)->GetNbins(); xbin++) {
      for (int xbin = 1; xbin <= h3->GetYaxis()->GetNbins(); xbin += xyRebin) {
        //if (xbin != (histogramFull->GetAxis(2)->GetNbins()/2 + 10)) continue;
        //if (xbin < h3->GetYaxis()->GetNbins()/2) continue;
        //if (xbin != 21) continue;
        //if (xbin < 25 || xbin > 125) continue;
        //if (xbin < 14 || xbin > 17) continue;
        //histogramFull->GetAxis(2)->SetRange(xbin, xbin);
        //std::cout << std::format("  column {}", xbin) << std::endl;

        if (histogramDeltaP) histogramDeltaP->GetAxis(2)->SetRange(xbin, xbin + xyRebin - 1);
        if (histogramDeltaQ) histogramDeltaQ->GetAxis(2)->SetRange(xbin, xbin + xyRebin - 1);

        //TH1* proj = histogramFull->Projection(0);
        TH1* proj = (TH1*)h3->ProjectionX(TString::Format("%s_%d_%d_%d_%d", h3->GetName(), nclusMin, nclusMax, xbin, ybin),
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
              if (ybin < (h3->GetZaxis()->GetNbins() / 2)) {
                histDCA->SetBinContent(xbin2, ybin2, fitResult[0] - bottomOffset);
              } else {
                histDCA->SetBinContent(xbin2, ybin2, fitResult[0] - topOffset);
              }
              histDCA->SetBinError(xbin2, ybin2, fitResult[1]);
            }
          }
        }

        TH1* proj2 = nullptr;
        if (h3MFT && histDCA2) {
          proj2 = (TH1*)h3MFT->ProjectionX(TString::Format("%s_%d_%d_%d_%d", h3MFT->GetName(), nclusMin, nclusMax, xbin, ybin),
              xbin, xbin + xyRebin - 1, ybin, ybin + xyRebin - 1);
          if (proj2) {
            auto fitResult = FitDCA(proj2, c);

            for (int ybin2 = ybin; ybin2 < ybin + xyRebin; ybin2++) {
              for (int xbin2 = xbin; xbin2 < xbin + xyRebin; xbin2++) {
                if (std::isnan(fitResult[0]) || std::isnan(fitResult[1])) {
                  histDCA2->SetBinContent(xbin2, ybin2, -100);
                  histDCA2->SetBinError(xbin2, ybin2, 0);
                } else {
                  if (ybin < (h3->GetZaxis()->GetNbins() / 2)) {
                    histDCA2->SetBinContent(xbin2, ybin2, fitResult[0] - bottomOffset);
                  } else {
                    histDCA2->SetBinContent(xbin2, ybin2, fitResult[0] - topOffset);
                  }
                  histDCA2->SetBinError(xbin2, ybin2, fitResult[1]);
                }
              }
            }
          }
        }

        if (true && proj->GetEntries() > 0 && !std::isnan(fitResult[0]) && !std::isnan(fitResult[1])) {
          std::string dcaFileName = std::format("mft_dca{}-clus_{}_{}-{}_{}.pdf",
              coordinate, nclusMin, nclusMax, xbin, ybin);
          c2.cd();
          proj->Draw();
          c2.SaveAs((dcaFileName + "(").c_str());
          proj->GetXaxis()->SetRangeUser(-0.1, 0.1);
          proj->Draw();
          c2.SaveAs(dcaFileName.c_str());
          if (proj2) {
            proj2->Draw();
            c2.SaveAs(dcaFileName.c_str());

            proj2->SetLineColor(kRed);
            proj->Draw();
            proj2->Draw("hist same");
            c2.SaveAs(dcaFileName.c_str());
          }

          std::cout << std::format("DCA stat: {} + {} (UF) + {} (OF)",
              proj->Integral(), proj->GetBinContent(0), proj->GetBinContent(proj->GetXaxis()->GetNbins()+1))
                    << std::endl;

          //if (std::isnan(fitResult[0]) || std::isnan(fitResult[1])) continue;
          //if (!histogramDeltaP && !histogramDeltaQ) continue;
          //if (histogramDeltaP && histogramDeltaP->GetEntries() < 1) continue;
          //if (histogramDeltaQ && histogramDeltaQ->GetEntries() < 1) continue;

          if (histogramDeltaP) {
            TH2* proj2 = histogramDeltaP->Projection(0, 1);
            proj2->SetName(TString::Format("%s_%d_%d_%d_%d", histogramDeltaP->GetName(), nclusMin, nclusMax, xbin, ybin));
            proj2->Draw("colz");
            c2.SaveAs(dcaFileName.c_str());
            std::cout << std::format("DeltaP stat: {}", proj2->Integral())
                      << std::endl;
            delete proj2;
          }

          if (histogramDeltaQ) {
            TH2* proj2 = histogramDeltaQ->Projection(0, 1);
            proj2->SetName(TString::Format("%s_%d_%d_%d_%d", histogramDeltaP->GetName(), nclusMin, nclusMax, xbin, ybin));
            proj2->Draw("colz");
            c2.SaveAs(dcaFileName.c_str());
            delete proj2;
          }

          c2.Clear();
          c2.SaveAs((dcaFileName + ")").c_str());
        }
      }
      std::cout << std::format("Row {} done.", ybin) << std::endl;
    }

    c.cd();
    histDCA->SetMinimum(-0.025);
    histDCA->SetMaximum(0.025);
    histDCA->Draw("colz");
    c.SaveAs(pdfFileName.c_str());

    c2.cd();
    c2.Clear();
    histDCA->Draw("colz");
    c2.SaveAs(std::format("DCA{}.root", coordinate).c_str());

    TGraph grVertical;
    // vertical strip
    for (int ybin = 1; ybin <= h3->GetZaxis()->GetNbins(); ybin++) {
      double total = 0;
      int counts = 0;
      int xbinMin = h3->GetYaxis()->FindBin(-1);
      int xbinMax = h3->GetYaxis()->FindBin(+1);
      std::cout << std::format("y={}", ybin) << std::endl;
      for (int xbin = xbinMin; xbin <= xbinMax; xbin++) {
        double value = histDCA->GetBinContent(xbin, ybin);
        if (value == -100) continue;
        total += value;
        counts += 1;
        std::cout << std::format("  x={} value={:0.3f}", xbin, value) << std::endl;
      }
      if (counts < 1) continue;
      std::cout << std::format("  average={:0.3f}", total / counts) << std::endl;

      grVertical.AddPoint(h3->GetZaxis()->GetBinCenter(ybin), total / counts);
    }
    c.cd();
    grVertical.Draw("AL*");
    grVertical.SetMinimum(-0.025);
    grVertical.SetMaximum(0.025);
    c.SaveAs(pdfFileName.c_str());

    // top/bottom averages
    double total[2]{ 0, 0 };
    int counts[2]{ 0, 0 };
    for (int ybin = 1; ybin <= h3->GetZaxis()->GetNbins(); ybin++) {
      for (int xbin = 1; xbin <= h3->GetYaxis()->GetNbins(); xbin++) {
        double value = histDCA->GetBinContent(xbin, ybin);
        if (value == -100) continue;

        int tb = 0;
        if (ybin < (h3->GetZaxis()->GetNbins() / 2)) {
          tb = 0;
        } else {
          tb = 1;
        }
        total[tb] += value;
        counts[tb] += 1;
      }
    }
    TH1F* htb = new TH1F("htb", TString::Format("DCA (%s) top/bottom", coordinate.c_str()), 2, 0, 2);
    htb->GetXaxis()->SetBinLabel(1, "B");
    htb->GetXaxis()->SetBinLabel(2, "T");
    htb->SetBinContent(1, (counts[0] > 0 ? total[0] / counts[0] : 0));
    htb->SetBinContent(2, (counts[1] > 0 ? total[1] / counts[1] : 0));

    std::cout << std::format("DCA ({}) top={:0.3f} bottom={:0.3f}", coordinate.c_str(),
        (counts[1] > 0 ? total[1] / counts[1] : 0), (counts[0] > 0 ? total[0] / counts[0] : 0)) << std::endl;

    c.cd();
    htb->Draw("hist");
    c.SaveAs(pdfFileName.c_str());

    if (histDCA2) {
      c.cd();
      histDCA2->SetMinimum(-0.025);
      histDCA2->SetMaximum(0.025);
      histDCA2->Draw("colz");
      c.SaveAs(pdfFileName.c_str());

      c2.cd();
      c2.Clear();
      histDCA2->Draw("colz");
      c2.SaveAs(std::format("DCAMFT{}.root", coordinate).c_str());
    }
  }
}

void muonGlobalAlignmentMftDCAoffs(const char* _rootFileName = "AnalysisResults.root", const char* _pdfFileName = "mftDCAoffs.pdf")
{
  //fAnalysisResults = new TFile("AnalysisResults.root");
  //fAnalysisResults = new TFile("AnalysisResults/AnalysisResultsFull.root");
  fAnalysisResults = new TFile(_rootFileName);
  pdfFileName = _pdfFileName;

  gStyle->SetOptStat(0);
  //gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(99);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  
  TCanvas c2("c2", "c2", 1200, 800);

  TCanvas c("c", "c", 1200, 800);
  c.SetGrid();
  c.SaveAs((pdfFileName + "(").c_str());

  std::string fullHistName = taskName + "/vertex_y_vs_x";
  auto h2 = GetTH2(fAnalysisResults, fullHistName);
  h2->Draw("col");
  c.SaveAs(pdfFileName.c_str());

  fullHistName = taskName + "/vertex_z";
  auto h1 = GetTH1(fAnalysisResults, fullHistName);
  h1->Draw("HIST");
  c.SaveAs(pdfFileName.c_str());

  TPaveText* title = new TPaveText(0.1, 0.4, 0.9, 0.6, "NDC");

  std::array<std::string, 2> coordinates{"x", "y"};
  //std::array<std::string, 2> coordinates{"x"};
  for (auto& coordinate : coordinates) {
    c.Clear();
    title->Clear();
    title->AddText(TString::Format("DCA(%s) offset vs. track XY", coordinate.c_str()));
    title->AddText("Integrated over vertex z");
    title->Draw();
    c.SaveAs(pdfFileName.c_str());
    ProcessDCA(coordinate, c, c2);
    ProcessDCAGlobalFwd(coordinate, c, c2);
    for (int l = 0; l < 1; l++) {
      //ProcessDCA(coordinate, l, false, c);
      //ProcessDCA(coordinate, l, true, c);
    }
    //break;
  }
  //ProcessDCA("x", c);
  //ProcessDCA("y", 9, false, c);

  c.Clear();
  c.SaveAs((pdfFileName + ")").c_str());
}
