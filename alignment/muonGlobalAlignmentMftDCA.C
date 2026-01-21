#include <MFTTracking/Constants.h>

TFile* fAnalysisResults;


TH2* GetTH2(TFile* f, TString histname)
{
  return (TH2*)f->Get(histname);
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

double DoubleSidedCBwithQuadBgd(double* x, double *par)
{
  return(par[0] * DoubleSidedCB2(x[0], par[1],par[2],par[3],par[4],par[5],par[6]) + par[7] + x[0] * par[8] + x[0] * x[0] * par[9]);
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
      std::cout << std::format("Fit status: {}", fitResult->Status()) << std::endl;
      //fitResult->Print();

      if (false && printFits) {
        proj->Draw("E");
        c.SaveAs("mft_dca_AO2D.pdf");
      }

      if (fitResult->Status() == 0) {
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
  //c.SaveAs("mft_dca_AO2D.pdf");
  //histogramSigma->Draw("E");
  //c.SaveAs("mft_dca_AO2D.pdf");

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

bool ProcessDCAvsZ(std::string histName, std::string coordinate, TCanvas& c, TCanvas& c2, bool printFits = false)
{
  std::string fullHistName = std::string("muon-global-alignment/") + histName;
  THnSparse* histogramFull = GetTHnSparse(fAnalysisResults, fullHistName);
  std::cout << fullHistName << " -> " << histogramFull << std::endl;
  if (!histogramFull) return false;

  // Axis assignments:
  // 0 -> DCA
  // 1 -> vertex z
  // 2 -> track x
  // 3 -> track y
  // 4 -> nClusters
  // 5 -> MFT track type


  c.cd();
  TH2* trackxy = histogramFull->Projection(3, 2);
  if (!trackxy) return false;
  trackxy->Draw("col");
  c.SaveAs("mft_dca_AO2D.pdf");

  for (int nclBin = 1; nclBin <= histogramFull->GetAxis(4)->GetNbins(); nclBin++) {
    histogramFull->GetAxis(4)->SetRange(nclBin, nclBin);
    for (int ttypeBin = 1; ttypeBin <= histogramFull->GetAxis(5)->GetNbins(); ttypeBin++) {
      histogramFull->GetAxis(5)->SetRange(ttypeBin, ttypeBin);
      trackxy = histogramFull->Projection(3, 2);
      if (!trackxy) continue;
      trackxy->SetTitle(std::format("Track (x,y) with nClus={} and trackType={}", nclBin + 4, (ttypeBin == 1 ? "LF" : "CA")).c_str());
      trackxy->Draw("col");
      c.SaveAs("mft_dca_AO2D.pdf");
    }
  }
  histogramFull->GetAxis(4)->SetRange(1, histogramFull->GetAxis(4)->GetNbins());
  histogramFull->GetAxis(5)->SetRange(1, histogramFull->GetAxis(5)->GetNbins());

  trackxy->Reset();
  trackxy->SetTitle(std::format("DCA({}) offset vs. MFT track position", coordinate).c_str());

  TH2* tracksxsy = (TH2*)trackxy->Clone();
  tracksxsy->SetTitle(std::format("DCA({}) slope vs. MFT track position", coordinate).c_str());

  TH2* dcaSigma = (TH2*)trackxy->Clone();
  dcaSigma->SetTitle(std::format("#sigma_{{DCA({})}} vs. MFT track position", coordinate).c_str());

  c2.cd();
  for (int ybin = 1; ybin <= histogramFull->GetAxis(3)->GetNbins(); ybin++) {
    histogramFull->GetAxis(3)->SetRange(ybin, ybin);
    for (int xbin = 1; xbin <= histogramFull->GetAxis(2)->GetNbins(); xbin++) {
      histogramFull->GetAxis(2)->SetRange(xbin, xbin);
      TH2* dcaVsZ = histogramFull->Projection(0, 1);
      if (!dcaVsZ) continue;
      dcaVsZ->SetName(TString::Format("%s_%d_%d", dcaVsZ->GetName(), xbin, ybin));
      //if (dcaVsZ->GetEntries() < 1000) continue;

      dcaVsZ->GetYaxis()->SetRangeUser(-0.1, 0.1);
      dcaVsZ->Draw("col");
      c2.SaveAs(std::format("mft_dca{}_vs_z_{}_{}.pdf(", coordinate, xbin, ybin).c_str());

      auto linfit = PlotDCAProjection(dcaVsZ, -0.1, 0.1, 1, c2, printFits);
      c2.SaveAs(std::format("mft_dca{}_vs_z_{}_{}.pdf", coordinate, xbin, ybin).c_str());

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
          dcaSigma->SetBinContent(xbin, ybin, fcb.GetParameter(2));
          dcaSigma->SetBinError(xbin, ybin, 0.1);
          std::cout << std::format("peak: {:0.3f} +/- {:0.3f}", fcb.GetParameter(1), fcb.GetParameter(2)) << std::endl;
        }
      }
      c2.SaveAs(std::format("mft_dca{}_vs_z_{}_{}.pdf", coordinate, xbin, ybin).c_str());

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
        trackxy->SetBinContent(xbin, ybin, std::get<0>(linfit));
        trackxy->SetBinError(xbin, ybin, 0.1);

        tracksxsy->SetBinContent(xbin, ybin, std::get<2>(linfit));
        tracksxsy->SetBinError(xbin, ybin, 0.1);
      }
      c2.SaveAs(std::format("mft_dca{}_vs_z_{}_{}.pdf)", coordinate, xbin, ybin).c_str());
    }
  }
  c.cd();

  float min = -0.02;
  float max = 0.02;
  for (int ybin = 1; ybin <= trackxy->GetYaxis()->GetNbins(); ybin++) {
    for (int xbin = 1; xbin <= trackxy->GetXaxis()->GetNbins(); xbin++) {
      if (trackxy->GetBinContent(xbin, ybin) <= min)
        trackxy->SetBinContent(xbin, ybin, min + 1.e-6);
    }
  }
  trackxy->SetMinimum(min);
  trackxy->SetMaximum(max);
  trackxy->Draw("colz1");
  c.SaveAs("mft_dca_AO2D.pdf");

  min = -0.002;
  max = 0.002;
  for (int ybin = 1; ybin <= tracksxsy->GetYaxis()->GetNbins(); ybin++) {
    for (int xbin = 1; xbin <= tracksxsy->GetXaxis()->GetNbins(); xbin++) {
      if (tracksxsy->GetBinContent(xbin, ybin) <= min)
        tracksxsy->SetBinContent(xbin, ybin, min + 1.e-6);
    }
  }
  tracksxsy->SetMinimum(min);
  tracksxsy->SetMaximum(max);
  tracksxsy->Draw("colz1");
  c.SaveAs("mft_dca_AO2D.pdf");

  min = 0;
  max = 0.05;
  /*for (int ybin = 1; ybin <= dcaSigma->GetYaxis()->GetNbins(); ybin++) {
    for (int xbin = 1; xbin <= dcaSigma->GetXaxis()->GetNbins(); xbin++) {
      if (dcaSigma->GetBinContent(xbin, ybin) >= min)
        dcaSigma->SetBinContent(xbin, ybin, min + 1.e-6);
    }
  }*/
  dcaSigma->SetMinimum(min);
  dcaSigma->SetMaximum(max);
  dcaSigma->Draw("colz1");
  c.SaveAs("mft_dca_AO2D.pdf");

  return true;
}

void muonGlobalAlignmentMftDCA()
{
  //fAnalysisResults = new TFile("AnalysisResults.root");
  fAnalysisResults = new TFile("AnalysisResults/AnalysisResultsFull.root");

  gStyle->SetOptStat(0);
  //gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  
  TCanvas c("c", "c", 1200, 800);
  c.SaveAs("mft_dca_AO2D.pdf(");

  TCanvas c2("c2", "c2", 1200, 800);

  c.cd();
  ProcessDCAvsZ("DCA/MFT/DCA_x", "x", c, c2, false);
  ProcessDCAvsZ("DCA/MFT/DCA_y", "y", c, c2, false);

  c.Clear();
  c.SaveAs("mft_dca_AO2D.pdf)");
}
