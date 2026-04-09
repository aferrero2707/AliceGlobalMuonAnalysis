TFile* fAnalysisResults;
std::string pdfFileName;


TH1* GetTH1(TFile* f, TString histname)
{
  return (TH1*)f->Get(histname);
}


TH2* GetTH2(TFile* f, TString histname)
{
  return (TH2*)f->Get(histname);
}


void PlotMftLayerEfficiencies(TCanvas& c)
{
  if (!GetTH2(fAnalysisResults, "muon-global-alignment/DCA/MFT/mftTrackEffNum_0")) return;

  TPaveText* title = new TPaveText(0.1, 0.4, 0.9, 0.6, "NDC");
  c.Clear();
  title->AddText("MFT layer efficiencies");
  title->Draw();
  c.SaveAs(pdfFileName.c_str());

  c.Clear();
  for (int i = 0; i < 10; i++) {
    std::string fullHistName = std::string("muon-global-alignment/DCA/MFT/mftTrackEffNum_") + std::to_string(i);
    TH2* histogramNum = GetTH2(fAnalysisResults, fullHistName);
    std::cout << fullHistName << " -> " << histogramNum << std::endl;
    if (!histogramNum) continue;

    fullHistName = std::string("muon-global-alignment/DCA/MFT/mftTrackEffDen_") + std::to_string(i);
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

void muonGlobalAlignmentMftPseudoeff(const char* _rootFileName = "AnalysisResults.root", const char* _pdfFileName = "mftPseudoeff.pdf")
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
  
  TCanvas c("c", "c", 1200, 800);
  c.SaveAs((pdfFileName + "(").c_str());

  PlotMftLayerEfficiencies(c);

  c.Clear();
  c.SaveAs((pdfFileName + ")").c_str());
}
