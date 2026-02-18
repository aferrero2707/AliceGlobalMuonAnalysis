#include <MFTTracking/Constants.h>

TCanvas c("c", "c", 1200, 800);
c.SetGridx(true);
c.SetGridy(true);

constexpr int nPoints = 3;

struct matchType
{
  std::vector<int> bins;
  std::string name;
};

std::vector<matchType> matchTypesVec = {
    {{1, 5}, "True matches"},
    {{2, 6}, "Wrong matches"},
    {{3, 7}, "Decay matches"},
    {{4, 8}, "Fake matches"}};

TH1* GetTH1(TFile* f, TString histname)
{
  TH1* result = (TH1*)f->Get(histname);
  //std::cout << "histname: " << histname << " -> " << result << std::endl;
  return result;
}

TH2* GetTH2(TFile* f, TString histname)
{
  TH2* result = (TH2*)f->Get(histname);
  //std::cout << "histname: " << histname << " -> " << result << std::endl;
  return result;
}

TH3* GetTH3(TFile* f, TString histname)
{
  TH3* result = (TH3*)f->Get(histname);
  //std::cout << "histname: " << histname << " -> " << result << std::endl;
  return result;
}

//----------------------------------------------------------------------

void PlotMatchEfficiancyAndPurity(TFile* rootFile)
{
  std::array<std::string, 3> matchingMethods{
    "Prod",
    "MatchXYPhiTanl",
    "MatchXYPhiTanlMom"};
  std::array<std::string, 4> variables{"p", "pt", "eta", "phi"};

  std::array<std::array<TH1*, variables.size()>, matchingMethods.size()> effPlots, purityPlots;

  TH1* hnum;
  TH1* hden;
  TPaveText* title = new TPaveText(0.1, 0.4, 0.9, 0.6, "NDC");
  TLegend* legend = new TLegend(0.15, 0.8, 0.85, 0.9);
  legend->SetNColumns(matchingMethods.size());

  //--
  // Matching efficiency
  //--
  c.Clear();
  title->Clear();
  title->AddText("Matching efficiency");
  title->Draw();
  c.SaveAs("matchingQA.pdf");

  for (int vi = 0; vi < variables.size(); vi++) {
    auto variable = variables[vi];

    for (int mi = 0, colorIndex = 1; mi < matchingMethods.size(); mi++, colorIndex++) {
      auto method = matchingMethods[mi];
      std::string path = std::string("qa-matching/matching/MC/") + method + "/matching-efficiency/";

      std::string histname = path + variable + "_num";
      hnum = GetTH1(rootFile, histname.c_str());
      if (!hnum) {
        std::cout << "Histogram \"" << histname << "\" not found" << std::endl;
        continue;
      }
      hnum->SetLineColor(colorIndex);

      histname = path + variable + "_den";
      hden = GetTH1(rootFile, histname.c_str());
      if (!hden) {
        std::cout << "Histogram \"" << histname << "\" not found" << std::endl;
        continue;
      }

      hnum->Divide(hden);
      hnum->SetTitle(Form("Matching efficiency vs. %s", variable.c_str()));

      effPlots[mi][vi] = hnum;
    }
  }

  for (int vi = 0; vi < variables.size(); vi++) {
    legend->Clear();
    for (int mi = 0; mi < matchingMethods.size(); mi++) {
      if (mi == 0) {
        effPlots[mi][vi]->SetMinimum(0);
        effPlots[mi][vi]->SetMaximum(1.2);
        effPlots[mi][vi]->Draw();
      } else {
        effPlots[mi][vi]->Draw("same");
      }
      legend->AddEntry(effPlots[mi][vi], matchingMethods[mi].c_str(), "l");
    }
    legend->Draw();
    c.SaveAs("matchingQA.pdf");
  }

  //--
  // Matching purity
  //--
  c.Clear();
  title->Clear();
  title->AddText("Matching purity");
  title->Draw();
  c.SaveAs("matchingQA.pdf");

  for (int vi = 0; vi < variables.size(); vi++) {
    auto variable = variables[vi];

    for (int mi = 0, colorIndex = 1; mi < matchingMethods.size(); mi++, colorIndex++) {
      auto method = matchingMethods[mi];
      std::string path = std::string("qa-matching/matching/MC/") + method + "/matching-purity/";

      std::string histname = path + variable + "_num";
      hnum = GetTH1(rootFile, histname.c_str());
      if (!hnum) {
        std::cout << "Histogram \"" << histname << "\" not found" << std::endl;
        continue;
      }
      hnum->SetLineColor(colorIndex);

      histname = path + variable + "_den";
      hden = GetTH1(rootFile, histname.c_str());
      if (!hden) {
        std::cout << "Histogram \"" << histname << "\" not found" << std::endl;
        continue;
      }

      hnum->Divide(hden);
      hnum->SetTitle(Form("Matching purity vs. %s", variable.c_str()));

      purityPlots[mi][vi] = hnum;
    }
  }

  for (int vi = 0; vi < variables.size(); vi++) {
    legend->Clear();
    for (int mi = 0; mi < matchingMethods.size(); mi++) {
      if (mi == 0) {
        purityPlots[mi][vi]->SetMinimum(0);
        purityPlots[mi][vi]->SetMaximum(1.2);
        purityPlots[mi][vi]->Draw();
      } else {
        purityPlots[mi][vi]->Draw("same");
      }
      legend->AddEntry(purityPlots[mi][vi], matchingMethods[mi].c_str(), "l");
    }
    legend->Draw();
    c.SaveAs("matchingQA.pdf");
  }
}

TH1F* GetRankingFraction(TH2* rankingHist, int binNum)
{
  TH1F* hist = new TH1F((std::string(rankingHist->GetName()) + "-ranking-fraction-" + std::to_string(binNum)).c_str(),
                        rankingHist->GetTitle(),
                        rankingHist->GetXaxis()->GetNbins(),
                        rankingHist->GetXaxis()->GetXmin(),
                        rankingHist->GetXaxis()->GetXmax());
  hist->GetXaxis()->SetTitle(rankingHist->GetXaxis()->GetTitle());

  for (int i = 1; i <= rankingHist->GetXaxis()->GetNbins(); i++) {
    float total = 0;
    // sum all bins, including overflow
    for (int j = 1; j <= rankingHist->GetYaxis()->GetNbins() + 1; j++) {
      total += rankingHist->GetBinContent(i, j);
    }

    // ratio between matches with correct ranking and all ranked matches
    float fraction = (total > 0) ? rankingHist->GetBinContent(i, binNum) / total : 0.f;
    hist->SetBinContent(i, fraction);
  }

  //hist->SetMinimum(0.5);

  return hist;
}

TH1F* GetGoodRankingFraction(TH2* rankingHist)
{
  TH1F* hist = GetRankingFraction(rankingHist, 2);
  hist->SetTitle("Fraction of correctly ranked matches");
  hist->SetMinimum(0.5);
  hist->SetMaximum(1.05);

  return hist;
}

TH1F* GetMissedRankingFraction(TH2* rankingHist)
{
  TH1F* hist = GetRankingFraction(rankingHist, 1);
  hist->SetTitle("Fraction of missed matches");
  hist->SetMinimum(0);
  hist->SetMaximum(0.1);

  return hist;
}

void MatchRankingCompareMatchingMethods(TFile* rootFile,
                                        const std::vector<std::string>& matchingMethods,
                                        std::string histName)
{
  TH1* h1;
  TH2* h2;
  TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);

  // same particle in MFT and MCH
  legend->Clear();
  c.Clear();
  c.SetLogy(kTRUE);
  int index = 0;
  for (auto method : matchingMethods) {
    std::string path = std::string("qa-matching/matching/MC/") + method + "/";
    h1 = GetTH1(rootFile, (path + histName).c_str());
    if (!h1) {
      std::cout << "Histogram \"" << (path + histName) << "\" not found" << std::endl;
      continue;
    }
    h1->SetLineColor(index + 1);
    // Normalize to unit area
    h1->Scale(1.0 / h1->Integral());
    if (index == 0) {
      h1->Draw("HIST");
      h1->SetMinimum(1.e-4);
      h1->SetMaximum(1.0);
    } else
      h1->Draw("same");
    legend->AddEntry(h1, method.c_str(), "l");
    index += 1;
  }
  legend->Draw();
  c.SaveAs("matchingQA.pdf");

  std::vector<std::string> variables{"P", "Pt", "McParticleDz", "MftTrackMult", "MftTrackType", "DeltaChi2"};
  for (auto variable : variables) {
    legend->Clear();
    c.Clear();
    c.SetLogy(kTRUE);

    legend->SetY1NDC(0.15);
    legend->SetY2NDC(0.35);
    int index = 0;
    auto method = matchingMethods[0];
    std::string path = std::string("qa-matching/matching/MC/") + method + "/";
    h2 = GetTH2(rootFile, (path + histName + "Vs" + variable).c_str());
    TH1* proj = (TH1*)h2->ProjectionX();
    proj->Draw("HIST");
    c.SaveAs("matchingQA.pdf");
    delete proj;
  }

  for (auto variable : variables) {
    legend->Clear();
    c.Clear();
    c.SetLogy(kFALSE);

    legend->SetY1NDC(0.15);
    legend->SetY2NDC(0.35);
    int index = 0;
    for (auto method : matchingMethods) {
      std::string path = std::string("qa-matching/matching/MC/") + method + "/";
      h2 = GetTH2(rootFile, (path + histName + "Vs" + variable).c_str());
      if (false && method == matchingMethods.front()) {
        h2->Draw("col");
        c.SetLogz(kTRUE);
        c.SaveAs("matchingQA.pdf");
        c.SetLogz(kFALSE);
      }
      h1 = GetGoodRankingFraction(h2);
      h1->SetLineColor(index + 1);
      if (index == 0)
        h1->Draw();
      else
        h1->Draw("same");
      legend->AddEntry(h1, method.c_str(), "l");
      index += 1;
    }
    legend->Draw();
    c.SaveAs("matchingQA.pdf");

    legend->SetY1NDC(0.65);
    legend->SetY2NDC(0.85);
    index = 0;
    for (auto method : matchingMethods) {
      std::string path = std::string("qa-matching/matching/MC/") + method + "/";
      h2 = GetTH2(rootFile, (path + histName + "Vs" + variable).c_str());
      h1 = GetMissedRankingFraction(h2);
      h1->SetLineColor(index + 1);
      if (index == 0)
        h1->Draw();
      else
        h1->Draw("same");
      //legend->AddEntry(h1, method.c_str(), "l");
      index += 1;
    }
    legend->Draw();
    c.SaveAs("matchingQA.pdf");
  }
}

void MatchingScoreCompareMatchingMethods(TFile* rootFile,
                                         const std::vector<std::string>& matchingMethods,
                                         std::string histName,
                                         std::vector<int> matchTypeBins,
                                         std::string histTitle)
{
  TH1* h1{nullptr};
  TH2* h2{nullptr};
  TLegend* legend = new TLegend(0.15, 0.8, 0.85, 0.9);
  legend->SetNColumns(3);

  legend->Clear();
  c.Clear();
  c.SetLogy(kFALSE);
  float max = 0;
  // matching score distributions
  std::vector<TH1*> h1vec;
  // cumulative distributions
  std::vector<TH1*> ch1vec;
  for (auto method : matchingMethods) {
    std::string path = std::string("qa-matching/matching/MC/") + method + "/";
    h2 = GetTH2(rootFile, (path + histName).c_str());
    h1 = nullptr;
    for (auto bin : matchTypeBins) {
      if (h1) {
        h1->Add((TH1*)h2->ProjectionY((path + histName + std::to_string(bin)).c_str(), bin, bin));
      } else {
        h1 = (TH1*)h2->ProjectionY((path + histName + std::to_string(bin)).c_str(), bin, bin);
        h1->SetTitle(histTitle.c_str());
      }
    }
    if (h1->GetMaximum() > max) {
      max = h1->GetMaximum();
    }
    h1vec.push_back(h1);

    TH1F* ch1 = (TH1F*)h1->Clone();
    ch1->SetTitle((histTitle + " (cumulative)").c_str());
    float integral = h1->Integral();
    int binMax = h1->GetXaxis()->GetNbins();
    for (int bin = 1; bin <= h1->GetXaxis()->GetNbins(); bin++) {
      float cumulant = h1->Integral(bin, binMax);
      ch1->SetBinContent(bin, cumulant / integral);
    }
    ch1vec.push_back(ch1);
  }

  for (int i = 0; i < h1vec.size(); i++) {
    auto* hist = h1vec[i];
    hist->SetLineColor(i + 1);
    if (i == 0) {
      hist->SetMaximum(2.f * max);
      hist->Draw();
    } else {
      hist->Draw("same");
    }
    legend->AddEntry(hist, matchingMethods[i].c_str(), "l");
  }
  legend->Draw();
  c.SetLogy(kTRUE);
  c.SaveAs("matchingQA.pdf");
  c.SetLogy(kFALSE);

  for (int i = 0; i < ch1vec.size(); i++) {
    auto* hist = ch1vec[i];
    hist->SetLineColor(i + 1);
    if (i == 0) {
      hist->SetMaximum(1.2f);
      hist->Draw();
    } else {
      hist->Draw("same");
    }
    // legend->AddEntry(hist, matchingMethods[i].c_str(), "l");
  }
  legend->Draw();
  c.SaveAs("matchingQA.pdf");
}

// Identical to the above, except splits into variable ranges for each matching method
// Each matching method will have its own histogram to compare the score distributions vs variable ranges
// at some point make a generic function that can plot X vs Y while varying Z
void MatchingScoreCompareVsVariable(TFile* rootFile,
                                    std::string matchingMethod,
                                    std::string histName,
                                    std::vector<int> matchTypeBins,
                                    std::string histTitle,
                                    std::string variable,
                                    const std::vector<std::pair<double, double>> ranges)
{
  TH1* h1;
  TH2* h2;
  TH3* h3;
  TLegend* legend = new TLegend(0.15, 0.8, 0.85, 0.9);
  legend->SetNColumns(3);

  legend->Clear();
  c.Clear();
  c.SetLogy(kFALSE);
  float max = 0;
  // matching score distributions
  std::vector<TH1*> h1vec;
  // cumulative distributions
  std::vector<TH1*> ch1vec;
  std::string path = std::string("qa-matching/matching/MC/") + matchingMethod + "/";
  h3 = GetTH3(rootFile, (path + histName + "Vs" + variable).c_str());

  for (const auto& range : ranges) {
    h1 = nullptr;
    for (auto bin : matchTypeBins) {
      if (h1) {
        h1->Add((TH1*)h3->ProjectionZ(TString::Format("%g<%s<%g_%d", range.first, variable.c_str(), range.second, bin).Data(),
            h3->GetXaxis()->FindBin(range.first + 0.5f * h3->GetXaxis()->GetBinWidth(1)),
            h3->GetXaxis()->FindBin(range.second - 0.5f * h3->GetXaxis()->GetBinWidth(1)),
            bin, bin));
      } else {
        h1 = (TH1*)h3->ProjectionZ(TString::Format("%g<%s<%g_%d", range.first, variable.c_str(), range.second, bin).Data(),
            h3->GetXaxis()->FindBin(range.first + 0.5f * h3->GetXaxis()->GetBinWidth(1)),
            h3->GetXaxis()->FindBin(range.second - 0.5f * h3->GetXaxis()->GetBinWidth(1)),
            bin, bin);
        h1->SetTitle(histTitle.c_str());
      }
    }

    h1->Scale(1.0 / h1->Integral()); // noramalize to 1
    if (h1->GetMaximum() > max) {
      max = h1->GetMaximum();
    }
    h1vec.push_back(h1);

    TH1F* ch1 = (TH1F*)h1->Clone();
    ch1->SetTitle((histTitle + " (cumulative)").c_str());
    float integral = h1->Integral();
    int binMax = h1->GetXaxis()->GetNbins();
    for (int bin = 1; bin <= h1->GetXaxis()->GetNbins(); bin++) {
      float cumulant = h1->Integral(bin, binMax);
      ch1->SetBinContent(bin, cumulant / integral);
    }
    ch1vec.push_back(ch1);
  }

  for (int i = 0; i < h1vec.size(); i++) {
    auto* hist = h1vec[i];
    hist->SetLineColor(i + 1);
    if (i == 0) {
      hist->SetMaximum(5.f * max);
      hist->Draw("HIST");
    } else {
      hist->Draw("HIST SAME");
    }
    legend->AddEntry(hist, hist->GetName(), "l");
  }
  legend->Draw();
  c.SetLogy(kTRUE);
  c.SaveAs("matchingQA.pdf");
  c.SetLogy(kFALSE);

  legend->Clear();
  for (int i = 0; i < ch1vec.size(); i++) {
    auto* hist = ch1vec[i];
    hist->SetLineColor(i + 1);
    if (i == 0) {
      hist->SetMaximum(1.2f);
      hist->Draw("HIST");
    } else {
      hist->Draw("HIST SAME");
    }
    legend->AddEntry(hist, hist->GetName(), "l");
  }
  legend->Draw();
  c.SaveAs("matchingQA.pdf");
}

void PlotMatchRanking(TFile* rootFile)
{
  std::vector<std::string> matchingMethods{
    "Prod",
    "MatchXYPhiTanl",
    "MatchXYPhiTanlMom"};

  TH1* h1;
  TH2* h2;
  TPaveText* title = new TPaveText(0.1, 0.4, 0.9, 0.6, "NDC");
  TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);
  const std::vector<std::pair<double, double>> ranges = {
    {0.f, 10.f},
    {10.f, 15.f},
    {15.f, 20.f},
    {20.f, 30.f},
    {30.f, 50.f},
    {50.f, 100.f}};

  c.Clear();
  title->Clear();
  title->AddText("True match ranking");
  title->AddText("Good MCH tracks, paired with MFT");
  title->Draw();
  c.SaveAs("matchingQA.pdf");
  MatchRankingCompareMatchingMethods(rootFile, matchingMethods, "matchRankingPairedGoodMCH");

  h2 = GetTH2(rootFile, std::string("qa-matching/matching/MC/Prod/matchScoreVsType").c_str());
  if (h2) {
    c.Clear();
    title->Clear();
    title->AddText("Matching score distribution");
    title->Draw();
    c.SaveAs("matchingQA.pdf");

    c.Clear();
    c.SetLogz(kTRUE);
    h2->Draw("col");
    c.SaveAs("matchingQA.pdf");
    c.SetLogz(kFALSE);
  }

  c.Clear();
  title->Clear();
  title->AddText("Matching score distribution");
  title->AddText("Good MCH tracks, true matches");
  title->Draw();
  c.SaveAs("matchingQA.pdf");
  MatchingScoreCompareMatchingMethods(rootFile, matchingMethods, "matchScoreVsType", {1, 5}, "Match score - true matches");

  c.Clear();
  title->Clear();
  title->AddText("Matching score distribution");
  title->AddText("Good MCH tracks, wrong matches");
  title->Draw();
  c.SaveAs("matchingQA.pdf");
  MatchingScoreCompareMatchingMethods(rootFile, matchingMethods, "matchScoreVsType", {2, 6}, "Match score - wrong matches");

  c.Clear();
  title->Clear();
  title->AddText("Matching score distribution");
  title->AddText("Good MCH tracks, decay matches");
  title->Draw();
  c.SaveAs("matchingQA.pdf");
  MatchingScoreCompareMatchingMethods(rootFile, matchingMethods, "matchScoreVsType", {3, 7}, "Match score - decay matches");

  c.Clear();
  title->Clear();
  title->AddText("Matching score distribution");
  title->AddText("Good MCH tracks, fake matches");
  title->Draw();
  c.SaveAs("matchingQA.pdf");
  MatchingScoreCompareMatchingMethods(rootFile, matchingMethods, "matchScoreVsType", {4, 8}, "Match score - fake matches");

  // beginning of my additions
  // opt for per method comparison of p ranges for simplicity, can of course swap toper p range comparison of methods if desired

  for (std::string entry : matchingMethods) {
    c.Clear();
    title->Clear();
    title->AddText(("Matching score distribution, for " + entry).c_str());
    title->AddText("Good MCH tracks, true matches");
    title->Draw();
    c.SaveAs("matchingQA.pdf");
    MatchingScoreCompareVsVariable(rootFile, entry, "matchScoreVsType", {1, 5}, "Match score - true matches", "P", ranges);
  }

  for (std::string entry : matchingMethods) {
    c.Clear();
    title->Clear();
    title->AddText(("Matching score distribution, P dependence, for " + entry).c_str());
    title->AddText("Good MCH tracks, wrong matches");
    title->Draw();
    c.SaveAs("matchingQA.pdf");
    MatchingScoreCompareVsVariable(rootFile, entry, "matchScoreVsType", {2, 6}, "Match score - wrong matches", "P", ranges);
  }

  for (std::string entry : matchingMethods) {
    c.Clear();
    title->Clear();
    title->AddText(("Matching score distribution, P dependence, for " + entry).c_str());
    title->AddText("Good MCH tracks, decay matches");
    title->Draw();
    c.SaveAs("matchingQA.pdf");
    MatchingScoreCompareVsVariable(rootFile, entry, "matchScoreVsType", {3, 7}, "Match score - decay matches", "P", ranges);
  }

  for (std::string entry : matchingMethods) {
    c.Clear();
    title->Clear();
    title->AddText(("Matching score distribution, P dependence, for " + entry).c_str());
    title->AddText("Good MCH tracks, fake matches");
    title->Draw();
    c.SaveAs("matchingQA.pdf");
    MatchingScoreCompareVsVariable(rootFile, entry, "matchScoreVsType", {4, 8}, "Match score - fake matches", "P", ranges);
  }
}

void PlotMatchType(TFile* rootFile)
{
  std::vector<std::string> matchingMethods{
    "Prod",
    "MatchXYPhiTanl",
    "MatchXYPhiTanlMom"};

  TH1* h1;
  TH2* h2;
  TPaveText* title = new TPaveText(0.1, 0.4, 0.9, 0.6, "NDC");
  TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);
  const std::vector<std::pair<double, double>> ranges = {
    {0.f, 10.f},
    {10.f, 15.f},
    {15.f, 20.f},
    {20.f, 30.f},
    {30.f, 50.f},
    {50.f, 100.f}};

  for (auto matchingMethod : matchingMethods) {
    // match type for good MCH tracks
    c.Clear();
    title->Clear();
    title->AddText((std::string("Match type for ") + matchingMethod).c_str());
    title->AddText("Good MCH tracks");
    title->Draw();
    c.SaveAs("matchingQA.pdf");

    std::string path = std::string("qa-matching/matching/MC/") + matchingMethod + "/";
    std::string histName = "matchType";
    std::string variable = "P";
    h2 = GetTH2(rootFile, (path + histName + "Vs" + variable).c_str());
    c.Clear();
    h2->Draw("col");
    c.SetLogz(kTRUE);
    c.SaveAs("matchingQA.pdf");
    c.SetLogz(kFALSE);
  }
}


TH1* GetInvmassForMatchTypes(TH2* hist2, std::vector<int> types)
{
  TH1* result{nullptr};
  for (auto type : types) {
    int ybin = hist2->GetYaxis()->FindBin(type);
    if (!result) {
      result = (TH1*)hist2->ProjectionX(TString::Format("%s_%d", hist2->GetName(), type), ybin, ybin);
    } else {
      TH1* proj = (TH1*)hist2->ProjectionX(TString::Format("%s_%d", hist2->GetName(), type), ybin, ybin);
      result->Add(proj);
    }
  }
  return result;
}


void PlotInvmassForMatchTypes(TFile* rootFile, std::vector<int> types)
{
  std::string typesstr;
  for (auto type : types) {
    if (typesstr.empty()) {
      typesstr = std::format("{:02d}", type);
    } else {
      typesstr += std::format("+{:02d}", type);
    }
  }

  TPaveText* title = new TPaveText(0.1, 0.4, 0.9, 0.6, "NDC");
  c.Clear();
  title->Clear();
  title->AddText(TString::Format("type = %s", typesstr.c_str()));
  title->Draw();
  c.SaveAs("matchingQA.pdf");

  c.Clear();

  TH2* h2 = GetTH2(rootFile, "qa-matching/dimuon/MC/invariantMass_MuonKine_GlobalMuonCuts_vs_match_type");
  TH1* h1 = GetInvmassForMatchTypes(h2, types);
  if (h1) {
    h1->SetTitle(TString::Format("%s, type = %s", h2->GetTitle(), typesstr.c_str()));
    h1->Draw();
    c.SaveAs("matchingQA.pdf");
  }

  h2 = GetTH2(rootFile, "qa-matching/dimuon/MC/invariantMass_ScaledMftKine_GlobalMuonCuts_vs_match_type");
  h1 = GetInvmassForMatchTypes(h2, types);
  if (h1) {
    h1->SetTitle(TString::Format("%s, type = %s", h2->GetTitle(), typesstr.c_str()));
    h1->Draw();
    c.SaveAs("matchingQA.pdf");
  }

  h2 = GetTH2(rootFile, "qa-matching/dimuon/MC/invariantMass_MuonKine_GlobalMuonCuts_GoodMatches_vs_match_type");
  h1 = GetInvmassForMatchTypes(h2, types);
  if (h1) {
    h1->SetTitle(TString::Format("%s, type = %s", h2->GetTitle(), typesstr.c_str()));
    h1->Draw();
    c.SaveAs("matchingQA.pdf");
  }

  h2 = GetTH2(rootFile, "qa-matching/dimuon/MC/invariantMass_ScaledMftKine_GlobalMuonCuts_GoodMatches_vs_match_type");
  h1 = GetInvmassForMatchTypes(h2, types);
  if (h1) {
    h1->SetTitle(TString::Format("%s, type = %s", h2->GetTitle(), typesstr.c_str()));
    h1->Draw();
    c.SaveAs("matchingQA.pdf");
  }
}

// Plots the matching score distributions as a breakdown for each of the match types in matchtypevec for each of the matching methods
// Should eventually be refactored into one function with MatchingScoreCompareVsVariable due to the high similarity
void MatchTypePlot(TFile* rootFile, std::string histName, std::string method)
{
  TH1* h1{nullptr};
  TH2* h2{nullptr};
  //TLegend* legend = new TLegend(0.7, 0.1, 0.9, 0.4);
  //legend->SetNColumns(1);
  TLegend* legend = new TLegend(0.15, 0.8, 0.85, 0.9);
  legend->SetNColumns(4);

  legend->Clear();
  c.Clear();
  c.SetLogy(kFALSE);
  float max = 0;
  // matching score distributions
  std::vector<TH1*> h1vec;
  // cumulative distributions
  std::vector<TH1*> ch1vec;
  std::string path = std::string("qa-matching/matching/MC/") + method + "/";

  for (matchType matchtype : matchTypesVec) {
    h2 = GetTH2(rootFile, (path + histName).c_str());
    h1 = nullptr;
    for (auto bin : matchtype.bins) {
      if (h1) {
        h1->Add((TH1 *)h2->ProjectionY((path + histName + std::to_string(bin)).c_str(), bin, bin));
      } else {
        h1 = (TH1 *)h2->ProjectionY((path + histName + std::to_string(bin)).c_str(), bin, bin);
        h1->SetTitle(("Match Score - "+ method).c_str());
      }
    }
    if (h1->GetMaximum() > max) {
      max = h1->GetMaximum();
    }
    h1vec.push_back(h1);

    TH1F* ch1 = (TH1F*)h1->Clone();
    ch1->SetTitle(("Match Score - "+ method + " (cumulative)").c_str());
    float integral = h1->Integral();
    int binMax = h1->GetXaxis()->GetNbins();
    for (int bin = 1; bin <= h1->GetXaxis()->GetNbins(); bin++) {
      float cumulant = h1->Integral(bin, binMax);
      ch1->SetBinContent(bin, cumulant / integral);
    }
    ch1vec.push_back(ch1);
  }

  for (int i = 0; i < h1vec.size(); i++) {
    auto* hist = h1vec[i];
    hist->SetLineColor(i + 1);
    if (i == 0) {
      hist->SetMaximum(10.f * max);
      hist->Draw();
    } else {
      hist->Draw("same");
    }
    legend->AddEntry(hist, matchTypesVec[i].name.c_str(), "l");
  }
  legend->Draw();
  c.SetLogy(kTRUE);
  c.SaveAs("matchingQA.pdf");
  c.SetLogy(kFALSE);

  for (int i = 0; i < ch1vec.size(); i++) {
    auto* hist = ch1vec[i];
    hist->SetLineColor(i + 1);
    if (i == 0) {
      hist->SetMinimum(0.f);
      hist->SetMaximum(1.2f);
      hist->Draw();
    } else {
      hist->Draw("same");
    }
    // legend->AddEntry(hist, matchingMethods[i].c_str(), "l");
  }
  //legend->SetX1NDC(0.1);
  //legend->SetX2NDC(0.3);
  legend->Draw();
  c.SaveAs("matchingQA.pdf");
}

void PlotMatchScoreVsMatchType(TFile* rootFile)
{
  std::vector<std::string> matchingMethods{
      "Prod",
      "MatchXYPhiTanl",
      "MatchXYPhiTanlMom"};

  for (auto matchingMethod : matchingMethods) {
    c.Clear();
    TPaveText* title = new TPaveText(0.1, 0.4, 0.9, 0.6, "NDC");
    title->Clear();
    title->AddText(("Matching score vs match type for " + matchingMethod).c_str());
    title->Draw();
    c.SaveAs("matchingQA.pdf");

    std::string path = std::string("qa-matching/matching/MC/") + matchingMethod + "/";
    MatchTypePlot(rootFile,"matchScoreVsType", matchingMethod);
  }
}

void PlotInvmass(TFile *rootFile)
{
  c.SetGridx(false);
  c.SetGridy(false);

  TPaveText* title = new TPaveText(0.1, 0.4, 0.9, 0.6, "NDC");
  c.Clear();
  title->Clear();
  title->AddText("Invariant mass distributions");
  title->Draw();
  c.SaveAs("matchingQA.pdf");

  // ============
  // Type 00
  // ============
  PlotInvmassForMatchTypes(rootFile, {0});

  // ============
  // Type 01+10
  // ============
  PlotInvmassForMatchTypes(rootFile, {1, 10});

  // ============
  // Type 11
  // ============
  PlotInvmassForMatchTypes(rootFile, {11});

  // ============
  // Type 02+20
  // ============
  PlotInvmassForMatchTypes(rootFile, {2, 20});

  // ============
  // Type 03+30
  // ============
  PlotInvmassForMatchTypes(rootFile, {3, 30});

  // ============
  // Type 33
  // ============
  PlotInvmassForMatchTypes(rootFile, {33});
}

void matchingQA()
{
  TFile* fAnalysisResults;

  //fAnalysisResults = new TFile("outputs/LHC25i4/AnalysisResults.root");
  fAnalysisResults = new TFile("AnalysisResults.root");

  gStyle->SetOptStat(0);
  // gStyle->SetOptStat(1111);
  // gStyle->SetOptFit(1111);

  c.SaveAs("matchingQA.pdf(");

  TH1* h1 = GetTH1(fAnalysisResults, "qa-matching/nTracksPerType");
  c.Clear();
  h1->Draw("HIST");
  c.SaveAs("matchingQA.pdf");

  h1 = GetTH1(fAnalysisResults, "qa-matching/tracksMultiplicityMFT");
  c.Clear();
  h1->GetXaxis()->SetRangeUser(1, h1->GetXaxis()->GetXmax());
  h1->Draw("HIST");
  c.SetLogx(kTRUE);
  c.SetLogy(kTRUE);
  c.SaveAs("matchingQA.pdf");
  c.SetLogx(kFALSE);
  c.SetLogy(kFALSE);

  h1 = GetTH1(fAnalysisResults, "qa-matching/tracksMultiplicityMCH");
  c.Clear();
  h1->GetXaxis()->SetRangeUser(1, h1->GetXaxis()->GetXmax());
  h1->Draw("HIST");
  c.SetLogx(kTRUE);
  c.SetLogy(kTRUE);
  c.SaveAs("matchingQA.pdf");
  c.SetLogx(kFALSE);
  c.SetLogy(kFALSE);

  c.Clear();

  PlotMatchEfficiancyAndPurity(fAnalysisResults);

  PlotMatchRanking(fAnalysisResults);

  PlotMatchType(fAnalysisResults);

  PlotMatchScoreVsMatchType(fAnalysisResults);

  PlotInvmass(fAnalysisResults);

  c.Clear();
  c.SaveAs("matchingQA.pdf)");
}
