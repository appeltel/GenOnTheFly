{

  const double PTMAX = 500.0;
  const double PTJMIN = 3.0;
  const int REBIN = 2;
  const double ZMIN = 1e-13;


    // pretty colz
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

  TFile * f = new TFile("AnaQCD_Z2_5020TeV_400k_Combined.root");

  TH1F * jet = (TH1F*) f->Get("qcdAna/jetspectrum");
  jet->GetYaxis()->SetTitle("Jet d#sigma/dp_{T}"); 
  TH2F * qsj = (TH2F*) f->Get("qcdAna/jet_qscale2D");
  qsj->GetYaxis()->SetTitle("p_{T} Jet [GeV/c]");
  qsj->GetXaxis()->SetTitle("#hat{p}_{T} [GeV/c]");


  TH1F * jet_trunc = jet->Clone("jet_trunc");
  TH1F * jetRE = jet->Clone("jetRE");
  TH1F * jet_truncRE = jet->Clone("jet_truncRE");
  TH2F * qsj_rel = qsj->Clone("qsj_rel");
  TH2F * qsj_relE = qsj->Clone("qsj_relE");

  TF1 * cut = new TF1("cut","[0]+[1]*x",0,1000);
  cut->SetParameters(5,1.25);

  for( int i=0; i<=qsj->GetNbinsY()+1; i++)
  {
    double pt = jet->GetBinCenter(i);
    double jetxsecCheck = 0.0;
    double jetxsecCheckErr2 = 0.0;
    double jetxsecTrunc = 0.0;
    double jetxsecTruncErr2 =0.0;
    double jetxsec = jet->GetBinContent(i);
    double jetxsecErr = jet->GetBinError(i);

    for( int j=0; j<=qsj->GetNbinsX()+1; j++)
    {
      double pthat = qsj->GetXaxis()->GetBinCenter(j);
      if( pt < cut->Eval(pthat) )
      {
        jetxsecTrunc += qsj->GetBinContent(j,i);
        jetxsecTruncErr2 += qsj->GetBinError(j,i)*qsj->GetBinError(j,i);
      }
      jetxsecCheck += qsj->GetBinContent(j,i);
      jetxsecCheckErr2 += qsj->GetBinError(j,i)*qsj->GetBinError(j,i);
      qsj_rel->SetBinContent(j,i, qsj->GetBinContent(j,i)/jetxsec ); 
      qsj_relE->SetBinContent(j,i, qsj->GetBinError(j,i)/jetxsecErr );
    }
 
    jet_trunc->SetBinContent(i,jetxsecTrunc);
    jet_trunc->SetBinError(i,sqrt(jetxsecTruncErr2));
    jetRE->SetBinContent(i,jetxsecErr/jetxsec);
    jetRE->SetBinError(i,0.0);
    jet_truncRE->SetBinContent(i,sqrt(jetxsecTruncErr2)/jetxsecTrunc);
    jet_truncRE->SetBinError(i,0.0);
    // check that I added right
    cout << pt << "   " << jetxsec << "  "  << jetxsecCheck 
               << "   " << jetxsecErr << "  "  << sqrt(jetxsecCheckErr2) << endl; 

  } 

  double pthatbins[9] = {20,30,50,80,120,170,230,300,380};
  TLine * lpthat[9];
  for( int i=0; i<9; i++)
  {
    lpthat[i] = new TLine(pthatbins[i],PTJMIN,pthatbins[i],PTMAX);
    lpthat[i]->SetLineWidth(1);
    lpthat[i]->SetLineStyle(2);
  }

  if ( REBIN > 1 )
  {
    qsj->Rebin2D(REBIN,REBIN);
    qsj_rel->Rebin2D(REBIN,REBIN);
    qsj_rel->Scale(1.0/(double)(REBIN*REBIN));
    qsj_relE->Rebin2D(REBIN,REBIN);
    qsj_relE->Scale(1.0/(double)(REBIN*REBIN));
   
  }


  TCanvas * c0 = new TCanvas("c0","Jet #sigma by #hat{p}_{T}",600,600);
  c0->cd();
  c0->SetLogz();

  qsj->GetXaxis()->SetRangeUser(0,PTMAX);
  qsj->GetYaxis()->SetRangeUser(PTJMIN,PTMAX);
  qsj->SetMinimum(ZMIN);
  qsj->Draw("colz");

  cut->SetLineWidth(2);
  cut->Draw("same");
  for( int i=0; i<9; i++) lpthat[i]->Draw();

  TCanvas * c1 = new TCanvas("c1","fraction of jet #sigma by #hat{p}_{T}",600,600);
  c1->cd();
  c1->SetLogz();

  qsj_rel->GetXaxis()->SetRangeUser(0,PTMAX);
  qsj_rel->GetYaxis()->SetRangeUser(PTJMIN,PTMAX);
  qsj_rel->SetMinimum(1e-4);
  qsj_rel->Draw("colz");

  cut->SetLineWidth(2);
  cut->Draw("same");
  for( int i=0; i<9; i++) lpthat[i]->Draw();

  TCanvas * c2 = new TCanvas("c2","fraction of jet #sigma error by #hat{p}_{T}",600,600);
  c2->cd();
  c2->SetLogz();

  qsj_relE->GetXaxis()->SetRangeUser(0,PTMAX);
  qsj_relE->GetYaxis()->SetRangeUser(PTJMIN,PTMAX);
  qsj_relE->SetMinimum(1e-4);
  qsj_relE->Draw("colz");
  cut->Draw("same");
  for( int i=0; i<9; i++) lpthat[i]->Draw();

  TCanvas * c3 = new TCanvas("c3","jet #sigma and truncated",600,600);
  c3->cd();
  c3->SetLogy();

  jet->SetMarkerStyle(20);
  jet->GetXaxis()->SetRangeUser(0,PTMAX);
  jet->Draw();
  jet_trunc->SetMarkerStyle(24);
  jet_trunc->SetMarkerColor(kRed);
  jet_trunc->SetLineColor(kRed);
  jet_trunc->Draw("same");

  TLegend * leg = new TLegend(0.6,0.6,0.9,0.8);
  leg->SetFillColor(kWhite);
  leg->AddEntry(jet,"All Jets","lp");
  leg->AddEntry(jet_trunc,"Truncated Jets","lp");
  leg->Draw(); 

  TCanvas * c4 = new TCanvas("c4","jet #sigma ration 4C-trunc / 4C",600,600);
  c4->cd();
  
  TH1F * rjet = jet_trunc->Clone("rjet");
  rjet->Divide(rjet,jet,1.0,1.0,"B");
  rjet->GetXaxis()->SetRangeUser(PTJMIN,PTMAX);
  rjet->Draw();

  TCanvas * c5 = new TCanvas("c5","jet #sigma rel error",600,600);
  c5->cd();
  c5->SetLogy();

  jetRE->SetMarkerStyle(20);
  jetRE->GetXaxis()->SetRangeUser(PTJMIN,PTMAX);
  jet_truncRE->GetXaxis()->SetRangeUser(PTJMIN,PTMAX);
  jetRE->GetYaxis()->SetTitle("#delta#sigma/#sigma");
  jetRE->Draw("P");
  jet_truncRE->SetMarkerStyle(24);
  jet_truncRE->SetMarkerColor(kRed);
  jet_truncRE->SetLineColor(kRed);
  jet_truncRE->Draw("same P");
  leg->Draw();

}
