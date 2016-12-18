R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/HWWLVJRooPdfs_cxx.so)
R__LOAD_LIBRARY(/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/PDFs/PdfDiagonalizer_cc.so)
#include "/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/readFitParam.h"
using namespace RooFit;

void rooFitTest(string channel, string catcut, bool pullTest=true){

  // Suppress all the INFO message

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(true);
  gROOT->ProcessLine("gErrorIgnoreLevel=kWarning;");

  // Input files and sum all backgrounds

  TChain* tree = new TChain("tree");

  tree->Add(Form("Zjets/DYJetsToLL_M-50_HT-100to200_13TeV_%sMiniTree.root", channel.data()));
  tree->Add(Form("Zjets/DYJetsToLL_M-50_HT-200to400_13TeV_%sMiniTree.root", channel.data()));
  tree->Add(Form("Zjets/DYJetsToLL_M-50_HT-400to600_13TeV_%sMiniTree.root", channel.data()));
  tree->Add(Form("Zjets/DYJetsToLL_M-50_HT-600toInf_13TeV_%sMiniTree.root", channel.data()));

  // Define all the variables from the trees

  RooRealVar cat     ("cat", "", 0, 2);
  RooRealVar mJet    ("prmass", "M_{jet}", 30., 300., "GeV");
  RooRealVar evWeight("evweight", "", 0., 1.e3);

  // Set the range in jet mass

  mJet.setRange("All",  30.,  300.);
  mJet.setRange("SB_l", 30.,  65.);
  mJet.setRange("SB_h", 135., 300.);
  mJet.setRange("SG",   105., 135.);

  RooBinning bin_mJet(54, 30, 300);

  TCut cut_bTag = Form("cat==%s", catcut.data());
  TCut cut_sb   = "prmass>30 && !(prmass>65 && prmass<135) && prmass<300";
  TCut cut_sg   = "prmass>105 && prmass<135";

  // Create a dataset from a tree -> to process an unbinned likelihood fitting

  RooDataSet set_Dom  ("set_Dom",   "set_Dom",   RooArgSet(cat, mJet, evWeight), Cut(cut_bTag),           WeightVar(evWeight), Import(*tree));
  RooDataSet set_sbDom("set_sbDom", "set_sbDom", RooArgSet(cat, mJet, evWeight), Cut(cut_bTag && cut_sb), WeightVar(evWeight), Import(*tree));

  // Total events number

  RooRealVar nEv_Dom  ("nEv_Dom",   "nEv_Dom",   set_Dom  .sumEntries(), set_Dom  .sumEntries()*0.5, set_Dom  .sumEntries()*1.5);
  RooRealVar nEv_sbDom("nEv_sbDom", "nEv_sbDom", set_sbDom.sumEntries(), set_sbDom.sumEntries()*0.5, set_sbDom.sumEntries()*1.5);

  // Set fit parameters for jet mass

  param myVal(channel.data(), catcut.data());

  RooRealVar j_mc  ("j_mc",   "j_mc",   myVal.value("j_mc"), myVal.value("j_mcMin"), myVal.value("j_mcMax"));
  RooRealVar j_sbmc("j_sbmc", "j_sbmc", myVal.value("j_mc"), myVal.value("j_mcMin"), myVal.value("j_mcMax"));

  // Create pdf for jet mass

  RooGenericPdf pdf_McJet  ("pdf_McJet",   "pdf_McJet",   "exp(-@0/@1)", RooArgSet(mJet, j_mc));
  RooGenericPdf pdf_sbMcJet("pdf_sbMcJet", "pdf_sbMcJet", "exp(-@0/@1)", RooArgSet(mJet, j_sbmc));

  // Extended pdf from RooGenericPdf

  RooExtendPdf ext_McJet  ("ext_McJet",   "ext_McJet",   pdf_McJet,   nEv_Dom);
  RooExtendPdf ext_sbMcJet("ext_sbMcJet", "ext_sbMcJet", pdf_sbMcJet, nEv_sbDom);

  // Fit jet mass

  RooLinkedList cmdJetList;

  cmdJetList.Add(Save(true).Clone());
  cmdJetList.Add(Minos(true).Clone());
  cmdJetList.Add(Offset(true).Clone());
  cmdJetList.Add(Extended(true).Clone());
  cmdJetList.Add(SumW2Error(false).Clone());
  cmdJetList.Add(NumCPU(8).Clone());
  cmdJetList.Add(Strategy(2).Clone());
  cmdJetList.Add(Range("All").Clone());
  cmdJetList.Add(Minimizer("Minuit2","migrad").Clone());

  RooLinkedList cmdsbJetList;

  cmdsbJetList.Add(Save(true).Clone());
  cmdsbJetList.Add(Minos(true).Clone());
  cmdsbJetList.Add(Offset(true).Clone());
  cmdsbJetList.Add(Extended(true).Clone());
  cmdsbJetList.Add(SumW2Error(false).Clone());
  cmdsbJetList.Add(NumCPU(8).Clone());
  cmdsbJetList.Add(Strategy(2).Clone());
  cmdsbJetList.Add(Range("SB_l,SB_h").Clone());
  cmdsbJetList.Add(Minimizer("Minuit2","migrad").Clone());

  RooFitResult* res_McJet   = ext_McJet  .fitTo(set_Dom,   cmdJetList);
  RooFitResult* res_sbMcJet = ext_sbMcJet.fitTo(set_sbDom, cmdsbJetList);

  fprintf(stdout, "j_mc   = %.3f +- %.3f\n", j_mc  .getVal(), j_mc  .getError());
  fprintf(stdout, "j_sbmc = %.3f +- %.3f\n", j_sbmc.getVal(), j_sbmc.getError());

  // Get real normalize factor from MC

  RooAbsReal* nFit_sg = ext_sbMcJet.createIntegral(mJet, Range("SG"));
  RooAbsReal* nFit_sb = ext_sbMcJet.createIntegral(mJet, Range("SB_l,SB_h"));

  float realNormFactor = nEv_sbDom.getVal()*(nFit_sg->getVal()/nFit_sb->getVal());

  // Produce n toyMCs to study fit bias and pull  
  // Properties of pull: mean is 0 if there is no bias; width is 1 if error is correct
  // Fit is converge: the fit really finds a set of parameter values that minimizes -log likelihood instead of finding a local minima

  TH1F* h_bias = new TH1F("h_bias", "", 19, -9.5, 9.5);
  TH1F* h_pull = new TH1F("h_pull", "", 19, -9.5, 9.5);

  // TFile toyOutput(Form("toyMCdataSet_%s_cat%s.root", channel.data(), catcut.data()), "recreate");

  for( int ntoy = 1000; ntoy > 0; --ntoy ){

    if( !pullTest ) break;

    RooDataSet* set_toyMc = ext_sbMcJet.generate(mJet, Extended(), ProtoData(set_Dom));
    RooDataSet  thisToyMc("thisToyMc", "thisToyMc", RooArgSet(mJet), Cut(cut_sb), Import(*set_toyMc));

    // thisToyMc.Write(Form("toy%04i", ntoy));

    RooRealVar nEv_toyMc("nEv_toyMc", "nEv_toyMc", thisToyMc.sumEntries(), thisToyMc.sumEntries()*0.5, thisToyMc.sumEntries()*1.5);
    RooRealVar j_toymc("j_toymc", "j_toymc", myVal.value("j_mc"), myVal.value("j_mcMin"), myVal.value("j_mcMax"));

    RooGenericPdf pdf_toyMcJet("pdf_toyMcJet", "pdf_toyMcJet", "exp(-@0/@1)", RooArgSet(mJet, j_toymc));
    RooExtendPdf  ext_toyMcJet("ext_toyMcJet", "ext_toyMcJet", pdf_toyMcJet, nEv_toyMc);

    RooLinkedList cmdtoyMcJetList;

    cmdtoyMcJetList.Add(Save(true).Clone());
    cmdtoyMcJetList.Add(Minos(true).Clone());
    cmdtoyMcJetList.Add(Offset(true).Clone());
    cmdtoyMcJetList.Add(Extended(true).Clone());
    cmdtoyMcJetList.Add(SumW2Error(false).Clone());
    cmdtoyMcJetList.Add(NumCPU(8).Clone());
    cmdtoyMcJetList.Add(Strategy(2).Clone());
    cmdtoyMcJetList.Add(Range("SB_l,SB_h").Clone());
    cmdtoyMcJetList.Add(Minimizer("Minuit2","migrad").Clone());

    RooFitResult* res_toyMcJet = ext_toyMcJet.fitTo(thisToyMc, cmdtoyMcJetList);

    fprintf(stdout, "nToy=%i\tj_mcToy=%f\tstatus=%i\n", ntoy, j_toymc.getVal(), res_toyMcJet->status());

    if( res_toyMcJet->status() != 0 ) continue;
   
    // calulate normalize factor

    RooAbsReal* nSIGFit = ext_toyMcJet.createIntegral(mJet, Range("SG"));
    RooAbsReal* nSBFit  = ext_toyMcJet.createIntegral(mJet, Range("SB_l,SB_h"));

    RooRealVar nSBHist("nSBHist", "nSBHist", 0, 1e4);

    nSBHist.setVal(thisToyMc.sumEntries());
    nSBHist.setConstant(true);

    float toyNormFactor = nSBHist.getVal()*(nSIGFit->getVal()/nSBFit->getVal());

    RooFormulaVar formula("formula", "formula", "@0*@1/@2", RooArgList(nSBHist, *nSIGFit, *nSBFit));

    h_bias->Fill((toyNormFactor - realNormFactor)/realNormFactor);
    h_pull->Fill((toyNormFactor - realNormFactor)/formula.getPropagatedError(*res_toyMcJet));

  } // End of ntoy loop

  // toyOutput.Close();

  RooRealVar bias("bias", "Bias", -9.5, 9.5);
  RooRealVar pull("pull", "Pull", -9.5, 9.5);

  RooDataHist hbias("hbias", "", bias, Import(*h_bias));
  RooDataHist hpull("hpull", "", pull, Import(*h_pull));

  RooRealVar gbmean("gbmean", "mean", 0.0, -5.0, 5.0);
  RooRealVar gbsigma("gbsigma", "sigma", 1.0, 0.1, 4.0);

  RooRealVar gpmean("gpmean", "mean", 0.0, -5.0, 5.0);
  RooRealVar gpsigma("gpsigma", "sigma", 1.0, 0.1, 4.0);

  RooGaussian gb("gb", "gauss", bias, gbmean, gbsigma);
  RooGaussian gp("gp", "gauss", pull, gpmean, gpsigma);

  gb.fitTo(hbias);
  gp.fitTo(hpull);

  // Plot the results on frame 

  RooPlot* mJetFrame        = mJet.frame();
  RooPlot* mJetSBFrame      = mJet.frame();
  RooPlot* biasFrame        = bias.frame();
  RooPlot* pullFrame        = pull.frame();
  RooPlot* mJetPullFrame    = mJet.frame();
  RooPlot* mJetSBPullFrame  = mJet.frame();

  set_Dom  .plotOn(mJetFrame, Binning(bin_mJet)); 
  ext_McJet.plotOn(mJetFrame, Range("All"), VisualizeError(*res_McJet,1,false), FillStyle(3002));
  set_Dom  .plotOn(mJetFrame, Binning(bin_mJet));
  ext_McJet.plotOn(mJetFrame, Range("All"));

  mJetPullFrame->addObject(mJetFrame->pullHist(), "P");

  set_sbDom  .plotOn(mJetSBFrame, Binning(bin_mJet));
  ext_sbMcJet.plotOn(mJetSBFrame, Range("All"), VisualizeError(*res_sbMcJet,1,false), FillStyle(3002));
  set_sbDom  .plotOn(mJetSBFrame, Binning(bin_mJet));
  ext_sbMcJet.plotOn(mJetSBFrame, Range("All"));

  mJetSBPullFrame->addObject(mJetSBFrame->pullHist(), "P");

  ext_McJet.plotOn(mJetSBFrame, Normalization(set_Dom.sumEntries(),RooAbsReal::NumEvent), Range("All"), LineStyle(7), LineColor(kRed));

  hbias.plotOn(biasFrame);
  gb.plotOn(biasFrame);
  gb.paramOn(biasFrame,Layout(0.65,0.9,0.8));

  hpull.plotOn(pullFrame);
  gp.plotOn(pullFrame);
  gp.paramOn(pullFrame,Layout(0.65,0.9,0.8));

  // Output results

  TLatex lar;

  lar.SetTextSize(0.035);
  lar.SetLineWidth(5);
  
  float up_height = 0.82;
  float dw_height = (1-up_height)*1.445;

  TCanvas c0("c0","",0,0,1000,800);
  TLegend leg0(0.60,0.67,0.85,0.80);
  
  c0.Divide(1,2);

  TPad* c0_up = (TPad*)c0.GetListOfPrimitives()->FindObject("c0_1");
  TPad* c0_dw = (TPad*)c0.GetListOfPrimitives()->FindObject("c0_2"); 

  c0_up->SetPad(0,1-up_height,1,1);
  c0_dw->SetPad(0,0,1,dw_height);
  c0_dw->SetBottomMargin(0.25);
  c0_up->cd()->SetLogy(1);

  leg0.AddEntry(mJetFrame->findObject(mJetFrame->nameOf(0)), "MC with statistical errors", "lep");
  leg0.AddEntry(mJetFrame->findObject(mJetFrame->nameOf(3)), "Fit curve with errors", "l");
  leg0.Draw();

  mJetFrame->addObject(&leg0);
  mJetFrame->SetTitle("");
  mJetFrame->SetMinimum(1e-3);
  mJetFrame->SetMaximum(catcut=="1"?100:10);
  mJetFrame->GetXaxis()->SetTitle("");
  mJetFrame->GetXaxis()->SetLabelOffset(999);
  mJetFrame->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.62, 0.92, "L = 2.51 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.83, Form("%s, %s btag", channel.data(), catcut.data()));

  c0_up->RedrawAxis();

  c0_dw->cd()->SetLogy(0);

  mJetPullFrame->SetTitle("");
  mJetPullFrame->GetYaxis()->SetTitle("Pulls");
  mJetPullFrame->GetYaxis()->SetTitleOffset(0.25);
  mJetPullFrame->GetXaxis()->SetLabelSize(0.125);
  mJetPullFrame->GetXaxis()->SetTitleSize(0.125);
  mJetPullFrame->GetYaxis()->SetLabelSize(0.125);
  mJetPullFrame->GetYaxis()->SetTitleSize(0.125);
  mJetPullFrame->GetYaxis()->SetNdivisions(505);
  mJetPullFrame->SetMinimum(-4);
  mJetPullFrame->SetMaximum(4);
  mJetPullFrame->Draw();

  c0.Draw();
  c0.Print(Form("rooFit_toyMC_%s_cat%s.pdf(", channel.data(), catcut.data()));

  TCanvas c1("c1","",0,0,1000,800);
  TLegend leg1(0.60,0.62,0.85,0.80);

  c1.Divide(1,2);

  TPad* c1_up = (TPad*)c1.GetListOfPrimitives()->FindObject("c1_1");
  TPad* c1_dw = (TPad*)c1.GetListOfPrimitives()->FindObject("c1_2");

  c1_up->SetPad(0,1-up_height,1,1);
  c1_dw->SetPad(0,0,1,dw_height);
  c1_dw->SetBottomMargin(0.25);
  c1_up->cd()->SetLogy(1);

  leg1.AddEntry(mJetSBFrame->findObject(mJetSBFrame->nameOf(0)), "MC with statistical errors", "lep");
  leg1.AddEntry(mJetSBFrame->findObject(mJetSBFrame->nameOf(3)), "Fit curve with errors", "l");
  leg1.AddEntry(mJetSBFrame->findObject(mJetSBFrame->nameOf(4)), "Fit curve of all range", "l");
  leg1.Draw();

  mJetSBFrame->addObject(&leg1);
  mJetSBFrame->SetTitle("");
  mJetSBFrame->SetMinimum(1e-3);
  mJetSBFrame->SetMaximum(catcut=="1"?100:10);
  mJetSBFrame->GetXaxis()->SetTitle("");
  mJetSBFrame->GetXaxis()->SetLabelOffset(999);
  mJetSBFrame->Draw();

  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.62, 0.92, "L = 2.51 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.15, 0.83, Form("%s, %s btag", channel.data(), catcut.data()));

  c1_up->RedrawAxis();
  c1_dw->cd()->SetLogy(0);

  mJetSBPullFrame->SetTitle("");
  mJetSBPullFrame->GetYaxis()->SetTitle("Pulls");
  mJetSBPullFrame->GetYaxis()->SetTitleOffset(0.25);
  mJetSBPullFrame->GetXaxis()->SetLabelSize(0.125);
  mJetSBPullFrame->GetXaxis()->SetTitleSize(0.125);
  mJetSBPullFrame->GetYaxis()->SetLabelSize(0.125);
  mJetSBPullFrame->GetYaxis()->SetTitleSize(0.125);
  mJetSBPullFrame->GetYaxis()->SetNdivisions(505);
  mJetSBPullFrame->SetMinimum(-4);
  mJetSBPullFrame->SetMaximum(4);
  mJetSBPullFrame->Draw();

  c1.Draw();
  c1.Print(Form("rooFit_toyMC_%s_cat%s.pdf",  channel.data(), catcut.data()));
  
  TCanvas c("c","",0,0,1000,800);

  c.Clear();
  c.cd();
  biasFrame->getAttText()->SetTextSize(0.025);
  biasFrame->SetTitle("");  
  biasFrame->Draw();
  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.55, 0.92, "L = 2.51 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.72, 0.83, Form("%s, %s btag", channel.data(), catcut.data()));
  c.Print(Form("rooFit_toyMC_%s_cat%s.pdf",  channel.data(), catcut.data()));

  c.Clear();
  c.cd();
  pullFrame->getAttText()->SetTextSize(0.025);
  pullFrame->SetTitle("");
  pullFrame->Draw();
  lar.DrawLatexNDC(0.12, 0.92, "CMS #it{#bf{Simulation}}");
  lar.DrawLatexNDC(0.55, 0.92, "L = 2.51 fb^{-1} at #sqrt{s} = 13 TeV");
  lar.DrawLatexNDC(0.72, 0.83, Form("%s, %s btag", channel.data(), catcut.data()));
  c.Print(Form("rooFit_toyMC_%s_cat%s.pdf)", channel.data(), catcut.data()));

}
