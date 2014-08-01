gROOT->Reset();

//TFile outFile("histograms/DecPi0_21012.root");
TFile inFile("histograms/DecPi0_20903.root");
TFile outFile("histograms/DecPi0_20898.root");
//TFile inFile("histograms/DecPi0_21013.root");
TFile compareFile("histograms/DecPi0_21052.root");

void main()
{
	c1 = new TCanvas("c1", "Data", 200, 10, 700, 500);
	c1->Divide(1,2);

	TH2D *plotcp = (TH2D*)compareFile.Get("THR_PhiCMP_v_MMissPi0P");
	TH2D *plotcr = (TH2D*)compareFile.Get("THR_PhiCMR_v_MMissPi0R");
	TH1D *scalersc = (TH1D*)compareFile.Get("SumScalers152to503");
	
	TH2D *plotp = (TH2D*)inFile.Get("THR_PhiCMP_v_MMissPi0P");
	TH2D *plotr = (TH2D*)inFile.Get("THR_PhiCMR_v_MMissPi0R");
	TH1D *scalers = (TH1D*)inFile.Get("SumScalers152to503");

	TH2D *plotp1 = (TH2D*)outFile.Get("THR_PhiCMP_v_MMissPi0P");
	TH2D *plotr1 = (TH2D*)outFile.Get("THR_PhiCMR_v_MMissPi0R");
	TH1D *scalers1 = (TH1D*)outFile.Get("SumScalers152to503");

	TH1D *gr_cp = plotcp->ProjectionY("five",100,200,"e");
	TH1D *gr_cr = plotcr->ProjectionY("six",100,200,"e");
	TH1D *gr_p = plotp->ProjectionY("one",100,200,"e");
	TH1D *gr1_p = plotp1->ProjectionY("three",100,200,"e");
	TH1D *gr_r = plotr->ProjectionY("two",100,200,"e");
	TH1D *gr1_r = plotr1->ProjectionY("four",100,200,"e");

	Double_t factor1 = scalers->Integral(229,245);
	Double_t factor2 = scalers1->Integral(229,245);
	Double_t factorrun = scalersc->Integral(229,245);
	Double_t scalepara = factorrun/factor1;
	Double_t scaleperp = factorrun/factor2;

	gr_cp->Add(gr_cr, -0.5);

	c1->cd(1);
	gr_p->Add(gr_r, -0.5);
	gr_p->Scale(scalepara);
	gr_p->SetTitle("The run is perp");
	TH1D *sum = (TH1D*)gr_p->Clone();
	TH1D *diff = (TH1D*)gr_p->Clone();
	sum->Add(gr_cp);
	diff->Add(gr_cp,-1);
	diff->Divide(sum);
	diff->Draw();

	c1->cd(2);
	gr1_p->Add(gr1_r,-0.5);
	gr1_p->Scale(scaleperp);
	gr1_p->SetTitle("The run is para");
	TH1D *sum1 = (TH1D*)gr1_p->Clone();
	TH1D *diff1 = (TH1D*)gr1_p->Clone();
	sum1->Add(gr_cp);
	diff1->Add(gr_cp,-1);
	diff1->Divide(sum1);
	diff1->Draw();

	
//	c1->Print("~/AcquRoot/4v3/2009.09.08/acqu/root/pdfFiles/perp_para_error.pdf");
}
