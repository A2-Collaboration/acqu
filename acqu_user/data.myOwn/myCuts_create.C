{
	TFile *f = new TFile("myCuts.root","RECREATE");

	TCutG *cutg = new TCutG("eCut",8);
	cutg->SetVarX("TA2CrystalBall");
	cutg->SetVarY("");
	cutg->SetTitle("AcquRoot 2D polygon electron cut");
	cutg->SetFillColor(1);
	cutg->SetLineColor(4);
	cutg->SetLineWidth(2);
	cutg->SetPoint(0,10,0.15);
	cutg->SetPoint(1,10,1);
	cutg->SetPoint(2,100,2);
	cutg->SetPoint(3,950,1.7);
	cutg->SetPoint(4,1070,1.1);
	cutg->SetPoint(5,1070,0.6);
	cutg->SetPoint(6,950,0.2);
	cutg->SetPoint(7,10,0.15);
	//cutg->Draw("SAME");
	//cutg->Print();
	cutg->Write();

	TCutG *cutg2 = new TCutG("pCut",8);
	cutg2->SetVarX("TA2CrystalBall");
	cutg2->SetVarY("");
	cutg2->SetTitle("AcquRoot 2D polygon proton cut");
	cutg2->SetFillColor(1);
	cutg2->SetLineColor(2);
	cutg2->SetLineWidth(2);
	cutg2->SetPoint(0,80,7.7);
	cutg2->SetPoint(1,250,7.4);
	cutg2->SetPoint(2,750,4.2);
	cutg2->SetPoint(3,750,2.1);
	cutg2->SetPoint(4,300,2.1);
	cutg2->SetPoint(5,200,2.4);
	cutg2->SetPoint(6,70,2.4);
	cutg2->SetPoint(7,80,7.7);
	//cutg2->Draw("SAME");
	cutg2->Write();

	TCutG *cutb = new TCutG("balanceCut",12);
	cutb->SetVarX("TA2CrystalBall");
	cutb->SetVarY("");
	cutb->SetTitle("AcquRoot 2D polygon balance cut");
	cutb->SetFillColor(1);
	cutb->SetLineColor(3);
	cutb->SetLineWidth(2);
	cutb->SetPoint(0,30,60);
	cutb->SetPoint(1,-40,3);
	cutb->SetPoint(2,-132,-111);
	cutb->SetPoint(3,-215,-264);
	cutb->SetPoint(4,-274,-378);
	cutb->SetPoint(5,-278,-429);
	cutb->SetPoint(6,-258,-422);
	cutb->SetPoint(7,-199,-353);
	cutb->SetPoint(8,-120,-238);
	cutb->SetPoint(9,-33,-124);
	cutb->SetPoint(10,42,-9);
	cutb->SetPoint(11,30,60);
	//cutg->Draw("SAME");
	cutb->Write();

	f->Close();
}
