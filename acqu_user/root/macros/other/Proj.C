gROOT->Reset();

TFile f( "histograms/ARH_test.root");

void Proj( Char_t orient[8], Char_t time[1], Char_t pr[1])
{
	Char_t orient[8], name[128];

	TCanvas *c1 = new TCanvas ( "c1", "Projection", 20, 20, 500, 700);
	c1->Divide(1,3);

	c1->cd(1);
	sprintf( name, "TaggHitsCut4%s%s_v_ThetaGG4%s%s_v_PhiPi0cm4%s%s", orient,
			time, orient, time, orient, time);
	h3 = (TH3D*)gROOT->FindObject( name);
//	h3->GetZaxis()->SetRange(21,284);
//	h3->GetYaxis()->SetRange(0,180);
//	h3->GetXaxis()->SetRange(-180,180);
	TH1D *h1 = h3->Project3D( pr);
	cout << h1->Integral() << endl;
	h1->Draw();

	c1->cd(2);
	if ( strcmp( pr, "x") == 0) sprintf( name,
			"THR_TaggHitsCut4%s%s_v_PhiPi0cm4%s%s", orient, time, orient, time);
	else if ( strcmp( pr, "y") == 0) sprintf( name,
			"THR_TaggHitsCut4%s%s_v_ThetaGG4%s%s", orient, time, orient, time);
	h2 = (TH2D*)gROOT->FindObject( name);
	TH1D *h4 = h2->ProjectionX();
//	TH1D *h4 = h2->ProjectionX("tmp", 21, 284);
	cout << h4->Integral() << endl;
	h4->Draw();

	c1->cd(3);
	if ( strcmp( pr, "x") == 0)
		sprintf( name, "THR_PhiPi0cm4%s%s", orient, time);
	else if ( strcmp( pr, "y") == 0)
		sprintf( name, "THR_ThetaGG4%s%s", orient, time);
	else if ( strcmp( pr, "z") == 0)
		sprintf( name, "THR_TaggHitsCut4%s%s", orient, time);
	h5 = (TH1D*)gROOT->FindObject( name);
	cout << h5->Integral() << endl;
	h5->Draw();
}

void Proj2( Char_t time[1], Char_t pr[1])
{
	Char_t name[128];

	TCanvas *c1 = new TCanvas ( "c1", "Projection", 20, 20, 500, 700);
	c1->Divide(1,3);

	c1->cd(1);
	sprintf( name, "TaggHitsCut%s_v_EmPi0%s_v_ThetaPi0cm%s", time, time, time);
	h3 = (TH3D*)gROOT->FindObject( name);
//	h3->GetZaxis()->SetRange(21,284);
//	h3->GetYaxis()->SetRange(0,180);
//	h3->GetXaxis()->SetRange(-180,180);
	TH1D *h1 = h3->Project3D( pr);
	cout << h1->Integral() << endl;
	h1->Draw();

	c1->cd(2);
	if ( strcmp( pr, "x") == 0) sprintf( name, "THR_TaggHitsCut%s_v_ThetaPi0cm%s", time, time);
	else if ( strcmp( pr, "y") == 0) sprintf( name, "THR_TaggHitsCut%s_v_EmPi0%s", time, time);
	h2 = (TH2D*)gROOT->FindObject( name);
	TH1D *h4 = h2->ProjectionX();
//	TH1D *h4 = h2->ProjectionX("tmp", 21, 284);
	cout << h4->Integral() << endl;
	h4->Draw();

	c1->cd(3);
	if ( strcmp( pr, "x") == 0) sprintf( name, "THR_ThetaPi0cm%s", time);
	else if ( strcmp( pr, "y") == 0) sprintf( name, "THR_EmPi0%s", time);
	else if ( strcmp( pr, "z") == 0) sprintf( name, "THR_TaggHitsCut%s", time);
	h5 = (TH1D*)gROOT->FindObject( name);
	cout << h5->Integral() << endl;
	h5->Draw();
}
