void graphme()
{
	TString filename;
	UInt_t i, j, k, count;
	Double_t junk;
	Double_t parax[67], parae[67], perpx[79], perpe[79], run[146], arrayx[146], arraye[146], parar[67], perpr[79], zero[79], paracheck[67];

	filename = "xs/Para.tmp";
	
	ifstream inFile(filename);
	if ( !inFile.is_open())
	{
		cout << "Error opening file ";
		cout << filename;
		cout << endl;
		break;
	}

	while( !inFile.eof()) 
	{
		inFile >> paracheck[i];
		i++;
	}
	inFile.close();


	filename = "xs/Alex.tmp";

	ifstream inFile(filename);
	if ( !inFile.is_open())
	{
		cout << "Error opening file ";
		cout << filename;
		cout << endl;
		break;
	}

	i = 0;
	j = 0;
	k = 0;
	count = 0;

	while( !inFile.eof()) 
	{
		inFile >> run[i] >> junk >> junk >> junk >> arrayx[i] >> arraye[i];
		i++;
	}
	inFile.close();

	for(i=0;i<146;i++)
	{
		flag = 0;

		for(count=0;count<67;count++)
		{
			if(run[i] == paracheck[count])
			{
				flag = 1;
				break;
			}
		}

		if(flag==0)
		{
			perpx[k] = arrayx[i];
			perpe[k] = arraye[i];
			perpr[k] = run[i];
			zero[k] = 0;
			k++;
		}

		if(flag==1)
		{
			parax[j] = arrayx[i];
			parae[j] = arraye[i];
			parar[j] = run[i];
			j++;		
		}
	}

	TGraphErrors *gr1 = new TGraphErrors(79, perpr, perpx, zero, perpe);
	TGraphErrors *gr2 = new TGraphErrors(67, parar, parax, zero, parae);

	c1 = new TCanvas("c1", "c1", 50, 10, 500, 700);

	gr2->SetMarkerStyle(21);
	gr2->SetMarkerColor(kRed);
	gr1->SetMarkerStyle(22);
	gr1->SetMarkerColor(kBlue);
	
	TMultiGraph *mg = new TMultiGraph();
	mg->SetTitle("Yield/Scalers*Factor with error");
	mg->Add(gr1);
	mg->Add(gr2);
	mg->Draw("AP");
	
	TLegend *tleg = new TLegend(0.6,0.25,0.8,0.35);
	tleg->SetFillColor(0);
	tleg->SetBorderSize(0);
	tleg->AddEntry(gr1,"Perp", "p");
	tleg->AddEntry(gr2,"Para", "p");
	tleg->Draw();

}

void TEpara()
{
	TString filename = "/home/alex/AcquRoot/4v4/2010.04.22/acqu/xs/TE_para_246.tmp";
	Double_t TE[10], error[10];
	Double_t run[10], zero[10];
	UInt_t i = 0;

	ifstream inFile(filename);
	if ( !inFile.is_open())
	{
		cout << "Error opening file ";
		cout << filename;
		cout << endl;
		break;
	}

	/*while( !inFile.eof()) 
	{
		inFile >> run[i] >> TE[i] >> error[i];
		zero[i] = 0;
		i++;
	}*/
	inFile.close();

//	filename = "xs/Alex.tmp";
	filename = "/home/alex/AcquRoot/4v4/2010.04.22/acqu/xs/Alex.tmp";

	i = 0;
	Double_t allrun, arrayx, arraye, junk;
	Double_t realrun[10], realTE[10], realerror[10];

	ifstream inFile(filename);
	if ( !inFile.is_open())
	{
		cout << "Error opening file ";
		cout << filename;
		cout << endl;
		break;
	}

	while( !inFile.eof()) 
	{
		inFile >> allrun >> junk >> junk >> junk >> arrayx >> arraye;
		if((allrun == 20880)||(allrun == 20904)||(allrun == 20948)||(allrun == 20962)||(allrun == 21052))
		{
			realrun[i] = allrun;
			realTE[i] = arrayx;
			realerror[i] = arraye;
			i++;
		}
	}
	inFile.close();

	/*c2 = new TCanvas("c2", "c2", 50, 10, 500, 700);
	gr3 = new TGraphErrors(5, run, TE, zero, error);
	gr4 = new TGraphErrors(5, realrun, realTE, zero, realerror);
	gr4->SetMarkerStyle(21);
	gr4->SetMarkerColor(kRed);
	gr3->SetMarkerStyle(22);
	gr3->SetMarkerColor(kBlue);
	
	TMultiGraph *mg = new TMultiGraph();
	mg->SetTitle("Yield/Scalers*Factor with error");
	mg->Add(gr3);
	mg->Add(gr4);
	mg->Draw("AP");*/
}
