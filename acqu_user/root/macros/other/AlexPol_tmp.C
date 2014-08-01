// This macro is designed to do the things we couldn't in TA2Tagger.
// Taking the output file from TA2Tagger (run#, dump#, tchan, peakchan, scaler)
// we will do a weighted polarization analysis.

void AlexPol( TString choice)
{
	Double_t Eg[352], Edge[352], PolChan[81][352];
	TString tag, string, filename;
	Double_t content, junk;
	Double_t convertChan[352], convertEg[352];
	UInt_t i, j, k;

	filename = "/home/alex/AcquRoot/4v4/2010.04.22/acqu/xs/tageng855.dat";

	ifstream inFile( filename);
	if ( !inFile.is_open())
	{
		cout << "Error opening file ";
		cout << filename;
		cout << endl;
		break;
	}

	i = 0;

	while( !inFile.eof()) 
	{
		inFile >> j >> content >> junk;
		convertChan[i] = j;
		convertEg[i] = content;
		i++;
	}
	inFile.close();

	if ( choice == "para") filename = "/home/alex/LinPolPara.dat";
	else if ( choice == "perp") filename = "/home/alex/LinPolPerp.dat";
	
	ifstream inFile( filename);
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

	Double_t x, newcontent;

	while( !inFile.eof()) 
	{
		inFile >> tag >> content;

		if (tag == "Channel:")
		{
			Eg[i] = content;
			//printf( " %3d  %lf\n", i, Eg[i]);
			i++;
		}
		if (tag == "Edge:")
		{
			for(Int_t index = 0; content < convertEg[index]; index++)
			{
				x = (content - convertEg[index+1])/(convertEg[index] - convertEg[index+1]);
				newcontent = x*convertChan[index] + (1-x)*convertChan[index+1];
			}
			Edge[j] = newcontent;
			//printf( "%3d  %lf\n", j, Edge[j]);
			j++;
			k = 0;
		}
		if (tag == "Pol:")
		{
			PolChan[j-1][k] = content;
			//printf( "%3d  %lf\n", k, PolChan[j-1][k]);
			k++;
		}
	}
	inFile.close();

	if ( choice == "para") filename = "/home/alex/ALEXPARA.dat";
	else if ( choice == "perp") filename = "/home/alex/ALEXPERP.dat";
	
	ifstream inFile( filename);
	if ( !inFile.is_open()) 
	{
		cout << "Error opening file ";
		cout << filename;
		cout << endl;
		break;
	}

	Double_t run, dump, edge, edgeenergy, polpeak, scalers, scallow, scalhigh,
				edgepos, polarization;
	UInt_t tchan;
	Double_t numerator[352], denominator[352], AlexPol[352];

//	c1 = new TCanvas("c1", "c1", 50, 10, 500, 700);
//	c1->Divide(3,4);
	h1 = new TH1D("h1", "h1" , 1000, 0, 1000);
	h2 = new TH1D("h2", "h2" , 1000, 0, 1000);
	h3 = new TH1D("h3", "h3" , 1000, 0, 1000);
	h4 = new TH1D("h4", "h4" , 1000, 0, 1000);
	h5 = new TH1D("h5", "h5" , 1000, 0, 1000);
	h6 = new TH1D("h6", "h6" , 1000, 0, 1000);
	h7 = new TH1D("h7", "h7" , 1000, 0, 1000);
	h8 = new TH1D("h8", "h8" , 1000, 0, 1000);
	h9 = new TH1D("h9", "h9" , 1000, 0, 1000);
	h10 = new TH1D("h10", "h10" , 1000, 0, 1000);
	h11 = new TH1D("h11", "h11" , 1000, 0, 1000);

	while( !inFile.eof()) 
	{
		inFile >> run >> dump >> tchan >> edge >> scalers;

		edgeenergy = 351 - edge;

		for( UInt_t cycle=0; cycle<80; cycle++)
		{
			if((edgeenergy<=Edge[cycle+1])&&(edgeenergy>Edge[cycle]))
			{
			   Double_t interp = edgeenergy - Edge[cycle]; //interpolate where the actual edge is relative to our 2 limits
			   Double_t diff = Edge[cycle+1] - Edge[cycle];
			   scalhigh = interp/diff; //scaling factor for the upper limit
			   scallow = 1 - scalhigh; //scaling factor for the lower limit
				edgepos = cycle;
			}
			else if(edgeenergy<Edge[0])
			printf("Bad fitting...check file!\n");
		}

		//unweightedPol[tchan-294] = polarization;
/*		if(tchan == 294)
		{
			polarization = scalhigh*PolChan[edgepos][294] + scallow*PolChan[edgepos+1][294];
			h1->Fill(dump,polarization);
			c1->cd(1);
			h1->Draw();
		}
		if(tchan == 295)
		{
			polarization = scalhigh*PolChan[edgepos][295] + scallow*PolChan[edgepos+1][295];
			h2->Fill(dump,polarization);
			c1->cd(2);
			h2->Draw();
		}
		if(tchan == 296)
		{
			polarization = scalhigh*PolChan[edgepos][296] + scallow*PolChan[edgepos+1][296];
			h3->Fill(dump,polarization);
			c1->cd(3);
			h3->Draw();
		}
		if(tchan == 297)
		{
			polarization = scalhigh*PolChan[edgepos][297] + scallow*PolChan[edgepos+1][297];
			h4->Fill(dump,polarization);
			c1->cd(4);
			h4->Draw();
		}
		if(tchan == 298)
		{
			polarization = scalhigh*PolChan[edgepos][298] + scallow*PolChan[edgepos+1][298];
			h5->Fill(dump,polarization);
			c1->cd(5);
			h5->Draw();
		}
		if(tchan == 299)
		{
			polarization = scalhigh*PolChan[edgepos][299] + scallow*PolChan[edgepos+1][299];
			h6->Fill(dump,polarization);
			c1->cd(6);
			h6->Draw();
		}
		if(tchan == 300)
		{
			polarization = scalhigh*PolChan[edgepos][300] + scallow*PolChan[edgepos+1][300];
			h7->Fill(dump,polarization);
			c1->cd(7);
			h7->Draw();
		}
		if(tchan == 301)
		{
			polarization = scalhigh*PolChan[edgepos][301] + scallow*PolChan[edgepos+1][301];
			h8->Fill(dump,polarization);
			c1->cd(8);
			h8->Draw();
		}
		if(tchan == 302)
		{
			polarization = scalhigh*PolChan[edgepos][302] + scallow*PolChan[edgepos+1][302];
			h9->Fill(dump,polarization);
			c1->cd(9);
			h9->Draw();
		}
		if(tchan == 303)
		{
			polarization = scalhigh*PolChan[edgepos][303] + scallow*PolChan[edgepos+1][303];
			h10->Fill(dump,polarization);
			c1->cd(10);
			h10->Draw();
		}
		if(tchan == 304)
		{
			polarization = scalhigh*PolChan[edgepos][304] + scallow*PolChan[edgepos+1][304];
			h11->Fill(dump,polarization);
			c1->cd(11);
			h11->Draw();
		}*/
		//cout << polarization << endl;

		//numerator[tchan] += (scalhigh*PolChan[edgepos+1][tchan] + scallow*PolChan[edgepos][tchan]);
		numerator[tchan] += (scalhigh*PolChan[edgepos][tchan] + scallow*PolChan[edgepos+1][tchan])*scalers;
		denominator[tchan] += scalers;
		AlexPol[tchan] = numerator[tchan]/denominator[tchan];

	}
	inFile.close();

	printf("%s %s\n","Channel","Weighted Polarization");
	for(tchan=245;tchan<261;tchan++)
	printf("%3d\t %3.3f\n",tchan, AlexPol[tchan]);
	//printf("%3d\t %3.3f\n",tchan, denominator[tchan]);

	for(tchan=290;tchan<309;tchan++)
	printf("%3d\t %3.3f\n",tchan, AlexPol[tchan]);
	//printf("%3d\t %3.3f\n",tchan, denominator[tchan]);

}
