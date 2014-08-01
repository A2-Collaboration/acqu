// Note that you must define the following structure with the proper
// variables in it before you can use this:
//
//	typedef struct {
//
//		Int_t egamma;
//		Double_t energy;
//		Double_t denergy;
//		Double_t etag;
//		Double_t d_etag;
//		Double_t edet[10];
//		Double_t d_edet[10];
//		Double_t pol[4];
//		Double_t pa_ratio;
//		Double_t xs;
//		Double_t dxs;
//	
//	} TData;
//	TData tcd[352];
//

//
// Read in channel tagger energies.
//
void ReadTagEng( TString file)
{
	UInt_t i;
	Double_t eff, deff;

	cout << "Reading in parameters from " << file << endl;

	ifstream inFile( file);
	if ( !inFile.is_open()) {
		cout << "ERROR: file ";
		cout << file;
		cout << " does not exist.";
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> i >> eff >> deff;
		tcd[i].energy = eff;
		tcd[i].denergy = deff/2;
		tcd[i].egamma = (int)(eff + 0.5);
	}
	inFile.close();
}

//
// Read in channel random subtraction factors.
//
void ReadSubt( TString file)
{
	UInt_t i;
	Double_t eff, deff;

	cout << "Reading in parameters from " << file << endl;

	ifstream inFile( file);
	if ( !inFile.is_open()) {
		cout << "ERROR: file ";
		cout << file;
		cout << " does not exist.";
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> i >> eff;
		tcd[i].pa_ratio = eff;
	}
	inFile.close();
}

//
// Read in channel detection efficiencies.
//
void ReadDetEff( TString exp, Bool_t th_flag = kTRUE)
{
	UInt_t i, j, n;
	Double_t eff, deff;
	TString file, cth, junk;

	if ( th_flag == kTRUE) {
		n = 10;
		cth = "";
	}
	else {
		n = 21;
		cth = "cth_";
	}

	cout << "Reading in detection efficiency parameters" << endl;

	for ( i = 0; i <= 351; i++) {
		for ( j = 0; j < n; j++) {
			if ( exp == "pi0") {
				tcd[i].edet[j] = 0.8;
				tcd[i].d_edet[j] = 0.01;
			}
			else if ( exp == "compton") {
				tcd[i].edet[j] = 0.2;
				tcd[i].d_edet[j] = 0.01;
			}
		}
		file = Form( "xs/%s/eff/DetEff_%schan_%d.out", exp.Data(), cth.Data(),
				i);
		ifstream inFile( file);
		if ( !inFile.is_open()) {
			cout << "Efficiency for channel ";
			cout << i;
			cout << " does not exist.  Using constant efficiency.";
			cout << endl;
		}
		else {
			for ( j = 0; j < n; j++) {
				inFile >> junk >> eff >> deff;
				tcd[i].edet[j] = eff;
				tcd[i].d_edet[j] = deff;
			}
		}
		inFile.close();
	}
}

//
// Read in weighted channel polarization.
//
// The columns are:
//
// 	chan / full-perp / full-para / empty-perp / empty-para
//
// so that the tags are:
//
// 	full-perp	0
// 	full-para	1
// 	empty-perp	2
// 	empty-para	3
//
void ReadPolData( TString file)
{
	UInt_t i, j;
	Double_t p[4];
	TString name;

	cout << "Reading in parameters from " << file << endl;

	// Initialize all channels
	for ( i = 0; i <= 351; i++)
		for ( j = 0; j < 4; tcd[i].pol[j++] = 0.0);

	ifstream inFile( file);
	if ( !inFile.is_open()) {
		cout << "ERROR: file ";
		cout << file;
		cout << " does not exist.";
		cout << endl;
		exit( -1);
	}
	while( !inFile.eof()) {
		name.ReadLine( inFile);
		if ( name[0] != '#') {
			sscanf( name.Data(), "%d%lf%lf%lf%lf", &i, &p[0], &p[1], &p[2], &p[3]);
			for ( j = 0; j < 4; j++) tcd[i].pol[j] = p[j];
		}
	}
	inFile.close();
}

//
// Read in channel tagging efficiencies.
//
void ReadTagEff( TString file, Bool_t errflag = kTRUE)
{
	UInt_t i;
	Double_t eff, deff;

	cout << "Reading in parameters from " << file << endl;

	ifstream inFile( file);
	if ( !inFile.is_open()) {
		cout << "ERROR: file ";
		cout << file;
		cout << " does not exist.";
		cout << endl;
		exit( -1);
	}
	while( !inFile.eof()) {
		if ( errflag == kTRUE) {
			inFile >> i >> eff >> deff;
			if ( eff != 0) {
				tcd[i].etag = eff;
				tcd[i].d_etag = deff;
			}
			else {
				tcd[i].etag = -1;
				tcd[i].d_etag = -1;
			}
		}
		else {
			inFile >> i >> eff;
			if ( eff != 0) tcd[i].etag = eff;
			else tcd[i].etag = -1;
			tcd[i].d_etag = -1;
		}
	}
	inFile.close();
}

//
// Read in Sergey Prakhov bin information.
//
void ReadSergeyBins( TString file)
{
	UInt_t bin, c_lo, c_hi;
	Double_t eg, deg;
	TString name;

	cout << "Reading in parameters from " << file << endl;

	ifstream inFile( file);
	if ( !inFile.is_open()) {
		cout << "ERROR: file ";
		cout << file;
		cout << " does not exist.";
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		name.ReadLine( inFile);
		if ( name[0] != '#') {
			sscanf( name.Data(), "%d%d%d%lf%lf", &bin, &c_lo, &c_hi, &eg, &deg);
			spbin[bin].ch_lo = c_lo;
			spbin[bin].ch_hi = c_hi;
			spbin[bin].eg = eg;
			spbin[bin].deg = deg;
		}
	}
	inFile.close();
}

//
// Read in weighted tagging efficiencies.
//
/*
void ReadTagEffW( TString file)
{
	UInt_t i, j;
	Double_t eff[4], deff[4];
	TString name;

	cout << "Reading in parameters from " << file << endl;

	// Initialize all channels
	for ( i = 0; i <= 351; i++)
		for ( j = 0; j < 4; j++) {
			tcd[i].etag[j] = 0.0;
			tcd[i].d_etag[j] = 0.0;
		}

	ifstream inFile( file);
	if ( !inFile.is_open()) {
		cout << "ERROR: file ";
		cout << file;
		cout << " does not exist.";
		cout << endl;
		exit( -1);
	}
	while( !inFile.eof()) {
		name.ReadLine( inFile);
		if ( name[0] != '#') {
			sscanf( name.Data(), "%d%lf%lf%lf%lf%lf%lf%lf%lf", &i,
					&eff[0], &deff[0], &eff[1], &deff[1], &eff[2], &deff[2],
					&eff[3], &deff[3]);
			for ( j = 0; j < 4; j++) {
				tcd[i].etag[j] = eff[j];
				tcd[i].d_etag[j] = deff[j];
			}
		}
	}
	inFile.close();
}
*/
