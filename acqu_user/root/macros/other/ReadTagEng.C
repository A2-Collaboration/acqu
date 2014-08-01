// Note that you must define the following structure with the proper
// variables in it before you can use this!
//
//	typedef struct {
//
//		Int_t egamma;
//		Double_t energy;
//		Double_t denergy;
//
//	} TData;
//
//	TData data[352];

void ReadTagEng( TString file)
{
	UInt_t i;
	Double_t eff, deff;

	ifstream inFile( file);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << file;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		File >> i >> eff >> deff;
		data[i].energy = eff;
		data[i].denergy = deff/2;
		data[i].egamma = (int)(eff + 0.5);
	}
	inFile.close();
}

void ReadSubt( TString file)
{
	UInt_t i;
	Double_t eff, deff;

	ifstream inFile( file);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << file;
		cout << endl;
		break;
	}
	while( !inFile.eof()) {
		inFile >> i >> eff;
		data[i].pa_ratio = eff;
	}
	inFile.close();
}

void ReadDetEff( TString file)
{
	UInt_t i;
	Double_t eff, deff;

	for ( i = 0; i <= 351; i++) {
		data[i].edet = 0.8;
		data[i].d_edet = 0.01;
	}
	ifstream inFile( file);
	if ( !inFile.is_open()) {
		cout << "Error opening file ";
		cout << file;
		cout << ".  Using constant efficiency!";
		cout << endl;
	}
	else {
		while( !inFile.eof()) {
			inFile >> i >> eff >> deff;
			data[i].edet = eff;
			data[i].d_edet = deff;
		}
	}
	inFile.close();
}

void ReadTagEff( TString file)
{
	UInt_t i;
	Double_t eff, deff;

	ifstream inFile( file);
	while( !inFile.eof()) {
		inFile >> i >> eff >> deff;
		if ( eff != 0) {
			data[i].etag = eff;
			data[i].d_etag = deff;
		}
		else {
			data[i].etag = -1;
			data[i].d_etag = -1;
		}
	}
	inFile.close();
}
