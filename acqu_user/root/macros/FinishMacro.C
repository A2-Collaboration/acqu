//#include <iostream>
//#include <string>
//#include <stdio.h>
#include <time.h>

// Get current date/time
/*const char currentDateTime()
{
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	// http://www.cplusplus.com/reference/clibrary/ctime/strftime/
	strftime(buf, sizeof(buf), "%d.%m.%Y", &tstruct);

	return buf;
}*/

void saveAnalysis()
{
	time_t now = time(0);
	struct tm tstruct;
	char file[80];
	tstruct = *localtime(&now);
	// http://www.cplusplus.com/reference/clibrary/ctime/strftime/
	strftime(file, sizeof(file), "Analysis_%Y-%m-%d_autosave.root", &tstruct);

	printf("Executing End-of-Run macro, saving analysis to %s\n", file);
	gUAN->SaveAll(file);
  	printf("Done.\n", file);
	//std::cout << "currentDateTime()=" << currentDateTime() << std::endl;  <-- doesn't work in separate method ?!
	//getchar();  // wait for keyboard input
}

void FinishMacro(Char_t* file = "ARHistograms.root")
{
	printf("\nEnd-of-Run macro executing:\n");
	TFile f(file, "recreate");
	if(!f) {
		printf("Open file %s for histogram save FAILED!!\n",file);
		return;
	}
	gROOT->GetList()->Write();
	f.Close();
  	printf("All histograms saved to %s\n\n", file);
}
