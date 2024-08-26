#include "call_libraries.h" // call libraries from ROOT and C++
#include "uiclogo.h"	    // call UIC logo and initialization
#include "CATree.h" 	    // call re-cluster for WTA axis: see https://github.com/FHead/PhysicsMiniProjects/tree/master/JetSmallSystem/24622_Recluster
#include "ntrkoff.h"        // get Ntrk offline
#define pimass 0.1396


std::map<unsigned long long, int> runLumiEvtToEntryMap;
unsigned long long keyFromRunLumiEvent(UInt_t run, UInt_t lumi, ULong64_t event);

/*
Main skim pPb data and MC

Written by Dener Lemos (dener.lemos@cern.ch)

--> Arguments
input_file: text file with a list of root input files: Forest or Skims
ouputfile: just a counting number to run on Condor
isMC: 0 for false --> data and > 0 for true --> MC
ntrkoff: 0 for no cut/selection or MC, 1 for MB [10,185], 2 for HM PD 1 to 6 [185,250] and 3 for HM PD 7 [250, inf]
largersample: 0 for p -> + eta and 2 for p -> - eta
input_ZDCfile: Input ZDC files
*/
void pPbSkim(TString input_file, TString ouputfile, int isMC, int ntrkoff, int largersample, TString input_ZDCfile){

	bool is_MC; if(isMC == 0){is_MC = false;}else{is_MC = true;}

	float jetrawptmin = 10.0;
	if(is_MC){jetrawptmin = 0.0;}
	bool storesoftdrop = false;
	bool storetracks = true;
	bool storepfcand = true;
	
	if(storepfcand) storesoftdrop = false; 

	TString outputFileName;
	outputFileName = Form("%s",ouputfile.Data());

	clock_t sec_start, sec_end;
	sec_start = clock(); // start timing measurement

	TDatime* date = new TDatime();

	printwelcome(true); // welcome message

	print_start(); // start timing print

	// Read the input file(s)
	fstream inputfile;
	inputfile.open(Form("%s",input_file.Data()), ios::in);
	if(!inputfile.is_open()){cout << "List of input files not founded!" << endl; return;}{cout << "List of input files founded! --> " << input_file.Data() << endl;}

	// Make a chain and a vector of file names
	std::vector<TString> file_name_vector;
	string file_chain;
	while(getline(inputfile, file_chain)){file_name_vector.push_back(Form("root://osg-se.sprace.org.br/%s",file_chain.c_str()));}
	inputfile.close();
	// Maximum size of arrays
	const Int_t nMaxJet = 500;				// Maximum number of jets in an event
	const Int_t nMaxTrack = 2000;		// Maximum number of tracks in an event
	
	// Define trees to be read from the files
	const int nJetTrees = 4;
	TChain *heavyIonTree = new TChain("hiEvtAnalyzer/HiTree");
	TChain *hltTree = new TChain("hltanalysis/HltTree");
	TChain *skimTree = new TChain("skimanalysis/HltTree");
	TChain *jetTree[nJetTrees];
	jetTree[0] = new TChain("ak4CaloJetAnalyzer/t");
	jetTree[1] = new TChain("ak4PFJetAnalyzer/t");
	jetTree[2] = new TChain("akCs4PFJetAnalyzer/t");
	jetTree[3] = new TChain("ak3PFJetAnalyzer/t");
	TChain *RhoTree = new TChain("hiFJRhoAnalyzer/t");
	TChain *trackTree = new TChain("ppTrack/trackTree");
	TChain *genTrackTree;
	if(is_MC){genTrackTree = new TChain("HiGenParticleAna/hi");}
	TChain *particleFlowCandidateTree = new TChain("pfcandAnalyzer/pfTree");
	TChain *checkFlatteningTree = new TChain("checkflattening/tree");

	// add all the trees to the chain
	for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++){
		//TFile testfile(*listIterator, "READ"); 
		//if(testfile.IsZombie() || testfile.TestBit(TFile::kRecovered)) continue;
		cout << "Adding file " << *listIterator << " to the chains" << endl;
		hltTree->Add(*listIterator);
		trackTree->Add(*listIterator);
		heavyIonTree->Add(*listIterator);
		for(int iJetType = 0; iJetType < nJetTrees; iJetType++){jetTree[iJetType]->Add(*listIterator);}
		skimTree->Add(*listIterator);
		RhoTree->Add(*listIterator);
		if(is_MC){genTrackTree->Add(*listIterator);}
		particleFlowCandidateTree->Add(*listIterator);
		checkFlatteningTree->Add(*listIterator);
	}

	// Read the input ZDC file(s)
	fstream inputfileZDC;
	inputfileZDC.open(Form("%s",input_ZDCfile.Data()), ios::in);
	if(!inputfileZDC.is_open()){cout << "List of ZDC input files not founded!" << endl; return;}{cout << "List of ZDC input files founded! --> " << input_ZDCfile.Data() << endl;}
	// Make a chain and a vector of file names
	std::vector<TString> file_name_vectorZDC;
	string file_chainZDC;
	while(getline(inputfileZDC, file_chainZDC)){file_name_vectorZDC.push_back(Form("%s",file_chainZDC.c_str()));}
	inputfileZDC.close();
	TChain *MainZDCTree = new TChain("demo/zdc_tree");
	// add all the trees to the chain
	for (std::vector<TString>::iterator listIterator = file_name_vectorZDC.begin(); listIterator != file_name_vectorZDC.end(); listIterator++){ cout << "Adding file " << *listIterator << " to the chains" << endl; MainZDCTree->Add(*listIterator);}

	// ZDC Branches
	TBranch *ZDC_runBranch;	 // Branch for run
	TBranch *ZDC_evtBranch;	 // Branch for event
	TBranch *ZDC_lumiBranch; // Branch for lumi
	TBranch *ZDC_SumNBranch; // Branch for Negative ZDC
	TBranch *ZDC_SumPBranch; // Branch for Positive ZDC
	UInt_t ZDC_run;			 // Run number
	ULong64_t ZDC_evt;			 // Event number
	UInt_t ZDC_lumi;			 // Luminosity block
	Float_t ZDC_SumN;		 // Negative ZDC
	Float_t ZDC_SumP;		 // Positive ZDC

	//Main tree (event info)
	MainZDCTree->SetBranchStatus("*",0);
	MainZDCTree->SetBranchStatus("runNumber",1);
	MainZDCTree->SetBranchAddress("runNumber",&ZDC_run,&ZDC_runBranch);
	MainZDCTree->SetBranchStatus("evNumber",1);
	MainZDCTree->SetBranchAddress("evNumber",&ZDC_evt,&ZDC_evtBranch);
	MainZDCTree->SetBranchStatus("LumiSection",1);
	MainZDCTree->SetBranchAddress("LumiSection",&ZDC_lumi,&ZDC_lumiBranch);
	MainZDCTree->SetBranchStatus("ZDCSumMinus",1);
	MainZDCTree->SetBranchAddress("ZDCSumMinus",&ZDC_SumN,&ZDC_SumNBranch);
	MainZDCTree->SetBranchStatus("ZDCSumPlus",1);
	MainZDCTree->SetBranchAddress("ZDCSumPlus",&ZDC_SumP,&ZDC_SumPBranch);

	
	// All the branches and leaves come in arrat of two, one for input and one for output
	
	// Branches for heavy ion tree
	TBranch *runBranch;						 	// Branch for run
	TBranch *eventBranch;					 	// Branch for event
	TBranch *lumiBranch;						// Branch for lumi
	TBranch *hiVzBranch;						// Branch for vertex z-position
	TBranch *hiVxBranch;						// Branch for vertex x-position
	TBranch *hiVyBranch;						// Branch for vertex y-position
	TBranch *hiHFplusBranch;					// Branch for HF+ energy deposity
	TBranch *hiHFminusBranch;					// Branch for HF- energy deposity
	TBranch *hiHFplusEta4Branch;				// Branch for HF+ energy deposity for |eta| > 4
	TBranch *hiHFminusEta4Branch;				// Branch for HF- energy deposity for |eta| > 4
	TBranch *hiZDCplusBranch;					// Branch for ZDC+ energy deposity
	TBranch *hiZDCminusBranch;					// Branch for ZDC- energy deposity

	TBranch *ptHatBranch;						// Branch for pT hat
	TBranch *eventWeightBranch;		 			// Branch for pthat weight for MC

	// Leaves for heavy ion tree
	UInt_t run;					 // Run number
	ULong64_t event;			 // Event number
	UInt_t lumi;				 // Luminosity block
	Float_t vertexZ;			 // Vertex z-position
	Float_t vertexX;			 // Vertex x-position
	Float_t vertexY;			 // Vertex y-position
	Float_t hiHFplus;			 // transverse energy sum of HF+ tower;
	Float_t hiHFminus;			 // transverse energy sum of HF- tower;
	Float_t hiHFplusEta4;		 // transverse energy sum of HF+ tower for |eta| > 4;
	Float_t hiHFminusEta4;		 // transverse energy sum of HF- tower for |eta| > 4;
	Int_t Ntroff;		 		 // Multiplicity --> Ntrkoffline

	Float_t hiZDCplus;			 // energy deposit in ZDC+;
	Float_t hiZDCminus;			 // energy deposit in ZDC-;

	Float_t hi_FRG;			 // FRG (Forward Rapidity Gap) variable for UPC;
	Float_t hi_FRG_noNsel;	 // FRG (Forward Rapidity Gap) variable for UPC without N select in 2.5<|eta|<3.0;
	Float_t hi_BRG;			 // BRG (Backward Rapidity Gap) variable for UPC
	Float_t hi_BRG_noNsel; 	 // BRG (Backward Rapidity Gap) variable for UPC without N select in 2.5<|eta|<3.0;
	double const pfE[20] = {13.4, 16.4, 15.3, 16.9, 13.4, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 13.4, 15.9, 31.7, 17.1, 13.6};

	Float_t ptHat;				 // pT hat
	Float_t eventWeight;			 // jet weight in the tree

	// Branches for EP tree
	TBranch *eventPlaneAngleBranch;						// Branch for event plane angles
	TBranch *eventPLaneQBranch;							// Branch for Q-vector magnitude in an event plane
	TBranch *eventPlaneQxBranch;						// Branch for Q-vector x-component in an event plane
	TBranch *eventPlaneQyBranch;						// Branch for Q-vector y-component in an event plane
	TBranch *eventPlaneMultiplicityBranch;	 			// Branch for event plane multiplicity
	
	const int numberofEPleaves = 182;				    		// Event plane leaves	
	// Name of EP in pPb
	TString EPNames = "HFm1/D:HFp1/D:trackmid1/D:trackm1/D:trackp1/D:trackm122/D:trackm118/D:trackm114/D:trackm110/D:trackm106/D:trackm102/D:trackp102/D:trackp106/D:trackp110/D:trackp114/D:trackp118/D:trackp122/D:trackmid1mc/D:trackm1mc/D:trackp1mc/D:trackm122mc/D:trackm118mc/D:trackm114mc/D:trackm110mc/D:trackm106mc/D:trackm102mc/D:trackp102mc/D:trackp106mc/D:trackp110mc/D:trackp114mc/D:trackp118mc/D:trackp122mc/D:HFm1a/D:HFm1b/D:HFm1c/D:HFm1d/D:HFm1e/D:HFm1f/D:HFp1a/D:HFp1b/D:HFp1c/D:HFp1d/D:HFp1e/D:HFp1f/D:HFm2/D:HFp2/D:trackmid2/D:trackm2/D:trackp2/D:trackm222/D:trackm218/D:trackm214/D:trackm210/D:trackm206/D:trackm202/D:trackp202/D:trackp206/D:trackp210/D:trackp214/D:trackp218/D:trackp222/D:HFm2a/D:HFm2b/D:HFm2c/D:HFm2d/D:HFm2e/D:HFm2f/D:HFp2a/D:HFp2b/D:HFp2c/D:HFp2d/D:HFp2e/D:HFp2f/D:HFm3/D:HFp3/D:trackmid3/D:trackm3/D:trackp3/D:trackm322/D:trackm318/D:trackm314/D:trackm310/D:trackm306/D:trackm302/D:trackp302/D:trackp306/D:trackp310/D:trackp314/D:trackp318/D:trackp322/D:HFm3a/D:HFm3b/D:HFm3c/D:HFm3d/D:HFm3e/D:HFm3f/D:HFp3a/D:HFp3b/D:HFp3c/D:HFp3d/D:HFp3e/D:HFp3f/D:HFm4/D:HFp4/D:trackmid4/D:trackm4/D:trackp4/D:trackm422/D:trackm418/D:trackm414/D:trackm410/D:trackm406/D:trackm402/D:trackp402/D:trackp406/D:trackp410/D:trackp414/D:trackp418/D:trackp422/D:HFm4a/D:HFm4b/D:HFm4c/D:HFm43d/D:HFm4e/D:HFm4f/D:HFp4a/D:HFp4b/D:HFp4c/D:HFp4d/D:HFp4e/D:HFp4f/D:HFm5/D:HFp5/D:trackmid5/D:trackm5/D:trackp5/D:trackm522/D:trackm518/D:trackm514/D:trackm510/D:trackm506/D:trackm502/D:trackp502/D:trackp506/D:trackp510/D:trackp514/D:trackp518/D:trackp522/D:HFm6/D:HFp6/D:trackmid6/D:trackm6/D:trackp6/D:trackm622/D:trackm618/D:trackm614/D:trackm610/D:trackm606/D:trackm602/D:trackp602/D:trackp606/D:trackp610/D:trackp614/D:trackp618/D:trackp622/D:HFm7/D:HFp7/D:trackmid7/D:trackm7/D:trackp7/D:trackm722/D:trackm718/D:trackm714/D:trackm710/D:trackm706/D:trackm702/D:trackp702/D:trackp706/D:trackp710/D:trackp714/D:trackp718/D:trackp722/D"; 
	Double_t eventPlaneAngle[numberofEPleaves] = {0};			// Event plane angles
	Double_t eventPlaneQ[numberofEPleaves] = {0};				// Magnitude of Q-vector in event plane
	Double_t eventPlaneQx[numberofEPleaves] = {0};				// x-component of the Q-vector
	Double_t eventPlaneQy[numberofEPleaves] = {0};				// y-component of the Q-vector
	Double_t eventPlaneMultiplicity[numberofEPleaves] = {0};	// Particle multiplicity in an event plane

	// Branches for HLT tree
	// HLT
	TBranch *MB_FirstCollisionAfterAbortGapBranch;	 			// Branch for MB trigger
	TBranch *MB_ForSkimBranch;						 			// Branch for MB trigger
	TBranch *MB_ForExpressBranch;	 							// Branch for MB trigger
	TBranch *MB_part1Branch;	 								// Branch for MB trigger
	TBranch *MB_part2Branch;	 								// Branch for MB trigger
	TBranch *MB_part3Branch;	 								// Branch for MB trigger
	TBranch *MB_part4Branch;	 								// Branch for MB trigger
	TBranch *MB_part5Branch;	 								// Branch for MB trigger
	TBranch *MB_part6Branch;	 								// Branch for MB trigger
	TBranch *MB_part7Branch;	 								// Branch for MB trigger
	TBranch *MB_part8Branch;	 								// Branch for MB trigger
	TBranch *MB_part9Branch;	 								// Branch for MB trigger
	TBranch *MB_part10Branch;	 								// Branch for MB trigger
	TBranch *MB_part11Branch;	 								// Branch for MB trigger
	TBranch *MB_part12Branch;	 								// Branch for MB trigger
	TBranch *MB_part13Branch;	 								// Branch for MB trigger
	TBranch *MB_part14Branch;	 								// Branch for MB trigger
	TBranch *MB_part15Branch;	 								// Branch for MB trigger
	TBranch *MB_part16Branch;	 								// Branch for MB trigger
	TBranch *MB_part17Branch;	 								// Branch for MB trigger
	TBranch *MB_part18Branch;	 								// Branch for MB trigger
	TBranch *MB_part19Branch;	 								// Branch for MB trigger
	TBranch *MB_part20Branch;	 								// Branch for MB trigger

	TBranch *caloJetFilterBranch60;				 		// Branch for calo jet 60 filter bit
	TBranch *caloJetFilterBranch80;				 		// Branch for calo jet 80 filter bit
	TBranch *caloJetFilterBranch100;			 		// Branch for calo jet 100 filter bit
	TBranch *pfJetFilterBranch60;				 		// Branch for PF jet 60 filter bit
	TBranch *pfJetFilterBranch80;				 		// Branch for PF jet 80 filter bit
	TBranch *pfJetFilterBranch100;						// Branch for PF jet 100 filter bit
	TBranch *pfJetFilterBranch120;						// Branch for PF jet 120 filter bit
	
	// Leaves for the HLT tree
	Int_t MB_FirstCollisionAfterAbortGap;	 			// Branch for MB trigger
	Int_t MB_ForSkim;						 			// Branch for MB trigger
	Int_t MB_ForExpress;	 							// Branch for MB trigger
	Int_t MB_part1;	 								// Branch for MB trigger
	Int_t MB_part2;	 								// Branch for MB trigger
	Int_t MB_part3;	 								// Branch for MB trigger
	Int_t MB_part4;	 								// Branch for MB trigger
	Int_t MB_part5;	 								// Branch for MB trigger
	Int_t MB_part6;	 								// Branch for MB trigger
	Int_t MB_part7;	 								// Branch for MB trigger
	Int_t MB_part8;	 								// Branch for MB trigger
	Int_t MB_part9;	 								// Branch for MB trigger
	Int_t MB_part10;	 								// Branch for MB trigger
	Int_t MB_part11;	 								// Branch for MB trigger
	Int_t MB_part12;	 								// Branch for MB trigger
	Int_t MB_part13;	 								// Branch for MB trigger
	Int_t MB_part14;	 								// Branch for MB trigger
	Int_t MB_part15;	 								// Branch for MB trigger
	Int_t MB_part16;	 								// Branch for MB trigger
	Int_t MB_part17;	 								// Branch for MB trigger
	Int_t MB_part18;	 								// Branch for MB trigger
	Int_t MB_part19;	 								// Branch for MB trigger
	Int_t MB_part20;	 								// Branch for MB trigger
	Int_t MB_all;	 									// Branch for MB trigger
	Int_t caloJetFilterBit60;							// Filter bit for calorimeter jets 60
	Int_t caloJetFilterBit80;							// Filter bit for calorimeter jets 80
	Int_t caloJetFilterBit100;							// Filter bit for calorimeter jets 100
	Int_t pfJetFilterBit60;								// Filter bit for particle flow flow jets 60
	Int_t pfJetFilterBit80;								// Filter bit for particle flow jets 80
	Int_t pfJetFilterBit100;							// Filter bit for particle flow jets 100
	Int_t pfJetFilterBit120;							// Filter bit for particle flow jets 100

	// Branches for skim tree
	TBranch *primaryVertexBranch;						// Branch for primary vertex filter bit
	TBranch *beamScrapingBranch;				// Branch for beam scraping filter bit
	TBranch *hBHENoiseBranchLoose;				// Branch for HB/HE noise filter bit loose
	TBranch *hBHENoiseBranchTight;				// Branch for HB/HE noise filter bit tight
	TBranch *hfCoincidenceBranch;				// Branch for energy recorded one HF tower above threshold on each side
	TBranch *pVertexFilterCutdz1p0Branch;		// Branch for PU Filter default
	TBranch *pVertexFilterCutGplusBranch;		// Branch for PU Filter GPlus
	TBranch *pVertexFilterCutVtx1Branch;		// Branch for PU Filter 1 vertex only

	// Leaves for the skim tree
	Int_t primaryVertexFilterBit;				// Filter bit for primary vertex
	Int_t beamScrapingFilterBit;				// Filter bit for beam scraping
	Int_t hBHENoiseFilterLooseBit;	 			// Filter bit for HB/HE noise loose
	Int_t hBHENoiseFilterTightBit; 				// Filter bit for HB/HE noise tight
	Int_t hfCoincidenceFilterBit;				// Filter bit or energy recorded one HF tower above threshold on each side
	Int_t pVertexFilterCutdz1p0Bit;				// Filter bit for PU Filter
	Int_t pVertexFilterCutGplusBit;				// Filter bit for PU Filter
	Int_t pVertexFilterCutVtx1Bit;					// Filter bit for PU Filter
	
	// Branches for jet tree
	TBranch *nJetsBranch[nJetTrees];				// Branch for number of jets in an event
	TBranch *jetRawPtBranch[nJetTrees];				// Branch for raw jet pT
	TBranch *jetMaxTrackPtBranch[nJetTrees];		// Maximum pT for a track inside a jet
	TBranch *jetPhiBranch[nJetTrees];				// Branch for jet phi
	TBranch *jetPhiBranchWTA[nJetTrees];			// Branch for jet phi with WTA axis
	TBranch *jetEtaBranch[nJetTrees];				// Branch for jet eta
	TBranch *jetEtaBranchWTA[nJetTrees];			// Branch for jet eta with WTA axis
	TBranch *jetMassBranch[nJetTrees];				// Branch for jet mass
	TBranch *jetPfNHFBranch[nJetTrees];				// Branch for PF neutral hadron energy fraction in jets in an event
	TBranch *jetPfNEFBranch[nJetTrees];				// Branch for PF neutral EM energy fraction in jets in an event
	TBranch *jetPfCHFBranch[nJetTrees];				// Branch for PF charged hadron energy fraction in jets in an event
	TBranch *jetPfMUFBranch[nJetTrees];				// Branch for PF muon energy fraction in jets in an event
	TBranch *jetPfCEFBranch[nJetTrees];				// Branch for PF charged EM energy fraction in jets in an event
	TBranch *jetPfCHMBranch[nJetTrees];				// Branch for PF charged hadron multiplicity in jets in an event
	TBranch *jetPfCEMBranch[nJetTrees];				// Branch for PF charged EM multiplicity in jets in an event
	TBranch *jetPfNHMBranch[nJetTrees];				// Branch for PF neutral hadron multiplicity in jets in an event
	TBranch *jetPfNEMBranch[nJetTrees];				// Branch for PF neutral EM multiplicity in jets in an event
	TBranch *jetPfMUMBranch[nJetTrees];				// Branch for PF muon multiplicity in jets in an event
	TBranch *jetHCALSUMBranch[nJetTrees];			// Branch for HCAL energy sum (calo) in jets in an event
	TBranch *jetECALSUMBranch[nJetTrees];			// Branch for ECAL energy sum (calo) in jets in an event
	TBranch *jetTRKSUMBranch[nJetTrees];				// Branch for Track energy sum (calo) in jets in an event
	TBranch *jetTRKNBranch[nJetTrees];				// Branch for Track number (calo) in jets in an event
	TBranch *jetMUNBranch[nJetTrees];				// Branch for Muon number (calo) in jets in an event

	TBranch *jetRefPtBranch[nJetTrees];				// Branch for reference generator level pT for a reconstructed jet
	TBranch *jetRefEtaBranch[nJetTrees];			// Branch for reference generator level eta for a reconstructed jet
	TBranch *jetRefPhiBranch[nJetTrees];			// Branch for reference generator level phi for a reconstructed jet
	TBranch *jetRefFlavorBranch[nJetTrees];			// Branch for flavor for the parton initiating the jet
	TBranch *jetRefFlavorForBBranch[nJetTrees];		// Branch for flavor for the parton initiating the jet
	TBranch *jetRefSubidBranch[nJetTrees];		    // Branch for jet subid
	TBranch *jetRefMassBranch[nJetTrees];			// Branch for reference generator level mass for a reconstructed jet

	TBranch *nGenJetsBranch[nJetTrees];				// Branch for the number of generator level jets in an event
	TBranch *genJetPtBranch[nJetTrees];				// Branch for the generator level jet pT
	TBranch *genJetEtaBranch[nJetTrees];			// Branch for the generetor level jet eta
	TBranch *genJetEtaBranchWTA[nJetTrees];			// Branch for the generetor level jet eta with WTA axis
	TBranch *genJetPhiBranch[nJetTrees];			// Branch for the generator level jet phi
	TBranch *genJetPhiBranchWTA[nJetTrees];			// Branch for the generator level jet phi with WTA axis
	TBranch *genJetSubidBranch[nJetTrees];          // Branch for the generator level jet subid
	TBranch *genJetMatchIndexBranch[nJetTrees];     // Branch for the generator level jet matched index
	TBranch *genJetMassBranch[nJetTrees];			// Branch for the generator level jet mass
	
	// Leaves for jet tree
	Int_t nJets[nJetTrees];									// number of jets in an event
	Float_t jetRawPtArray[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the jets in an event
	Float_t jetMaxTrackPtArray[nJetTrees][nMaxJet] = {{0}}; // maximum track pT inside a jet for all the jets in an event
	Float_t jetPhiArray[nJetTrees][nMaxJet] = {{0}};		// phi of all the jets in an event
	Float_t jetPhiArrayWTA[nJetTrees][nMaxJet] = {{0}};		// phi of all the jets in an event	with WTA axis
	Float_t jetEtaArray[nJetTrees][nMaxJet] = {{0}};		// eta of all the jets in an event
	Float_t jetEtaArrayWTA[nJetTrees][nMaxJet] = {{0}};		// eta of all the jets in an event	with WTA axis
	Float_t jetMassArray[nJetTrees][nMaxJet] = {{0}};		// Mass of all the jets in an event
	Float_t jetPfNHFArray[nJetTrees][nMaxJet] = {{0}};		// PF neutral hadron energy fraction in jets in an event
	Float_t jetPfNEFArray[nJetTrees][nMaxJet] = {{0}};		// PF neutral EM energy fraction in jets in an event
	Float_t jetPfCHFArray[nJetTrees][nMaxJet] = {{0}};		// PF charged hadron energy fraction in jets in an event
	Float_t jetPfMUFArray[nJetTrees][nMaxJet] = {{0}};		// PF muon energy fraction in jets in an event
	Float_t jetPfCEFArray[nJetTrees][nMaxJet] = {{0}};		// PF charged EM energy fraction in jets in an event
	Int_t jetPfCHMArray[nJetTrees][nMaxJet] = {{0}};		// PF charged hadron multiplicity in jets in an event
	Int_t jetPfCEMArray[nJetTrees][nMaxJet] = {{0}};		// PF charged EM multiplicity in jets in an event
	Int_t jetPfNHMArray[nJetTrees][nMaxJet] = {{0}};		// PF neutral hadron multiplicity in jets in an event
	Int_t jetPfNEMArray[nJetTrees][nMaxJet] = {{0}};		// PF neutral EM multiplicity in jets in an event
	Int_t jetPfMUMArray[nJetTrees][nMaxJet] = {{0}};		// PF muon multiplicity in jets in an event
	Float_t jetHCALSUMArray[nJetTrees][nMaxJet] = {{0}};	// HCAL energy sum (calo) in jets in an event
	Float_t jetECALSUMArray[nJetTrees][nMaxJet] = {{0}};	// ECAL energy sum (calo) in jets in an event
	Float_t jetTRKSUMArray[nJetTrees][nMaxJet] = {{0}};		// Track energy sum (calo) in jets in an event
	Int_t jetTRKNArray[nJetTrees][nMaxJet] = {{0}};			// Track number (calo) in jets in an event
	Int_t jetMUNArray[nJetTrees][nMaxJet] = {{0}};			// Muon number (calo) in jets in an event

	Float_t jetRefPtArray[nJetTrees][nMaxJet] = {{0}};		// reference generator level pT for a reconstructed jet
	Float_t jetRefEtaArray[nJetTrees][nMaxJet] = {{0}};		// reference generator level pT for a reconstructed jet
	Float_t jetRefPhiArray[nJetTrees][nMaxJet] = {{0}};		// reference generator level pT for a reconstructed jet
	Int_t jetRefFlavorArray[nJetTrees][nMaxJet] = {{0}};	// flavor for initiating parton for the reference gen jet
	Int_t jetRefFlavorForBArray[nJetTrees][nMaxJet] = {{0}};// heavy flavor for initiating parton for the reference gen jet
	Int_t jetRefSubidArray[nJetTrees][nMaxJet] = {{0}};     // jet subid
	Float_t jetRefMassArray[nJetTrees][nMaxJet] = {{0}};	// reference generator level mass for a reconstructed jet

	Int_t nGenJets[nJetTrees];								// number of generator level jets in an event
	Float_t genJetPtArray[nJetTrees][nMaxJet] = {{0}};		// pT of all the generator level jets in an event
	Float_t genJetPhiArray[nJetTrees][nMaxJet] = {{0}};		// phi of all the generator level jets in an event
	Float_t genJetPhiArrayWTA[nJetTrees][nMaxJet] = {{0}};	// phi of all the generator level jets in an event with WTA axis
	Float_t genJetEtaArray[nJetTrees][nMaxJet] = {{0}};		// eta of all the generator level jets in an event
	Float_t genJetEtaArrayWTA[nJetTrees][nMaxJet] = {{0}};	// eta of all the generator level jets in an event with WTA axis
	Int_t genJetSubidArray[nJetTrees][nMaxJet] = {{0}};     // subid of all the generator level jets in an event
	Int_t genJetMatchIndexArray[nJetTrees][nMaxJet] = {{0}};// matched index of all the generator level jets in an event
	Float_t genJetMassArray[nJetTrees][nMaxJet] = {{0}};	// mass of all the generator level jets in an event

	// Branches for track tree
	TBranch *nTracksBranch;									// Branch for number of tracks
	TBranch *trackPtBranch;									// Branch for track pT
	TBranch *trackPtErrorBranch;							// Branch for track pT error
	TBranch *trackPhiBranch;								// Branch for track phi
	TBranch *trackEtaBranch;								// Branch for track eta
	TBranch *trackHighPurityBranch;							// Branch for high purity of the track
	TBranch *trackVertexDistanceZBranch;			 		// Branch for track distance from primary vertex in z-direction
	TBranch *trackVertexDistanceZErrorBranch;				// Branch for error for track distance from primary vertex in z-direction
	TBranch *trackVertexDistanceXYBranch;					// Branch for track distance from primary vertex in xy-direction
	TBranch *trackVertexDistanceXYErrorBranch; 				// Branch for error for track distance from primary vertex in xy-direction
	TBranch *trackEnergyEcalBranch;							// Branch for track energy in ECal
	TBranch *trackEnergyHcalBranch;							// Branch for track energy in HCal
	TBranch *trackChargeBranch;								// Branch for track charge
	TBranch *PixelnHitsTrackBranch;							// Branch for number of valid pixel hits for the track
	
	// Leaves for the track tree
	Int_t nTracks;														// Number of tracks
	Float_t trackPtArray[nMaxTrack] = {0};								// Array for track pT
	Float_t trackPtErrorArray[nMaxTrack] = {0};							// Array for track pT errors
	Float_t trackPhiArray[nMaxTrack] = {0};								// Array for track phis
	Float_t trackEtaArray[nMaxTrack] = {0};								// Array for track etas
	Bool_t trackHighPurityArray[nMaxTrack] = {0};						// Array for the high purity of tracks
	Float_t trackVertexDistanceZArray[nMaxTrack] = {0};			 		// Array for track distance from primary vertex in z-direction
	Float_t trackVertexDistanceZErrorArray[nMaxTrack] = {0};			// Array for error for track distance from primary vertex in z-direction
	Float_t trackVertexDistanceXYArray[nMaxTrack] = {0};				// Array for track distance from primary vertex in xy-direction
	Float_t trackVertexDistanceXYErrorArray[nMaxTrack] = {0}; 			// Array for error for track distance from primary vertex in xy-direction
	Float_t trackEnergyEcalArray[nMaxTrack] = {0};						// Array for track energy in ECal
	Float_t trackEnergyHcalArray[nMaxTrack] = {0};						// Array for track energy in HCal
	Int_t trackChargeArray[nMaxTrack] = {0}; 										// Array for track charge
	UChar_t PixelnHitsTrackArray[nMaxTrack] = {0}; 							// Array for number of valid pixel hits for the track

	// Branches for generator level track tree
	TBranch *genTrackPtBranch;				 // Branch for generator level track pT:s
	TBranch *genTrackPhiBranch;				 // Branch for generator level track phis
	TBranch *genTrackEtaBranch;				 // Branch for generator level track etas
	TBranch *genTrackPdgBranch;				 // Branch for generator level track PDG code
	TBranch *genTrackChargeBranch;				 // Branch for generator level track charges
	TBranch *genTrackSubeventBranch;			 // Branch for generator level track subevent indices (0 = PYTHIA, (>0) = other MC)
	TBranch *genTrackMassBranch;				 // Branch for generator level track masses
	
	// Leaves for generator level track tree
	vector<float> *genTrackPtArray;			 // Vector for generator level track pT
	vector<float> *genTrackPhiArray;		 // Vector for generator level track phi
	vector<float> *genTrackEtaArray;		 // Vector for generator level track eta
	vector<int>	 *genTrackPdgArray;			 // Vector for generator level track PDG code
	vector<int>	 *genTrackChargeArray;		 // Vector for generator level track charges
	vector<int>	 *genTrackSubeventArray;	 // Vector for generator level track subevent indices (0 = PYTHIA, (>0) = other MC)
	vector<float> *genTrackMassArray;		 // Vector for generator level track mass
	
	// Branches for particle flow candidate ID tree
	TBranch *particleFlowCandidatePtBranch;		// Branch for particle flow candidate pT
	TBranch *particleFlowCandidateEtaBranch;	 // Branch for particle flow candidate eta
	TBranch *particleFlowCandidatePhiBranch;	 // Branch for particle flow candidate phi
	TBranch *particleFlowCandidateEnergyBranch;	 // Branch for particle flow candidate energy
	TBranch *particleFlowCandidateIDBranch;	 // Branch for particle flow candidate ID --> See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideParticleFlow#Particle_Identification_and_Reco
	TBranch *particleFlowCandidateMassBranch;	 // Branch for particle flow candidate mass
	
	// Leaves for particle flow candidate tree
	vector<float> *particleFlowCandidatePtVector;		// Vector for particle flow candidate pT
	vector<float> *particleFlowCandidateEtaVector;		// Vector for particle flow candidate eta
	vector<float> *particleFlowCandidatePhiVector;		// Vector for particle flow candidate phi
	vector<float> *particleFlowCandidateEnergyVector;	// Vector for particle flow candidate energy
	vector<int> *particleFlowCandidateIDVector;			// Vector for particle flow candidate ID --> See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideParticleFlow#Particle_Identification_and_Reco
	vector<float> *particleFlowCandidateMassVector;	    // Vector for particle flow candidate mass

	// Branches for rho
	TBranch *RhoetamaxBranch;				 // Branch for rho etamax
	TBranch *RhoetaminBranch;				 // Branch for rho etamin
	TBranch *RhoBranch;				 		 // Branch for rho
	TBranch *RhomBranch;					 // Branch for rhom
	// Leaves for rho	
	vector<double> *RhoetamaxVector;		// Vector for particle flow candidate pT
	vector<double> *RhoetaminVector;		// Vector for particle flow candidate pT
	vector<double> *RhoVector;		// Vector for particle flow candidate pT
	vector<double> *RhomVector;		// Vector for particle flow candidate pT


	// ========================================== //
	// Read all the branches from the input trees //
	// ========================================== //
	
	// Connect the branches of the heavy ion tree
	heavyIonTree->SetBranchStatus("*",0); // remove all branchs to read it fast
	heavyIonTree->SetBranchStatus("run",1);
	heavyIonTree->SetBranchAddress("run",&run,&runBranch);
	heavyIonTree->SetBranchStatus("evt",1);
	heavyIonTree->SetBranchAddress("evt",&event,&eventBranch);
	heavyIonTree->SetBranchStatus("lumi",1);
	heavyIonTree->SetBranchAddress("lumi",&lumi,&lumiBranch);
	heavyIonTree->SetBranchStatus("vz",1);
	heavyIonTree->SetBranchAddress("vz",&vertexZ,&hiVzBranch);
	heavyIonTree->SetBranchStatus("vx",1);
	heavyIonTree->SetBranchAddress("vx",&vertexX,&hiVxBranch);
	heavyIonTree->SetBranchStatus("vy",1);
	heavyIonTree->SetBranchAddress("vy",&vertexY,&hiVyBranch);
	heavyIonTree->SetBranchStatus("hiHFplus",1);
	heavyIonTree->SetBranchAddress("hiHFplus",&hiHFplus,&hiHFplusBranch);
	heavyIonTree->SetBranchStatus("hiHFminus",1);
	heavyIonTree->SetBranchAddress("hiHFminus",&hiHFminus,&hiHFminusBranch);
	heavyIonTree->SetBranchStatus("hiHFplusEta4",1);
	heavyIonTree->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4,&hiHFplusEta4Branch);
	heavyIonTree->SetBranchStatus("hiHFminusEta4",1);
	heavyIonTree->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4,&hiHFminusEta4Branch);
	heavyIonTree->SetBranchStatus("hiZDCplus",1);
	heavyIonTree->SetBranchAddress("hiZDCplus",&hiZDCplus,&hiZDCplusBranch);
	heavyIonTree->SetBranchStatus("hiZDCminus",1);
	heavyIonTree->SetBranchAddress("hiZDCminus",&hiZDCminus,&hiZDCminusBranch);
	
	// Event plane
	checkFlatteningTree->SetBranchStatus("epang",1);
	checkFlatteningTree->SetBranchAddress("epang",&eventPlaneAngle,&eventPlaneAngleBranch);
	checkFlatteningTree->SetBranchStatus("q",1);
	checkFlatteningTree->SetBranchAddress("q",&eventPlaneQ,&eventPLaneQBranch);
	checkFlatteningTree->SetBranchStatus("qx",1);
	checkFlatteningTree->SetBranchAddress("qx",&eventPlaneQx,&eventPlaneQxBranch);
	checkFlatteningTree->SetBranchStatus("qy",1);
	checkFlatteningTree->SetBranchAddress("qy",&eventPlaneQy,&eventPlaneQyBranch);
	checkFlatteningTree->SetBranchStatus("mult",1);
	checkFlatteningTree->SetBranchAddress("mult",&eventPlaneMultiplicity,&eventPlaneMultiplicityBranch);

	// ptHat and event weight only for MC
	if(is_MC){
		heavyIonTree->SetBranchStatus("pthat",1);
		heavyIonTree->SetBranchAddress("pthat",&ptHat,&ptHatBranch);
		heavyIonTree->SetBranchStatus("weight",1);
		heavyIonTree->SetBranchAddress("weight",&eventWeight,&eventWeightBranch);
	}
	
	// Connect the branches to the HLT tree
	hltTree->SetBranchStatus("*",0);

	hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_FirstCollisionAfterAbortGap_v1",1);
	hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_FirstCollisionAfterAbortGap_v1",&MB_FirstCollisionAfterAbortGap,&MB_FirstCollisionAfterAbortGapBranch);
	hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_ForSkim_v1",1);
	hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_ForSkim_v1",&MB_ForSkim,&MB_ForSkimBranch);
	hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_ForExpress_v1",1);
	hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_ForExpress_v1",&MB_ForExpress,&MB_ForExpressBranch);
	if(largersample == 1){
		hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1",1);
		hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v1",&MB_part1,&MB_part1Branch);
		hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1",1);
		hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v1",&MB_part2,&MB_part2Branch);
		hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1",1);
		hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v1",&MB_part3,&MB_part3Branch);
		hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1",1);
		hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v1",&MB_part4,&MB_part4Branch);
		hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1",1);
		hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v1",&MB_part5,&MB_part5Branch);
		hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1",1);
		hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v1",&MB_part6,&MB_part6Branch);
		hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1",1);
		hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v1",&MB_part7,&MB_part7Branch);
		hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1",1);
		hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v1",&MB_part8,&MB_part8Branch);
	}else if(largersample == 0){
		hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v2",1);
		hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part1_v2",&MB_part1,&MB_part1Branch);
		hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v2",1);
		hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part2_v2",&MB_part2,&MB_part2Branch);
		hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v2",1);
		hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part3_v2",&MB_part3,&MB_part3Branch);
		hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v2",1);
		hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part4_v2",&MB_part4,&MB_part4Branch);
		hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v2",1);
		hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part5_v2",&MB_part5,&MB_part5Branch);
		hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v2",1);
		hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part6_v2",&MB_part6,&MB_part6Branch);
		hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v2",1);
		hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part7_v2",&MB_part7,&MB_part7Branch);
		hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v2",1);
		hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part8_v2",&MB_part8,&MB_part8Branch);
		if(!is_MC){
			hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part9_v1",1);
			hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part9_v1",&MB_part9,&MB_part9Branch);
			hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part10_v1",1);
			hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part10_v1",&MB_part10,&MB_part10Branch);
			hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part11_v1",1);
			hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part11_v1",&MB_part11,&MB_part11Branch);
			hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part12_v1",1);
			hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part12_v1",&MB_part12,&MB_part12Branch);
			hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part13_v1",1);
			hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part13_v1",&MB_part13,&MB_part13Branch);
			hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part14_v1",1);
			hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part14_v1",&MB_part14,&MB_part14Branch);
			hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part15_v1",1);
			hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part15_v1",&MB_part15,&MB_part15Branch);
			hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part16_v1",1);
			hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part16_v1",&MB_part16,&MB_part16Branch);	
			hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part17_v1",1);
			hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part17_v1",&MB_part17,&MB_part17Branch);
			hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part18_v1",1);
			hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part18_v1",&MB_part18,&MB_part18Branch);
			hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part19_v1",1);
			hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part19_v1",&MB_part19,&MB_part19Branch);
			hltTree->SetBranchStatus("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part20_v1",1);
			hltTree->SetBranchAddress("HLT_PAL1MinimumBiasHF_OR_SinglePixelTrack_part20_v1",&MB_part20,&MB_part20Branch);
		}
	}

	
	hltTree->SetBranchStatus("HLT_PAAK4CaloJet60_Eta5p1_v3",1);
	hltTree->SetBranchAddress("HLT_PAAK4CaloJet60_Eta5p1_v3",&caloJetFilterBit60,&caloJetFilterBranch60);
	hltTree->SetBranchStatus("HLT_PAAK4CaloJet80_Eta5p1_v3",1);
	hltTree->SetBranchAddress("HLT_PAAK4CaloJet80_Eta5p1_v3",&caloJetFilterBit80,&caloJetFilterBranch80);
	hltTree->SetBranchStatus("HLT_PAAK4CaloJet100_Eta5p1_v3",1);
	hltTree->SetBranchAddress("HLT_PAAK4CaloJet100_Eta5p1_v3",&caloJetFilterBit100,&caloJetFilterBranch100);
	hltTree->SetBranchStatus("HLT_PAAK4PFJet60_Eta5p1_v4",1);

	hltTree->SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v4",&pfJetFilterBit60,&pfJetFilterBranch60);
	hltTree->SetBranchStatus("HLT_PAAK4PFJet80_Eta5p1_v3",1);
	hltTree->SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v3",&pfJetFilterBit80,&pfJetFilterBranch80);
	hltTree->SetBranchStatus("HLT_PAAK4PFJet100_Eta5p1_v3",1);
	hltTree->SetBranchAddress("HLT_PAAK4PFJet100_Eta5p1_v3",&pfJetFilterBit100,&pfJetFilterBranch100);
	hltTree->SetBranchStatus("HLT_PAAK4PFJet120_Eta5p1_v2",1);
	hltTree->SetBranchAddress("HLT_PAAK4PFJet120_Eta5p1_v2",&pfJetFilterBit120,&pfJetFilterBranch120);


	// Connect the branches to the skim tree
	skimTree->SetBranchStatus("*",0);
	skimTree->SetBranchStatus("pPAprimaryVertexFilter",1);
	skimTree->SetBranchAddress("pPAprimaryVertexFilter",&primaryVertexFilterBit,&primaryVertexBranch);
	skimTree->SetBranchStatus("pBeamScrapingFilter",1);
	skimTree->SetBranchAddress("pBeamScrapingFilter",&beamScrapingFilterBit,&beamScrapingBranch);
	skimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
	skimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&hBHENoiseFilterLooseBit,&hBHENoiseBranchLoose);
	skimTree->SetBranchStatus("HBHENoiseFilterResultRun2Tight",1);
	skimTree->SetBranchAddress("HBHENoiseFilterResultRun2Tight",&hBHENoiseFilterTightBit,&hBHENoiseBranchTight);
	skimTree->SetBranchStatus("phfCoincFilter",1);
	skimTree->SetBranchAddress("phfCoincFilter", &hfCoincidenceFilterBit, &hfCoincidenceBranch);
	skimTree->SetBranchStatus("pVertexFilterCutdz1p0",1);
	skimTree->SetBranchAddress("pVertexFilterCutdz1p0", &pVertexFilterCutdz1p0Bit, &pVertexFilterCutdz1p0Branch);
	skimTree->SetBranchStatus("pVertexFilterCutGplus",1);
	skimTree->SetBranchAddress("pVertexFilterCutGplus", &pVertexFilterCutGplusBit, &pVertexFilterCutGplusBranch);
	skimTree->SetBranchStatus("pVertexFilterCutVtx1",1);
	skimTree->SetBranchAddress("pVertexFilterCutVtx1", &pVertexFilterCutVtx1Bit, &pVertexFilterCutVtx1Branch);

	RhoTree->SetBranchStatus("*",0);
	RhoTree->SetBranchStatus("etaMin",1);
	RhoTree->SetBranchAddress("etaMin",&RhoetaminVector,&RhoetaminBranch);
	RhoTree->SetBranchStatus("etaMax",1);
	RhoTree->SetBranchAddress("etaMax",&RhoetamaxVector,&RhoetamaxBranch);
	RhoTree->SetBranchStatus("rho",1);
	RhoTree->SetBranchAddress("rho",&RhoVector,&RhoBranch);
	RhoTree->SetBranchStatus("rhom",1);
	RhoTree->SetBranchAddress("rhom",&RhomVector,&RhomBranch);

	// Same branch names for all jet collections
	for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
		
		// Connect the branches to the jet tree
		jetTree[iJetType]->SetBranchStatus("*",0);

		jetTree[iJetType]->SetBranchStatus("nref",1);
		jetTree[iJetType]->SetBranchAddress("nref",&nJets[iJetType],&nJetsBranch[iJetType]);
		jetTree[iJetType]->SetBranchStatus("rawpt",1);
		jetTree[iJetType]->SetBranchAddress("rawpt",&jetRawPtArray[iJetType],&jetRawPtBranch[iJetType]);
		jetTree[iJetType]->SetBranchStatus("trackMax",1);
		jetTree[iJetType]->SetBranchAddress("trackMax",&jetMaxTrackPtArray[iJetType],&jetMaxTrackPtBranch[iJetType]);

		// Jet phi with E-scheme, WTA axes are calculated later
		jetTree[iJetType]->SetBranchStatus("jtphi",1);
		jetTree[iJetType]->SetBranchAddress("jtphi",&jetPhiArray[iJetType],&jetPhiBranch[iJetType]);
		
		// Jet eta with E-scheme, WTA axes are calculated later
		jetTree[iJetType]->SetBranchStatus("jteta",1);
		jetTree[iJetType]->SetBranchAddress("jteta",&jetEtaArray[iJetType],&jetEtaBranch[iJetType]);
		
		// Jet mass
		jetTree[iJetType]->SetBranchStatus("jtm",1);
		jetTree[iJetType]->SetBranchAddress("jtm",&jetMassArray[iJetType],&jetMassBranch[iJetType]);	

		// Jet ID stuff
		jetTree[iJetType]->SetBranchStatus("jtPfNHF",1);
		jetTree[iJetType]->SetBranchAddress("jtPfNHF",&jetPfNHFArray[iJetType],&jetPfNHFBranch[iJetType]);	
		jetTree[iJetType]->SetBranchStatus("jtPfNEF",1);
		jetTree[iJetType]->SetBranchAddress("jtPfNEF",&jetPfNEFArray[iJetType],&jetPfNEFBranch[iJetType]);	
		jetTree[iJetType]->SetBranchStatus("jtPfCHF",1);
		jetTree[iJetType]->SetBranchAddress("jtPfCHF",&jetPfCHFArray[iJetType],&jetPfCHFBranch[iJetType]);	
		jetTree[iJetType]->SetBranchStatus("jtPfMUF",1);
		jetTree[iJetType]->SetBranchAddress("jtPfMUF",&jetPfMUFArray[iJetType],&jetPfMUFBranch[iJetType]);	
		jetTree[iJetType]->SetBranchStatus("jtPfCEF",1);
		jetTree[iJetType]->SetBranchAddress("jtPfCEF",&jetPfCEFArray[iJetType],&jetPfCEFBranch[iJetType]);	
		jetTree[iJetType]->SetBranchStatus("jtPfCHM",1);
		jetTree[iJetType]->SetBranchAddress("jtPfCHM",&jetPfCHMArray[iJetType],&jetPfCHMBranch[iJetType]);	
		jetTree[iJetType]->SetBranchStatus("jtPfCEM",1);
		jetTree[iJetType]->SetBranchAddress("jtPfCEM",&jetPfCEMArray[iJetType],&jetPfCEMBranch[iJetType]);	
		jetTree[iJetType]->SetBranchStatus("jtPfNHM",1);
		jetTree[iJetType]->SetBranchAddress("jtPfNHM",&jetPfNHMArray[iJetType],&jetPfNHMBranch[iJetType]);	
		jetTree[iJetType]->SetBranchStatus("jtPfNEM",1);
		jetTree[iJetType]->SetBranchAddress("jtPfNEM",&jetPfNEMArray[iJetType],&jetPfNEMBranch[iJetType]);	
		jetTree[iJetType]->SetBranchStatus("jtPfMUM",1);
		jetTree[iJetType]->SetBranchAddress("jtPfMUM",&jetPfMUMArray[iJetType],&jetPfMUMBranch[iJetType]);	
		jetTree[iJetType]->SetBranchStatus("hcalSum",1);
		jetTree[iJetType]->SetBranchAddress("hcalSum",&jetHCALSUMArray[iJetType],&jetHCALSUMBranch[iJetType]);
		jetTree[iJetType]->SetBranchStatus("ecalSum",1);
		jetTree[iJetType]->SetBranchAddress("ecalSum",&jetECALSUMArray[iJetType],&jetECALSUMBranch[iJetType]);
		jetTree[iJetType]->SetBranchStatus("trackSum",1);
		jetTree[iJetType]->SetBranchAddress("trackSum",&jetTRKSUMArray[iJetType],&jetTRKSUMBranch[iJetType]);
		jetTree[iJetType]->SetBranchStatus("trackN",1);
		jetTree[iJetType]->SetBranchAddress("trackN",&jetTRKNArray[iJetType],&jetTRKNBranch[iJetType]);
		jetTree[iJetType]->SetBranchStatus("muN",1);
		jetTree[iJetType]->SetBranchAddress("muN",&jetMUNArray[iJetType],&jetMUNBranch[iJetType]);

		
		// If we are looking at Monte Carlo, connect the reference pT and parton arrays
		if(is_MC){
			// Matched jet variables
			jetTree[iJetType]->SetBranchStatus("refpt",1);
			jetTree[iJetType]->SetBranchAddress("refpt",&jetRefPtArray[iJetType],&jetRefPtBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refeta",1);
			jetTree[iJetType]->SetBranchAddress("refeta",&jetRefEtaArray[iJetType],&jetRefEtaBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refphi",1);
			jetTree[iJetType]->SetBranchAddress("refphi",&jetRefPhiArray[iJetType],&jetRefPhiBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refparton_flavor",1);
			jetTree[iJetType]->SetBranchAddress("refparton_flavor",&jetRefFlavorArray[iJetType],&jetRefFlavorBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refparton_flavorForB",1);
			jetTree[iJetType]->SetBranchAddress("refparton_flavorForB", &jetRefFlavorForBArray[iJetType], &jetRefFlavorForBBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("subid",1);
			jetTree[iJetType]->SetBranchAddress("subid", &jetRefSubidArray[iJetType], &jetRefSubidBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("refm",1);
			jetTree[iJetType]->SetBranchAddress("refm",&jetRefMassArray[iJetType],&jetRefMassBranch[iJetType]);	
			
			// Gen jet variables
			jetTree[iJetType]->SetBranchStatus("ngen",1);
			jetTree[iJetType]->SetBranchAddress("ngen",&nGenJets[iJetType],&nGenJetsBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("genpt",1);
			jetTree[iJetType]->SetBranchAddress("genpt",&genJetPtArray[iJetType],&genJetPtBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("genphi",1);
			jetTree[iJetType]->SetBranchAddress("genphi",&genJetPhiArray[iJetType],&genJetPhiBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("geneta",1);
			jetTree[iJetType]->SetBranchAddress("geneta",&genJetEtaArray[iJetType],&genJetEtaBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("genmatchindex",1);
			jetTree[iJetType]->SetBranchAddress("genmatchindex",&genJetMatchIndexArray[iJetType],&genJetMatchIndexBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("gensubid",1);
			jetTree[iJetType]->SetBranchAddress("gensubid",&genJetSubidArray[iJetType],&genJetSubidBranch[iJetType]);
			jetTree[iJetType]->SetBranchStatus("genm",1);
			jetTree[iJetType]->SetBranchAddress("genm",&genJetMassArray[iJetType],&genJetMassBranch[iJetType]);	
			
		}
		
	} // Loop over different jet collections


	// Connect the branches to the track tree
	trackTree->SetBranchStatus("*",0);

	trackTree->SetBranchStatus("nTrk",1);
	trackTree->SetBranchAddress("nTrk",&nTracks,&nTracksBranch);
	trackTree->SetBranchStatus("highPurity",1);
	trackTree->SetBranchAddress("highPurity",&trackHighPurityArray,&trackHighPurityBranch);
	trackTree->SetBranchStatus("trkPt",1);
	trackTree->SetBranchAddress("trkPt",&trackPtArray,&trackPtBranch);
	trackTree->SetBranchStatus("trkPtError",1);
	trackTree->SetBranchAddress("trkPtError",&trackPtErrorArray,&trackPtErrorBranch);
	trackTree->SetBranchStatus("trkPhi",1);
	trackTree->SetBranchAddress("trkPhi",&trackPhiArray,&trackPhiBranch);
	trackTree->SetBranchStatus("trkEta",1);
	trackTree->SetBranchAddress("trkEta",&trackEtaArray,&trackEtaBranch);
	trackTree->SetBranchStatus("trkDz1",1);
	trackTree->SetBranchAddress("trkDz1",&trackVertexDistanceZArray,&trackVertexDistanceZBranch);
	trackTree->SetBranchStatus("trkDzError1",1);
	trackTree->SetBranchAddress("trkDzError1",&trackVertexDistanceZErrorArray,&trackVertexDistanceZErrorBranch);
	trackTree->SetBranchStatus("trkDxy1",1);
	trackTree->SetBranchAddress("trkDxy1",&trackVertexDistanceXYArray,&trackVertexDistanceXYBranch);
	trackTree->SetBranchStatus("trkDxyError1",1);
	trackTree->SetBranchAddress("trkDxyError1",&trackVertexDistanceXYErrorArray,&trackVertexDistanceXYErrorBranch);
	trackTree->SetBranchStatus("trkNPixelHit",1);
	trackTree->SetBranchAddress("trkNPixelHit",&PixelnHitsTrackArray,&PixelnHitsTrackBranch);
	trackTree->SetBranchStatus("pfEcal",1);
	trackTree->SetBranchAddress("pfEcal",&trackEnergyEcalArray,&trackEnergyEcalBranch);
	trackTree->SetBranchStatus("pfHcal",1);
	trackTree->SetBranchAddress("pfHcal",&trackEnergyHcalArray,&trackEnergyHcalBranch);
	trackTree->SetBranchStatus("trkCharge",1);
	trackTree->SetBranchAddress("trkCharge",&trackChargeArray,&trackChargeBranch);
	
	// Generator level tracks only in Monte Carlo
	if(is_MC){
		
		// Connect the branches to generator level track tree
		genTrackTree->SetBranchStatus("*",0);
		genTrackTree->SetBranchStatus("pt",1);
		genTrackTree->SetBranchAddress("pt",&genTrackPtArray,&genTrackPtBranch);
		genTrackTree->SetBranchStatus("phi",1);
		genTrackTree->SetBranchAddress("phi",&genTrackPhiArray,&genTrackPhiBranch);
		genTrackTree->SetBranchStatus("eta",1);
		genTrackTree->SetBranchAddress("eta",&genTrackEtaArray,&genTrackEtaBranch);
		genTrackTree->SetBranchStatus("pdg",1);
		genTrackTree->SetBranchAddress("pdg",&genTrackPdgArray,&genTrackPdgBranch);
		genTrackTree->SetBranchStatus("chg",1);
		genTrackTree->SetBranchAddress("chg",&genTrackChargeArray,&genTrackChargeBranch);
		genTrackTree->SetBranchStatus("sube",1);
		genTrackTree->SetBranchAddress("sube",&genTrackSubeventArray,&genTrackSubeventBranch);
		genTrackTree->SetBranchStatus("mass",1);
		genTrackTree->SetBranchAddress("mass",&genTrackMassArray,&genTrackMassBranch);
		
	}
	
	// Connect the branches to the particle flow candidate tree if requested
	particleFlowCandidateTree->SetBranchStatus("*",0);
	particleFlowCandidateTree->SetBranchStatus("pfPt",1);
	particleFlowCandidateTree->SetBranchAddress("pfPt",&particleFlowCandidatePtVector,&particleFlowCandidatePtBranch);
	particleFlowCandidateTree->SetBranchStatus("pfPhi",1);
	particleFlowCandidateTree->SetBranchAddress("pfPhi",&particleFlowCandidatePhiVector,&particleFlowCandidatePhiBranch);
	particleFlowCandidateTree->SetBranchStatus("pfEta",1);
	particleFlowCandidateTree->SetBranchAddress("pfEta",&particleFlowCandidateEtaVector,&particleFlowCandidateEtaBranch);
	particleFlowCandidateTree->SetBranchStatus("pfEnergy",1);
	particleFlowCandidateTree->SetBranchAddress("pfEnergy",&particleFlowCandidateEnergyVector,&particleFlowCandidateEnergyBranch);
	particleFlowCandidateTree->SetBranchStatus("pfId",1);
	particleFlowCandidateTree->SetBranchAddress("pfId",&particleFlowCandidateIDVector,&particleFlowCandidateIDBranch);
	particleFlowCandidateTree->SetBranchStatus("pfM",1);
	particleFlowCandidateTree->SetBranchAddress("pfM",&particleFlowCandidateMassVector,&particleFlowCandidateMassBranch);

	// ========================================== //
	//			 Define output trees
	// ========================================== //
	
	// Copy the heavy ion tree to the output
	TTree *heavyIonTreeOutput = new TTree("HiTree","");
	// Connect the branches of the heavy ion tree
	heavyIonTreeOutput->Branch("run",&run,"run/i");
	heavyIonTreeOutput->Branch("evt",&event,"evt/l");
	heavyIonTreeOutput->Branch("lumi",&lumi,"lumi/i");
	heavyIonTreeOutput->Branch("vz",&vertexZ,"vz/F");
	heavyIonTreeOutput->Branch("vx",&vertexX,"vx/F");
	heavyIonTreeOutput->Branch("vy",&vertexY,"vy/F");
	heavyIonTreeOutput->Branch("hiHFplus",&hiHFplus,"hiHFplus/F");
	heavyIonTreeOutput->Branch("hiHFminus",&hiHFminus,"hiHFminus/F");
	heavyIonTreeOutput->Branch("hiHFplusEta4",&hiHFplusEta4,"hiHFplusEta4/F");
	heavyIonTreeOutput->Branch("hiHFminusEta4",&hiHFminusEta4,"hiHFminusEta4/F");	
	if(!storetracks) heavyIonTreeOutput->Branch("Ntroff",&Ntroff,"Ntroff/I");	
	heavyIonTreeOutput->Branch("hiZDCplus",&hiZDCplus,"hiZDCplus/F");
	heavyIonTreeOutput->Branch("hiZDCminus",&hiZDCminus,"hiZDCminus/F");
	heavyIonTreeOutput->Branch("hi_FRG",&hi_FRG,"hi_FRG/F");	
	heavyIonTreeOutput->Branch("hi_FRG_noNsel",&hi_FRG_noNsel,"hi_FRG_noNsel/F");
	heavyIonTreeOutput->Branch("hi_BRG",&hi_BRG,"hi_BRG/F");
	heavyIonTreeOutput->Branch("hi_BRG_noNsel",&hi_BRG_noNsel,"hi_BRG_noNsel/F");
	
	
	// Event plane
	TTree *checkFlatteningTreeOutput = new TTree("tree","");
	Float_t epang_HFm2, epang_HFp2,epang_HFm3, epang_HFp3, epang_HFm4, epang_HFp4;	
	Float_t q_HFm2, q_HFp2,q_HFm3, q_HFp3, q_HFm4, q_HFp4;
	Float_t qx_HFm2, qx_HFp2,qx_HFm3, qx_HFp3, qx_HFm4, qx_HFp4;
	Float_t qy_HFm2, qy_HFp2,qy_HFm3, qy_HFp3, qy_HFm4, qy_HFp4;
	Float_t mult_HFm2, mult_HFp2,mult_HFm3, mult_HFp3, mult_HFm4, mult_HFp4;

	checkFlatteningTreeOutput->Branch("epang_HFm2",&epang_HFm2,"epang_HFm2/F");
	checkFlatteningTreeOutput->Branch("epang_HFp2",&epang_HFp2,"epang_HFp2/F");
	checkFlatteningTreeOutput->Branch("epang_HFm3",&epang_HFm3,"epang_HFm3/F");
	checkFlatteningTreeOutput->Branch("epang_HFp3",&epang_HFp3,"epang_HFp3/F");
	checkFlatteningTreeOutput->Branch("epang_HFm4",&epang_HFm4,"epang_HFm4/F");
	checkFlatteningTreeOutput->Branch("epang_HFp4",&epang_HFp4,"epang_HFp4/F");

	checkFlatteningTreeOutput->Branch("q_HFm2",&q_HFm2,"q_HFm2/F");
	checkFlatteningTreeOutput->Branch("q_HFp2",&q_HFp2,"q_HFp2/F");
	checkFlatteningTreeOutput->Branch("q_HFm3",&q_HFm3,"q_HFm3/F");
	checkFlatteningTreeOutput->Branch("q_HFp3",&q_HFp3,"q_HFp3/F");
	checkFlatteningTreeOutput->Branch("q_HFm4",&q_HFm4,"q_HFm4/F");
	checkFlatteningTreeOutput->Branch("q_HFp4",&q_HFp4,"q_HFp4/F");

	checkFlatteningTreeOutput->Branch("qx_HFm2",&qx_HFm2,"qx_HFm2/F");
	checkFlatteningTreeOutput->Branch("qx_HFp2",&qx_HFp2,"qx_HFp2/F");
	checkFlatteningTreeOutput->Branch("qx_HFm3",&qx_HFm3,"qx_HFm3/F");
	checkFlatteningTreeOutput->Branch("qx_HFp3",&qx_HFp3,"qx_HFp3/F");
	checkFlatteningTreeOutput->Branch("qx_HFm4",&qx_HFm4,"qx_HFm4/F");
	checkFlatteningTreeOutput->Branch("qx_HFp4",&qx_HFp4,"qx_HFp4/F");

	checkFlatteningTreeOutput->Branch("qy_HFm2",&qy_HFm2,"qy_HFm2/F");
	checkFlatteningTreeOutput->Branch("qy_HFp2",&qy_HFp2,"qy_HFp2/F");
	checkFlatteningTreeOutput->Branch("qy_HFm3",&qy_HFm3,"qy_HFm3/F");
	checkFlatteningTreeOutput->Branch("qy_HFp3",&qy_HFp3,"qy_HFp3/F");
	checkFlatteningTreeOutput->Branch("qy_HFm4",&qy_HFm4,"qy_HFm4/F");
	checkFlatteningTreeOutput->Branch("qy_HFp4",&qy_HFp4,"qy_HFp4/F");

	checkFlatteningTreeOutput->Branch("mult_HFm2",&mult_HFm2,"mult_HFm2/F");
	checkFlatteningTreeOutput->Branch("mult_HFp2",&mult_HFp2,"mult_HFp2/F");
	checkFlatteningTreeOutput->Branch("mult_HFm3",&mult_HFm3,"mult_HFm3/F");
	checkFlatteningTreeOutput->Branch("mult_HFp3",&mult_HFp3,"mult_HFp3/F");
	checkFlatteningTreeOutput->Branch("mult_HFm4",&mult_HFm4,"mult_HFm4/F");
	checkFlatteningTreeOutput->Branch("mult_HFp4",&mult_HFp4,"mult_HFp4/F");
	
	// ptHat and event weight only for MC
	if(is_MC){
		heavyIonTreeOutput->Branch("pthat",&ptHat,"pthat/F");
		heavyIonTreeOutput->Branch("weight",&eventWeight,"weight/F");
	}
	
	// Copy the HLT tree to the output
	TTree *hltTreeOutput = new TTree("HltTree","");	
	// Connect the branches of the HLT tree
	hltTreeOutput->Branch("HLT_MinimumBiasAll",&MB_all,"HLT_MinimumBiasAll/I");
	hltTreeOutput->Branch("HLT_PAAK4CaloJet60_Eta5p1_v3",&caloJetFilterBit60,"HLT_PAAK4CaloJet60_Eta5p1_v3/I");
	hltTreeOutput->Branch("HLT_PAAK4CaloJet80_Eta5p1_v3",&caloJetFilterBit80,"HLT_PAAK4CaloJet80_Eta5p1_v3/I");
	hltTreeOutput->Branch("HLT_PAAK4CaloJet100_Eta5p1_v3",&caloJetFilterBit100,"HLT_PAAK4CaloJet100_Eta5p1_v3/I");
	hltTreeOutput->Branch("HLT_PAAK4PFJet60_Eta5p1_v4",&pfJetFilterBit60,"HLT_PAAK4PFJet60_Eta5p1_v4/I");
	hltTreeOutput->Branch("HLT_PAAK4PFJet80_Eta5p1_v3",&pfJetFilterBit80,"HLT_PAAK4PFJet80_Eta5p1_v3/I");
	hltTreeOutput->Branch("HLT_PAAK4PFJet100_Eta5p1_v3",&pfJetFilterBit100,"HLT_PAAK4PFJet100_Eta5p1_v3/I");
	hltTreeOutput->Branch("HLT_PAAK4PFJet120_Eta5p1_v2",&pfJetFilterBit120,"HLT_PAAK4PFJet120_Eta5p1_v2/I");

	// Copy the skim tree to the output
	TTree *skimTreeOutput = new TTree("HltTree","");
	skimTreeOutput->Branch("pPAprimaryVertexFilter",&primaryVertexFilterBit,"pPAprimaryVertexFilter/I");
	skimTreeOutput->Branch("pBeamScrapingFilter",&beamScrapingFilterBit,"pBeamScrapingFilter/I");
	skimTreeOutput->Branch("HBHENoiseFilterResultRun2Loose",&hBHENoiseFilterLooseBit,"HBHENoiseFilterResultRun2Loose/I");
	skimTreeOutput->Branch("HBHENoiseFilterResultRun2Tight",&hBHENoiseFilterTightBit,"HBHENoiseFilterResultRun2Tight/I");
	skimTreeOutput->Branch("phfCoincFilter", &hfCoincidenceFilterBit, "phfCoincFilter/I");
	skimTreeOutput->Branch("pVertexFilterCutdz1p0", &pVertexFilterCutdz1p0Bit, "pVertexFilterCutdz1p0/I");
	skimTreeOutput->Branch("pVertexFilterCutGplus",&pVertexFilterCutGplusBit,"pVertexFilterCutGplus/I");
	skimTreeOutput->Branch("pVertexFilterCutVtx1",&pVertexFilterCutVtx1Bit,"pVertexFilterCutVtx1/I");


 	// Copy the jet trees to the output
	TTree *jetTreeOutput[nJetTrees];
	
	// Leaves for jet tree
	Int_t nJetsOutput[nJetTrees];										// number of jets in an event
	Float_t jetPhiArrayOutput[nJetTrees][nMaxJet] = {{0}};				// phi of all the jets in an event
	Float_t jetPhiArrayWTAOutput[nJetTrees][nMaxJet] = {{0}};			// phi of all the jets in an event	with WTA axis
	Float_t jetEtaArrayOutput[nJetTrees][nMaxJet] = {{0}};				// eta of all the jets in an event
	Float_t jetEtaArrayWTAOutput[nJetTrees][nMaxJet] = {{0}};			// eta of all the jets in an event	with WTA axis
	Float_t jetRawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};			// raw jet pT for all the jets in an event
	Float_t jetMaxTrackPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// maximum track pT inside a jet for all the jets in an event
	Float_t jetMassArrayOutput[nJetTrees][nMaxJet] = {{0}};				// jet mass for all the jets in an event
	Float_t jetMassCalcArrayOutput[nJetTrees][nMaxJet] = {{0}};			// jet mass calculated for all the jets in an event
	Float_t jetPfNHFArrayOutput[nJetTrees][nMaxJet] = {{0}};			// PF neutral hadron energy fraction in jets in an event
	Float_t jetPfNEFArrayOutput[nJetTrees][nMaxJet] = {{0}};			// PF neutral EM energy fraction in jets in an event
	Float_t jetPfCHFArrayOutput[nJetTrees][nMaxJet] = {{0}};			// PF charged hadron energy fraction in jets in an event
	Float_t jetPfMUFArrayOutput[nJetTrees][nMaxJet] = {{0}};			// PF muon energy fraction in jets in an event
	Float_t jetPfCEFArrayOutput[nJetTrees][nMaxJet] = {{0}};			// PF charged EM energy fraction in jets in an event
	Int_t jetPfCHMArrayOutput[nJetTrees][nMaxJet] = {{0}};				// PF charged hadron multiplicity in jets in an event
	Int_t jetPfCEMArrayOutput[nJetTrees][nMaxJet] = {{0}};				// PF charged EM multiplicity in jets in an event
	Int_t jetPfNHMArrayOutput[nJetTrees][nMaxJet] = {{0}};				// PF neutral hadron multiplicity in jets in an event
	Int_t jetPfNEMArrayOutput[nJetTrees][nMaxJet] = {{0}};				// PF neutral EM multiplicity in jets in an event
	Int_t jetPfMUMArrayOutput[nJetTrees][nMaxJet] = {{0}};				// PF muon multiplicity in jets in an event
	Float_t jetHCALSUMArrayOutput[nJetTrees][nMaxJet] = {{0}};			// HCAL energy sum (calo) in jets in an event
	Float_t jetECALSUMArrayOutput[nJetTrees][nMaxJet] = {{0}};			// ECAL energy sum (calo) in jets in an event
	Float_t jetTRKSUMArrayOutput[nJetTrees][nMaxJet] = {{0}};			// Track energy sum (calo) in jets in an event
	Int_t jetTRKNArrayOutput[nJetTrees][nMaxJet] = {{0}};				// Track number (calo) in jets in an event
	Int_t jetMUNArrayOutput[nJetTrees][nMaxJet] = {{0}};				// Muon number (calo) in jets in an event

	Float_t jetRefPtArrayOutput[nJetTrees][nMaxJet] = {{0}};			// reference generator level pT for a reconstructed jet
	Float_t jetRefEtaArrayOutput[nJetTrees][nMaxJet] = {{0}};			// reference generator level eta for a reconstructed jet
	Float_t jetRefPhiArrayOutput[nJetTrees][nMaxJet] = {{0}};			// reference generator level phi for a reconstructed jet
	Int_t jetRefFlavorArrayOutput[nJetTrees][nMaxJet] = {{0}};			// flavor for initiating parton for the reference gen jet
	Int_t jetRefFlavorForBArrayOutput[nJetTrees][nMaxJet] = {{0}};		// heavy flavor for initiating parton for the reference gen jet
	Int_t jetRefSubidArrayOutput[nJetTrees][nMaxJet] = {{0}};           // jet subid
	Float_t jetRefMassArrayOutput[nJetTrees][nMaxJet] = {{0}};			// reference generator level mass for a reconstructed jet

	// soft drop
	// zcut 0.1 and beta 0.0
	Float_t Subjet1_z0p1_b0p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the leading subjets in an event after SD
	Float_t Subjet1_z0p1_b0p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the leading subjets in an event after SD
	Float_t Subjet1_z0p1_b0p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the leading subjets in an event after SD
	Float_t Subjet1_z0p1_b0p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t Subjet1_z0p1_b0p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t Subjet1_z0p1_b0p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the leading subjets in an event after SD
	Float_t Subjet2_z0p1_b0p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the subleading subjets in an event after SD
	Float_t Subjet2_z0p1_b0p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the subleading subjets in an event after SD
	Float_t Subjet2_z0p1_b0p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the subleading subjets in an event after SD
	Float_t Subjet2_z0p1_b0p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t Subjet2_z0p1_b0p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t Subjet2_z0p1_b0p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the subleading subjets in an event after SD

	// zcut 0.1 and beta 1.0
	Float_t Subjet1_z0p1_b1p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the leading subjets in an event after SD
	Float_t Subjet1_z0p1_b1p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the leading subjets in an event after SD
	Float_t Subjet1_z0p1_b1p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the leading subjets in an event after SD
	Float_t Subjet1_z0p1_b1p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t Subjet1_z0p1_b1p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t Subjet1_z0p1_b1p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the leading subjets in an event after SD
	Float_t Subjet2_z0p1_b1p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the subleading subjets in an event after SD
	Float_t Subjet2_z0p1_b1p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the subleading subjets in an event after SD
	Float_t Subjet2_z0p1_b1p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the subleading subjets in an event after SD
	Float_t Subjet2_z0p1_b1p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t Subjet2_z0p1_b1p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t Subjet2_z0p1_b1p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the subleading subjets in an event after SD

	// zcut 0.1 and beta 2.0
	Float_t Subjet1_z0p1_b2p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the leading subjets in an event after SD
	Float_t Subjet1_z0p1_b2p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the leading subjets in an event after SD
	Float_t Subjet1_z0p1_b2p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the leading subjets in an event after SD
	Float_t Subjet1_z0p1_b2p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t Subjet1_z0p1_b2p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t Subjet1_z0p1_b2p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the leading subjets in an event after SD
	Float_t Subjet2_z0p1_b2p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the subleading subjets in an event after SD
	Float_t Subjet2_z0p1_b2p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the subleading subjets in an event after SD
	Float_t Subjet2_z0p1_b2p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the subleading subjets in an event after SD
	Float_t Subjet2_z0p1_b2p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t Subjet2_z0p1_b2p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t Subjet2_z0p1_b2p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the subleading subjets in an event after SD

	// zcut 0.2 and beta 0.0
	Float_t Subjet1_z0p2_b0p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the leading subjets in an event after SD
	Float_t Subjet1_z0p2_b0p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the leading subjets in an event after SD
	Float_t Subjet1_z0p2_b0p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the leading subjets in an event after SD
	Float_t Subjet1_z0p2_b0p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t Subjet1_z0p2_b0p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t Subjet1_z0p2_b0p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the leading subjets in an event after SD
	Float_t Subjet2_z0p2_b0p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the subleading subjets in an event after SD
	Float_t Subjet2_z0p2_b0p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the subleading subjets in an event after SD
	Float_t Subjet2_z0p2_b0p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the subleading subjets in an event after SD
	Float_t Subjet2_z0p2_b0p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t Subjet2_z0p2_b0p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t Subjet2_z0p2_b0p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the subleading subjets in an event after SD

	// zcut 0.4 and beta 0.0
	Float_t Subjet1_z0p4_b0p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the leading subjets in an event after SD
	Float_t Subjet1_z0p4_b0p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the leading subjets in an event after SD
	Float_t Subjet1_z0p4_b0p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the leading subjets in an event after SD
	Float_t Subjet1_z0p4_b0p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t Subjet1_z0p4_b0p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t Subjet1_z0p4_b0p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the leading subjets in an event after SD
	Float_t Subjet2_z0p4_b0p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the subleading subjets in an event after SD
	Float_t Subjet2_z0p4_b0p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the subleading subjets in an event after SD
	Float_t Subjet2_z0p4_b0p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the subleading subjets in an event after SD
	Float_t Subjet2_z0p4_b0p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t Subjet2_z0p4_b0p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t Subjet2_z0p4_b0p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the subleading subjets in an event after SD

	// zcut 0.5 and beta 1.0
	Float_t Subjet1_z0p5_b1p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the leading subjets in an event after SD
	Float_t Subjet1_z0p5_b1p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the leading subjets in an event after SD
	Float_t Subjet1_z0p5_b1p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the leading subjets in an event after SD
	Float_t Subjet1_z0p5_b1p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t Subjet1_z0p5_b1p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t Subjet1_z0p5_b1p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the leading subjets in an event after SD
	Float_t Subjet2_z0p5_b1p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the subleading subjets in an event after SD
	Float_t Subjet2_z0p5_b1p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the subleading subjets in an event after SD
	Float_t Subjet2_z0p5_b1p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the subleading subjets in an event after SD
	Float_t Subjet2_z0p5_b1p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t Subjet2_z0p5_b1p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t Subjet2_z0p5_b1p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the subleading subjets in an event after SD

	// zcut 0.5 and beta 1.5
	Float_t Subjet1_z0p5_b1p5_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the leading subjets in an event after SD
	Float_t Subjet1_z0p5_b1p5_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the leading subjets in an event after SD
	Float_t Subjet1_z0p5_b1p5_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the leading subjets in an event after SD
	Float_t Subjet1_z0p5_b1p5_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t Subjet1_z0p5_b1p5_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t Subjet1_z0p5_b1p5_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the leading subjets in an event after SD
	Float_t Subjet2_z0p5_b1p5_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the subleading subjets in an event after SD
	Float_t Subjet2_z0p5_b1p5_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the subleading subjets in an event after SD
	Float_t Subjet2_z0p5_b1p5_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the subleading subjets in an event after SD
	Float_t Subjet2_z0p5_b1p5_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t Subjet2_z0p5_b1p5_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t Subjet2_z0p5_b1p5_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the subleading subjets in an event after SD

	Int_t nGenJetsOutput[nJetTrees];								 	// number of generator level jets in an event
	Float_t genJetPtArrayOutput[nJetTrees][nMaxJet] = {{0}};			// pT of all the generator level jets in an event
	Float_t genJetPhiArrayOutput[nJetTrees][nMaxJet] = {{0}};			// phi of all the generator level jets in an event
	Float_t genJetPhiArrayWTAOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the generator level jets in an event with WTA axis
	Float_t genJetEtaArrayOutput[nJetTrees][nMaxJet] = {{0}};			// eta of all the generator level jets in an event
	Float_t genJetEtaArrayWTAOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the generator level jets in an event with WTA axis
	Int_t genJetSubidArrayOutput[nJetTrees][nMaxJet] = {{0}};     		// subid of all the generator level jets in an event
	Int_t genJetMatchIndexArrayOutput[nJetTrees][nMaxJet] = {{0}};		// matched index of all the generator level jets in an event
	Float_t genJetMassArrayOutput[nJetTrees][nMaxJet] = {{0}};			// jet mass for all the generator level jets in an event
	Float_t genJetMassCalcArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass calculated for all the generator level jets in an event

	// soft drop
	// zcut 0.1 and beta 0.0
	Float_t SubjetGen1_z0p1_b0p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p1_b0p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p1_b0p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p1_b0p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p1_b0p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p1_b0p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b0p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b0p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b0p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b0p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b0p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b0p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the subleading subjets in an event after SD

	// zcut 0.1 and beta 1.0
	Float_t SubjetGen1_z0p1_b1p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p1_b1p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p1_b1p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p1_b1p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p1_b1p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p1_b1p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b1p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b1p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b1p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b1p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA phi of all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b1p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};	// WTA eta of all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b1p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the subleading subjets in an event after SD

	// zcut 0.1 and beta 2.0
	Float_t SubjetGen1_z0p1_b2p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p1_b2p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p1_b2p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p1_b2p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA phi of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p1_b2p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA eta of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p1_b2p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b2p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b2p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b2p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b2p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA phi of all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b2p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA eta of all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p1_b2p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the subleading subjets in an event after SD

	// zcut 0.2 and beta 0.0
	Float_t SubjetGen1_z0p2_b0p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p2_b0p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p2_b0p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p2_b0p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA phi of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p2_b0p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA eta of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p2_b0p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p2_b0p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p2_b0p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p2_b0p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p2_b0p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA phi of all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p2_b0p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA eta of all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p2_b0p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the subleading subjets in an event after SD

	// zcut 0.4 and beta 0.0
	Float_t SubjetGen1_z0p4_b0p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p4_b0p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p4_b0p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p4_b0p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA phi of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p4_b0p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA eta of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p4_b0p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p4_b0p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p4_b0p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p4_b0p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p4_b0p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA phi of all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p4_b0p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA eta of all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p4_b0p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the subleading subjets in an event after SD

	// zcut 0.5 and beta 1.0
	Float_t SubjetGen1_z0p5_b1p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p5_b1p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p5_b1p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p5_b1p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA phi of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p5_b1p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA eta of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p5_b1p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p5_b1p0_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p5_b1p0_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p5_b1p0_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p5_b1p0_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA phi of all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p5_b1p0_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA eta of all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p5_b1p0_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the subleading subjets in an event after SD

	// zcut 0.5 and beta 1.5
	Float_t SubjetGen1_z0p5_b1p5_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p5_b1p5_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p5_b1p5_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p5_b1p5_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA phi of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p5_b1p5_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA eta of all the leading subjets in an event after SD
	Float_t SubjetGen1_z0p5_b1p5_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p5_b1p5_RawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};		// raw jet pT for all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p5_b1p5_PhiArrayOutput[nJetTrees][nMaxJet] = {{0}};		// phi of all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p5_b1p5_EtaArrayOutput[nJetTrees][nMaxJet] = {{0}};		// eta of all the subleading subjets in an event after SD
	Float_t SubjetGen2_z0p5_b1p5_PhiWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA phi of all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p5_b1p5_EtaWTAArrayOutput[nJetTrees][nMaxJet] = {{0}};		// WTA eta of all the leading subjets in an event after SD
	Float_t SubjetGen2_z0p5_b1p5_MassArrayOutput[nJetTrees][nMaxJet] = {{0}};		// jet mass for all the subleading subjets in an event after SD

	
	for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
		
		jetTreeOutput[iJetType] = new TTree("t","");		
		jetTreeOutput[iJetType]->Branch("nref",&nJetsOutput[iJetType],"nref/I");
		jetTreeOutput[iJetType]->Branch("rawpt",&jetRawPtArrayOutput[iJetType],"rawpt[nref]/F");
		jetTreeOutput[iJetType]->Branch("trackMax",&jetMaxTrackPtArrayOutput[iJetType],"trackMax[nref]/F");
		jetTreeOutput[iJetType]->Branch("jtm",&jetMassArrayOutput[iJetType],"jtm[nref]/F");
		jetTreeOutput[iJetType]->Branch("jtmcalc",&jetMassCalcArrayOutput[iJetType],"jtmcalc[nref]/F");

		// Jet eta with E-scheme and WTA axes
		jetTreeOutput[iJetType]->Branch("jtphi",&jetPhiArrayOutput[iJetType],"jtphi[nref]/F");
		jetTreeOutput[iJetType]->Branch("WTAphi",&jetPhiArrayWTAOutput[iJetType],"WTAphi[nref]/F");
		// Jet phi with E-scheme and WTA axes
		jetTreeOutput[iJetType]->Branch("jteta",&jetEtaArrayOutput[iJetType],"jteta[nref]/F");
		jetTreeOutput[iJetType]->Branch("WTAeta",&jetEtaArrayWTAOutput[iJetType],"WTAeta[nref]/F");
		
		// Jet ID stuff
		jetTreeOutput[iJetType]->Branch("jtPfNHF",&jetPfNHFArrayOutput[iJetType],"jtPfNHF[nref]/F");	
		jetTreeOutput[iJetType]->Branch("jtPfNEF",&jetPfNEFArrayOutput[iJetType],"jtPfNEF[nref]/F");	
		jetTreeOutput[iJetType]->Branch("jtPfCHF",&jetPfCHFArrayOutput[iJetType],"jtPfCHF[nref]/F");	
		jetTreeOutput[iJetType]->Branch("jtPfMUF",&jetPfMUFArrayOutput[iJetType],"jtPfMUF[nref]/F");	
		jetTreeOutput[iJetType]->Branch("jtPfCEF",&jetPfCEFArrayOutput[iJetType],"jtPfCEF[nref]/F");	
		jetTreeOutput[iJetType]->Branch("jtPfCHM",&jetPfCHMArrayOutput[iJetType],"jtPfCHM[nref]/I");	
		jetTreeOutput[iJetType]->Branch("jtPfCEM",&jetPfCEMArrayOutput[iJetType],"jtPfCEM[nref]/I");	
		jetTreeOutput[iJetType]->Branch("jtPfNHM",&jetPfNHMArrayOutput[iJetType],"jtPfNHM[nref]/I");	
		jetTreeOutput[iJetType]->Branch("jtPfNEM",&jetPfNEMArrayOutput[iJetType],"jtPfNEM[nref]/I");	
		jetTreeOutput[iJetType]->Branch("jtPfMUM",&jetPfMUMArrayOutput[iJetType],"jtPfMUM[nref]/I");
		jetTreeOutput[iJetType]->Branch("hcalSum",&jetHCALSUMArrayOutput[iJetType],"hcalSum[nref]/F");
		jetTreeOutput[iJetType]->Branch("ecalSum",&jetECALSUMArrayOutput[iJetType],"ecalSum[nref]/F");
		jetTreeOutput[iJetType]->Branch("trackSum",&jetTRKSUMArrayOutput[iJetType],"trackSum[nref]/F");
		jetTreeOutput[iJetType]->Branch("trackN",&jetTRKNArrayOutput[iJetType],"trackN[nref]/I");
		jetTreeOutput[iJetType]->Branch("muN",&jetMUNArrayOutput[iJetType],"muN[nref]/I");		
		
		if(storesoftdrop && (iJetType == 1 || iJetType == 2)){ // SoftDrop for ak4PF and akCs4PF

			jetTreeOutput[iJetType]->Branch("Subjet1_rawpt_SD_z0p1_b0p0",&Subjet1_z0p1_b0p0_RawPtArrayOutput[iJetType],"Subjet1_rawpt_SD_z0p1_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_eta_SD_z0p1_b0p0",&Subjet1_z0p1_b0p0_EtaArrayOutput[iJetType],"Subjet1_eta_SD_z0p1_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_phi_SD_z0p1_b0p0",&Subjet1_z0p1_b0p0_PhiArrayOutput[iJetType],"Subjet1_phi_SD_z0p1_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_etaWTA_SD_z0p1_b0p0",&Subjet1_z0p1_b0p0_EtaWTAArrayOutput[iJetType],"Subjet1_etaWTA_SD_z0p1_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_phiWTA_SD_z0p1_b0p0",&Subjet1_z0p1_b0p0_PhiWTAArrayOutput[iJetType],"Subjet1_phiWTA_SD_z0p1_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_mass_SD_z0p1_b0p0",&Subjet1_z0p1_b0p0_MassArrayOutput[iJetType],"Subjet1_mass_SD_z0p1_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_rawpt_SD_z0p1_b0p0",&Subjet2_z0p1_b0p0_RawPtArrayOutput[iJetType],"Subjet2_rawpt_SD_z0p1_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_eta_SD_z0p1_b0p0",&Subjet2_z0p1_b0p0_EtaArrayOutput[iJetType],"Subjet2_eta_SD_z0p1_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_phi_SD_z0p1_b0p0",&Subjet2_z0p1_b0p0_PhiArrayOutput[iJetType],"Subjet2_phi_SD_z0p1_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_etaWTA_SD_z0p1_b0p0",&Subjet2_z0p1_b0p0_EtaWTAArrayOutput[iJetType],"Subjet2_etaWTA_SD_z0p1_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_phiWTA_SD_z0p1_b0p0",&Subjet2_z0p1_b0p0_PhiWTAArrayOutput[iJetType],"Subjet2_phiWTA_SD_z0p1_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_mass_SD_z0p1_b0p0",&Subjet2_z0p1_b0p0_MassArrayOutput[iJetType],"Subjet2_mass_SD_z0p1_b0p0[nref]/F");

			jetTreeOutput[iJetType]->Branch("Subjet1_rawpt_SD_z0p1_b1p0",&Subjet1_z0p1_b1p0_RawPtArrayOutput[iJetType],"Subjet1_rawpt_SD_z0p1_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_eta_SD_z0p1_b1p0",&Subjet1_z0p1_b1p0_EtaArrayOutput[iJetType],"Subjet1_eta_SD_z0p1_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_phi_SD_z0p1_b1p0",&Subjet1_z0p1_b1p0_PhiArrayOutput[iJetType],"Subjet1_phi_SD_z0p1_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_etaWTA_SD_z0p1_b1p0",&Subjet1_z0p1_b1p0_EtaWTAArrayOutput[iJetType],"Subjet1_etaWTA_SD_z0p1_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_phiWTA_SD_z0p1_b1p0",&Subjet1_z0p1_b1p0_PhiWTAArrayOutput[iJetType],"Subjet1_phiWTA_SD_z0p1_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_mass_SD_z0p1_b1p0",&Subjet1_z0p1_b1p0_MassArrayOutput[iJetType],"Subjet1_mass_SD_z0p1_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_rawpt_SD_z0p1_b1p0",&Subjet2_z0p1_b1p0_RawPtArrayOutput[iJetType],"Subjet2_rawpt_SD_z0p1_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_eta_SD_z0p1_b1p0",&Subjet2_z0p1_b1p0_EtaArrayOutput[iJetType],"Subjet2_eta_SD_z0p1_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_phi_SD_z0p1_b1p0",&Subjet2_z0p1_b1p0_PhiArrayOutput[iJetType],"Subjet2_phi_SD_z0p1_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_etaWTA_SD_z0p1_b1p0",&Subjet2_z0p1_b1p0_EtaWTAArrayOutput[iJetType],"Subjet2_etaWTA_SD_z0p1_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_phiWTA_SD_z0p1_b1p0",&Subjet2_z0p1_b1p0_PhiWTAArrayOutput[iJetType],"Subjet2_phiWTA_SD_z0p1_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_mass_SD_z0p1_b1p0",&Subjet2_z0p1_b1p0_MassArrayOutput[iJetType],"Subjet2_mass_SD_z0p1_b1p0[nref]/F");

			jetTreeOutput[iJetType]->Branch("Subjet1_rawpt_SD_z0p1_b2p0",&Subjet1_z0p1_b2p0_RawPtArrayOutput[iJetType],"Subjet1_rawpt_SD_z0p1_b2p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_eta_SD_z0p1_b2p0",&Subjet1_z0p1_b2p0_EtaArrayOutput[iJetType],"Subjet1_eta_SD_z0p1_b2p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_phi_SD_z0p1_b2p0",&Subjet1_z0p1_b2p0_PhiArrayOutput[iJetType],"Subjet1_phi_SD_z0p1_b2p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_etaWTA_SD_z0p1_b2p0",&Subjet1_z0p1_b2p0_EtaWTAArrayOutput[iJetType],"Subjet1_etaWTA_SD_z0p1_b2p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_phiWTA_SD_z0p1_b2p0",&Subjet1_z0p1_b2p0_PhiWTAArrayOutput[iJetType],"Subjet1_phiWTA_SD_z0p1_b2p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_mass_SD_z0p1_b2p0",&Subjet1_z0p1_b2p0_MassArrayOutput[iJetType],"Subjet1_mass_SD_z0p1_b2p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_rawpt_SD_z0p1_b2p0",&Subjet2_z0p1_b2p0_RawPtArrayOutput[iJetType],"Subjet2_rawpt_SD_z0p1_b2p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_eta_SD_z0p1_b2p0",&Subjet2_z0p1_b2p0_EtaArrayOutput[iJetType],"Subjet2_eta_SD_z0p1_b2p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_phi_SD_z0p1_b2p0",&Subjet2_z0p1_b2p0_PhiArrayOutput[iJetType],"Subjet2_phi_SD_z0p1_b2p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_etaWTA_SD_z0p1_b2p0",&Subjet2_z0p1_b2p0_EtaWTAArrayOutput[iJetType],"Subjet2_etaWTA_SD_z0p1_b2p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_phiWTA_SD_z0p1_b2p0",&Subjet2_z0p1_b2p0_PhiWTAArrayOutput[iJetType],"Subjet2_phiWTA_SD_z0p1_b2p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_mass_SD_z0p1_b2p0",&Subjet2_z0p1_b2p0_MassArrayOutput[iJetType],"Subjet2_mass_SD_z0p1_b2p0[nref]/F");

			jetTreeOutput[iJetType]->Branch("Subjet1_rawpt_SD_z0p2_b0p0",&Subjet1_z0p2_b0p0_RawPtArrayOutput[iJetType],"Subjet1_rawpt_SD_z0p2_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_eta_SD_z0p2_b0p0",&Subjet1_z0p2_b0p0_EtaArrayOutput[iJetType],"Subjet1_eta_SD_z0p2_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_phi_SD_z0p2_b0p0",&Subjet1_z0p2_b0p0_PhiArrayOutput[iJetType],"Subjet1_phi_SD_z0p2_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_etaWTA_SD_z0p2_b0p0",&Subjet1_z0p2_b0p0_EtaWTAArrayOutput[iJetType],"Subjet1_etaWTA_SD_z0p2_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_phiWTA_SD_z0p2_b0p0",&Subjet1_z0p2_b0p0_PhiWTAArrayOutput[iJetType],"Subjet1_phiWTA_SD_z0p2_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_mass_SD_z0p2_b0p0",&Subjet1_z0p2_b0p0_MassArrayOutput[iJetType],"Subjet1_mass_SD_z0p2_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_rawpt_SD_z0p2_b0p0",&Subjet2_z0p2_b0p0_RawPtArrayOutput[iJetType],"Subjet2_rawpt_SD_z0p2_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_eta_SD_z0p2_b0p0",&Subjet2_z0p2_b0p0_EtaArrayOutput[iJetType],"Subjet2_eta_SD_z0p2_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_phi_SD_z0p2_b0p0",&Subjet2_z0p2_b0p0_PhiArrayOutput[iJetType],"Subjet2_phi_SD_z0p2_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_etaWTA_SD_z0p2_b0p0",&Subjet2_z0p2_b0p0_EtaWTAArrayOutput[iJetType],"Subjet2_etaWTA_SD_z0p2_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_phiWTA_SD_z0p2_b0p0",&Subjet2_z0p2_b0p0_PhiWTAArrayOutput[iJetType],"Subjet2_phiWTA_SD_z0p2_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_mass_SD_z0p2_b0p0",&Subjet2_z0p2_b0p0_MassArrayOutput[iJetType],"Subjet2_mass_SD_z0p2_b0p0[nref]/F");

			jetTreeOutput[iJetType]->Branch("Subjet1_rawpt_SD_z0p4_b0p0",&Subjet1_z0p4_b0p0_RawPtArrayOutput[iJetType],"Subjet1_rawpt_SD_z0p4_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_eta_SD_z0p4_b0p0",&Subjet1_z0p4_b0p0_EtaArrayOutput[iJetType],"Subjet1_eta_SD_z0p4_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_phi_SD_z0p4_b0p0",&Subjet1_z0p4_b0p0_PhiArrayOutput[iJetType],"Subjet1_phi_SD_z0p4_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_etaWTA_SD_z0p4_b0p0",&Subjet1_z0p4_b0p0_EtaWTAArrayOutput[iJetType],"Subjet1_etaWTA_SD_z0p4_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_phiWTA_SD_z0p4_b0p0",&Subjet1_z0p4_b0p0_PhiWTAArrayOutput[iJetType],"Subjet1_phiWTA_SD_z0p4_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_mass_SD_z0p4_b0p0",&Subjet1_z0p4_b0p0_MassArrayOutput[iJetType],"Subjet1_mass_SD_z0p4_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_rawpt_SD_z0p4_b0p0",&Subjet2_z0p4_b0p0_RawPtArrayOutput[iJetType],"Subjet2_rawpt_SD_z0p4_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_eta_SD_z0p4_b0p0",&Subjet2_z0p4_b0p0_EtaArrayOutput[iJetType],"Subjet2_eta_SD_z0p4_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_phi_SD_z0p4_b0p0",&Subjet2_z0p4_b0p0_PhiArrayOutput[iJetType],"Subjet2_phi_SD_z0p4_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_etaWTA_SD_z0p4_b0p0",&Subjet2_z0p4_b0p0_EtaWTAArrayOutput[iJetType],"Subjet2_etaWTA_SD_z0p4_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_phiWTA_SD_z0p4_b0p0",&Subjet2_z0p4_b0p0_PhiWTAArrayOutput[iJetType],"Subjet2_phiWTA_SD_z0p4_b0p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_mass_SD_z0p4_b0p0",&Subjet2_z0p4_b0p0_MassArrayOutput[iJetType],"Subjet2_mass_SD_z0p4_b0p0[nref]/F");

			jetTreeOutput[iJetType]->Branch("Subjet1_rawpt_SD_z0p5_b1p0",&Subjet1_z0p5_b1p0_RawPtArrayOutput[iJetType],"Subjet1_rawpt_SD_z0p5_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_eta_SD_z0p5_b1p0",&Subjet1_z0p5_b1p0_EtaArrayOutput[iJetType],"Subjet1_eta_SD_z0p5_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_phi_SD_z0p5_b1p0",&Subjet1_z0p5_b1p0_PhiArrayOutput[iJetType],"Subjet1_phi_SD_z0p5_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_etaWTA_SD_z0p5_b1p0",&Subjet1_z0p5_b1p0_EtaWTAArrayOutput[iJetType],"Subjet1_etaWTA_SD_z0p5_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_phiWTA_SD_z0p5_b1p0",&Subjet1_z0p5_b1p0_PhiWTAArrayOutput[iJetType],"Subjet1_phiWTA_SD_z0p5_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_mass_SD_z0p5_b1p0",&Subjet1_z0p5_b1p0_MassArrayOutput[iJetType],"Subjet1_mass_SD_z0p5_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_rawpt_SD_z0p5_b1p0",&Subjet2_z0p5_b1p0_RawPtArrayOutput[iJetType],"Subjet2_rawpt_SD_z0p5_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_eta_SD_z0p5_b1p0",&Subjet2_z0p5_b1p0_EtaArrayOutput[iJetType],"Subjet2_eta_SD_z0p5_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_phi_SD_z0p5_b1p0",&Subjet2_z0p5_b1p0_PhiArrayOutput[iJetType],"Subjet2_phi_SD_z0p5_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_etaWTA_SD_z0p5_b1p0",&Subjet2_z0p5_b1p0_EtaWTAArrayOutput[iJetType],"Subjet2_etaWTA_SD_z0p5_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_phiWTA_SD_z0p5_b1p0",&Subjet2_z0p5_b1p0_PhiWTAArrayOutput[iJetType],"Subjet2_phiWTA_SD_z0p5_b1p0[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_mass_SD_z0p5_b1p0",&Subjet2_z0p5_b1p0_MassArrayOutput[iJetType],"Subjet2_mass_SD_z0p5_b1p0[nref]/F");

			jetTreeOutput[iJetType]->Branch("Subjet1_rawpt_SD_z0p5_b1p5",&Subjet1_z0p5_b1p5_RawPtArrayOutput[iJetType],"Subjet1_rawpt_SD_z0p5_b1p5[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_eta_SD_z0p5_b1p5",&Subjet1_z0p5_b1p5_EtaArrayOutput[iJetType],"Subjet1_eta_SD_z0p5_b1p5[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_phi_SD_z0p5_b1p5",&Subjet1_z0p5_b1p5_PhiArrayOutput[iJetType],"Subjet1_phi_SD_z0p5_b1p5[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_etaWTA_SD_z0p5_b1p5",&Subjet1_z0p5_b1p5_EtaWTAArrayOutput[iJetType],"Subjet1_etaWTA_SD_z0p5_b1p5[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_phiWTA_SD_z0p5_b1p5",&Subjet1_z0p5_b1p5_PhiWTAArrayOutput[iJetType],"Subjet1_phiWTA_SD_z0p5_b1p5[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet1_mass_SD_z0p5_b1p5",&Subjet1_z0p5_b1p5_MassArrayOutput[iJetType],"Subjet1_mass_SD_z0p5_b1p5[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_rawpt_SD_z0p5_b1p5",&Subjet2_z0p5_b1p5_RawPtArrayOutput[iJetType],"Subjet2_rawpt_SD_z0p5_b1p5[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_eta_SD_z0p5_b1p5",&Subjet2_z0p5_b1p5_EtaArrayOutput[iJetType],"Subjet2_eta_SD_z0p5_b1p5[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_phi_SD_z0p5_b1p5",&Subjet2_z0p5_b1p5_PhiArrayOutput[iJetType],"Subjet2_phi_SD_z0p5_b1p5[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_etaWTA_SD_z0p5_b1p5",&Subjet2_z0p5_b1p5_EtaWTAArrayOutput[iJetType],"Subjet2_etaWTA_SD_z0p5_b1p5[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_phiWTA_SD_z0p5_b1p5",&Subjet2_z0p5_b1p5_PhiWTAArrayOutput[iJetType],"Subjet2_phiWTA_SD_z0p5_b1p5[nref]/F");
			jetTreeOutput[iJetType]->Branch("Subjet2_mass_SD_z0p5_b1p5",&Subjet2_z0p5_b1p5_MassArrayOutput[iJetType],"Subjet2_mass_SD_z0p5_b1p5[nref]/F");
		
		}
		
		// If we are looking at Monte Carlo, connect the reference pT and parton arrays
		if(is_MC){
			jetTreeOutput[iJetType]->Branch("refpt",&jetRefPtArrayOutput[iJetType],"refpt[nref]/F");
			jetTreeOutput[iJetType]->Branch("refeta",&jetRefEtaArrayOutput[iJetType],"refeta[nref]/F");
			jetTreeOutput[iJetType]->Branch("refphi",&jetRefPhiArrayOutput[iJetType],"refphi[nref]/F");
			jetTreeOutput[iJetType]->Branch("refparton_flavor", &jetRefFlavorArrayOutput[iJetType], "refparton_flavor[nref]/I");
			jetTreeOutput[iJetType]->Branch("refparton_flavorForB", &jetRefFlavorForBArrayOutput[iJetType], "refparton_flavorForB[nref]/I");
			jetTreeOutput[iJetType]->Branch("subid", &jetRefSubidArrayOutput[iJetType], "subid[nref]/I");
			jetTreeOutput[iJetType]->Branch("refm",&jetRefMassArrayOutput[iJetType],"refm[nref]/F");
	
			jetTreeOutput[iJetType]->Branch("ngen",&nGenJetsOutput[iJetType],"ngen/I");
			jetTreeOutput[iJetType]->Branch("genpt",&genJetPtArrayOutput[iJetType],"genpt[ngen]/F");		
			// Gen jet phi for e-scheme and WTA axes
			jetTreeOutput[iJetType]->Branch("genphi",&genJetPhiArrayOutput[iJetType],"genphi[ngen]/F");
			jetTreeOutput[iJetType]->Branch("WTAgenphi",&genJetPhiArrayWTAOutput[iJetType],"WTAgenphi[ngen]/F");
			// Gen jet eta for e-scheme and WTA axes
			jetTreeOutput[iJetType]->Branch("geneta",&genJetEtaArrayOutput[iJetType],"geneta[ngen]/F");
			jetTreeOutput[iJetType]->Branch("WTAgeneta",&genJetEtaArrayWTAOutput[iJetType],"WTAgeneta[ngen]/F");
			// Gen match and subid
			jetTreeOutput[iJetType]->Branch("genmatchindex",&genJetMatchIndexArrayOutput[iJetType],"genmatchindex[ngen]/F");
			jetTreeOutput[iJetType]->Branch("gensubid",&genJetSubidArrayOutput[iJetType],"gensubid[ngen]/F");
			jetTreeOutput[iJetType]->Branch("genm",&genJetMassArrayOutput[iJetType],"genm[ngen]/F");
			jetTreeOutput[iJetType]->Branch("genmcalc",&genJetMassCalcArrayOutput[iJetType],"genmcalc[ngen]/F");
		
			if(storesoftdrop && (iJetType == 1 || iJetType == 2)){

				jetTreeOutput[iJetType]->Branch("SubjetGen1_pt_SD_z0p1_b0p0",&SubjetGen1_z0p1_b0p0_RawPtArrayOutput[iJetType],"SubjetGen1_pt_SD_z0p1_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_eta_SD_z0p1_b0p0",&SubjetGen1_z0p1_b0p0_EtaArrayOutput[iJetType],"SubjetGen1_eta_SD_z0p1_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_phi_SD_z0p1_b0p0",&SubjetGen1_z0p1_b0p0_PhiArrayOutput[iJetType],"SubjetGen1_phi_SD_z0p1_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_etaWTA_SD_z0p1_b0p0",&SubjetGen1_z0p1_b0p0_EtaWTAArrayOutput[iJetType],"SubjetGen1_etaWTA_SD_z0p1_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_phiWTA_SD_z0p1_b0p0",&SubjetGen1_z0p1_b0p0_PhiWTAArrayOutput[iJetType],"SubjetGen1_phiWTA_SD_z0p1_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_mass_SD_z0p1_b0p0",&SubjetGen1_z0p1_b0p0_MassArrayOutput[iJetType],"SubjetGen1_mass_SD_z0p1_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_pt_SD_z0p1_b0p0",&SubjetGen2_z0p1_b0p0_RawPtArrayOutput[iJetType],"SubjetGen2_pt_SD_z0p1_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_eta_SD_z0p1_b0p0",&SubjetGen2_z0p1_b0p0_EtaArrayOutput[iJetType],"SubjetGen2_eta_SD_z0p1_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_phi_SD_z0p1_b0p0",&SubjetGen2_z0p1_b0p0_PhiArrayOutput[iJetType],"SubjetGen2_phi_SD_z0p1_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_etaWTA_SD_z0p1_b0p0",&SubjetGen2_z0p1_b0p0_EtaWTAArrayOutput[iJetType],"SubjetGen2_etaWTA_SD_z0p1_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_phiWTA_SD_z0p1_b0p0",&SubjetGen2_z0p1_b0p0_PhiWTAArrayOutput[iJetType],"SubjetGen2_phiWTA_SD_z0p1_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_mass_SD_z0p1_b0p0",&SubjetGen2_z0p1_b0p0_MassArrayOutput[iJetType],"SubjetGen2_mass_SD_z0p1_b0p0[ngen]/F");

				jetTreeOutput[iJetType]->Branch("SubjetGen1_pt_SD_z0p1_b1p0",&SubjetGen1_z0p1_b1p0_RawPtArrayOutput[iJetType],"SubjetGen1_pt_SD_z0p1_b1p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_eta_SD_z0p1_b1p0",&SubjetGen1_z0p1_b1p0_EtaArrayOutput[iJetType],"SubjetGen1_eta_SD_z0p1_b1p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_phi_SD_z0p1_b1p0",&SubjetGen1_z0p1_b1p0_PhiArrayOutput[iJetType],"SubjetGen1_phi_SD_z0p1_b1p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_etaWTA_SD_z0p1_b1p0",&SubjetGen1_z0p1_b1p0_EtaWTAArrayOutput[iJetType],"SubjetGen1_etaWTA_SD_z0p1_b1p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_phiWTA_SD_z0p1_b1p0",&SubjetGen1_z0p1_b1p0_PhiWTAArrayOutput[iJetType],"SubjetGen1_phiWTA_SD_z0p1_b1p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_mass_SD_z0p1_b1p0",&SubjetGen1_z0p1_b1p0_MassArrayOutput[iJetType],"SubjetGen1_mass_SD_z0p1_b1p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_pt_SD_z0p1_b1p0",&SubjetGen2_z0p1_b1p0_RawPtArrayOutput[iJetType],"SubjetGen2_pt_SD_z0p1_b1p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_eta_SD_z0p1_b1p0",&SubjetGen2_z0p1_b1p0_EtaArrayOutput[iJetType],"SubjetGen2_eta_SD_z0p1_b1p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_phi_SD_z0p1_b1p0",&SubjetGen2_z0p1_b1p0_PhiArrayOutput[iJetType],"SubjetGen2_phi_SD_z0p1_b1p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_etaWTA_SD_z0p1_b1p0",&SubjetGen2_z0p1_b1p0_EtaWTAArrayOutput[iJetType],"SubjetGen2_etaWTA_SD_z0p1_b1p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_phiWTA_SD_z0p1_b1p0",&SubjetGen2_z0p1_b1p0_PhiWTAArrayOutput[iJetType],"SubjetGen2_phiWTA_SD_z0p1_b1p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_mass_SD_z0p1_b1p0",&SubjetGen2_z0p1_b1p0_MassArrayOutput[iJetType],"SubjetGen2_mass_SD_z0p1_b1p0[ngen]/F");

				jetTreeOutput[iJetType]->Branch("SubjetGen1_pt_SD_z0p1_b2p0",&SubjetGen1_z0p1_b2p0_RawPtArrayOutput[iJetType],"SubjetGen1_pt_SD_z0p1_b2p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_eta_SD_z0p1_b2p0",&SubjetGen1_z0p1_b2p0_EtaArrayOutput[iJetType],"SubjetGen1_eta_SD_z0p1_b2p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_phi_SD_z0p1_b2p0",&SubjetGen1_z0p1_b2p0_PhiArrayOutput[iJetType],"SubjetGen1_phi_SD_z0p1_b2p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_etaWTA_SD_z0p1_b2p0",&SubjetGen1_z0p1_b2p0_EtaWTAArrayOutput[iJetType],"SubjetGen1_etaWTA_SD_z0p1_b2p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_phiWTA_SD_z0p1_b2p0",&SubjetGen1_z0p1_b2p0_PhiWTAArrayOutput[iJetType],"SubjetGen1_phiWTA_SD_z0p1_b2p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_mass_SD_z0p1_b2p0",&SubjetGen1_z0p1_b2p0_MassArrayOutput[iJetType],"SubjetGen1_mass_SD_z0p1_b2p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_pt_SD_z0p1_b2p0",&SubjetGen2_z0p1_b2p0_RawPtArrayOutput[iJetType],"SubjetGen2_pt_SD_z0p1_b2p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_eta_SD_z0p1_b2p0",&SubjetGen2_z0p1_b2p0_EtaArrayOutput[iJetType],"SubjetGen2_eta_SD_z0p1_b2p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_phi_SD_z0p1_b2p0",&SubjetGen2_z0p1_b2p0_PhiArrayOutput[iJetType],"SubjetGen2_phi_SD_z0p1_b2p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_etaWTA_SD_z0p1_b2p0",&SubjetGen2_z0p1_b2p0_EtaWTAArrayOutput[iJetType],"SubjetGen2_etaWTA_SD_z0p1_b2p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_phiWTA_SD_z0p1_b2p0",&SubjetGen2_z0p1_b2p0_PhiWTAArrayOutput[iJetType],"SubjetGen2_phiWTA_SD_z0p1_b2p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_mass_SD_z0p1_b2p0",&SubjetGen2_z0p1_b2p0_MassArrayOutput[iJetType],"SubjetGen2_mass_SD_z0p1_b2p0[ngen]/F");

				jetTreeOutput[iJetType]->Branch("SubjetGen1_pt_SD_z0p2_b0p0",&SubjetGen1_z0p2_b0p0_RawPtArrayOutput[iJetType],"SubjetGen1_pt_SD_z0p2_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_eta_SD_z0p2_b0p0",&SubjetGen1_z0p2_b0p0_EtaArrayOutput[iJetType],"SubjetGen1_eta_SD_z0p2_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_phi_SD_z0p2_b0p0",&SubjetGen1_z0p2_b0p0_PhiArrayOutput[iJetType],"SubjetGen1_phi_SD_z0p2_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_etaWTA_SD_z0p2_b0p0",&SubjetGen1_z0p2_b0p0_EtaWTAArrayOutput[iJetType],"SubjetGen1_etaWTA_SD_z0p2_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_phiWTA_SD_z0p2_b0p0",&SubjetGen1_z0p2_b0p0_PhiWTAArrayOutput[iJetType],"SubjetGen1_phiWTA_SD_z0p2_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_mass_SD_z0p2_b0p0",&SubjetGen1_z0p2_b0p0_MassArrayOutput[iJetType],"SubjetGen1_mass_SD_z0p2_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_pt_SD_z0p2_b0p0",&SubjetGen2_z0p2_b0p0_RawPtArrayOutput[iJetType],"SubjetGen2_pt_SD_z0p2_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_eta_SD_z0p2_b0p0",&SubjetGen2_z0p2_b0p0_EtaArrayOutput[iJetType],"SubjetGen2_eta_SD_z0p2_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_phi_SD_z0p2_b0p0",&SubjetGen2_z0p2_b0p0_PhiArrayOutput[iJetType],"SubjetGen2_phi_SD_z0p2_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_etaWTA_SD_z0p2_b0p0",&SubjetGen2_z0p2_b0p0_EtaWTAArrayOutput[iJetType],"SubjetGen2_etaWTA_SD_z0p2_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_phiWTA_SD_z0p2_b0p0",&SubjetGen2_z0p2_b0p0_PhiWTAArrayOutput[iJetType],"SubjetGen2_phiWTA_SD_z0p2_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_mass_SD_z0p2_b0p0",&SubjetGen2_z0p2_b0p0_MassArrayOutput[iJetType],"SubjetGen2_mass_SD_z0p2_b0p0[ngen]/F");

				jetTreeOutput[iJetType]->Branch("SubjetGen1_pt_SD_z0p4_b0p0",&SubjetGen1_z0p4_b0p0_RawPtArrayOutput[iJetType],"SubjetGen1_pt_SD_z0p4_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_eta_SD_z0p4_b0p0",&SubjetGen1_z0p4_b0p0_EtaArrayOutput[iJetType],"SubjetGen1_eta_SD_z0p4_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_phi_SD_z0p4_b0p0",&SubjetGen1_z0p4_b0p0_PhiArrayOutput[iJetType],"SubjetGen1_phi_SD_z0p4_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_etaWTA_SD_z0p4_b0p0",&SubjetGen1_z0p4_b0p0_EtaWTAArrayOutput[iJetType],"SubjetGen1_etaWTA_SD_z0p4_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_phiWTA_SD_z0p4_b0p0",&SubjetGen1_z0p4_b0p0_PhiWTAArrayOutput[iJetType],"SubjetGen1_phiWTA_SD_z0p4_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_mass_SD_z0p4_b0p0",&SubjetGen1_z0p4_b0p0_MassArrayOutput[iJetType],"SubjetGen1_mass_SD_z0p4_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_pt_SD_z0p4_b0p0",&SubjetGen2_z0p4_b0p0_RawPtArrayOutput[iJetType],"SubjetGen2_pt_SD_z0p4_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_eta_SD_z0p4_b0p0",&SubjetGen2_z0p4_b0p0_EtaArrayOutput[iJetType],"SubjetGen2_eta_SD_z0p4_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_phi_SD_z0p4_b0p0",&SubjetGen2_z0p4_b0p0_PhiArrayOutput[iJetType],"SubjetGen2_phi_SD_z0p4_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_etaWTA_SD_z0p4_b0p0",&SubjetGen2_z0p4_b0p0_EtaWTAArrayOutput[iJetType],"SubjetGen2_etaWTA_SD_z0p4_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_phiWTA_SD_z0p4_b0p0",&SubjetGen2_z0p4_b0p0_PhiWTAArrayOutput[iJetType],"SubjetGen2_phiWTA_SD_z0p4_b0p0[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_mass_SD_z0p4_b0p0",&SubjetGen2_z0p4_b0p0_MassArrayOutput[iJetType],"SubjetGen2_mass_SD_z0p4_b0p0[ngen]/F");

				jetTreeOutput[iJetType]->Branch("SubjetGen1_rawpt_SD_z0p5_b1p0",&SubjetGen1_z0p5_b1p0_RawPtArrayOutput[iJetType],"SubjetGen1_rawpt_SD_z0p5_b1p0[nref]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_eta_SD_z0p5_b1p0",&SubjetGen1_z0p5_b1p0_EtaArrayOutput[iJetType],"SubjetGen1_eta_SD_z0p5_b1p0[nref]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_phi_SD_z0p5_b1p0",&SubjetGen1_z0p5_b1p0_PhiArrayOutput[iJetType],"SubjetGen1_phi_SD_z0p5_b1p0[nref]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_etaWTA_SD_z0p5_b1p0",&SubjetGen1_z0p5_b1p0_EtaWTAArrayOutput[iJetType],"SubjetGen1_etaWTA_SD_z0p5_b1p0[nref]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_phiWTA_SD_z0p5_b1p0",&SubjetGen1_z0p5_b1p0_PhiWTAArrayOutput[iJetType],"SubjetGen1_phiWTA_SD_z0p5_b1p0[nref]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_mass_SD_z0p5_b1p0",&SubjetGen1_z0p5_b1p0_MassArrayOutput[iJetType],"SubjetGen1_mass_SD_z0p5_b1p0[nref]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_rawpt_SD_z0p5_b1p0",&SubjetGen2_z0p5_b1p0_RawPtArrayOutput[iJetType],"SubjetGen2_rawpt_SD_z0p5_b1p0[nref]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_eta_SD_z0p5_b1p0",&SubjetGen2_z0p5_b1p0_EtaArrayOutput[iJetType],"SubjetGen2_eta_SD_z0p5_b1p0[nref]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_phi_SD_z0p5_b1p0",&SubjetGen2_z0p5_b1p0_PhiArrayOutput[iJetType],"SubjetGen2_phi_SD_z0p5_b1p0[nref]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_etaWTA_SD_z0p5_b1p0",&SubjetGen2_z0p5_b1p0_EtaWTAArrayOutput[iJetType],"SubjetGen2_etaWTA_SD_z0p5_b1p0[nref]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_phiWTA_SD_z0p5_b1p0",&SubjetGen2_z0p5_b1p0_PhiWTAArrayOutput[iJetType],"SubjetGen2_phiWTA_SD_z0p5_b1p0[nref]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_mass_SD_z0p5_b1p0",&SubjetGen2_z0p5_b1p0_MassArrayOutput[iJetType],"SubjetGen2_mass_SD_z0p5_b1p0[nref]/F");

				jetTreeOutput[iJetType]->Branch("SubjetGen1_pt_SD_z0p5_b1p5",&SubjetGen1_z0p5_b1p5_RawPtArrayOutput[iJetType],"SubjetGen1_pt_SD_z0p5_b1p5[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_eta_SD_z0p5_b1p5",&SubjetGen1_z0p5_b1p5_EtaArrayOutput[iJetType],"SubjetGen1_eta_SD_z0p5_b1p5[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_phi_SD_z0p5_b1p5",&SubjetGen1_z0p5_b1p5_PhiArrayOutput[iJetType],"SubjetGen1_phi_SD_z0p5_b1p5[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_etaWTA_SD_z0p5_b1p5",&SubjetGen1_z0p5_b1p5_EtaWTAArrayOutput[iJetType],"SubjetGen1_etaWTA_SD_z0p5_b1p5[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_phiWTA_SD_z0p5_b1p5",&SubjetGen1_z0p5_b1p5_PhiWTAArrayOutput[iJetType],"SubjetGen1_phiWTA_SD_z0p5_b1p5[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen1_mass_SD_z0p5_b1p5",&SubjetGen1_z0p5_b1p5_MassArrayOutput[iJetType],"SubjetGen1_mass_SD_z0p5_b1p5[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_pt_SD_z0p5_b1p5",&SubjetGen2_z0p5_b1p5_RawPtArrayOutput[iJetType],"SubjetGen2_pt_SD_z0p5_b1p5[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_eta_SD_z0p5_b1p5",&SubjetGen2_z0p5_b1p5_EtaArrayOutput[iJetType],"SubjetGen2_eta_SD_z0p5_b1p5[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_phi_SD_z0p5_b1p5",&SubjetGen2_z0p5_b1p5_PhiArrayOutput[iJetType],"SubjetGen2_phi_SD_z0p5_b1p5[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_etaWTA_SD_z0p5_b1p5",&SubjetGen2_z0p5_b1p5_EtaWTAArrayOutput[iJetType],"SubjetGen2_etaWTA_SD_z0p5_b1p5[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_phiWTA_SD_z0p5_b1p5",&SubjetGen2_z0p5_b1p5_PhiWTAArrayOutput[iJetType],"SubjetGen2_phiWTA_SD_z0p5_b1p5[ngen]/F");
				jetTreeOutput[iJetType]->Branch("SubjetGen2_mass_SD_z0p5_b1p5",&SubjetGen2_z0p5_b1p5_MassArrayOutput[iJetType],"SubjetGen2_mass_SD_z0p5_b1p5[ngen]/F");
		
			}
		} // Branches only for MC

	} // Jet type loop

	
	// Copy the track trees to the output
	TTree *trackTreeOutput = new TTree("trackTree","");
	Int_t nTracksOutput;										 // Number of tracks
	Float_t trackPtOutput[nMaxTrack] = {0};		 				 // Array for track pT:s
	Float_t trackPtErrorOutput[nMaxTrack] = {0};				 // Array for track pT errors
	Float_t trackPhiOutput[nMaxTrack] = {0};					 // Array for track phis
	Float_t trackEtaOutput[nMaxTrack] = {0};					 // Array for track etas
	Bool_t trackHighPurityOutput[nMaxTrack] = {0};				 // Array for the high purity of tracks
	Float_t trackVertexDistanceZOutput[nMaxTrack] = {0};		 // Array for track distance from primary vertex in z-direction
	Float_t trackVertexDistanceZErrorOutput[nMaxTrack] = {0};	 // Array for error for track distance from primary vertex in z-direction
	Float_t trackVertexDistanceXYOutput[nMaxTrack] = {0};		 // Array for track distance from primary vertex in xy-direction
	Float_t trackVertexDistanceXYErrorOutput[nMaxTrack] = {0};	 // Array for error for track distance from primary vertex in xy-direction
	Float_t trackEnergyEcalOutput[nMaxTrack] = {0};							// Array for track energy in ECal
	Float_t trackEnergyHcalOutput[nMaxTrack] = {0};							// Array for track energy in HCal
	UChar_t PixelnHitsTrackOutput[nMaxTrack] = {0};				 // Array for number of hits for the track
	Int_t trackChargeOutput[nMaxTrack] = {0};					 // Array for track charge
	
	trackTreeOutput->Branch("nTrk",&nTracksOutput,"nTrk/I");
	trackTreeOutput->Branch("trkPt",&trackPtOutput,"trkPt[nTrk]/F");
	trackTreeOutput->Branch("trkPtError",&trackPtErrorOutput,"trkPtError[nTrk]/F");
	trackTreeOutput->Branch("trkPhi",&trackPhiOutput,"trkPhi[nTrk]/F");
	trackTreeOutput->Branch("trkEta",&trackEtaOutput,"trkEta[nTrk]/F");
	trackTreeOutput->Branch("highPurity",&trackHighPurityOutput,"highPurity[nTrk]/O");
	trackTreeOutput->Branch("trkDz1",&trackVertexDistanceZOutput,"trkDz1[nTrk]/F");
	trackTreeOutput->Branch("trkDzError1",&trackVertexDistanceZErrorOutput,"trkDzError1[nTrk]/F");
	trackTreeOutput->Branch("trkDxy1",&trackVertexDistanceXYOutput,"trkDxy1[nTrk]/F");
	trackTreeOutput->Branch("trkDxyError1",&trackVertexDistanceXYErrorOutput,"trkDxyError1[nTrk]/F");
	trackTreeOutput->Branch("trkNPixelHit",&PixelnHitsTrackOutput,"trkNPixelHit[nTrk]/b");
	trackTreeOutput->Branch("pfEcal",&trackEnergyEcalOutput,"pfEcal[nTrk]/F");
	trackTreeOutput->Branch("pfHcal",&trackEnergyHcalOutput,"pfHcal[nTrk]/F");
	trackTreeOutput->Branch("trkCharge",&trackChargeOutput,"trkCharge[nTrk]/I");
	
	// Generator level tracks only in Monte Carlo
	TTree *genTrackTreeOutput = new TTree("hi","");
	std::vector<float> *genTrackPtVector = new std::vector<float>(); genTrackPtVector->clear();
	std::vector<float> *genTrackPhiVector = new std::vector<float>(); genTrackPhiVector->clear();
	std::vector<float> *genTrackEtaVector = new std::vector<float>(); genTrackEtaVector->clear();
	std::vector<int> *genTrackPdgVector = new std::vector<int>(); genTrackPdgVector->clear();
	std::vector<int> *genTrackChargeVector = new std::vector<int>(); genTrackChargeVector->clear();
	std::vector<int> *genTrackSubeventVector = new std::vector<int>(); genTrackSubeventVector->clear();
	// Connect the branches to generator level track tree
	if(is_MC){
		genTrackTreeOutput->Branch("pt","vector<float>", &genTrackPtVector);
		genTrackTreeOutput->Branch("phi","vector<float>", &genTrackPhiVector);
		genTrackTreeOutput->Branch("eta","vector<float>", &genTrackEtaVector);
		genTrackTreeOutput->Branch("pdg","vector<int>", &genTrackPdgVector);
		genTrackTreeOutput->Branch("chg","vector<int>", &genTrackChargeVector);
		genTrackTreeOutput->Branch("sube","vector<int>", &genTrackSubeventVector);
	}
	
	// Copy the track trees to the output
	TTree *pfTrackTreeOutput = new TTree("pfTree","");
	// pf tracks only in Monte Carlo
	std::vector<float> *pfTrackPtVector = new std::vector<float>(); pfTrackPtVector->clear();
	std::vector<float> *pfTrackPhiVector = new std::vector<float>(); pfTrackPhiVector->clear();
	std::vector<float> *pfTrackEtaVector = new std::vector<float>(); pfTrackEtaVector->clear();
	std::vector<float> *pfTrackEnergyVector = new std::vector<float>(); pfTrackEnergyVector->clear();
	std::vector<float> *pfTrackMassVector = new std::vector<float>(); pfTrackMassVector->clear();
	std::vector<int> *pfTrackIDVector = new std::vector<int>(); pfTrackIDVector->clear();

	pfTrackTreeOutput->Branch("pfPt","vector<float>", &pfTrackPtVector);
	pfTrackTreeOutput->Branch("pfPhi","vector<float>", &pfTrackPhiVector);
	pfTrackTreeOutput->Branch("pfEta","vector<float>", &pfTrackEtaVector);
	pfTrackTreeOutput->Branch("pfEnergy","vector<float>", &pfTrackEnergyVector);
	pfTrackTreeOutput->Branch("pfM","vector<float>", &pfTrackMassVector);
	pfTrackTreeOutput->Branch("pfId","vector<int>", &pfTrackIDVector);

	TTree *RhoTreeOutput = new TTree("rhotree","");
	RhoTreeOutput->Branch("etaMin","vector<double>", &RhoetaminVector);
	RhoTreeOutput->Branch("etaMax","vector<double>", &RhoetamaxVector);
	RhoTreeOutput->Branch("rho","vector<double>", &RhoVector);
	RhoTreeOutput->Branch("rhom","vector<double>", &RhomVector);

	// ========================================== //
	//			Starting matching events	      //
	// ========================================== //

	Int_t jet_events = heavyIonTree->GetEntries(); // number of events
    // loop through jets and create a key for each event
    for(int ii_entry = 0; ii_entry < jet_events; ii_entry++){
       heavyIonTree->GetEntry(ii_entry);
       unsigned long long key = keyFromRunLumiEvent(run, lumi, event);
       runLumiEvtToEntryMap[key] = ii_entry;
    }

	// ========================================== //
	//				Loop over all events 		  //
	// ========================================== //
	
	bool passTrackCuts;
	bool passJetCuts;
	int iTrackOutput;
	int iJetOutput;

	int nEvents = MainZDCTree->GetEntries();
	cout << "There are " << nEvents << " events" << endl;

	for(int iEvent = 0; iEvent < nEvents; iEvent++) {
		
		if( iEvent % 1000 == 0 )	std::cout << "iEvent: " << iEvent <<	" of " << nEvents << std::endl;


		// ========================================== //
		//			Start with the ZDCs	              //
		// ========================================== //
		MainZDCTree->GetEntry(iEvent);

		//Find matching ZDC event
		if (ZDC_evt < 0) continue;
		unsigned long long key = keyFromRunLumiEvent((UInt_t)ZDC_run,(UInt_t)ZDC_lumi,(ULong64_t)ZDC_evt);
       	long long i_entry = -1;
       	if(runLumiEvtToEntryMap.count(key) == 0){ continue; // skip reco event if there is no event match
       	}else{ i_entry = runLumiEvtToEntryMap.at(key); }

		// ========================================== //
		//	Read the event to input trees	      //
		// ========================================== //
		
		heavyIonTree->GetEntry(i_entry);
		hltTree->GetEntry(i_entry);
		skimTree->GetEntry(i_entry);
		trackTree->GetEntry(i_entry);
		RhoTree->GetEntry(i_entry);
		if(is_MC) genTrackTree->GetEntry(i_entry);
		particleFlowCandidateTree->GetEntry(i_entry);
		checkFlatteningTree->GetEntry(i_entry);

		Ntroff = 0;
		MB_all = 0;

		if ( MB_FirstCollisionAfterAbortGap != 0 || MB_ForSkim != 0 || MB_ForExpress != 0 || MB_part1 != 0 ||  MB_part2 != 0 ||  MB_part3 != 0 ||  MB_part4 != 0 ||  MB_part5 != 0 ||  MB_part6 != 0 ||  MB_part7 != 0 ||  MB_part8 != 0 ||  MB_part9 != 0 ||  MB_part10 != 0 ||  MB_part11 != 0 ||  MB_part12 != 0 ||  MB_part13 != 0 ||  MB_part14 != 0 ||  MB_part15 != 0 ||  MB_part16 != 0 ||  MB_part17 != 0 ||  MB_part18 != 0 ||  MB_part19 != 0 ||  MB_part20 != 0 ) MB_all = 1;

		for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
			jetTree[iJetType]->GetEntry(i_entry);
		}

		int multiplicity = get_Ntrkoff(nTracks, trackEtaArray, trackPtArray, trackChargeArray, trackHighPurityArray, trackPtErrorArray, trackVertexDistanceXYArray, trackVertexDistanceXYErrorArray, trackVertexDistanceZArray, trackVertexDistanceZErrorArray);
		bool multsel = true;
		if(multiplicity > 40) multsel = false;
		if(multsel == false) continue;		
		Ntroff = multiplicity;

		//moving zdc to calibrated one
		hiZDCplus = ZDC_SumP;
		hiZDCminus = ZDC_SumN;	

		double track_gap = 10000.0;

		hltTreeOutput->Fill();
		skimTreeOutput->Fill();
		RhoTreeOutput->Fill();
		
		//Event plane (just what we want EP from 2 to 4)
		epang_HFm2 = (float) eventPlaneAngle[44];
		epang_HFp2 = (float) eventPlaneAngle[45];
		epang_HFm3 = (float) eventPlaneAngle[73];
		epang_HFp3 = (float) eventPlaneAngle[74];
		epang_HFm4 = (float) eventPlaneAngle[102];
		epang_HFp4 = (float) eventPlaneAngle[103];

		q_HFm2 = (float) eventPlaneQ[44];
		q_HFp2 = (float) eventPlaneQ[45];
		q_HFm3 = (float) eventPlaneQ[73];
		q_HFp3 = (float) eventPlaneQ[74];
		q_HFm4 = (float) eventPlaneQ[102];
		q_HFp4 = (float) eventPlaneQ[103];

		qx_HFm2 = (float) eventPlaneQx[44];
		qx_HFp2 = (float) eventPlaneQx[45];
		qx_HFm3 = (float) eventPlaneQx[73];
		qx_HFp3 = (float) eventPlaneQx[74];
		qx_HFm4 = (float) eventPlaneQx[102];
		qx_HFp4 = (float) eventPlaneQx[103];
		
		qy_HFm2 = (float) eventPlaneQy[44];
		qy_HFp2 = (float) eventPlaneQy[45];
		qy_HFm3 = (float) eventPlaneQy[73];
		qy_HFp3 = (float) eventPlaneQy[74];
		qy_HFm4 = (float) eventPlaneQy[102];
		qy_HFp4 = (float) eventPlaneQy[103];

		mult_HFm2 = (float) eventPlaneMultiplicity[44];
		mult_HFp2 = (float) eventPlaneMultiplicity[45];
		mult_HFm3 = (float) eventPlaneMultiplicity[73];
		mult_HFp3 = (float) eventPlaneMultiplicity[74];
		mult_HFm4 = (float) eventPlaneMultiplicity[102];
		mult_HFp4 = (float) eventPlaneMultiplicity[103];

		checkFlatteningTreeOutput->Fill();

    	// Fill jet histograms using basic jet cuts
		for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
		
			iJetOutput = 0;
			nJetsOutput[iJetType] = nJets[iJetType];
		
			for(int iJet = 0; iJet < nJets[iJetType]; iJet++){ // loop over reco jets

				passJetCuts = true;
				
				vector<Node *> NodesWTAScheme;
				NodesWTAScheme.clear();
				
				vector<Node *> NodesCAWTATree;
				NodesCAWTATree.clear();				
				
				vector<Node *> NodesCATree;
				NodesCATree.clear();

				double jetR = 0.4;
				if( iJetType == 3 ) jetR = 0.3;
				Float_t jetPhiWTA = -999;
				Float_t jetEtaWTA = -999;
				Float_t jetmassC = -999;
				
				// zcut 0.1 and beta 0.0
				Float_t SD_subjet1_pt_0p1_0p0 = -999;
				Float_t SD_subjet1_eta_0p1_0p0 = -999;
				Float_t SD_subjet1_phi_0p1_0p0 = -999;
				Float_t SD_subjet1_etaWTA_0p1_0p0 = -999;
				Float_t SD_subjet1_phiWTA_0p1_0p0 = -999;
				Float_t SD_subjet1_mass_0p1_0p0 = -999;
				Float_t SD_subjet2_pt_0p1_0p0 = -999;
				Float_t SD_subjet2_eta_0p1_0p0 = -999;
				Float_t SD_subjet2_phi_0p1_0p0 = -999;
				Float_t SD_subjet2_etaWTA_0p1_0p0 = -999;
				Float_t SD_subjet2_phiWTA_0p1_0p0 = -999;
				Float_t SD_subjet2_mass_0p1_0p0 = -999;

				// zcut 0.1 and beta 1.0
				Float_t SD_subjet1_pt_0p1_1p0 = -999;
				Float_t SD_subjet1_eta_0p1_1p0 = -999;
				Float_t SD_subjet1_phi_0p1_1p0 = -999;
				Float_t SD_subjet1_etaWTA_0p1_1p0 = -999;
				Float_t SD_subjet1_phiWTA_0p1_1p0 = -999;
				Float_t SD_subjet1_mass_0p1_1p0 = -999;
				Float_t SD_subjet2_pt_0p1_1p0 = -999;
				Float_t SD_subjet2_eta_0p1_1p0 = -999;
				Float_t SD_subjet2_phi_0p1_1p0 = -999;
				Float_t SD_subjet2_etaWTA_0p1_1p0 = -999;
				Float_t SD_subjet2_phiWTA_0p1_1p0 = -999;
				Float_t SD_subjet2_mass_0p1_1p0 = -999;

				// zcut 0.1 and beta 2.0
				Float_t SD_subjet1_pt_0p1_2p0 = -999;
				Float_t SD_subjet1_eta_0p1_2p0 = -999;
				Float_t SD_subjet1_phi_0p1_2p0 = -999;
				Float_t SD_subjet1_etaWTA_0p1_2p0 = -999;
				Float_t SD_subjet1_phiWTA_0p1_2p0 = -999;
				Float_t SD_subjet1_mass_0p1_2p0 = -999;
				Float_t SD_subjet2_pt_0p1_2p0 = -999;
				Float_t SD_subjet2_eta_0p1_2p0 = -999;
				Float_t SD_subjet2_phi_0p1_2p0 = -999;
				Float_t SD_subjet2_etaWTA_0p1_2p0 = -999;
				Float_t SD_subjet2_phiWTA_0p1_2p0 = -999;
				Float_t SD_subjet2_mass_0p1_2p0 = -999;

				// zcut 0.2 and beta 0.0
				Float_t SD_subjet1_pt_0p2_0p0 = -999;
				Float_t SD_subjet1_eta_0p2_0p0 = -999;
				Float_t SD_subjet1_phi_0p2_0p0 = -999;
				Float_t SD_subjet1_etaWTA_0p2_0p0 = -999;
				Float_t SD_subjet1_phiWTA_0p2_0p0 = -999;
				Float_t SD_subjet1_mass_0p2_0p0 = -999;
				Float_t SD_subjet2_pt_0p2_0p0 = -999;
				Float_t SD_subjet2_eta_0p2_0p0 = -999;
				Float_t SD_subjet2_phi_0p2_0p0 = -999;
				Float_t SD_subjet2_etaWTA_0p2_0p0 = -999;
				Float_t SD_subjet2_phiWTA_0p2_0p0 = -999;
				Float_t SD_subjet2_mass_0p2_0p0 = -999;
				
				// zcut 0.4 and beta 0.0
				Float_t SD_subjet1_pt_0p4_0p0 = -999;
				Float_t SD_subjet1_eta_0p4_0p0 = -999;
				Float_t SD_subjet1_phi_0p4_0p0 = -999;
				Float_t SD_subjet1_etaWTA_0p4_0p0 = -999;
				Float_t SD_subjet1_phiWTA_0p4_0p0 = -999;
				Float_t SD_subjet1_mass_0p4_0p0 = -999;
				Float_t SD_subjet2_pt_0p4_0p0 = -999;
				Float_t SD_subjet2_eta_0p4_0p0 = -999;
				Float_t SD_subjet2_phi_0p4_0p0 = -999;
				Float_t SD_subjet2_etaWTA_0p4_0p0 = -999;
				Float_t SD_subjet2_phiWTA_0p4_0p0 = -999;
				Float_t SD_subjet2_mass_0p4_0p0 = -999;
				
				// zcut 0.5 and beta 1.0
				Float_t SD_subjet1_pt_0p5_1p0 = -999;
				Float_t SD_subjet1_eta_0p5_1p0 = -999;
				Float_t SD_subjet1_phi_0p5_1p0 = -999;
				Float_t SD_subjet1_etaWTA_0p5_1p0 = -999;
				Float_t SD_subjet1_phiWTA_0p5_1p0 = -999;
				Float_t SD_subjet1_mass_0p5_1p0 = -999;
				Float_t SD_subjet2_pt_0p5_1p0 = -999;
				Float_t SD_subjet2_eta_0p5_1p0 = -999;
				Float_t SD_subjet2_phi_0p5_1p0 = -999;
				Float_t SD_subjet2_etaWTA_0p5_1p0 = -999;
				Float_t SD_subjet2_phiWTA_0p5_1p0 = -999;
				Float_t SD_subjet2_mass_0p5_1p0 = -999;

				// zcut 0.5 and beta 1.5
				Float_t SD_subjet1_pt_0p5_1p5 = -999;
				Float_t SD_subjet1_eta_0p5_1p5 = -999;
				Float_t SD_subjet1_phi_0p5_1p5 = -999;
				Float_t SD_subjet1_etaWTA_0p5_1p5 = -999;
				Float_t SD_subjet1_phiWTA_0p5_1p5 = -999;
				Float_t SD_subjet1_mass_0p5_1p5 = -999;
				Float_t SD_subjet2_pt_0p5_1p5 = -999;
				Float_t SD_subjet2_eta_0p5_1p5 = -999;
				Float_t SD_subjet2_phi_0p5_1p5 = -999;
				Float_t SD_subjet2_etaWTA_0p5_1p5 = -999;
				Float_t SD_subjet2_phiWTA_0p5_1p5 = -999;
				Float_t SD_subjet2_mass_0p5_1p5 = -999;
				
				// recluster and find WTA axis
				if( iJetType == 0){
					for(int itrk = 0; itrk < nTracks; itrk++) {
     					// Set particle kinematics
     					double deltaphi = jetPhiArray[iJetType][iJet] - trackPhiArray[itrk];
     					if(deltaphi > (TMath::Pi())){deltaphi += -2*TMath::Pi();}
     					if(deltaphi < (-1.0*TMath::Pi())){deltaphi += 2*TMath::Pi();}
     					double deltaeta = jetEtaArray[iJetType][iJet] - trackEtaArray[itrk];
     					double deltaR = sqrt(pow(deltaphi,2) + pow(deltaeta,2));
     					if(deltaR >= jetR) continue;
      					FourVector P;
      					P.SetPtEtaPhiMass(trackPtArray[itrk], trackEtaArray[itrk], trackPhiArray[itrk], pimass);
      					// Add into the node object vector
      					NodesWTAScheme.push_back(new Node(P));
      					NodesCAWTATree.push_back(new Node(P));
      					NodesCATree.push_back(new Node(P));
					}					
				}else{
					for(int pfi = 0; pfi < particleFlowCandidatePtVector->size(); pfi++) {
     					// Set particle kinematics
     					double deltaphi = jetPhiArray[iJetType][iJet] - particleFlowCandidatePhiVector->at(pfi);
     					if(deltaphi > (TMath::Pi())){deltaphi += -2*TMath::Pi();}
     					if(deltaphi < (-1.0*TMath::Pi())){deltaphi += 2*TMath::Pi();}
     					double deltaeta = jetEtaArray[iJetType][iJet] - particleFlowCandidateEtaVector->at(pfi);
     					double deltaR = sqrt(pow(deltaphi,2) + pow(deltaeta,2));
     					if(deltaR >= jetR) continue;
      					FourVector P;
      					P.SetPtEtaPhiMass(particleFlowCandidatePtVector->at(pfi), particleFlowCandidateEtaVector->at(pfi), particleFlowCandidatePhiVector->at(pfi), particleFlowCandidateMassVector->at(pfi));
      					// Add into the node object vector
      					NodesWTAScheme.push_back(new Node(P));
      					NodesCAWTATree.push_back(new Node(P));
      					NodesCATree.push_back(new Node(P));
   					}
				}	

				// Do the reclustering!
   				if(NodesWTAScheme.size()>0){
   					BuildCATree(NodesWTAScheme, -1, WTAScheme);
  					jetPhiWTA = NodesWTAScheme[0]->P.GetPhi();
  					jetEtaWTA = NodesWTAScheme[0]->P.GetEta();
  					jetmassC = NodesWTAScheme[0]->P.GetMass();
	  				delete NodesWTAScheme[0];
	  				NodesWTAScheme.clear();
				}

   				if(NodesCATree.size()>0 && storesoftdrop && (iJetType == 1 || iJetType == 2)){
					std::pair<float, float> zcut_beta(0.1,0.0);
  					BuildCATree(NodesCATree, 0, EScheme);	
					Node *SDNode = FindSDNode(NodesCATree[0], zcut_beta.first, zcut_beta.second, jetR);
  					if(SDNode->N > 1){
  						FourVector Subjet1 = SDNode->Child1->P;
						FourVector Subjet2 = SDNode->Child2->P;
      					SD_subjet1_pt_0p1_0p0 =  Subjet1.GetPT();
      					SD_subjet1_eta_0p1_0p0 =  Subjet1.GetEta();
      					SD_subjet1_phi_0p1_0p0 =  Subjet1.GetPhi();
      					SD_subjet1_mass_0p1_0p0 =  Subjet1.GetMass();
      					SD_subjet2_pt_0p1_0p0 =  Subjet2.GetPT();
      					SD_subjet2_eta_0p1_0p0 =  Subjet2.GetEta();
      					SD_subjet2_phi_0p1_0p0 =  Subjet2.GetPhi();
      					SD_subjet2_mass_0p1_0p0 =  Subjet2.GetMass();
  					}
  					zcut_beta = std::make_pair(0.1,1.0);
					SDNode = FindSDNode(NodesCATree[0], zcut_beta.first, zcut_beta.second, jetR);
  					if(SDNode->N > 1){
  						FourVector Subjet1 = SDNode->Child1->P;
						FourVector Subjet2 = SDNode->Child2->P;
      					SD_subjet1_pt_0p1_1p0 =  Subjet1.GetPT();
      					SD_subjet1_eta_0p1_1p0 =  Subjet1.GetEta();
      					SD_subjet1_phi_0p1_1p0 =  Subjet1.GetPhi();
      					SD_subjet1_mass_0p1_1p0 =  Subjet1.GetMass();
      					SD_subjet2_pt_0p1_1p0 =  Subjet2.GetPT();
      					SD_subjet2_eta_0p1_1p0 =  Subjet2.GetEta();
      					SD_subjet2_phi_0p1_1p0 =  Subjet2.GetPhi();
      					SD_subjet2_mass_0p1_1p0 =  Subjet2.GetMass();
  					}
  					zcut_beta = std::make_pair(0.1,2.0);
					SDNode = FindSDNode(NodesCATree[0], zcut_beta.first, zcut_beta.second, jetR);
  					if(SDNode->N > 1){
  						FourVector Subjet1 = SDNode->Child1->P;
						FourVector Subjet2 = SDNode->Child2->P;
      					SD_subjet1_pt_0p1_2p0 =  Subjet1.GetPT();
      					SD_subjet1_eta_0p1_2p0 =  Subjet1.GetEta();
      					SD_subjet1_phi_0p1_2p0 =  Subjet1.GetPhi();
      					SD_subjet1_mass_0p1_2p0 =  Subjet1.GetMass();
      					SD_subjet2_pt_0p1_2p0 =  Subjet2.GetPT();
      					SD_subjet2_eta_0p1_2p0 =  Subjet2.GetEta();
      					SD_subjet2_phi_0p1_2p0 =  Subjet2.GetPhi();
      					SD_subjet2_mass_0p1_2p0 =  Subjet2.GetMass();
  					}
  					zcut_beta = std::make_pair(0.2,0.0);
					SDNode = FindSDNode(NodesCATree[0], zcut_beta.first, zcut_beta.second, jetR);
  					if(SDNode->N > 1){
  						FourVector Subjet1 = SDNode->Child1->P;
						FourVector Subjet2 = SDNode->Child2->P;
      					SD_subjet1_pt_0p2_0p0 =  Subjet1.GetPT();
      					SD_subjet1_eta_0p2_0p0 =  Subjet1.GetEta();
      					SD_subjet1_phi_0p2_0p0 =  Subjet1.GetPhi();
      					SD_subjet1_mass_0p2_0p0 =  Subjet1.GetMass();
      					SD_subjet2_pt_0p2_0p0 =  Subjet2.GetPT();
      					SD_subjet2_eta_0p2_0p0 =  Subjet2.GetEta();
      					SD_subjet2_phi_0p2_0p0 =  Subjet2.GetPhi();
      					SD_subjet2_mass_0p2_0p0 =  Subjet2.GetMass();
  					}
  					zcut_beta = std::make_pair(0.4,0.0);
					SDNode = FindSDNode(NodesCATree[0], zcut_beta.first, zcut_beta.second, jetR);
  					if(SDNode->N > 1){
  						FourVector Subjet1 = SDNode->Child1->P;
						FourVector Subjet2 = SDNode->Child2->P;
      					SD_subjet1_pt_0p4_0p0 =  Subjet1.GetPT();
      					SD_subjet1_eta_0p4_0p0 =  Subjet1.GetEta();
      					SD_subjet1_phi_0p4_0p0 =  Subjet1.GetPhi();
      					SD_subjet1_mass_0p4_0p0 =  Subjet1.GetMass();
      					SD_subjet2_pt_0p4_0p0 =  Subjet2.GetPT();
      					SD_subjet2_eta_0p4_0p0 =  Subjet2.GetEta();
      					SD_subjet2_phi_0p4_0p0 =  Subjet2.GetPhi();
      					SD_subjet2_mass_0p4_0p0 =  Subjet2.GetMass();
  					}
  					zcut_beta = std::make_pair(0.5,1.0);
					SDNode = FindSDNode(NodesCATree[0], zcut_beta.first, zcut_beta.second, jetR);
  					if(SDNode->N > 1){
  						FourVector Subjet1 = SDNode->Child1->P;
						FourVector Subjet2 = SDNode->Child2->P;
      					SD_subjet1_pt_0p5_1p0 =  Subjet1.GetPT();
      					SD_subjet1_eta_0p5_1p0 =  Subjet1.GetEta();
      					SD_subjet1_phi_0p5_1p0 =  Subjet1.GetPhi();
      					SD_subjet1_mass_0p5_1p0 =  Subjet1.GetMass();
      					SD_subjet2_pt_0p5_1p0 =  Subjet2.GetPT();
      					SD_subjet2_eta_0p5_1p0 =  Subjet2.GetEta();
      					SD_subjet2_phi_0p5_1p0 =  Subjet2.GetPhi();
      					SD_subjet2_mass_0p5_1p0 =  Subjet2.GetMass();
  					}
  					zcut_beta = std::make_pair(0.5,1.5);
					SDNode = FindSDNode(NodesCATree[0], zcut_beta.first, zcut_beta.second, jetR);
  					if(SDNode->N > 1){
  						FourVector Subjet1 = SDNode->Child1->P;
						FourVector Subjet2 = SDNode->Child2->P;
      					SD_subjet1_pt_0p5_1p5 =  Subjet1.GetPT();
      					SD_subjet1_eta_0p5_1p5 =  Subjet1.GetEta();
      					SD_subjet1_phi_0p5_1p5 =  Subjet1.GetPhi();
      					SD_subjet1_mass_0p5_1p5 =  Subjet1.GetMass();
      					SD_subjet2_pt_0p5_1p5 =  Subjet2.GetPT();
      					SD_subjet2_eta_0p5_1p5 =  Subjet2.GetEta();
      					SD_subjet2_phi_0p5_1p5 =  Subjet2.GetPhi();
      					SD_subjet2_mass_0p5_1p5 =  Subjet2.GetMass();
  					}
  					delete NodesCATree[0];
  					NodesCATree.clear();
  				}
  				
  				if(NodesCAWTATree.size()>0 && storesoftdrop && (iJetType == 1 || iJetType == 2)){
					std::pair<float, float> zcut_beta(0.1,0.0);
  				  	BuildCATree(NodesCAWTATree, 0, WTAScheme);
					Node *SDNodeWTA = FindSDNode(NodesCAWTATree[0], zcut_beta.first, zcut_beta.second, jetR);
  					if(SDNodeWTA->N > 1){
  						FourVector Subjet1 = SDNodeWTA->Child1->P;
						FourVector Subjet2 = SDNodeWTA->Child2->P;
						SD_subjet1_etaWTA_0p1_0p0 =  Subjet1.GetEta();
      					SD_subjet1_phiWTA_0p1_0p0 =  Subjet1.GetPhi();
						SD_subjet2_etaWTA_0p1_0p0 =  Subjet2.GetEta();
      					SD_subjet2_phiWTA_0p1_0p0 =  Subjet2.GetPhi();						
  					}
					zcut_beta = std::make_pair(0.1,1.0);
					SDNodeWTA = FindSDNode(NodesCAWTATree[0], zcut_beta.first, zcut_beta.second, jetR);
  					if(SDNodeWTA->N > 1){
  						FourVector Subjet1 = SDNodeWTA->Child1->P;
						FourVector Subjet2 = SDNodeWTA->Child2->P;
						SD_subjet1_etaWTA_0p1_1p0 =  Subjet1.GetEta();
      					SD_subjet1_phiWTA_0p1_1p0 =  Subjet1.GetPhi();
						SD_subjet2_etaWTA_0p1_1p0 =  Subjet2.GetEta();
      					SD_subjet2_phiWTA_0p1_1p0 =  Subjet2.GetPhi();						
  					}
					zcut_beta = std::make_pair(0.1,2.0);
					SDNodeWTA = FindSDNode(NodesCAWTATree[0], zcut_beta.first, zcut_beta.second, jetR);
  					if(SDNodeWTA->N > 1){
  						FourVector Subjet1 = SDNodeWTA->Child1->P;
						FourVector Subjet2 = SDNodeWTA->Child2->P;
						SD_subjet1_etaWTA_0p1_2p0 =  Subjet1.GetEta();
      					SD_subjet1_phiWTA_0p1_2p0 =  Subjet1.GetPhi();
						SD_subjet2_etaWTA_0p1_2p0 =  Subjet2.GetEta();
      					SD_subjet2_phiWTA_0p1_2p0 =  Subjet2.GetPhi();						
  					}
					zcut_beta = std::make_pair(0.2,0.0);
					SDNodeWTA = FindSDNode(NodesCAWTATree[0], zcut_beta.first, zcut_beta.second, jetR);
  					if(SDNodeWTA->N > 1){
  						FourVector Subjet1 = SDNodeWTA->Child1->P;
						FourVector Subjet2 = SDNodeWTA->Child2->P;
						SD_subjet1_etaWTA_0p2_0p0 =  Subjet1.GetEta();
      					SD_subjet1_phiWTA_0p2_0p0 =  Subjet1.GetPhi();
						SD_subjet2_etaWTA_0p2_0p0 =  Subjet2.GetEta();
      					SD_subjet2_phiWTA_0p2_0p0 =  Subjet2.GetPhi();						
  					}
					zcut_beta = std::make_pair(0.4,0.0);
					SDNodeWTA = FindSDNode(NodesCAWTATree[0], zcut_beta.first, zcut_beta.second, jetR);
  					if(SDNodeWTA->N > 1){
  						FourVector Subjet1 = SDNodeWTA->Child1->P;
						FourVector Subjet2 = SDNodeWTA->Child2->P;
						SD_subjet1_etaWTA_0p4_0p0 =  Subjet1.GetEta();
      					SD_subjet1_phiWTA_0p4_0p0 =  Subjet1.GetPhi();
						SD_subjet2_etaWTA_0p4_0p0 =  Subjet2.GetEta();
      					SD_subjet2_phiWTA_0p4_0p0 =  Subjet2.GetPhi();						
  					}
					zcut_beta = std::make_pair(0.5,1.0);
					SDNodeWTA = FindSDNode(NodesCAWTATree[0], zcut_beta.first, zcut_beta.second, jetR);
  					if(SDNodeWTA->N > 1){
  						FourVector Subjet1 = SDNodeWTA->Child1->P;
						FourVector Subjet2 = SDNodeWTA->Child2->P;
						SD_subjet1_etaWTA_0p5_1p0 =  Subjet1.GetEta();
      					SD_subjet1_phiWTA_0p5_1p0 =  Subjet1.GetPhi();
						SD_subjet2_etaWTA_0p5_1p0 =  Subjet2.GetEta();
      					SD_subjet2_phiWTA_0p5_1p0 =  Subjet2.GetPhi();						
  					}
					zcut_beta = std::make_pair(0.5,1.5);
					SDNodeWTA = FindSDNode(NodesCAWTATree[0], zcut_beta.first, zcut_beta.second, jetR);
  					if(SDNodeWTA->N > 1){
  						FourVector Subjet1 = SDNodeWTA->Child1->P;
						FourVector Subjet2 = SDNodeWTA->Child2->P;
						SD_subjet1_etaWTA_0p5_1p5 =  Subjet1.GetEta();
      					SD_subjet1_phiWTA_0p5_1p5 =  Subjet1.GetPhi();
						SD_subjet2_etaWTA_0p5_1p5 =  Subjet2.GetEta();
      					SD_subjet2_phiWTA_0p5_1p5 =  Subjet2.GetPhi();						
  					}
  					delete NodesCAWTATree[0];
  					NodesCAWTATree.clear();
  				}
  				  				
				// Apply very basic jet cuts
				if(jetRawPtArray[iJetType][iJet] < jetrawptmin) passJetCuts = false;    // Minumum pT cut of 30 GeV
			
				// Fill the jet arrays with reconstructed jets
				if(passJetCuts){
					jetRawPtArrayOutput[iJetType][iJetOutput] = jetRawPtArray[iJetType][iJet];
					jetPhiArrayOutput[iJetType][iJetOutput] = jetPhiArray[iJetType][iJet];
					jetPhiArrayWTAOutput[iJetType][iJetOutput] = jetPhiWTA;
					jetEtaArrayOutput[iJetType][iJetOutput] = jetEtaArray[iJetType][iJet];
					jetEtaArrayWTAOutput[iJetType][iJetOutput] = jetEtaWTA;
					jetMaxTrackPtArrayOutput[iJetType][iJetOutput] = jetMaxTrackPtArray[iJetType][iJet];
					jetMassArrayOutput[iJetType][iJetOutput] = jetMassArray[iJetType][iJet];
					jetMassCalcArrayOutput[iJetType][iJetOutput] = jetmassC;
					
					jetPfNHFArrayOutput[iJetType][iJetOutput] = jetPfNHFArray[iJetType][iJet];
					jetPfNEFArrayOutput[iJetType][iJetOutput] = jetPfNEFArray[iJetType][iJet];
					jetPfCHFArrayOutput[iJetType][iJetOutput] = jetPfCHFArray[iJetType][iJet];
					jetPfMUFArrayOutput[iJetType][iJetOutput] = jetPfMUFArray[iJetType][iJet];
					jetPfCEFArrayOutput[iJetType][iJetOutput] = jetPfCEFArray[iJetType][iJet];
					jetPfCHMArrayOutput[iJetType][iJetOutput] = jetPfCHMArray[iJetType][iJet];
					jetPfCEMArrayOutput[iJetType][iJetOutput] = jetPfCEMArray[iJetType][iJet];
					jetPfNHMArrayOutput[iJetType][iJetOutput] = jetPfNHMArray[iJetType][iJet];	
					jetPfNEMArrayOutput[iJetType][iJetOutput] = jetPfNEMArray[iJetType][iJet];
					jetPfMUMArrayOutput[iJetType][iJetOutput] = jetPfMUMArray[iJetType][iJet];				

					jetHCALSUMArrayOutput[iJetType][iJetOutput]	= jetHCALSUMArray[iJetType][iJet];
					jetECALSUMArrayOutput[iJetType][iJetOutput]	= jetECALSUMArray[iJetType][iJet];			
					jetTRKSUMArrayOutput[iJetType][iJetOutput]	= jetTRKSUMArray[iJetType][iJet];
					jetTRKNArrayOutput[iJetType][iJetOutput]	= jetTRKNArray[iJetType][iJet];
					jetMUNArrayOutput[iJetType][iJetOutput]		= jetMUNArray[iJetType][iJet];
					
					Subjet1_z0p1_b0p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjet1_pt_0p1_0p0;
					Subjet1_z0p1_b0p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjet1_phi_0p1_0p0;
					Subjet1_z0p1_b0p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjet1_eta_0p1_0p0;
					Subjet1_z0p1_b0p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjet1_phiWTA_0p1_0p0;
					Subjet1_z0p1_b0p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjet1_etaWTA_0p1_0p0;
					Subjet1_z0p1_b0p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjet1_mass_0p1_0p0;
					Subjet2_z0p1_b0p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjet2_pt_0p1_0p0;
					Subjet2_z0p1_b0p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjet2_phi_0p1_0p0;
					Subjet2_z0p1_b0p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjet2_eta_0p1_0p0;
					Subjet2_z0p1_b0p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjet2_phiWTA_0p1_0p0;
					Subjet2_z0p1_b0p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjet2_etaWTA_0p1_0p0;
					Subjet2_z0p1_b0p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjet2_mass_0p1_0p0;

					Subjet1_z0p1_b1p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjet1_pt_0p1_1p0;
					Subjet1_z0p1_b1p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjet1_phi_0p1_1p0;
					Subjet1_z0p1_b1p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjet1_eta_0p1_1p0;
					Subjet1_z0p1_b1p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjet1_phiWTA_0p1_1p0;
					Subjet1_z0p1_b1p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjet1_etaWTA_0p1_1p0;
					Subjet1_z0p1_b1p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjet1_mass_0p1_1p0;
					Subjet2_z0p1_b1p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjet2_pt_0p1_1p0;
					Subjet2_z0p1_b1p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjet2_phi_0p1_1p0;
					Subjet2_z0p1_b1p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjet2_eta_0p1_1p0;
					Subjet2_z0p1_b1p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjet2_phiWTA_0p1_1p0;
					Subjet2_z0p1_b1p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjet2_etaWTA_0p1_1p0;
					Subjet2_z0p1_b1p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjet2_mass_0p1_1p0;

					Subjet1_z0p1_b2p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjet1_pt_0p1_2p0;
					Subjet1_z0p1_b2p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjet1_phi_0p1_2p0;
					Subjet1_z0p1_b2p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjet1_eta_0p1_2p0;
					Subjet1_z0p1_b2p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjet1_phiWTA_0p1_2p0;
					Subjet1_z0p1_b2p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjet1_etaWTA_0p1_2p0;
					Subjet1_z0p1_b2p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjet1_mass_0p1_2p0;
					Subjet2_z0p1_b2p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjet2_pt_0p1_2p0;
					Subjet2_z0p1_b2p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjet2_phi_0p1_2p0;
					Subjet2_z0p1_b2p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjet2_eta_0p1_2p0;
					Subjet2_z0p1_b2p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjet2_phiWTA_0p1_2p0;
					Subjet2_z0p1_b2p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjet2_etaWTA_0p1_2p0;
					Subjet2_z0p1_b2p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjet2_mass_0p1_2p0;

					Subjet1_z0p2_b0p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjet1_pt_0p2_0p0;
					Subjet1_z0p2_b0p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjet1_phi_0p2_0p0;
					Subjet1_z0p2_b0p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjet1_eta_0p2_0p0;
					Subjet1_z0p2_b0p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjet1_phiWTA_0p2_0p0;
					Subjet1_z0p2_b0p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjet1_etaWTA_0p2_0p0;
					Subjet1_z0p2_b0p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjet1_mass_0p2_0p0;
					Subjet2_z0p2_b0p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjet2_pt_0p2_0p0;
					Subjet2_z0p2_b0p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjet2_phi_0p2_0p0;
					Subjet2_z0p2_b0p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjet2_eta_0p2_0p0;
					Subjet2_z0p2_b0p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjet2_phiWTA_0p2_0p0;
					Subjet2_z0p2_b0p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjet2_etaWTA_0p2_0p0;
					Subjet2_z0p2_b0p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjet2_mass_0p2_0p0;

					Subjet1_z0p4_b0p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjet1_pt_0p4_0p0;
					Subjet1_z0p4_b0p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjet1_phi_0p4_0p0;
					Subjet1_z0p4_b0p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjet1_eta_0p4_0p0;
					Subjet1_z0p4_b0p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjet1_phiWTA_0p4_0p0;
					Subjet1_z0p4_b0p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjet1_etaWTA_0p4_0p0;
					Subjet1_z0p4_b0p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjet1_mass_0p4_0p0;
					Subjet2_z0p4_b0p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjet2_pt_0p4_0p0;
					Subjet2_z0p4_b0p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjet2_phi_0p4_0p0;
					Subjet2_z0p4_b0p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjet2_eta_0p4_0p0;
					Subjet2_z0p4_b0p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjet2_phiWTA_0p4_0p0;
					Subjet2_z0p4_b0p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjet2_etaWTA_0p4_0p0;
					Subjet2_z0p4_b0p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjet2_mass_0p4_0p0;

					Subjet1_z0p5_b1p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjet1_pt_0p5_1p0;
					Subjet1_z0p5_b1p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjet1_phi_0p5_1p0;
					Subjet1_z0p5_b1p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjet1_eta_0p5_1p0;
					Subjet1_z0p5_b1p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjet1_phiWTA_0p5_1p0;
					Subjet1_z0p5_b1p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjet1_etaWTA_0p5_1p0;
					Subjet1_z0p5_b1p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjet1_mass_0p5_1p0;
					Subjet2_z0p5_b1p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjet2_pt_0p5_1p0;
					Subjet2_z0p5_b1p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjet2_phi_0p5_1p0;
					Subjet2_z0p5_b1p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjet2_eta_0p5_1p0;
					Subjet2_z0p5_b1p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjet2_phiWTA_0p5_1p0;
					Subjet2_z0p5_b1p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjet2_etaWTA_0p5_1p0;
					Subjet2_z0p5_b1p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjet2_mass_0p5_1p0;

					Subjet1_z0p5_b1p5_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjet1_pt_0p5_1p5;
					Subjet1_z0p5_b1p5_PhiArrayOutput[iJetType][iJetOutput] = SD_subjet1_phi_0p5_1p5;
					Subjet1_z0p5_b1p5_EtaArrayOutput[iJetType][iJetOutput] = SD_subjet1_eta_0p5_1p5;
					Subjet1_z0p5_b1p5_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjet1_phiWTA_0p5_1p5;
					Subjet1_z0p5_b1p5_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjet1_etaWTA_0p5_1p5;
					Subjet1_z0p5_b1p5_MassArrayOutput[iJetType][iJetOutput] = SD_subjet1_mass_0p5_1p5;
					Subjet2_z0p5_b1p5_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjet2_pt_0p5_1p5;
					Subjet2_z0p5_b1p5_PhiArrayOutput[iJetType][iJetOutput] = SD_subjet2_phi_0p5_1p5;
					Subjet2_z0p5_b1p5_EtaArrayOutput[iJetType][iJetOutput] = SD_subjet2_eta_0p5_1p5;
					Subjet2_z0p5_b1p5_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjet2_phiWTA_0p5_1p5;
					Subjet2_z0p5_b1p5_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjet2_etaWTA_0p5_1p5;
					Subjet2_z0p5_b1p5_MassArrayOutput[iJetType][iJetOutput] = SD_subjet2_mass_0p5_1p5;

					if(is_MC){
						jetRefPtArrayOutput[iJetType][iJetOutput] = jetRefPtArray[iJetType][iJet];
						jetRefEtaArrayOutput[iJetType][iJetOutput] = jetRefEtaArray[iJetType][iJet];
						jetRefPhiArrayOutput[iJetType][iJetOutput] = jetRefPhiArray[iJetType][iJet];
						jetRefFlavorArrayOutput[iJetType][iJetOutput] = jetRefFlavorArray[iJetType][iJet];
						jetRefFlavorForBArrayOutput[iJetType][iJetOutput] = jetRefFlavorForBArray[iJetType][iJet];
						jetRefSubidArrayOutput[iJetType][iJetOutput] = jetRefSubidArray[iJetType][iJet];
						jetRefMassArrayOutput[iJetType][iJetOutput] = jetRefMassArray[iJetType][iJet];						
					}
					
					iJetOutput++;
				} else {nJetsOutput[iJetType]--;}
			} // Reconstructed jet loop
		
			if(is_MC){
			
				iJetOutput = 0;
				nGenJetsOutput[iJetType] = nGenJets[iJetType];
			
				for(int iJet = 0; iJet < nGenJets[iJetType]; iJet++){
				
					passJetCuts = true;

					vector<Node *> NodesWTASchemeGen;
					NodesWTASchemeGen.clear();
					
					vector<Node *> NodesCAWTATreeGen;
					NodesCAWTATreeGen.clear();				
				
					vector<Node *> NodesCATreeGen;
					NodesCATreeGen.clear();

					double jetRGen = 0.4;
					if( iJetType == 3 ) jetRGen = 0.3;
					Float_t jetPhiWTAGen = -999;
					Float_t jetEtaWTAGen = -999;
					Float_t jetmassCGen = -999;
					
					// zcut 0.1 and beta 0.0
					Float_t SD_subjetgen1_pt_0p1_0p0 = -999;
					Float_t SD_subjetgen1_eta_0p1_0p0 = -999;
					Float_t SD_subjetgen1_phi_0p1_0p0 = -999;
					Float_t SD_subjetgen1_etaWTA_0p1_0p0 = -999;
					Float_t SD_subjetgen1_phiWTA_0p1_0p0 = -999;
					Float_t SD_subjetgen1_mass_0p1_0p0 = -999;
					Float_t SD_subjetgen2_pt_0p1_0p0 = -999;
					Float_t SD_subjetgen2_eta_0p1_0p0 = -999;
					Float_t SD_subjetgen2_phi_0p1_0p0 = -999;
					Float_t SD_subjetgen2_etaWTA_0p1_0p0 = -999;
					Float_t SD_subjetgen2_phiWTA_0p1_0p0 = -999;
					Float_t SD_subjetgen2_mass_0p1_0p0 = -999;

					// zcut 0.1 and beta 1.0
					Float_t SD_subjetgen1_pt_0p1_1p0 = -999;
					Float_t SD_subjetgen1_eta_0p1_1p0 = -999;
					Float_t SD_subjetgen1_phi_0p1_1p0 = -999;
					Float_t SD_subjetgen1_etaWTA_0p1_1p0 = -999;
					Float_t SD_subjetgen1_phiWTA_0p1_1p0 = -999;
					Float_t SD_subjetgen1_mass_0p1_1p0 = -999;
					Float_t SD_subjetgen2_pt_0p1_1p0 = -999;
					Float_t SD_subjetgen2_eta_0p1_1p0 = -999;
					Float_t SD_subjetgen2_phi_0p1_1p0 = -999;
					Float_t SD_subjetgen2_etaWTA_0p1_1p0 = -999;
					Float_t SD_subjetgen2_phiWTA_0p1_1p0 = -999;
					Float_t SD_subjetgen2_mass_0p1_1p0 = -999;

					// zcut 0.1 and beta 2.0
					Float_t SD_subjetgen1_pt_0p1_2p0 = -999;
					Float_t SD_subjetgen1_eta_0p1_2p0 = -999;
					Float_t SD_subjetgen1_phi_0p1_2p0 = -999;
					Float_t SD_subjetgen1_etaWTA_0p1_2p0 = -999;
					Float_t SD_subjetgen1_phiWTA_0p1_2p0 = -999;
					Float_t SD_subjetgen1_mass_0p1_2p0 = -999;
					Float_t SD_subjetgen2_pt_0p1_2p0 = -999;
					Float_t SD_subjetgen2_eta_0p1_2p0 = -999;
					Float_t SD_subjetgen2_phi_0p1_2p0 = -999;
					Float_t SD_subjetgen2_etaWTA_0p1_2p0 = -999;
					Float_t SD_subjetgen2_phiWTA_0p1_2p0 = -999;
					Float_t SD_subjetgen2_mass_0p1_2p0 = -999;

					// zcut 0.2 and beta 0.0
					Float_t SD_subjetgen1_pt_0p2_0p0 = -999;
					Float_t SD_subjetgen1_eta_0p2_0p0 = -999;
					Float_t SD_subjetgen1_phi_0p2_0p0 = -999;
					Float_t SD_subjetgen1_etaWTA_0p2_0p0 = -999;
					Float_t SD_subjetgen1_phiWTA_0p2_0p0 = -999;
					Float_t SD_subjetgen1_mass_0p2_0p0 = -999;
					Float_t SD_subjetgen2_pt_0p2_0p0 = -999;
					Float_t SD_subjetgen2_eta_0p2_0p0 = -999;
					Float_t SD_subjetgen2_phi_0p2_0p0 = -999;
					Float_t SD_subjetgen2_etaWTA_0p2_0p0 = -999;
					Float_t SD_subjetgen2_phiWTA_0p2_0p0 = -999;
					Float_t SD_subjetgen2_mass_0p2_0p0 = -999;
				
					// zcut 0.4 and beta 0.0
					Float_t SD_subjetgen1_pt_0p4_0p0 = -999;
					Float_t SD_subjetgen1_eta_0p4_0p0 = -999;
					Float_t SD_subjetgen1_phi_0p4_0p0 = -999;
					Float_t SD_subjetgen1_etaWTA_0p4_0p0 = -999;
					Float_t SD_subjetgen1_phiWTA_0p4_0p0 = -999;
					Float_t SD_subjetgen1_mass_0p4_0p0 = -999;
					Float_t SD_subjetgen2_pt_0p4_0p0 = -999;
					Float_t SD_subjetgen2_eta_0p4_0p0 = -999;
					Float_t SD_subjetgen2_phi_0p4_0p0 = -999;
					Float_t SD_subjetgen2_etaWTA_0p4_0p0 = -999;
					Float_t SD_subjetgen2_phiWTA_0p4_0p0 = -999;
					Float_t SD_subjetgen2_mass_0p4_0p0 = -999;
				
					// zcut 0.5 and beta 1.0
					Float_t SD_subjetgen1_pt_0p5_1p0 = -999;
					Float_t SD_subjetgen1_eta_0p5_1p0 = -999;
					Float_t SD_subjetgen1_phi_0p5_1p0 = -999;
					Float_t SD_subjetgen1_etaWTA_0p5_1p0 = -999;
					Float_t SD_subjetgen1_phiWTA_0p5_1p0 = -999;
					Float_t SD_subjetgen1_mass_0p5_1p0 = -999;
					Float_t SD_subjetgen2_pt_0p5_1p0 = -999;
					Float_t SD_subjetgen2_eta_0p5_1p0 = -999;
					Float_t SD_subjetgen2_phi_0p5_1p0 = -999;
					Float_t SD_subjetgen2_etaWTA_0p5_1p0 = -999;
					Float_t SD_subjetgen2_phiWTA_0p5_1p0 = -999;
					Float_t SD_subjetgen2_mass_0p5_1p0 = -999;

					// zcut 0.5 and beta 1.5
					Float_t SD_subjetgen1_pt_0p5_1p5 = -999;
					Float_t SD_subjetgen1_eta_0p5_1p5 = -999;
					Float_t SD_subjetgen1_phi_0p5_1p5 = -999;
					Float_t SD_subjetgen1_etaWTA_0p5_1p5 = -999;
					Float_t SD_subjetgen1_phiWTA_0p5_1p5 = -999;
					Float_t SD_subjetgen1_mass_0p5_1p5 = -999;
					Float_t SD_subjetgen2_pt_0p5_1p5 = -999;
					Float_t SD_subjetgen2_eta_0p5_1p5 = -999;
					Float_t SD_subjetgen2_phi_0p5_1p5 = -999;
					Float_t SD_subjetgen2_etaWTA_0p5_1p5 = -999;
					Float_t SD_subjetgen2_phiWTA_0p5_1p5 = -999;
					Float_t SD_subjetgen2_mass_0p5_1p5 = -999;

					for(int gpi = 0; gpi < genTrackPtArray->size(); gpi++) {
     					// Set particle kinematics
     					double deltaphi = genJetPhiArray[iJetType][iJet] - genTrackPhiArray->at(gpi);
     					if(deltaphi > (TMath::Pi())){deltaphi += -2*TMath::Pi();}
     					if(deltaphi < (-1.0*TMath::Pi())){deltaphi += 2*TMath::Pi();}
     					double deltaeta = genJetEtaArray[iJetType][iJet] - genTrackEtaArray->at(gpi);
     					double deltaR = sqrt(pow(deltaphi,2) + pow(deltaeta,2));
     					if(deltaR >= jetRGen) continue;
      					FourVector P;
      					P.SetPtEtaPhiMass(genTrackPtArray->at(gpi), genTrackEtaArray->at(gpi), genTrackPhiArray->at(gpi), genTrackMassArray->at(gpi));
      					// Add into the node object vector
      					NodesWTASchemeGen.push_back(new Node(P));
      					NodesCAWTATreeGen.push_back(new Node(P));
      					NodesCATreeGen.push_back(new Node(P));
   					}

   					// Do the reclustering!
   					if(NodesWTASchemeGen.size()>0){
   						BuildCATree(NodesWTASchemeGen, -1, WTAScheme);
  						jetPhiWTAGen = NodesWTASchemeGen[0]->P.GetPhi();
  						jetEtaWTAGen = NodesWTASchemeGen[0]->P.GetEta();
  						jetmassCGen = NodesWTASchemeGen[0]->P.GetMass();
  						delete NodesWTASchemeGen[0];
  						NodesWTASchemeGen.clear();
					}
					
   					if(NodesCATreeGen.size()>0 && storesoftdrop && (iJetType == 1 || iJetType == 2)){
						std::pair<float, float> zcut_beta(0.1,0.0);
  						BuildCATree(NodesCATreeGen, 0, EScheme);	
						Node *SDNode = FindSDNode(NodesCATreeGen[0], zcut_beta.first, zcut_beta.second, jetRGen);
  						if(SDNode->N > 1){
  							FourVector Subjet1 = SDNode->Child1->P;
							FourVector Subjet2 = SDNode->Child2->P;
      						SD_subjetgen1_pt_0p1_0p0 =  Subjet1.GetPT();
      						SD_subjetgen1_eta_0p1_0p0 =  Subjet1.GetEta();
      						SD_subjetgen1_phi_0p1_0p0 =  Subjet1.GetPhi();
      						SD_subjetgen1_mass_0p1_0p0 =  Subjet1.GetMass();
      						SD_subjetgen2_pt_0p1_0p0 =  Subjet2.GetPT();
      						SD_subjetgen2_eta_0p1_0p0 =  Subjet2.GetEta();
      						SD_subjetgen2_phi_0p1_0p0 =  Subjet2.GetPhi();
      						SD_subjetgen2_mass_0p1_0p0 =  Subjet2.GetMass();
  						}
  						zcut_beta = std::make_pair(0.1,1.0);
						SDNode = FindSDNode(NodesCATreeGen[0], zcut_beta.first, zcut_beta.second, jetRGen);
  						if(SDNode->N > 1){
  							FourVector Subjet1 = SDNode->Child1->P;
							FourVector Subjet2 = SDNode->Child2->P;
      						SD_subjetgen1_pt_0p1_1p0 =  Subjet1.GetPT();
      						SD_subjetgen1_eta_0p1_1p0 =  Subjet1.GetEta();
      						SD_subjetgen1_phi_0p1_1p0 =  Subjet1.GetPhi();
      						SD_subjetgen1_mass_0p1_1p0 =  Subjet1.GetMass();
      						SD_subjetgen2_pt_0p1_1p0 =  Subjet2.GetPT();
      						SD_subjetgen2_eta_0p1_1p0 =  Subjet2.GetEta();
      						SD_subjetgen2_phi_0p1_1p0 =  Subjet2.GetPhi();
      						SD_subjetgen2_mass_0p1_1p0 =  Subjet2.GetMass();
  						}
  						zcut_beta = std::make_pair(0.1,2.0);
						SDNode = FindSDNode(NodesCATreeGen[0], zcut_beta.first, zcut_beta.second, jetRGen);
  						if(SDNode->N > 1){
  							FourVector Subjet1 = SDNode->Child1->P;
							FourVector Subjet2 = SDNode->Child2->P;
      						SD_subjetgen1_pt_0p1_2p0 =  Subjet1.GetPT();
      						SD_subjetgen1_eta_0p1_2p0 =  Subjet1.GetEta();
      						SD_subjetgen1_phi_0p1_2p0 =  Subjet1.GetPhi();
      						SD_subjetgen1_mass_0p1_2p0 =  Subjet1.GetMass();
      						SD_subjetgen2_pt_0p1_2p0 =  Subjet2.GetPT();
      						SD_subjetgen2_eta_0p1_2p0 =  Subjet2.GetEta();
      						SD_subjetgen2_phi_0p1_2p0 =  Subjet2.GetPhi();
      						SD_subjetgen2_mass_0p1_2p0 =  Subjet2.GetMass();
  						}
  						zcut_beta = std::make_pair(0.2,0.0);
						SDNode = FindSDNode(NodesCATreeGen[0], zcut_beta.first, zcut_beta.second, jetRGen);
  						if(SDNode->N > 1){
  							FourVector Subjet1 = SDNode->Child1->P;
							FourVector Subjet2 = SDNode->Child2->P;
      						SD_subjetgen1_pt_0p2_0p0 =  Subjet1.GetPT();
      						SD_subjetgen1_eta_0p2_0p0 =  Subjet1.GetEta();
      						SD_subjetgen1_phi_0p2_0p0 =  Subjet1.GetPhi();
      						SD_subjetgen1_mass_0p2_0p0 =  Subjet1.GetMass();
      						SD_subjetgen2_pt_0p2_0p0 =  Subjet2.GetPT();
      						SD_subjetgen2_eta_0p2_0p0 =  Subjet2.GetEta();
      						SD_subjetgen2_phi_0p2_0p0 =  Subjet2.GetPhi();
      						SD_subjetgen2_mass_0p2_0p0 =  Subjet2.GetMass();
  						}
  						zcut_beta = std::make_pair(0.4,0.0);
						SDNode = FindSDNode(NodesCATreeGen[0], zcut_beta.first, zcut_beta.second, jetRGen);
  						if(SDNode->N > 1){
  							FourVector Subjet1 = SDNode->Child1->P;
							FourVector Subjet2 = SDNode->Child2->P;
      						SD_subjetgen1_pt_0p4_0p0 =  Subjet1.GetPT();
      						SD_subjetgen1_eta_0p4_0p0 =  Subjet1.GetEta();
      						SD_subjetgen1_phi_0p4_0p0 =  Subjet1.GetPhi();
      						SD_subjetgen1_mass_0p4_0p0 =  Subjet1.GetMass();
      						SD_subjetgen2_pt_0p4_0p0 =  Subjet2.GetPT();
      						SD_subjetgen2_eta_0p4_0p0 =  Subjet2.GetEta();
      						SD_subjetgen2_phi_0p4_0p0 =  Subjet2.GetPhi();
      						SD_subjetgen2_mass_0p4_0p0 =  Subjet2.GetMass();
  						}
  						zcut_beta = std::make_pair(0.5,1.0);
						SDNode = FindSDNode(NodesCATreeGen[0], zcut_beta.first, zcut_beta.second, jetRGen);
  						if(SDNode->N > 1){
  							FourVector Subjet1 = SDNode->Child1->P;
							FourVector Subjet2 = SDNode->Child2->P;
      						SD_subjetgen1_pt_0p5_1p0 =  Subjet1.GetPT();
      						SD_subjetgen1_eta_0p5_1p0 =  Subjet1.GetEta();
      						SD_subjetgen1_phi_0p5_1p0 =  Subjet1.GetPhi();
      						SD_subjetgen1_mass_0p5_1p0 =  Subjet1.GetMass();
      						SD_subjetgen2_pt_0p5_1p0 =  Subjet2.GetPT();
      						SD_subjetgen2_eta_0p5_1p0 =  Subjet2.GetEta();
      						SD_subjetgen2_phi_0p5_1p0 =  Subjet2.GetPhi();
      						SD_subjetgen2_mass_0p5_1p0 =  Subjet2.GetMass();
  						}
  						zcut_beta = std::make_pair(0.5,1.5);
						SDNode = FindSDNode(NodesCATreeGen[0], zcut_beta.first, zcut_beta.second, jetRGen);
  						if(SDNode->N > 1){
  							FourVector Subjet1 = SDNode->Child1->P;
							FourVector Subjet2 = SDNode->Child2->P;
      						SD_subjetgen1_pt_0p5_1p5 =  Subjet1.GetPT();
      						SD_subjetgen1_eta_0p5_1p5 =  Subjet1.GetEta();
      						SD_subjetgen1_phi_0p5_1p5 =  Subjet1.GetPhi();
      						SD_subjetgen1_mass_0p5_1p5 =  Subjet1.GetMass();
      						SD_subjetgen2_pt_0p5_1p5 =  Subjet2.GetPT();
      						SD_subjetgen2_eta_0p5_1p5 =  Subjet2.GetEta();
      						SD_subjetgen2_phi_0p5_1p5 =  Subjet2.GetPhi();
      						SD_subjetgen2_mass_0p5_1p5 =  Subjet2.GetMass();
  						}
  						delete NodesCATreeGen[0];
  						NodesCATreeGen.clear();
  					}
  				
  					if(NodesCAWTATreeGen.size()>0 && storesoftdrop && (iJetType == 1 || iJetType == 2)){
						std::pair<float, float> zcut_beta(0.1,0.0);
  				  		BuildCATree(NodesCAWTATreeGen, 0, WTAScheme);
						Node *SDNodeWTA = FindSDNode(NodesCAWTATreeGen[0], zcut_beta.first, zcut_beta.second, jetRGen);
  						if(SDNodeWTA->N > 1){
  							FourVector Subjet1 = SDNodeWTA->Child1->P;
							FourVector Subjet2 = SDNodeWTA->Child2->P;
							SD_subjetgen1_etaWTA_0p1_0p0 =  Subjet1.GetEta();
      						SD_subjetgen1_phiWTA_0p1_0p0 =  Subjet1.GetPhi();
							SD_subjetgen2_etaWTA_0p1_0p0 =  Subjet2.GetEta();
      						SD_subjetgen2_phiWTA_0p1_0p0 =  Subjet2.GetPhi();						
  						}
						zcut_beta = std::make_pair(0.1,1.0);
						SDNodeWTA = FindSDNode(NodesCAWTATreeGen[0], zcut_beta.first, zcut_beta.second, jetRGen);
  						if(SDNodeWTA->N > 1){
  							FourVector Subjet1 = SDNodeWTA->Child1->P;
							FourVector Subjet2 = SDNodeWTA->Child2->P;
							SD_subjetgen1_etaWTA_0p1_1p0 =  Subjet1.GetEta();
   		   					SD_subjetgen1_phiWTA_0p1_1p0 =  Subjet1.GetPhi();
							SD_subjetgen2_etaWTA_0p1_1p0 =  Subjet2.GetEta();
      						SD_subjetgen2_phiWTA_0p1_1p0 =  Subjet2.GetPhi();						
  						}
						zcut_beta = std::make_pair(0.1,2.0);
						SDNodeWTA = FindSDNode(NodesCAWTATreeGen[0], zcut_beta.first, zcut_beta.second, jetRGen);
  						if(SDNodeWTA->N > 1){
  							FourVector Subjet1 = SDNodeWTA->Child1->P;
							FourVector Subjet2 = SDNodeWTA->Child2->P;
							SD_subjetgen1_etaWTA_0p1_2p0 =  Subjet1.GetEta();
      						SD_subjetgen1_phiWTA_0p1_2p0 =  Subjet1.GetPhi();
							SD_subjetgen2_etaWTA_0p1_2p0 =  Subjet2.GetEta();
      						SD_subjetgen2_phiWTA_0p1_2p0 =  Subjet2.GetPhi();						
  						}
						zcut_beta = std::make_pair(0.2,0.0);
						SDNodeWTA = FindSDNode(NodesCAWTATreeGen[0], zcut_beta.first, zcut_beta.second, jetRGen);
  						if(SDNodeWTA->N > 1){
  							FourVector Subjet1 = SDNodeWTA->Child1->P;
							FourVector Subjet2 = SDNodeWTA->Child2->P;
							SD_subjetgen1_etaWTA_0p2_0p0 =  Subjet1.GetEta();
      						SD_subjetgen1_phiWTA_0p2_0p0 =  Subjet1.GetPhi();
							SD_subjetgen2_etaWTA_0p2_0p0 =  Subjet2.GetEta();
      						SD_subjetgen2_phiWTA_0p2_0p0 =  Subjet2.GetPhi();						
  						}
						zcut_beta = std::make_pair(0.4,0.0);
						SDNodeWTA = FindSDNode(NodesCAWTATreeGen[0], zcut_beta.first, zcut_beta.second, jetRGen);
  						if(SDNodeWTA->N > 1){
  							FourVector Subjet1 = SDNodeWTA->Child1->P;
							FourVector Subjet2 = SDNodeWTA->Child2->P;
							SD_subjetgen1_etaWTA_0p4_0p0 =  Subjet1.GetEta();
      						SD_subjetgen1_phiWTA_0p4_0p0 =  Subjet1.GetPhi();
							SD_subjetgen2_etaWTA_0p4_0p0 =  Subjet2.GetEta();
      						SD_subjetgen2_phiWTA_0p4_0p0 =  Subjet2.GetPhi();						
  						}
						zcut_beta = std::make_pair(0.5,1.0);
						SDNodeWTA = FindSDNode(NodesCAWTATreeGen[0], zcut_beta.first, zcut_beta.second, jetRGen);
  						if(SDNodeWTA->N > 1){
  							FourVector Subjet1 = SDNodeWTA->Child1->P;
							FourVector Subjet2 = SDNodeWTA->Child2->P;
							SD_subjetgen1_etaWTA_0p5_1p0 =  Subjet1.GetEta();
      						SD_subjetgen1_phiWTA_0p5_1p0 =  Subjet1.GetPhi();
							SD_subjetgen2_etaWTA_0p5_1p0 =  Subjet2.GetEta();
      						SD_subjetgen2_phiWTA_0p5_1p0 =  Subjet2.GetPhi();						
  						}
						zcut_beta = std::make_pair(0.5,1.5);
						SDNodeWTA = FindSDNode(NodesCAWTATreeGen[0], zcut_beta.first, zcut_beta.second, jetRGen);
  						if(SDNodeWTA->N > 1){
  							FourVector Subjet1 = SDNodeWTA->Child1->P;
							FourVector Subjet2 = SDNodeWTA->Child2->P;
							SD_subjetgen1_etaWTA_0p5_1p5 =  Subjet1.GetEta();
      						SD_subjetgen1_phiWTA_0p5_1p5 =  Subjet1.GetPhi();
							SD_subjetgen2_etaWTA_0p5_1p5 =  Subjet2.GetEta();
      						SD_subjetgen2_phiWTA_0p5_1p5 =  Subjet2.GetPhi();						
  						}
  						delete NodesCAWTATreeGen[0];
  						NodesCAWTATreeGen.clear();
  					}
  	
					// Apply very basic jet cuts
					if(genJetPtArray[iJetType][iJet] < jetrawptmin) passJetCuts = false;    // Minimum pT cut of 30 GeV
				
					// Fill the jet arrays with generated jets
					if(passJetCuts){
				
						genJetPtArrayOutput[iJetType][iJetOutput] = genJetPtArray[iJetType][iJet];
						genJetPhiArrayOutput[iJetType][iJetOutput] = genJetPhiArray[iJetType][iJet];
						genJetPhiArrayWTAOutput[iJetType][iJetOutput] = jetPhiWTAGen;
						genJetEtaArrayOutput[iJetType][iJetOutput] = genJetEtaArray[iJetType][iJet];
						genJetEtaArrayWTAOutput[iJetType][iJetOutput] = jetEtaWTAGen;
						genJetSubidArrayOutput[iJetType][iJetOutput] = genJetSubidArray[iJetType][iJet];
						genJetMatchIndexArrayOutput[iJetType][iJetOutput] = genJetMatchIndexArray[iJetType][iJet];
						genJetMassArrayOutput[iJetType][iJetOutput] = genJetMassArray[iJetType][iJet];
						genJetMassCalcArrayOutput[iJetType][iJetOutput] = jetmassCGen;
					
						SubjetGen1_z0p1_b0p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_pt_0p1_0p0;
						SubjetGen1_z0p1_b0p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_phi_0p1_0p0;
						SubjetGen1_z0p1_b0p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_eta_0p1_0p0;
						SubjetGen1_z0p1_b0p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_phiWTA_0p1_0p0;
						SubjetGen1_z0p1_b0p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_etaWTA_0p1_0p0;
						SubjetGen1_z0p1_b0p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_mass_0p1_0p0;
						SubjetGen2_z0p1_b0p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_pt_0p1_0p0;
						SubjetGen2_z0p1_b0p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_phi_0p1_0p0;
						SubjetGen2_z0p1_b0p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_eta_0p1_0p0;
						SubjetGen2_z0p1_b0p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_phiWTA_0p1_0p0;
						SubjetGen2_z0p1_b0p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_etaWTA_0p1_0p0;
						SubjetGen2_z0p1_b0p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_mass_0p1_0p0;

						SubjetGen1_z0p1_b1p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_pt_0p1_1p0;
						SubjetGen1_z0p1_b1p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_phi_0p1_1p0;
						SubjetGen1_z0p1_b1p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_eta_0p1_1p0;
						SubjetGen1_z0p1_b1p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_phiWTA_0p1_1p0;
						SubjetGen1_z0p1_b1p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_etaWTA_0p1_1p0;
						SubjetGen1_z0p1_b1p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_mass_0p1_1p0;
						SubjetGen2_z0p1_b1p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_pt_0p1_1p0;
						SubjetGen2_z0p1_b1p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_phi_0p1_1p0;
						SubjetGen2_z0p1_b1p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_eta_0p1_1p0;
						SubjetGen2_z0p1_b1p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_phiWTA_0p1_1p0;
						SubjetGen2_z0p1_b1p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_etaWTA_0p1_1p0;
						SubjetGen2_z0p1_b1p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_mass_0p1_1p0;

						SubjetGen1_z0p1_b2p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_pt_0p1_2p0;
						SubjetGen1_z0p1_b2p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_phi_0p1_2p0;
						SubjetGen1_z0p1_b2p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_eta_0p1_2p0;
						SubjetGen1_z0p1_b2p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_phiWTA_0p1_2p0;
						SubjetGen1_z0p1_b2p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_etaWTA_0p1_2p0;
						SubjetGen1_z0p1_b2p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_mass_0p1_2p0;
						SubjetGen2_z0p1_b2p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_pt_0p1_2p0;
						SubjetGen2_z0p1_b2p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_phi_0p1_2p0;
						SubjetGen2_z0p1_b2p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_eta_0p1_2p0;
						SubjetGen2_z0p1_b2p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_phiWTA_0p1_2p0;
						SubjetGen2_z0p1_b2p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_etaWTA_0p1_2p0;
						SubjetGen2_z0p1_b2p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_mass_0p1_2p0;

						SubjetGen1_z0p2_b0p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_pt_0p2_0p0;
						SubjetGen1_z0p2_b0p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_phi_0p2_0p0;
						SubjetGen1_z0p2_b0p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_eta_0p2_0p0;
						SubjetGen1_z0p2_b0p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_phiWTA_0p2_0p0;
						SubjetGen1_z0p2_b0p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_etaWTA_0p2_0p0;
						SubjetGen1_z0p2_b0p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_mass_0p2_0p0;
						SubjetGen2_z0p2_b0p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_pt_0p2_0p0;
						SubjetGen2_z0p2_b0p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_phi_0p2_0p0;
						SubjetGen2_z0p2_b0p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_eta_0p2_0p0;
						SubjetGen2_z0p2_b0p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_phiWTA_0p2_0p0;
						SubjetGen2_z0p2_b0p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_etaWTA_0p2_0p0;
						SubjetGen2_z0p2_b0p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_mass_0p2_0p0;

						SubjetGen1_z0p4_b0p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_pt_0p4_0p0;
						SubjetGen1_z0p4_b0p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_phi_0p4_0p0;
						SubjetGen1_z0p4_b0p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_eta_0p4_0p0;
						SubjetGen1_z0p4_b0p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_phiWTA_0p4_0p0;
						SubjetGen1_z0p4_b0p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_etaWTA_0p4_0p0;
						SubjetGen1_z0p4_b0p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_mass_0p4_0p0;
						SubjetGen2_z0p4_b0p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_pt_0p4_0p0;
						SubjetGen2_z0p4_b0p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_phi_0p4_0p0;
						SubjetGen2_z0p4_b0p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_eta_0p4_0p0;
						SubjetGen2_z0p4_b0p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_phiWTA_0p4_0p0;
						SubjetGen2_z0p4_b0p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_etaWTA_0p4_0p0;
						SubjetGen2_z0p4_b0p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_mass_0p4_0p0;

						SubjetGen1_z0p5_b1p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_pt_0p5_1p0;
						SubjetGen1_z0p5_b1p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_phi_0p5_1p0;
						SubjetGen1_z0p5_b1p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_eta_0p5_1p0;
						SubjetGen1_z0p5_b1p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_phiWTA_0p5_1p0;
						SubjetGen1_z0p5_b1p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_etaWTA_0p5_1p0;
						SubjetGen1_z0p5_b1p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_mass_0p5_1p0;
						SubjetGen2_z0p5_b1p0_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_pt_0p5_1p0;
						SubjetGen2_z0p5_b1p0_PhiArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_phi_0p5_1p0;
						SubjetGen2_z0p5_b1p0_EtaArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_eta_0p5_1p0;
						SubjetGen2_z0p5_b1p0_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_phiWTA_0p5_1p0;
						SubjetGen2_z0p5_b1p0_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_etaWTA_0p5_1p0;
						SubjetGen2_z0p5_b1p0_MassArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_mass_0p5_1p0;

						SubjetGen1_z0p5_b1p5_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_pt_0p5_1p5;
						SubjetGen1_z0p5_b1p5_PhiArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_phi_0p5_1p5;
						SubjetGen1_z0p5_b1p5_EtaArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_eta_0p5_1p5;
						SubjetGen1_z0p5_b1p5_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_phiWTA_0p5_1p5;
						SubjetGen1_z0p5_b1p5_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_etaWTA_0p5_1p5;
						SubjetGen1_z0p5_b1p5_MassArrayOutput[iJetType][iJetOutput] = SD_subjetgen1_mass_0p5_1p5;
						SubjetGen2_z0p5_b1p5_RawPtArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_pt_0p5_1p5;
						SubjetGen2_z0p5_b1p5_PhiArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_phi_0p5_1p5;
						SubjetGen2_z0p5_b1p5_EtaArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_eta_0p5_1p5;
						SubjetGen2_z0p5_b1p5_PhiWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_phiWTA_0p5_1p5;
						SubjetGen2_z0p5_b1p5_EtaWTAArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_etaWTA_0p5_1p5;
						SubjetGen2_z0p5_b1p5_MassArrayOutput[iJetType][iJetOutput] = SD_subjetgen2_mass_0p5_1p5;

						iJetOutput++;

					} else {nGenJetsOutput[iJetType]--;}
				} // Generator level jet loop
			} // If for filling generator jet loop
			jetTreeOutput[iJetType]->Fill();
    	} // Loop over jet collections

	    // Reco track loop
    	nTracksOutput = nTracks;
    	iTrackOutput = 0;
    	for(int iTrack = 0; iTrack < nTracks; iTrack++){

    		// Difference with respect eta = -5.0
    		if((trackEtaArray[iTrack]+5.0 < track_gap) && (trackPtArray[iTrack] > 0.2) && (trackHighPurityArray[iTrack] == 1)) {track_gap = trackEtaArray[iTrack]+5.0;}

    		passTrackCuts = true;
      
    		// Do basic track cuts for the reconstructed tracks
    		if(trackHighPurityArray[iTrack] != 1) passTrackCuts = false;
    		if(fabs(trackPtErrorArray[iTrack]/trackPtArray[iTrack]) > 0.15) passTrackCuts = false;
    		if(fabs(trackVertexDistanceZArray[iTrack]/trackVertexDistanceZErrorArray[iTrack]) > 5.0) passTrackCuts = false;
    		if(fabs(trackVertexDistanceXYArray[iTrack]/trackVertexDistanceXYErrorArray[iTrack]) > 5.0) passTrackCuts = false;
    		if(fabs(trackEtaArray[iTrack]) >= 2.4) passTrackCuts = false;  //acceptance of the tracker
    		if(trackPtArray[iTrack] <= 0.2) passTrackCuts = false;   // Minimum track pT
      
    		if(passTrackCuts){
    			trackPtOutput[iTrackOutput] = trackPtArray[iTrack];
    			trackPtErrorOutput[iTrackOutput] = trackPtErrorArray[iTrack];
    			trackPhiOutput[iTrackOutput] = trackPhiArray[iTrack];
        		trackEtaOutput[iTrackOutput] = trackEtaArray[iTrack];
        		trackHighPurityOutput[iTrackOutput] = trackHighPurityArray[iTrack];
        		trackVertexDistanceZOutput[iTrackOutput] = trackVertexDistanceZArray[iTrack];
        		trackVertexDistanceZErrorOutput[iTrackOutput] = trackVertexDistanceZErrorArray[iTrack];
        		trackVertexDistanceXYOutput[iTrackOutput] = trackVertexDistanceXYArray[iTrack];
        		trackVertexDistanceXYErrorOutput[iTrackOutput] = trackVertexDistanceXYErrorArray[iTrack];
        		trackEnergyEcalOutput[iTrackOutput] = trackEnergyEcalArray[iTrack];
        		trackEnergyHcalOutput[iTrackOutput] = trackEnergyHcalArray[iTrack];
        		trackChargeOutput[iTrackOutput] = trackChargeArray[iTrack];
        		PixelnHitsTrackOutput[iTrackOutput] = PixelnHitsTrackArray[iTrack];
        		iTrackOutput++;
      		} else {nTracksOutput--;} 
    	}
    	trackTreeOutput->Fill();

		if(storepfcand){
	   		for(int ipfTrack = 0; ipfTrack < particleFlowCandidatePtVector->size(); ipfTrack++){
   				// Cut away low pT tracks and tracks with eta outside of tracker acceptance
		        // Fill the output vectors with gen particles surviving the cuts
	       		pfTrackPtVector->push_back(particleFlowCandidatePtVector->at(ipfTrack));
    	   		pfTrackPhiVector->push_back(particleFlowCandidatePhiVector->at(ipfTrack));
       			pfTrackEtaVector->push_back(particleFlowCandidateEtaVector->at(ipfTrack));
    			pfTrackEnergyVector->push_back(particleFlowCandidateEnergyVector->at(ipfTrack));
 	      		pfTrackMassVector->push_back(particleFlowCandidateMassVector->at(ipfTrack));
 	     		pfTrackIDVector->push_back(particleFlowCandidateIDVector->at(ipfTrack));
 	     	}
      
      		pfTrackTreeOutput->Fill();
    		// Clear the vectors before the next event! Otherwise all the tracks pile up cumulatively
			pfTrackPtVector->clear();
     		pfTrackPhiVector->clear();
      		pfTrackEtaVector->clear();
      		pfTrackEnergyVector->clear();
      		pfTrackMassVector->clear();
      		pfTrackIDVector->clear();
	    }

    	// Gen track loop
    	if(is_MC){
    		for(int igTrack = 0; igTrack < genTrackPtArray->size(); igTrack++){
    			// Cut away low pT tracks and tracks with eta outside of tracker acceptance
				if(TMath::Abs(genTrackEtaArray->at(igTrack)) >= 2.4) continue; //acceptance of the tracker
				if(genTrackPtArray->at(igTrack) <= 0.2) continue;   // Minimum track pT
		        // Fill the output vectors with gen particles surviving the cuts
        		genTrackPtVector->push_back(genTrackPtArray->at(igTrack));
        		genTrackPhiVector->push_back(genTrackPhiArray->at(igTrack));
        		genTrackEtaVector->push_back(genTrackEtaArray->at(igTrack));
       			genTrackPdgVector->push_back(genTrackPdgArray->at(igTrack));
        		genTrackChargeVector->push_back(genTrackChargeArray->at(igTrack));
        		genTrackSubeventVector->push_back(genTrackSubeventArray->at(igTrack));
      		}
      
      		genTrackTreeOutput->Fill();
    		// Clear the vectors before the next event! Otherwise all the tracks pile up cumulatively
       		genTrackPtVector->clear();
      		genTrackPhiVector->clear();
      		genTrackEtaVector->clear();
      		genTrackPdgVector->clear();
      		genTrackChargeVector->clear();
      		genTrackSubeventVector->clear();

    	} // Filling gen tracks for MC

        TH1D htempRG("htempRG", "htempRG", 20, -5., 5.);
        TH1D htempRG_noNsel("htempRG_noNsel", "htempRG_noNsel", 20, -5., 5.);

    	for(int iPF = 0; iPF < particleFlowCandidateEnergyVector->size(); iPF++){

    		double pfeta = particleFlowCandidateEtaVector->at(iPF);
            double fabspfeta = fabs(pfeta);
            double pfenergy = particleFlowCandidateEnergyVector->at(iPF);
            int pfid = particleFlowCandidateIDVector->at(iPF);
                
            if(fabspfeta < 2.5){htempRG.Fill(pfeta, pfenergy); htempRG_noNsel.Fill(pfeta, pfenergy);
            }else if(fabspfeta >= 2.5 && fabspfeta < 3.0){htempRG_noNsel.Fill(pfeta, pfenergy); if(pfid == 5){htempRG.Fill(pfeta, pfenergy);}
            }else if(fabspfeta >= 3.0){htempRG.Fill(pfeta, pfenergy); htempRG_noNsel.Fill(pfeta, pfenergy);}

    	}

    	float FRG = -99.; 
    	float FRG_noNsel = -99.; 
    	float BRG = -99.; 
    	float BRG_noNsel = -99.;

    	// FRG
        for(int i = 0; i < 15; i++){ if( htempRG.GetBinContent(i+1) > pfE[i] || (track_gap < 0.5*(i+1) ) ){ FRG = 0.5 * i; break;} }
        for(int i = 0; i < 15; i++){ if( htempRG_noNsel.GetBinContent(i+1) > pfE[i] || (track_gap < 0.5*(i+1) ) ){ FRG_noNsel = 0.5 * i; break;} }
        // BRG
		for(int i = 0; i < 20; i++){if ( htempRG.GetBinContent(20-i) > pfE[19-i] ) { BRG = 0.5 * i; break; } }
		for(int i = 0; i < 20; i++){if ( htempRG_noNsel.GetBinContent(20-i) > pfE[19-i] ) { BRG_noNsel = 0.5 * i; break; } }

		hi_FRG = FRG;
		hi_FRG_noNsel = FRG_noNsel;
		hi_BRG = BRG;
		hi_BRG_noNsel = BRG_noNsel;

		heavyIonTreeOutput->Fill();


	} // End loop over events

	// Write the skimmed trees to the output file
  	TFile *outputFile = new TFile(outputFileName, "RECREATE");
  	outputFile->SetCompressionLevel(1);

	gDirectory->mkdir("hiEvtAnalyzer");
	gDirectory->cd("hiEvtAnalyzer");
	heavyIonTreeOutput->Write();

	gDirectory->cd("../");
	gDirectory->mkdir("hltanalysis");
	gDirectory->cd("hltanalysis");
	hltTreeOutput->Write();
  
	gDirectory->cd("../");
	gDirectory->mkdir("skimanalysis");
	gDirectory->cd("skimanalysis");
	skimTreeOutput->Write();

	gDirectory->cd("../");
	gDirectory->mkdir("checkflattening");
	gDirectory->cd("checkflattening");
	checkFlatteningTreeOutput->Write();

	gDirectory->cd("../");
	const char *jetDirectories[] = {"ak4CaloJetAnalyzer","ak4PFJetAnalyzer","akCs4PFJetAnalyzer","ak3PFJetAnalyzer"};
	for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
		gDirectory->mkdir(jetDirectories[iJetType]);
		gDirectory->cd(jetDirectories[iJetType]);
		jetTreeOutput[iJetType]->Write();
		gDirectory->cd("../");
	} // Loop over jet types

  	gDirectory->mkdir("hiFJRhoAnalyzer");
  	gDirectory->cd("hiFJRhoAnalyzer");
	RhoTreeOutput->Write();
	gDirectory->cd("../");	

  	if(storetracks) {
  		gDirectory->mkdir("ppTrack");
  		gDirectory->cd("ppTrack");
		trackTreeOutput->Write();
		gDirectory->cd("../");
	
	}
		
  	if(storepfcand) {
  		gDirectory->mkdir("pfcandAnalyzer");
  		gDirectory->cd("pfcandAnalyzer");
		pfTrackTreeOutput->Write();
		gDirectory->cd("../");
	
	}
	// Generator particles only present in MC
	if((storetracks && is_MC) || (storepfcand && is_MC)){
		gDirectory->mkdir("HiGenParticleAna");
		gDirectory->cd("HiGenParticleAna");
		genTrackTreeOutput->Write();
		gDirectory->cd("../");
	}	

	outputFile->Close();

	cout << endl;
	cout << "------------------------------------- SKIMMER DONE --------------------------------------" << endl;
	cout << endl;


	sec_end = clock(); // stop time counting
	cout << "========================================" << endl;
	cout << "Total running time: " << (double)(sec_end - sec_start) / CLOCKS_PER_SEC << " [s]" << endl;
	cout << "========================================" << endl;

	print_stop(); // Print time, date and hour when it stops

}

unsigned long long keyFromRunLumiEvent(UInt_t run, UInt_t lumi, ULong64_t event){

  const unsigned long long runMult = 1;
  const unsigned long long lumiMult = 1000000;
  const unsigned long long evtMult = 10000000000;
  const unsigned long long evtLimit = 10000000000;

  unsigned long long key = 0;
  if(event >= evtLimit){
    std::cout << "RUN-LUMI-EVENTKEY WARNING : \'" << event << "\' is greated that event limit 10^10. returning key 0" << std::endl;
    return key;
  }

  key += runMult*static_cast<unsigned long long>(run);
  key += lumiMult*static_cast<unsigned long long>(lumi);
  key += evtMult*event;
  
  return key;
  
}


int main(int argc, char** argv){
				TString firstArgument(argv[1]);
				TString outfile(argv[2]);
				int mc = atoi(argv[3]);
				int ntrkoffline = atoi(argv[4]);
				int psidearg = atoi(argv[5]);
				TString secArgument(argv[6]);	
				pPbSkim(firstArgument,outfile,mc,ntrkoffline,psidearg,secArgument);
}
