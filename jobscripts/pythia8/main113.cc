// main113.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This test program will generate Pb-Pb collisions at
// sqrt(S_NN)=2.76TeV using the Angantyr model for Heavy Ion
// collisions. The analysis will divide the event in centrality
// classes using the same observable as was used for p-Pb in the ATLAS
// analysis in arXiv:1508.00848 [hep-ex] (see main112.cc). The
// centrality classes are same as in the ALICE analysis in
// arXiv:1012.1657 [nucl-ex] although the actual observable used is
// not the same. Histograms of multiplicity distributions are measured
// for each centrality percentile.

// Note that heavy ion collisions are computationally quite CPU
// intensive and generating a single event will take around a second
// on a reasonable desktop. To get reasonable statistics, this program
// will take a couple of hours to run.

//This is the old main113.cc file with Antonio Silva's changes to file

#include "Pythia8/Pythia.h"

// You need to include this to get access to the HIInfo object for
// HeavyIons.
#include "Pythia8/HeavyIons.h"
#include "Pythia8Plugins/Pythia8Rivet.h"
//#include "HepMCInterface.h"
//#include "Pythia8Plugins/HepMC2.h"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
//Heavy ion for HepMC
#include "HepMC/HeavyIon.h"
// Following line to be used with HepMC 2.04 onwards.
#ifdef HEPMC_HAS_UNITS
#include "HepMC/Units.h"
#endif
#include <cstdio>

using namespace Pythia8;

int main(int argc, char *argv[]) {
//int main() {

    //float pTHatMin = atof(argv[1]);
    //float pTHatMax = atof(argv[2]);


  // Interface for conversion from Pythia8::Event to HepMC one.
  HepMC::Pythia8ToHepMC ToHepMC;

  // Specify file where HepMC events will be stored.
  //const string outFile = "hepmc_PbPb_pTHatInfo_CentLim_" + std::to_string(pTHatMin) + "_" + std::to_string(pTHatMax) + ".hepmc";

  std::string directory = argv[1];
  std::string filename = argv[2];
  std::string sNevents = argv[3];
  std::string sbeam1 = argv[4];
  std::string sbeam2 = argv[5];
  std::string cmsEnergy = "Beams:eCM = " + string(argv[6]);
  std::string seed = "Random:seed = " + string(argv[7]);

  //string fnumber = argv[1];
  //const string outFile = "hepmc_PbPb_MB_" + fnumber + ".hepmc";
  HepMC::IO_GenEvent ascii_io(directory + "/" + filename, std::ios::out);

  std::string beam1 = "Beams:idA = ";
  std::string beam2 = "Beams:idB = ";

  if(sbeam1.compare("p") == 0) beam1 += "2212";
  else if(sbeam1.compare("Au") == 0) beam1 += "1000791970";
  else if(sbeam1.compare("Cu") == 0) beam1 += "1000290630";
  else if(sbeam1.compare("U") == 0) beam1 += "1000922380";
  else if(sbeam1.compare("He") == 0) beam1 += "1000020040";
  else if(sbeam1.compare("d") == 0) beam1 += "1000010020";
  else if(sbeam1.compare("Pb") == 0) beam1 += "1000822080";

  if(sbeam2.compare("p") == 0) beam2 += "2212";
  else if(sbeam2.compare("Au") == 0) beam2 += "1000791970";
  else if(sbeam2.compare("Cu") == 0) beam2 += "1000290630";
  else if(sbeam2.compare("U") == 0) beam2 += "1000922380";
  else if(sbeam2.compare("He") == 0) beam2 += "1000020040";
  else if(sbeam2.compare("d") == 0) beam2 += "1000010020";
  else if(sbeam2.compare("Pb") == 0) beam2 += "1000822080";

  Pythia pythia;

  cout << "beam1: " << beam1 << endl;
  cout << "beam2: " << beam2 << endl;
  cout << "energy: " << cmsEnergy << endl;

  pythia.readString(beam1);
  pythia.readString(beam2); // The ions.
  pythia.readString(cmsEnergy);
  pythia.readString("Beams:frameType = 1");
  pythia.readString("PartonLevel:all = on");
  pythia.readString("ProcessLevel:all = on");
  if(sbeam1.compare("p") == 0 && sbeam2.compare("p") == 0) pythia.readString("SoftQCD:inelastic = on");
  //if(sbeam1.compare("p") == 0 && sbeam2.compare("p") == 0) pythia.readString("HardQCD:all = on");
  // Pick new random number seed for each run, based on clock.
  pythia.readString("Random:setSeed = on");
  pythia.readString(seed);

  bool isHeavyIon = false;

  if((sbeam1.compare("p") != 0) || (sbeam2.compare("p") != 0))
  {
      isHeavyIon = true;

      pythia.readString("HeavyIon:mode = 1");
      pythia.readString("HeavyIon:showInit = on");

      // Initialize the Angantyr model to fit the total and semi-includive
      // cross sections in Pythia within some tolerance.
      pythia.readString("HeavyIon:SigFitErr = 0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0");
      // These parameters are typicall suitable for sqrt(S_NN)=5TeV
      pythia.readString("HeavyIon:SigFitDefPar = 17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0");
      // A simple genetic algorithm is run for 20 generations to fit the
      // parameters.
      pythia.readString("HeavyIon:SigFitNGen = 20");
  }

  const string spTHatMin = "PhaseSpace:pTHatMin = " + std::string(argv[8]);
  const string spTHatMax = "PhaseSpace:pTHatMax = " + std::string(argv[9]);
  if(!((atof(argv[8])==0) && (atof(argv[9])==-1)))
  {
      cout << "Using pT-hat bins" << endl;
      pythia.readString("HardQCD:all = on");
      pythia.readString(spTHatMin);
      pythia.readString(spTHatMax);
  }

  // Initialise Pythia.
  pythia.init();

  // Loop over events.
  int nEvents = std::stoi(sNevents);
  while(nEvents > 0)
  {
    if ( !pythia.next() ) continue;

    //HeavyIons *hipythia = pythia.getHeavyIonsPtr();

    //cout << "Event #" << 100-nEvents << endl;
    //Set HeavyIon with arbitrary values
    HepMC::HeavyIon ion;

    if(isHeavyIon)
    {
        ion.set_Ncoll_hard(pythia.info.hiinfo->nCollNDTot());
        ion.set_Ncoll(pythia.info.hiinfo->nAbsProj() +
                    pythia.info.hiinfo->nDiffProj() +
                    pythia.info.hiinfo->nAbsTarg() +
                    pythia.info.hiinfo->nDiffTarg() -
                    pythia.info.hiinfo->nCollND() -
                    pythia.info.hiinfo->nCollDD());
        ion.set_Npart_proj(pythia.info.hiinfo->nAbsProj() +
                        pythia.info.hiinfo->nDiffProj());
        ion.set_Npart_targ(pythia.info.hiinfo->nAbsTarg() +
                        pythia.info.hiinfo->nDiffTarg());

        //ion.set_impact_parameter(pythia.info.pTHat()); //This is a trick I use if I want to store the pt-hat of the event on HepMC
        ion.set_impact_parameter(pythia.info.hiinfo->b());
    }

    double x1=0., x2=0., q=0., xf1=0., xf2=0.;
    HepMC::PdfInfo pdf( 2, 3, x1, x2, q, xf1, xf2, 230, 230);

    // Construct new empty HepMC event.
#ifdef HEPMC_HAS_UNITS
    // This form with arguments is only meaningful for HepMC 2.04 onwards,
    // and even then unnecessary if HepMC was built with GeV and mm as units..
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
    //HepMC::GenEvent* hepmcevt = new HepMC::GenEvent( HepMC::Units::GEV, HepMC::Units::MM, 0, 0, 0, std::vector<double>(), std::vector<long>(), ion, pdf);
#else
    // This form is needed for backwards compatibility.
    // In HepMCInterface.cc a conversion from GeV to MeV will be done.
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
#endif

    if(isHeavyIon) hepmcevt->set_heavy_ion(ion);

    // Fill HepMC event, including PDF info.
    ToHepMC.fill_next_event( pythia, hepmcevt );
    // This alternative older method fills event, without PDF info.
    // ToHepMC.fill_next_event( pythia.event, hepmcevt );

    // Write the HepMC event to file. Done with it.
    ascii_io << hepmcevt;
    delete hepmcevt;

    nEvents--;
  }

  // The run is over, so we write out some statistics.

  // And we're done!
  return 0;
}