// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PHENIX_2011_I886590 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PHENIX_2011_I886590);

    // Function to get bin centers for scaling
    bool getDeltaPt(YODA::Histo1D hist, double pT, double &deltaPt)
    {
        if(pT > hist.xMin() && pT < hist.xMax())
        {
        	deltaPt = hist.bin(hist.binIndexAt(pT)).xMid();
                return true;
        }
        else return false;
    }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const ALICE::PrimaryParticles fsPI0(Cuts::abseta < 0.35 && Cuts::pT > 0.3*GeV && Cuts::pT < 5*GeV);
      declare(fsPI0, "fsPI0");
      
      const ALICE::PrimaryParticles fsPI(Cuts::abseta < 0.35 && Cuts::pT > 0.3*GeV && Cuts::pT < 3*GeV);
      declare(fsPI, "fsPI");

      const ALICE::PrimaryParticles fsK(Cuts::abseta < 0.35 && Cuts::pT > 0.4*GeV && Cuts::pT < 2*GeV);
      declare(fsK, "fsK");

      const ALICE::PrimaryParticles fsP(Cuts::abseta < 0.35 && Cuts::pT > 0.5*GeV && Cuts::pT < 4.5*GeV);
      declare(fsP, "fsP");

      // Beam options
      beamOpt = getOption<string>("beam", "NONE");
      if (beamOpt == "PP200") collsys = pp200;
      if (beamOpt == "PP62") collsys = pp62;

      // Book histograms
      // Histos from HEPdata at 200GeV
      book(_h["xsec_piplus_200"], 1, 1, 1);
      book(_h["xsec_piminus_200"], 1, 1, 2);

      book(_h["xsec_kplus_200"], 2, 1, 1);
      book(_h["xsec_kminus_200"], 2, 1, 2);

      book(_h["xsec_p_noFD_200_1"], 3, 1, 1);
      book(_h["xsec_p_noFD_200_2"], 9, 1, 1);
      book(_h["xsec_pbar_noFD_200_1"], 3, 1, 2);
      book(_h["xsec_pbar_noFD_200_2"], 9, 1, 2);

      book(_h["xsec_p_withFD_200_1"], 4, 1, 1);
      book(_h["xsec_p_withFD_200_2"], 10, 1, 1);
      book(_h["xsec_pbar_withFD_200_1"], 4, 1, 2);
      book(_h["xsec_pbar_withFD_200_2"], 10, 1, 2);

      // Histos from HEPdata at 62.4GeV
      book(_h["xsec_piplus_62"], 5, 1, 1);
      book(_h["xsec_piminus_62"], 5, 1, 2);

      book(_h["xsec_kplus_62"], 6, 1, 1);
      book(_h["xsec_kminus_62"], 6, 1, 2);

      book(_h["xsec_p_noFD_62"], 7, 1, 1);
      book(_h["xsec_pbar_noFD_62"], 7, 1, 2);

      book(_h["xsec_p_withFD_62"], 8, 1, 1);
      book(_h["xsec_pbar_withFD_62"], 8, 1, 2);

      // Counter histos for event weights
      book(_c["sow_pp200"], "_sow_pp200");
      book(_c["sow_pp62"], "_sow_pp62");

      //Cross sections
      book(_c["xsec_pp200"], "xsec_pp200");
      book(_c["xsec_pp62"], "xsec_pp62");

      // New Histos added by Shannon Harris

      std::cout << "BOOKED ORIGINAL HISTOS" << std::endl;

      //Figure 5 (paper) Table 10 (HEPData)

      //pi0 ratio 200 GeV/62.4 GeV
      // book(_h["pi0_200A"], 11, 1, 1);
      // book(_h["pi0_62A"], 11, 1, 1);
      // book(_s["pi0_200/62A"], 11, 1, 1);

      // //pi0 ratio 200 GeV/62.4 GeV
      string refname1 = mkAxisCode(11, 1, 1);
      const Scatter2D& refdata1 = refData(refname1);
      book(_h["pi0_200A"], refname1 + "pi0_200A", refdata1);
      book(_h["pi0_62A"], refname1 + "pi0_62A", refdata1);
      book(_s["pi0_200/62A"], refname1, refdata1);

      //std::cout << "booked 11,1,1" << std::endl;

      //pi+ ratio 200 GeV/62.4 GeV
      // book(_h["piplus_200B"], 12, 1, 1);
      // book(_h["piplus_62B"], 12, 1, 1);
      // book(_s["piplus_200/62B"], 12, 1, 1);

      //pi+ ratio 200 GeV/62.4 GeV
      string refname2 = mkAxisCode(12, 1, 1);
      const Scatter2D& refdata2 = refData(refname2);
      book(_h["piplus_200B"], refname2 + "piplus_200B", refdata2);
      book(_h["piplus_62B"], refname2 + "piplus_62B", refdata2);
      book(_s["piplus_200/62B"], refname2, refdata2);

      //std::cout << "booked 12,1,1" << std::endl;

      //pi- ratio 200 GeV/62.4 GeV
      // book(_h["piminus_200C"], 13, 1, 1);
      // book(_h["piminus_62C"], 13, 1, 1);
      // book(_s["piminus_200/62C"], 13, 1, 1);

      //pi- ratio 200 GeV/62.4 GeV
      string refname3 = mkAxisCode(13, 1, 1);
      const Scatter2D& refdata3 = refData(refname3);
      book(_h["piminus_200C"], refname3 + "piminus_200C", refdata3);
      book(_h["piminus_62C"], refname3 + "piminus_62C", refdata3);
      book(_s["piminus_200/62C"], refname3, refdata3);

      //std::cout << "booked 13,1,1" << std::endl;

      //K+ ratio 200 GeV/62.4 GeV
      // book(_h["kplus_200D"], 14, 1, 1);
      // book(_h["kplus_62D"], 14, 1, 1);
      // book(_s["kplus_200/62D"], 14, 1, 1);

      //K+ ratio 200 GeV/62.4 GeV
      string refname4 = mkAxisCode(14, 1, 1);
      const Scatter2D& refdata4 = refData(refname4);
      book(_h["kplus_200D"], refname4 + "kplus_200D", refdata4);
      book(_h["kplus_62D"], refname4 + "kplus_62D", refdata4);
      book(_s["kplus_200/62D"], refname4, refdata4);

      //std::cout << "booked 14,1,1" << std::endl;

      //K- ratio 200 GeV/62.4 GeV
      // book(_h["kminus_200E"], 15, 1, 1);
      // book(_h["kminus_62E"], 15, 1, 1);
      // book(_s["kminus_200/62E"], 15, 1, 1);

      //K- ratio 200 GeV/62.4 GeV
      string refname5 = mkAxisCode(15, 1, 1);
      const Scatter2D& refdata5 = refData(refname5);
      book(_h["kminus_200E"], refname5 + "kminus_200E", refdata5);
      book(_h["kminus_62E"], refname5 + "kminus_62E", refdata5);
      book(_s["kminus_200/62E"], refname5, refdata5);

      //std::cout << "booked 15,1,1" << std::endl;

      //p ratio 200 GeV/62.4 GeV
      // book(_h["p_200F"], 16, 1, 1);
      // book(_h["p_62F"], 16, 1, 1);
      // book(_s["p_200/62F"], 16, 1, 1);

      //p ratio 200 GeV/62.4 GeV
      string refname6 = mkAxisCode(16, 1, 1);
      const Scatter2D& refdata6 = refData(refname6);
      book(_h["p_200F"], refname6 + "p_200F", refdata6);
      book(_h["p_62F"], refname6 + "p_62F", refdata6);
      book(_s["p_200/62F"], refname6, refdata6);

      //std::cout << "booked 16,1,1" << std::endl;

      //pbar ratio 200 GeV/62.4 GeV
      // book(_h["pbar_200G"], 17, 1, 1);
      // book(_h["pbar_62G"], 17, 1, 1);
      // book(_s["pbar_200/62G"], 17, 1, 1);

      //pbar ratio 200 GeV/62.4 GeV
      string refname7 = mkAxisCode(17, 1, 1);
      const Scatter2D& refdata7 = refData(refname7);
      book(_h["pbar_200G"], refname7 + "pbar_200G", refdata7);
      book(_h["pbar_62G"], refname7 + "pbar_62G", refdata7);
      book(_s["pbar_200/62G"], refname7, refdata7);

      //std::cout << "booked 17,1,1" << std::endl;

      // //Figure 11 (paper) Table 11 (HEPData)

      //pi-/pi+ ratio 200 GeV
      // book(_h["piminus_200H"], 18, 1, 1);
      // book(_h["piplus_200H"], 18, 1, 1);
      // book(_s["piminus/piplus_200H"], 18, 1, 1);

      //pi-/pi+ ratio 200 GeV
      string refname8 = mkAxisCode(18, 1, 1);
      const Scatter2D& refdata8 = refData(refname8);
      book(_h["piminus_200H"], refname8 + "piminus_200H", refdata8);
      book(_h["piplus_200H"], refname8 + "piplus_200H", refdata8);
      book(_s["piminus/piplus_200H"], refname8, refdata8);

      //std::cout << "booked 18,1,1" << std::endl;

      //pi-/pi+ ratio 62.4 GeV
      // book(_h["piminus_62I"], 19, 1, 1);
      // book(_h["piplus_62I"], 19, 1, 1);
      // book(_s["piminus/piplus_62I"], 19, 1, 1);

      //pi-/pi+ ratio 62.4 GeV
      string refname9 = mkAxisCode(19, 1, 1);
      const Scatter2D& refdata9 = refData(refname9);
      book(_h["piminus_62I"], refname9 + "piminus_62I", refdata9);
      book(_h["piplus_62I"], refname9 + "piplus_62I", refdata9);
      book(_s["piminus/piplus_62I"], refname9, refdata9);

      //std::cout << "booked 19,1,1" << std::endl;

      //K-/K+ ratio 200 GeV
      // book(_h["kminus_200J"], 20, 1, 1);
      // book(_h["kplus_200J"], 20, 1, 1);
      // book(_s["kminus/kplus_200J"], 20, 1, 1);

      //K-/K+ ratio 200 GeV
      string refname10 = mkAxisCode(20, 1, 1);
      const Scatter2D& refdata10 = refData(refname10);
      book(_h["kminus_200J"], refname10 + "kminus_200J", refdata10);
      book(_h["kplus_200J"], refname10 + "kplus_200J", refdata10);
      book(_s["kminus/kplus_200J"], refname10, refdata10);

      //std::cout << "booked 20,1,1" << std::endl;

      //K-/K+ ratio 62.4 GeV
      // book(_h["kminus_62K"], 21, 1, 1);
      // book(_h["kplus_62K"], 21, 1, 1);
      // book(_s["kminus/kplus_62K"], 21, 1, 1);

      //K-/K+ ratio 62.4 GeV
      string refname11 = mkAxisCode(21, 1, 1);
      const Scatter2D& refdata11 = refData(refname11);
      book(_h["kminus_62K"], refname11 + "kminus_62K", refdata11);
      book(_h["kplus_62K"], refname11 + "kplus_62K", refdata11);
      book(_s["kminus/kplus_62K"], refname11, refdata11);

      //std::cout << "booked 21,1,1" << std::endl;

      //pbar/p ratio 200 GeV
      // book(_h["pbar_200L"], 22, 1, 1);
      // book(_h["p_200L"], 22, 1, 1);
      // book(_s["pbar/p_200L"], 22, 1, 1);



//
      // //pbar/p ratio 200 GeV
      // string refname12 = mkAxisCode(22, 1, 1);
      // const Scatter2D& refdata12 = refData(refname12);
      // book(_h["pbar_200L"], refname12 + "pbar_200L", refdata12);
      // book(_h["p_200L"], refname12 + "p_200L", refdata12);
      // book(_s["pbar/p_200L"], refname12, refdata12);

      // std::cout << "booked 22,1,1" << std::endl;
//



      //pbar/p ratio 62.4 GeV
      // book(_h["pbar_62M"], 23, 1, 1);
      // book(_h["p_62M"], 23, 1, 1);
      // book(_s["pbar/p_62M"], 23, 1, 1);

      //pbar/p ratio 62.4 GeV
      string refname13 = mkAxisCode(23, 1, 1);
      const Scatter2D& refdata13 = refData(refname13);
      book(_h["pbar_62M"], refname13 + "pbar_62M", refdata13);
      book(_h["p_62M"], refname13 + "p_62M", refdata13);
      book(_s["pbar/p_62M"], refname13, refdata13);

      //std::cout << "booked 23,1,1" << std::endl;

      // //Figure 12 (paper) Table 12 (HEPData)

      //K+/pi+ ratio 200 GeV
      // book(_h["kplus_200N"], 24, 1, 1);
      // book(_h["piplus_200N"], 24, 1, 1);
      // book(_s["kplus/piplus_200N"], 24, 1, 1);

      //K+/pi+ ratio 200 GeV
      string refname14 = mkAxisCode(24, 1, 1);
      const Scatter2D& refdata14 = refData(refname14);
      book(_h["kplus_200N"], refname14 + "kplus_200N", refdata14);
      book(_h["piplus_200N"], refname14 + "piplus_200N", refdata14);
      book(_s["kplus/piplus_200N"], refname14, refdata14);

      //std::cout << "booked 24,1,1" << std::endl;

      //K+/pi+ ratio 62.4 GeV
      // book(_h["kplus_62O"], 24, 1, 2);
      // book(_h["piplus_62O"], 24, 1, 2);
      // book(_s["kplus/piplus_62O"], 24, 1, 2);

      //K+/pi+ ratio 62.4 GeV
      string refname15 = mkAxisCode(24, 1, 2);
      const Scatter2D& refdata15 = refData(refname15);
      book(_h["kplus_62O"], refname15 + "kplus_62O", refdata15);
      book(_h["piplus_62O"], refname15 + "piplus_62O", refdata15);
      book(_s["kplus/piplus_62O"], refname15, refdata15);

      //K-/pi- ratio 200 GeV
      // book(_h["kminus_200P"], 25, 1, 1);
      // book(_h["piminus_200P"], 25, 1, 1);
      // book(_s["kminus/piminus_200P"], 25, 1, 1);

      //K-/pi- ratio 200 GeV
      string refname16 = mkAxisCode(25, 1, 1);
      const Scatter2D& refdata16 = refData(refname16);
      book(_h["kminus_200P"], refname16 + "kminus_200P", refdata16);
      book(_h["piminus_200P"], refname16 + "piminus_200P", refdata16);
      book(_s["kminus/piminus_200P"], refname16, refdata16);

      //K-/pi- ratio 62.4 GeV
      // book(_h["kminus_62Q"], 25, 1, 2);
      // book(_h["piminus_62Q"], 25, 1, 2);
      // book(_s["kminus/piminus_62Q"], 25, 1, 2);

      //K-/pi- ratio 62.4 GeV
      string refname17 = mkAxisCode(25, 1, 2);
      const Scatter2D& refdata17 = refData(refname17);
      book(_h["kminus_62Q"], refname17 + "kminus_62Q", refdata17);
      book(_h["piminus_62Q"], refname17 + "piminus_62Q", refdata17);
      book(_s["kminus/piminus_62Q"], refname17, refdata17);

      //p/pi+ ratio 200 GeV
      // book(_h["p_200R"], 26, 1, 1);
      // book(_h["piplus_200R"], 26, 1, 1);
      // book(_s["p/piplus_200R"], 26, 1, 1);

      //p/pi+ ratio 200 GeV
      string refname18 = mkAxisCode(26, 1, 1);
      const Scatter2D& refdata18 = refData(refname18);
      book(_h["p_200R"], refname18 + "p_200R", refdata18);
      book(_h["piplus_200R"], refname18 + "piplus_200R", refdata18);
      book(_s["p/piplus_200R"], refname18, refdata18);

      //p/pi+ ratio 62.4 GeV
      // book(_h["p_62S"], 27, 1, 1);
      // book(_h["piplus_62S"], 27, 1, 1);
      // book(_s["p/piplus_62S"], 27, 1, 1);

      //p/pi+ ratio 62.4 GeV
      string refname19 = mkAxisCode(27, 1, 1);
      const Scatter2D& refdata19 = refData(refname19);
      book(_h["p_62S"], refname19 + "p_62S", refdata19);
      book(_h["piplus_62S"], refname19 + "piplus_62S", refdata19);
      book(_s["p/piplus_62S"], refname19, refdata19);

      //pbar/pi- ratio 200 GeV
      // book(_h["pbar_200T"], 28, 1, 1);
      // book(_h["piminus_200T"], 28, 1, 1);
      // book(_s["pbar/piminus_200T"], 28, 1, 1);

      //pbar/pi- ratio 200 GeV
      string refname20 = mkAxisCode(28, 1, 1);
      const Scatter2D& refdata20 = refData(refname20);
      book(_h["pbar_200T"], refname20 + "pbar_200T", refdata20);
      book(_h["piminus_200T"], refname20 + "piminus_200T", refdata20);
      book(_s["pbar/piminus_200T"], refname20, refdata20);

      //pbar/pi- ratio 62.4 GeV
      // book(_h["pbar_62U"], 29, 1, 1);
      // book(_h["piminus_62U"], 29, 1, 1);
      // book(_s["pbar/piminus_62U"], 29, 1, 1);

      //pbar/pi- ratio 62.4 GeV
      string refname21 = mkAxisCode(29, 1, 1);
      const Scatter2D& refdata21 = refData(refname21);
      book(_h["pbar_62U"], refname21 + "pbar_62U", refdata21);
      book(_h["piminus_62U"], refname21 + "piminus_62U", refdata21);
      book(_s["pbar/piminus_62U"], refname21, refdata21);

      //p/pi0 ratio 200 GeV
      // book(_h["p_200V"], 30, 1, 1);
      // book(_h["pi0_200V"], 30, 1, 1);
      // book(_s["p/pi0_200V"], 30, 1, 1);

      //p/pi0 ratio 200 GeV
      string refname22 = mkAxisCode(30, 1, 1);
      const Scatter2D& refdata22 = refData(refname22);
      book(_h["p_200V"], refname22 + "p_200V", refdata22);
      book(_h["pi0_200V"], refname22 + "pi0_200V", refdata22);
      book(_s["p/pi0_200V"], refname22, refdata22);

      //p/pi0 ratio 62.4 GeV
      // book(_h["p_62W"], 31, 1, 1);
      // book(_h["pi0_62W"], 31, 1, 1);
      // book(_s["p/pi0_62W"], 31, 1, 1);

      //p/pi0 ratio 62.4 GeV
      string refname23 = mkAxisCode(31, 1, 1);
      const Scatter2D& refdata23 = refData(refname23);
      book(_h["p_62W"], refname23 + "p_62W", refdata23);
      book(_h["pi0_62W"], refname23 + "pi0_62W", refdata23);
      book(_s["p/pi0_62W"], refname23, refdata23);

      //pbar/pi0 ratio 200 GeV
      // book(_h["pbar_200X"], 32, 1, 1);
      // book(_h["pi0_200X"], 32, 1, 1);
      // book(_s["pbar/pi0_200X"], 32, 1, 1);

      //pbar/pi0 ratio 200 GeV
      string refname24 = mkAxisCode(32, 1, 1);
      const Scatter2D& refdata24 = refData(refname24);
      book(_h["pbar_200X"], refname24 + "pbar_200X", refdata24);
      book(_h["pi0_200X"], refname24 + "pi0_200X", refdata24);
      book(_s["pbar/pi0_200X"], refname24, refdata24);

      //pbar/pi0 ratio 62.4 GeV
      // book(_h["pbar_62Y"], 33, 1, 1);
      // book(_h["pi0_62Y"], 33, 1, 1);
      // book(_s["pbar/pi0_62Y"], 33, 1, 1);

      //pbar/pi0 ratio 62.4 GeV
      string refname25 = mkAxisCode(33, 1, 1);
      const Scatter2D& refdata25 = refData(refname25);
      book(_h["pbar_62Y"], refname25 + "pbar_62Y", refdata25);
      book(_h["pi0_62Y"], refname25 + "pi0_62Y", refdata25);
      book(_s["pbar/pi0_62Y"], refname25, refdata25);

      std::cout << "BOOKED NEW HISTOS/FINISHED INIT" << std::endl;
    }

    std::cout << "FINISHED INIT; STARTING ANALYZE" << std::endl;

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Particles fsPI0Particles = applyProjection<ALICE::PrimaryParticles>(event,"fsPI0").particles();
      Particles fsPIParticles = applyProjection<ALICE::PrimaryParticles>(event,"fsPI").particles();
      Particles fsKParticles = applyProjection<ALICE::PrimaryParticles>(event,"fsK").particles();
      Particles fsPParticles = applyProjection<ALICE::PrimaryParticles>(event,"fsP").particles();

      pair<double,double> cs = HepMCUtils::crossSection(*event.genEvent());

// std::cout << "starting beamOpt" << std::endl;

      if(beamOpt == "NONE")
        		{

  			int NN = 1.;
  			if (fuzzyEquals(sqrtS()/GeV, 200*NN, 1E-3))
        {
          collsys = pp200;
          _c["sow_pp200"]->fill();
          _c["xsec_pp200"]->fill(cs.first);
        }
  			if (fuzzyEquals(sqrtS()/GeV, 62.4*NN, 1E-3))
        {
          collsys = pp62;
          _c["sow_pp62"]->fill();
          _c["xsec_pp62"]->fill(cs.first);
        }
  		}

// std::cout << "finished beamOpt" << std::endl;


      //Pi zeroes
      for( const Particle& pPI0 : fsPI0Particles)
      {
                double PtPI0 = pPI0.pT()/GeV;
                double PtPI0_weight = 1./(2.*M_PI*0.7);
                double deltaPtPI0 = 0;

    // Fill histos 200GeV

    if (collsys == pp200)
    {

        if(pPI0.pid() == 111)
        {
          _h["pi0_200A"]->fill(pPI0.pT()/GeV);
          _h["pi0_200V"]->fill(pPI0.pT()/GeV);
          _h["pi0_200X"]->fill(pPI0.pT()/GeV);
        }

    }

    //Fill histos 62.4GeV
    if (collsys == pp62)
    {
    
        if(pPI0.pid() == 111)
        {
          _h["pi0_62A"]->fill(pPI0.pT()/GeV);
          _h["pi0_62W"]->fill(pPI0.pT()/GeV);
          _h["pi0_62Y"]->fill(pPI0.pT()/GeV);

    }

      }


      // (charged) Pions
      for( const Particle& pPI : fsPIParticles)
      {
	      	double PtPI = pPI.pT()/GeV;
                double PtPI_weight = 1./(2.*M_PI*0.7);
                double deltaPtPI = 0;

		// Fill histos 200GeV



		if (collsys == pp200)
		{


			if(getDeltaPt(*_h["xsec_piplus_200"],PtPI,deltaPtPI))
			{
      				if(pPI.pid() == 211)
              {
                PtPI_weight /= deltaPtPI;
                _h["xsec_piplus_200"]->fill(pPI.pT()/GeV, PtPI_weight);
              }
			}

			if(getDeltaPt(*_h["xsec_piminus_200"],PtPI,deltaPtPI))
			{
				if(pPI.pid() == -211)
        {
          PtPI_weight /= deltaPtPI;
          _h["xsec_piminus_200"]->fill(pPI.pT()/GeV, PtPI_weight);
        }
			}
		}

    //Shannon
        if(pPI.pid() == 211)
        {
          _h["piplus_200B"]->fill(pPI.pT()/GeV);
          _h["piplus_200H"]->fill(pPI.pT()/GeV);
          _h["piplus_200N"]->fill(pPI.pT()/GeV);
          _h["piplus_200R"]->fill(pPI.pT()/GeV);
        }

        if(pPI.pid() == -211)
        {
          _h["piminus_200C"]->fill(pPI.pT()/GeV);
          _h["piminus_200H"]->fill(pPI.pT()/GeV);
          _h["piminus_200P"]->fill(pPI.pT()/GeV);
          _h["piminus_200T"]->fill(pPI.pT()/GeV);
        }

		// Fill histos 62.4GeV
		if (collsys == pp62)
                {

			if(getDeltaPt(*_h["xsec_piplus_62"],PtPI,deltaPtPI))
			{
                		if(pPI.pid() == 211)
                    {
                      PtPI_weight /= deltaPtPI;
                      _h["xsec_piplus_62"]->fill(pPI.pT()/GeV, PtPI_weight);
                    }
			}

			if(getDeltaPt(*_h["xsec_piminus_62"],PtPI,deltaPtPI))
			{

                		if(pPI.pid() == -211)
                    {
                      PtPI_weight /= deltaPtPI;
                      _h["xsec_piminus_62"]->fill(pPI.pT()/GeV, PtPI_weight);
                    }
			}
		}

      //Shannon
        if(pPI.pid() == 211)
        {
          _h["piplus_62B"]->fill(pPI.pT()/GeV);
          _h["piplus_62I"]->fill(pPI.pT()/GeV);
          _h["piplus_62O"]->fill(pPI.pT()/GeV);
          _h["piplus_62S"]->fill(pPI.pT()/GeV);
        }

        if(pPI.pid() == -211)
        {
          _h["piminus_62C"]->fill(pPI.pT()/GeV);
          _h["piminus_62I"]->fill(pPI.pT()/GeV);
          _h["piminus_62Q"]->fill(pPI.pT()/GeV);
          _h["piminus_62U"]->fill(pPI.pT()/GeV);
        }



      }

      // Kaons
      for( const Particle& pK : fsKParticles)
      {
                double PtK = pK.pT()/GeV;
                double PtK_weight = 1./(2.*M_PI*0.7);
                double deltaPtK = 0;

		// Fill histos 200GeV
		if (collsys == pp200)
		{

                        if(getDeltaPt(*_h["xsec_kplus_200"],PtK,deltaPtK))
                        {

                		if(pK.pid() == 321)
                    {
                      PtK_weight /= deltaPtK;
                      _h["xsec_kplus_200"]->fill(pK.pT()/GeV, PtK_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_kminus_200"],PtK,deltaPtK))
                        {

                		if(pK.pid() == -321)
                    {
                      PtK_weight /= deltaPtK;
                      _h["xsec_kminus_200"]->fill(pK.pT()/GeV, PtK_weight);
                    }
			}

        //Shannon
                if(pK.pid() == 321)
                {
                  _h["kplus_200D"]->fill(pK.pT()/GeV);
                  _h["kplus_200J"]->fill(pK.pT()/GeV);
                  _h["kplus_200N"]->fill(pK.pT()/GeV);
                  _h["kplus_200R"]->fill(pK.pT()/GeV);
                }

                if(pK.pid() == -321)
                {
                  _h["kminus_200E"]->fill(pK.pT()/GeV);
                  _h["kminus_200J"]->fill(pK.pT()/GeV);
                  _h["kminus_200P"]->fill(pK.pT()/GeV);
                  _h["kminus_200T"]->fill(pK.pT()/GeV);
                }

		}

		// Fill histos 62.4GeV
		if (collsys == pp62)
		{

                        if(getDeltaPt(*_h["xsec_kplus_62"],PtK,deltaPtK))
                        {

                		if(pK.pid() == 321)
                    {
                      PtK_weight /= deltaPtK;
                      _h["xsec_kplus_62"]->fill(pK.pT()/GeV, PtK_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_kminus_62"],PtK,deltaPtK))
                        {

                		if(pK.pid() == -321)
                    {
                      PtK_weight /= deltaPtK;
                      _h["xsec_kminus_62"]->fill(pK.pT()/GeV, PtK_weight);
                    }
			}

        //Shannon
                if(pK.pid() == 321)
                {
                  _h["kplus_62D"]->fill(pK.pT()/GeV);
                  _h["kplus_62K"]->fill(pK.pT()/GeV);
                  _h["kplus_62O"]->fill(pK.pT()/GeV);
                }

                if(pK.pid() == -321)
                {
                  _h["kminus_62E"]->fill(pK.pT()/GeV);
                  _h["kminus_62K"]->fill(pK.pT()/GeV);
                  _h["kminus_62Q"]->fill(pK.pT()/GeV);
                }

		}

      }

      // Protons and antiprotons
      for( const Particle& pP : fsPParticles)
      {
                double PtP = pP.pT()/GeV;
                double PtP_weight = 1./(2.*M_PI*0.7);
                double deltaPtP = 0;

		// Fill histos 200GeV
		if (collsys == pp200)
		{

                        if(getDeltaPt(*_h["xsec_p_noFD_200_1"],PtP,deltaPtP))
                        {

                		if(pP.pid() == 2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_p_noFD_200_1"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_p_noFD_200_2"],PtP,deltaPtP))
                        {

                		if(pP.pid() == 2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_p_noFD_200_2"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_pbar_noFD_200_1"],PtP,deltaPtP))
                        {

                		if(pP.pid() == -2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_pbar_noFD_200_1"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_pbar_noFD_200_2"],PtP,deltaPtP))
                        {

                		if(pP.pid() == -2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_pbar_noFD_200_2"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_p_withFD_200_1"],PtP,deltaPtP))
                        {

                		if(pP.pid() == 2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_p_withFD_200_1"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_p_withFD_200_2"],PtP,deltaPtP))
                        {

                		if(pP.pid() == 2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_p_withFD_200_2"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_pbar_withFD_200_1"],PtP,deltaPtP))
                        {

                		if(pP.pid() == -2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_pbar_withFD_200_1"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_pbar_withFD_200_2"],PtP,deltaPtP))
                        {

                		if(pP.pid() == -2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_pbar_withFD_200_2"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

        //Shannon
                if(pP.pid() == 2212)
                {
                  _h["p_200F"]->fill(pP.pT()/GeV);
                  _h["p_200L"]->fill(pP.pT()/GeV);
                  _h["p_200R"]->fill(pP.pT()/GeV);
                  _h["p_200V"]->fill(pP.pT()/GeV);
                }

                if(pP.pid() == -2212)
                {
                  _h["pbar_200G"]->fill(pP.pT()/GeV);
                  // _h["pbar_200L"]->fill(pP.pT()/GeV);
                  _h["pbar_200T"]->fill(pP.pT()/GeV);
                  _h["pbar_200X"]->fill(pP.pT()/GeV);
                }

		}

		// Fill histos 62.4GeV
		if (collsys == pp62)
		{

                        if(getDeltaPt(*_h["xsec_p_noFD_62"],PtP,deltaPtP))
                        {

                		if(pP.pid() == 2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_p_noFD_62"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_pbar_noFD_62"],PtP,deltaPtP))
                        {

                		if(pP.pid() == -2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_pbar_noFD_62"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_p_withFD_62"],PtP,deltaPtP))
                        {

                		if(pP.pid() == 2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_p_withFD_62"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

                        if(getDeltaPt(*_h["xsec_pbar_withFD_62"],PtP,deltaPtP))
                        {

                		if(pP.pid() == -2212)
                    {
                      PtP_weight /= deltaPtP;
                      _h["xsec_pbar_withFD_62"]->fill(pP.pT()/GeV, PtP_weight);
                    }
			}

        //Shannon
                if(pP.pid() == 2212)
                {
                  _h["p_62F"]->fill(pP.pT()/GeV);
                  _h["p_62M"]->fill(pP.pT()/GeV);
                  _h["p_62S"]->fill(pP.pT()/GeV);
                  _h["p_62W"]->fill(pP.pT()/GeV);
                }

                if(pP.pid() == -2212)
                {
                  _h["pbar_62G"]->fill(pP.pT()/GeV);
                  _h["pbar_62M"]->fill(pP.pT()/GeV);
                  _h["pbar_62U"]->fill(pP.pT()/GeV);
                  _h["pbar_62Y"]->fill(pP.pT()/GeV);
                }

		}

      }

    }

    }

    std::cout << "FINISHED ANALYZE; STARTING FINALIZE" << std::endl;

    /// Normalise histograms etc., after the run
    void finalize() {

      double xsec200 = (_c["xsec_pp200"]->sumW()/millibarn)/_c["sow_pp200"]->sumW();
      double xsec62 = (_c["xsec_pp62"]->sumW()/millibarn)/_c["sow_pp62"]->sumW();

		_h["xsec_piplus_200"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
                _h["xsec_piminus_200"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
                _h["xsec_kplus_200"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
                _h["xsec_kminus_200"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
                _h["xsec_p_noFD_200_1"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
                _h["xsec_p_noFD_200_2"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
                _h["xsec_pbar_noFD_200_1"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
                _h["xsec_pbar_noFD_200_2"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
                _h["xsec_p_withFD_200_1"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
                _h["xsec_p_withFD_200_2"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
                _h["xsec_pbar_withFD_200_1"]->scaleW(xsec200/_c["sow_pp200"]->sumW());
                _h["xsec_pbar_withFD_200_2"]->scaleW(xsec200/_c["sow_pp200"]->sumW());


		_h["xsec_piplus_62"]->scaleW(xsec62/_c["sow_pp200"]->sumW());
                _h["xsec_piminus_62"]->scaleW(xsec62/_c["sow_pp200"]->sumW());
                _h["xsec_kplus_62"]->scaleW(xsec62/_c["sow_pp200"]->sumW());
                _h["xsec_kminus_62"]->scaleW(xsec62/_c["sow_pp200"]->sumW());
                _h["xsec_p_noFD_62"]->scaleW(xsec62/_c["sow_pp200"]->sumW());
                _h["xsec_pbar_noFD_62"]->scaleW(xsec62/_c["sow_pp200"]->sumW());
                _h["xsec_p_withFD_62"]->scaleW(xsec62/_c["sow_pp200"]->sumW());
                _h["xsec_pbar_withFD_62"]->scaleW(xsec62/_c["sow_pp200"]->sumW());


      // New Histos added by Shannon Harris

      //Figure 5 (paper) Table 10 (HEPData)

      //pi0 ratio 200 GeV/62.4 GeV
      //_h["pi0_200/62A"]->divide(_h["200"], _h["62"]);
      //alternatively
      divide(_h["pi0_200A"], _h["pi0_62A"], _s["pi0_200/62A"]);

      //pi+ ratio 200 GeV/62.4 GeV
      //_h["piplus_200/62B"]->divide(_h["200"], _h["62"]);
      divide(_h["piplus_200B"], _h["piplus_62B"], _s["piplus_200/62B"]);

      //pi- ratio 200 GeV/62.4 GeV
      //_h["piminus_200/62C"]->divide(_h["200"], _h["62"]);
      divide(_h["piminus_200C"], _h["piminus_62C"], _s["piminus_200/62C"]);
      
      //K+ ratio 200 GeV/62.4 GeV
      //_h["kplus_200/62D"]->divide(_h["200"], _h["62"]);
      divide(_h["kplus_200D"], _h["kplus_62D"], _s["kplus_200/62D"]);

      //K- ratio 200 GeV/62.4 GeV
      //_h["kminus_200/62E"]->divide(_h["200"], _h["62"]);
      divide(_h["kminus_200E"], _h["kminus_62E"], _s["kminus_200/62E"]);

      //p ratio 200 GeV/62.4 GeV
      //_h["p_200/62F"]->divide(_h["200"], _h["62"]);
      divide(_h["p_200F"], _h["p_62F"], _s["p_200/62F"]);

      //pbar ratio 200 GeV/62.4 GeV
      //_h["pbar_200/62G"]->divide(_h["200"], _h["62"]);
      divide(_h["pbar_200G"], _h["pbar_62G"], _s["pbar_200/62G"]);

      //Figure 11 (paper) Table 11 (HEPData)

      //pi-/pi+ ratio 200 GeV
      //hPionNegPosPt["piminus/piplus_200H"]->divide(_h["200"], _h["200"]);
      divide(_h["piminus_200H"], _h["piplus_200H"], _s["piminus/piplus_200H"]);
      
      //pi-/pi+ ratio 62.4 GeV
      //hPionNegPosPt["piminus/piplus_62I"]->divide(_h["62"], _h["62"]);
      divide(_h["piminus_62I"], _h["piplus_62I"], _s["piminus/piplus_62I"]);

      //K-/K+ ratio 200 GeV
      //hKaonNegPosPt["kminus/kplus_200J"]->divide(_h["200"], _h["200"]);
      divide(_h["kminus_200J"], _h["kplus_200J"], _s["kminus/kplus_200J"]);

      //K-/K+ ratio 62.4 GeV
      //hKaonNegPosPt["kminus/kplus_62K"]->divide(_h["62"], _h["62"]);
      divide(_h["kminus_62K"], _h["kplus_62K"], _s["kminus/kplus_62K"]);

      //pbar/p ratio 200 GeV
      //hProtonNegPosPt["pbar/p_200L"]->divide(_h["200"], _h["200"]);
      // divide(_h["pbar_200L"], _h["p_200L"], _s["pbar/p_200L"]);

      //pbar/p ratio 62.4 GeV
      //hProtonNegPosPt["pbar/p_62M"]->divide(_h["62"], _h["62"]);
      divide(_h["pbar_62M"], _h["p_62M"], _s["pbar/p_62M"]);

      //Figure 12 (paper) Table 12 (HEPData)

      //K+/pi+ ratio 200 GeV
      //hKaonPosPionPosPt["kplus/piplus_200N"]->divide(_h["200"], _h["200"]);
      divide(_h["kplus_200N"], _h["piplus_200N"], _s["kplus/piplus_200N"]);

      //K+/pi+ ratio 62.4 GeV
      //hKaonPosPionPosPt["kplus/piplus_62O"]->divide(_h["62"], _h["62"]);
      divide(_h["kplus_62O"], _h["piplus_62O"], _s["kplus/piplus_62O"]);

      //K-/pi- ratio 200 GeV
      //hKaonNegPionNegPt["kminus/piminus_200P"]->divide(_h["200"], _h["200"]);
      divide(_h["kminus_200P"], _h["piminus_200P"], _s["kminus/piminus_200P"]);

      //K-/pi- ratio 62.4 GeV
      //hKaonNegPionNegPt["kminus/piminus_62Q"]->divide(_h["62"], _h["62"]);
      divide(_h["kminus_62Q"], _h["piminus_62Q"], _s["kminus/piminus_62Q"]);

      //p/pi+ ratio 200 GeV
      //hProtonPosPionPosPt["p/piplus_200R"]->divide(_h["200"], _h["200"]);
      divide(_h["p_200R"], _h["piplus_200R"], _s["p/piplus_200R"]);

      //p/pi+ ratio 62.4 GeV
      //hProtonPosPionPosPt["p/piplus_62S"]->divide(_h["62"], _h["62"]);
      divide(_h["p_62S"], _h["piplus_62S"], _s["p/piplus_62S"]);

      //pbar/pi- ratio 200 GeV
      //hProtonNegPionNegPt["pbar/piminus_200T"]->divide(_h["200"], _h["200"]);
      divide(_h["pbar_200T"], _h["piminus_200T"], _s["pbar/piminus_200T"]);

      //pbar/pi- ratio 62.4 GeV
      //hProtonNegPionNegPt["pbar/piminus_62U"]->divide(_h["62"], _h["62"]);
      divide(_h["pbar_62U"], _h["piminus_62U"], _s["pbar/piminus_62U"]);

      //p/pi0 ratio 200 GeV
      //hProtonPosPionZeroPt["p/pi0_200V"]->divide(_h["200"], _h["200"]);
      divide(_h["p_200V"], _h["pi0_200V"], _s["p/pi0_200V"]);

      //p/pi0 ratio 62.4 GeV
      //hProtonPosPionZeroPt["p/pi0_62W"]->divide(_h["62"], _h["62"]);
      divide(_h["p_62W"], _h["pi0_62W"], _s["p/pi0_62W"]);

      //pbar/pi0 ratio 200 GeV
      //hProtonNegPionZeroPt["pbar/pi0_200X"]->divide(_h["200"], _h["200"]);
      divide(_h["pbar_200X"], _h["pi0_200X"], _s["pbar/pi0_200X"]);

      //pbar/pi0 ratio 62.4 GeV
      //hProtonNegPionZeroPt["pbar/pi0_62Y"]->divide(_h["62"], _h["62"]);
      divide(_h["pbar_62Y"], _h["pi0_62Y"], _s["pbar/pi0_62Y"]);

    }

    std::cout << "FINISHED FINALIZE" << std::endl;

    //@}

    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    map<string, Scatter2DPtr> _s;

    // New Histos added by Shannon Harris
    // map<string, Histo1DPtr> _h;
    // map<string, Histo1DPtr> hPionPosPt;
    // map<string, Histo1DPtr> _h;
    // map<string, Histo1DPtr> _h;
    // map<string, Histo1DPtr> _h;
    // map<string, Histo1DPtr> _h;
    // map<string, Histo1DPtr> _h;

    // map<string, Scatter2DPtr> RatioPionZero;
    // map<string, Scatter2DPtr> RatioPion;
    // map<string, Scatter2DPtr> _s;
    // map<string, Scatter2DPtr> _s;
    // map<string, Scatter2DPtr> _s;
    // map<string, Scatter2DPtr> _s;
    // map<string, Scatter2DPtr> _s;
    // map<string, Scatter2DPtr> _s;

    string beamOpt = "NONE";
    enum CollisionSystem {pp200, pp62};
    CollisionSystem collsys;
    //@}

  };


  DECLARE_RIVET_PLUGIN(PHENIX_2011_I886590);


}