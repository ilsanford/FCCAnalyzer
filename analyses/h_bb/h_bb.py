
import functions
import helpers
import ROOT
import argparse
import logging

logger = logging.getLogger("fcclogger")

parser = functions.make_def_argparser()
args = parser.parse_args()
functions.set_threads(args)

functions.add_include_file("analyses/higgs_mass_xsec/functions.h")
functions.add_include_file("analyses/higgs_mass_xsec/functions_gen.h")


# define histograms

bins_m = (250, 0, 250)
bins_p = (200, 0, 200)
bins_m_zoom = (200, 110, 130) # 100 MeV


bins_theta = (500, 0, 5)
bins_phi = (400, -4, 4)

bins_count = (100, 0, 100)
bins_pdgid = (60, -30, 30)
bins_charge = (10, -5, 5)

bins_resolution = (10000, 0.95, 1.05)
bins_resolution_1 = (20000, 0, 2)

jet_energy = (1000, 0, 100) # 100 MeV bins
dijet_m = (2000, 0, 200) # 100 MeV bins
visMass = (2000, 0, 200) # 100 MeV bins
missEnergy  = (2000, 0, 200) # 100 MeV bins

dijet_m_final = (500, 50, 100) # 100 MeV bins

bins_cos = (100, -1, 1)
bins_aco = (1000,-360,360)
bins_cosThetaMiss = (10000, 0, 1)



def build_graph(df, dataset):

    logging.info(f"build graph {dataset.name}")
    results, cols = [], []

    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")
    df = helpers.defineCutFlowVars(df) # make the cutX=X variables
    
    # define collections
    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")


    # muons
    df = df.Alias("Muon0", "Muon#0.index")
    df = df.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
    df = df.Define("muons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_all)")
    df = df.Define("muons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons_all)")
    df = df.Define("muons_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons_all)")
    df = df.Define("muons_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons_all)")
    df = df.Define("muons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(muons_all)")

    df = df.Define("muons", "FCCAnalyses::ReconstructedParticle::sel_p(25)(muons_all)")
    df = df.Define("muons_p", "FCCAnalyses::ReconstructedParticle::get_p(muons)")
    df = df.Define("muons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons)")
    df = df.Define("muons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons)")
    df = df.Define("muons_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons)")
    df = df.Define("muons_no", "FCCAnalyses::ReconstructedParticle::get_n(muons)")

    
    # electrons
    df = df.Alias("Electron0", "Electron#0.index")
    df = df.Define("electrons_all", "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)")
    df = df.Define("electrons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons_all)")
    df = df.Define("electrons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(electrons_all)")
    df = df.Define("electrons_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(electrons_all)")
    df = df.Define("electrons_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(electrons_all)")
    df = df.Define("electrons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(electrons_all)")

    df = df.Define("electrons", "FCCAnalyses::ReconstructedParticle::sel_p(25)(electrons_all)")
    df = df.Define("electrons_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons)")
    df = df.Define("electrons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(electrons)")
    df = df.Define("electrons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(electrons)")
    df = df.Define("electrons_q", "FCCAnalyses::ReconstructedParticle::get_charge(electrons)")
    df = df.Define("electrons_no", "FCCAnalyses::ReconstructedParticle::get_n(electrons)")


    # jets (bugged)
    #df = df.Alias("Jet0", "Jet#0.index")
    #df = df.Define("jets_all", "FCCAnalyses::ReconstructedParticle::get(Jet0, ReconstructedParticles)")
    #df = df.Define("jets_all_p", "FCCAnalyses::ReconstructedParticle::get_p(jets_all)")
    #df = df.Define("jets_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(jets_all)")
    #df = df.Define("jets_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(jets_all)")
    #df = df.Define("jets_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(jets_all)")
    #df = df.Define("jets_all_no", "FCCAnalyses::ReconstructedParticle::get_n(jets_all)")

    #df = df.Define("jets1", "FCCAnalyses::ReconstructedParticle::sel_p(25)(jets_all)")
    #df = df.Define("jets_p", "FCCAnalyses::ReconstructedParticle::get_p(jets1)")
    #df = df.Define("jets_theta", "FCCAnalyses::ReconstructedParticle::get_theta(jets1)")
    #df = df.Define("jets_phi", "FCCAnalyses::ReconstructedParticle::get_phi(jets1)")
    #df = df.Define("jets_q", "FCCAnalyses::ReconstructedParticle::get_charge(jets1)")
    #df = df.Define("jets_no", "FCCAnalyses::ReconstructedParticle::get_n(jets1)")


    # lepton kinematic histograms
    results.append(df.Histo1D(("muons_all_p_cut0", "", *bins_p), "muons_all_p"))
    results.append(df.Histo1D(("muons_all_theta_cut0", "", *bins_theta), "muons_all_theta"))
    results.append(df.Histo1D(("muons_all_phi_cut0", "", *bins_phi), "muons_all_phi"))
    results.append(df.Histo1D(("muons_all_q_cut0", "", *bins_charge), "muons_all_q"))
    results.append(df.Histo1D(("muons_all_no_cut0", "", *bins_count), "muons_all_no"))

    results.append(df.Histo1D(("electrons_all_p_cut0", "", *bins_p), "electrons_all_p"))
    results.append(df.Histo1D(("electrons_all_theta_cut0", "", *bins_theta), "electrons_all_theta"))
    results.append(df.Histo1D(("electrons_all_phi_cut0", "", *bins_phi), "electrons_all_phi"))
    results.append(df.Histo1D(("electrons_all_q_cut0", "", *bins_charge), "electrons_all_q"))
    results.append(df.Histo1D(("electrons_all_no_cut0", "", *bins_count), "electrons_all_no"))


    #########
    ### CUT 0: all events
    #########
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0"))
    
    
    #########
    ### CUT 1: cos theta(miss)
    #########
    df = df.Define("missingEnergy_rp", "FCCAnalyses::missingEnergy(240., ReconstructedParticles)")
    df = df.Define("missingEnergy", "missingEnergy_rp[0].energy")
    df = df.Define("cosTheta_miss", "FCCAnalyses::get_cosTheta_miss(missingEnergy_rp)")
    results.append(df.Histo1D(("cosThetaMiss_nOne", "", *bins_cosThetaMiss), "cosTheta_miss"))
    df = df.Filter("cosTheta_miss < 0.98")
    results.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1"))
    

    #########
    ### CUT 2: select Z decay product
    #########
    select_mumu = "muons_no == 2 && electrons_no == 0 && missingEnergy < 30"
    select_ee   = "muons_no == 0 && electrons_no == 2 && missingEnergy < 30"
    select_nunu = "muons_no == 0 && electrons_no == 0 && missingEnergy > 95"
    select_qq = f"! ( ({select_mumu}) || ({select_ee}) || ({select_nunu}) )"
    
    df_mumu = df.Filter(select_mumu)
    df_ee   = df.Filter(select_ee)
    df_nunu = df.Filter(select_nunu)
    df_qq   = df.Filter(select_qq)
    
    results.append(df_mumu.Histo1D(("cutFlow", "", *bins_count), "cut2"))
    #results.append(df_ee.Histo1D(("cutFlow", "", *bins_count), "cut2"))
    #results.append(df_nunu.Histo1D(("cutFlow", "", *bins_count), "cut2"))
    #results.append(df_qq.Histo1D(("cutFlow", "", *bins_count), "cut2"))


    #########
    ### CUT 3: we want to detect Z->mumu/ee (so we don't have Z->qq interfering with measurement of H->bb)
    ###        Z->qq and Z->nunu do not get the next few cuts
    #########

    # build the Z resonance based on the available leptons. Returns the best lepton pair compatible with the Z mass and recoil at 125 GeV
    # technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system (Z), index and 2 the leptons of the pair
    
    # muons
    df_mumu = df_mumu.Define("hbuilder_result", "FCCAnalyses::resonanceBuilder_mass_recoil(91.2, 125, 0, 240, false)(muons, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    
    df_mumu = df_mumu.Filter("hbuilder_result.size() > 0")
    results.append(df_mumu.Histo1D(("cutFlow", "", *bins_count), "cut3"))

    df_mumu = df_mumu.Define("zmumu", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{hbuilder_result[0]}") # the Z
    df_mumu = df_mumu.Define("zmumu_tlv", "FCCAnalyses::makeLorentzVectors(zmumu)") # the muons
    df_mumu = df_mumu.Define("zmumu_leps", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{hbuilder_result[1],hbuilder_result[2]}")
    df_mumu = df_mumu.Define("zmumu_leps_tlv", "FCCAnalyses::makeLorentzVectors(zmumu_leps)")

    df_mumu = df_mumu.Define("zmumu_m", "FCCAnalyses::ReconstructedParticle::get_mass(zmumu)[0]")
    df_mumu = df_mumu.Define("zmumu_p", "FCCAnalyses::ReconstructedParticle::get_p(zmumu)[0]")
    df_mumu = df_mumu.Define("zmumu_recoil", "FCCAnalyses::ReconstructedParticle::recoilBuilder(240)(zmumu)")
    df_mumu = df_mumu.Define("zmumu_recoil_m", "FCCAnalyses::ReconstructedParticle::get_mass(zmumu_recoil)[0]")

    # electrons
    df_ee = df_ee.Define("hbuilder_result", "FCCAnalyses::resonanceBuilder_mass_recoil(91.2, 125, 0, 240, false)(electrons, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    
    df_ee = df_ee.Filter("hbuilder_result.size() > 0")
    #results.append(df_ee.Histo1D(("cutFlow", "", *bins_count), "cut3"))

    df_ee = df_ee.Define("zee", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{hbuilder_result[0]}") # the Z
    df_ee = df_ee.Define("zee_tlv", "FCCAnalyses::makeLorentzVectors(zee)") # the electrons
    df_ee = df_ee.Define("zee_leps", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{hbuilder_result[1],hbuilder_result[2]}")
    df_ee = df_ee.Define("zee_leps_tlv", "FCCAnalyses::makeLorentzVectors(zee_leps)")

    df_ee = df_ee.Define("zee_m", "FCCAnalyses::ReconstructedParticle::get_mass(zee)[0]")
    df_ee = df_ee.Define("zee_p", "FCCAnalyses::ReconstructedParticle::get_p(zee)[0]")
    df_ee = df_ee.Define("zee_recoil", "FCCAnalyses::ReconstructedParticle::recoilBuilder(240)(zee)")
    df_ee = df_ee.Define("zee_recoil_m", "FCCAnalyses::ReconstructedParticle::get_mass(zee_recoil)[0]")



    #########
    ### CUT 4: recoil cut (H mass)
    #########  
    results.append(df_mumu.Histo1D(("mumu_recoil_m_nOne", "", *bins_m), "zmumu_recoil_m"))
    df_mumu = df_mumu.Filter("zmumu_recoil_m > 100 && zmumu_recoil_m < 150")
    results.append(df_mumu.Histo1D(("cutFlow", "", *bins_count), "cut4"))
    
    results.append(df_ee.Histo1D(("ee_recoil_m_nOne", "", *bins_m), "zee_recoil_m"))
    df_ee = df_ee.Filter("zee_recoil_m > 100 && zee_recoil_m < 150")
    #results.append(df_ee.Histo1D(("cutFlow", "", *bins_count), "cut4"))


    #########
    ### CUT 5: momentum
    #########
    results.append(df_mumu.Histo1D(("mumu_p_nOne", "", *bins_p), "zmumu_p"))
    df_mumu = df_mumu.Filter("zmumu_p > 35 && zmumu_p < 65")
    results.append(df_mumu.Histo1D(("cutFlow", "", *bins_count), "cut5"))

    results.append(df_ee.Histo1D(("ee_p_nOne", "", *bins_p), "zee_p"))
    df_ee = df_ee.Filter("zee_p > 35 && zee_p < 65")
    #results.append(df_ee.Histo1D(("cutFlow", "", *bins_count), "cut5"))

    #########
    ### CUT 6: acolinearity
    #########
    df_mumu = df_mumu.Define("acoplanarity", "FCCAnalyses::acoplanarity(zmumu_leps)")
    df_mumu = df_mumu.Define("acolinearity", "FCCAnalyses::acolinearity(zmumu_leps)")
    results.append(df_mumu.Histo1D(("acolinearity_mumu", "", *bins_aco), "acolinearity"))

    df_mumu = df_mumu.Filter("acolinearity > 0.05")
    results.append(df_mumu.Histo1D(("cutFlow", "", *bins_count), "cut6"))


    df_ee = df_ee.Define("acoplanarity", "FCCAnalyses::acoplanarity(zee_leps)")
    df_ee = df_ee.Define("acolinearity", "FCCAnalyses::acolinearity(zee_leps)")
    results.append(df_ee.Histo1D(("acolinearity_ee", "", *bins_aco), "acolinearity"))
    
    df_ee = df_ee.Filter("acolinearity > 0.05")
    #results.append(df_ee.Histo1D(("cutFlow", "", *bins_count), "cut6"))


    #########
    ### CUT 7: cut on Z mass
    #########
    results.append(df_mumu.Histo1D(("zmumu_m_nOne", "", *bins_m), "zmumu_m"))
    df_mumu = df_mumu.Filter("zmumu_m > 80 && zmumu_m < 100")
    results.append(df_mumu.Histo1D(("cutFlow", "", *bins_count), "cut7"))

    results.append(df_ee.Histo1D(("zee_m_nOne", "", *bins_m), "zee_m"))
    df_ee = df_ee.Filter("zee_m > 80 && zee_m < 100")
    #results.append(df_ee.Histo1D(("cutFlow", "", *bins_count), "cut7"))

    results.append(df_mumu.Histo1D(("zmumu_m_nocat", "", *bins_m_zoom), "zmumu_m"))


    # jet clustering
    for leps, df in [("muons", df_mumu), ("electrons", df_ee)]:
        # define PF candidates collection by removing the leptons
        df = df.Define("rps_no_leps", f"FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, {leps})")
        df = df.Define("RP_px", "FCCAnalyses::ReconstructedParticle::get_px(rps_no_leps)")
        df = df.Define("RP_py", "FCCAnalyses::ReconstructedParticle::get_py(rps_no_leps)")
        df = df.Define("RP_pz","FCCAnalyses::ReconstructedParticle::get_pz(rps_no_leps)")
        df = df.Define("RP_e", "FCCAnalyses::ReconstructedParticle::get_e(rps_no_leps)")
        df = df.Define("RP_m", "FCCAnalyses::ReconstructedParticle::get_mass(rps_no_leps)")
        df = df.Define("RP_q", "FCCAnalyses::ReconstructedParticle::get_charge(rps_no_leps)")
        df = df.Define("pseudo_jets", "FCCAnalyses::JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
        
        df = df.Define("clustered_jets", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets)")
        df = df.Define("jets", "FCCAnalyses::JetClusteringUtils::get_pseudoJets(clustered_jets)")
        df = df.Define("jetconstituents", "FCCAnalyses::JetClusteringUtils::get_constituents(clustered_jets)")
        df = df.Define("jets_e", "FCCAnalyses::JetClusteringUtils::get_e(jets)")
        df = df.Define("jets_px", "FCCAnalyses::JetClusteringUtils::get_px(jets)")
        df = df.Define("jets_py", "FCCAnalyses::JetClusteringUtils::get_py(jets)")
        df = df.Define("jets_pz", "FCCAnalyses::JetClusteringUtils::get_pz(jets)")
        df = df.Define("jets_m", "FCCAnalyses::JetClusteringUtils::get_m(jets)")
        
        df = df.Define("jet1", "ROOT::Math::PxPyPzEVector(jets_px[0], jets_py[0], jets_pz[0], jets_e[0])")
        df = df.Define("jet2", "ROOT::Math::PxPyPzEVector(jets_px[1], jets_py[1], jets_pz[1], jets_e[1])")
        df = df.Define("jet1_p", "jet1.P()")
        df = df.Define("jet2_p", "jet2.P()")
        df = df.Define("dijet", "jet1+jet2")
        df = df.Define("dijet_m", "dijet.M()")
        df = df.Define("dijet_p", "dijet.P()")
        
        results.append(df.Histo1D((f"z{leps}_m", "", *bins_m), "dijet_m"))


    return results, weightsum


if __name__ == "__main__":

    datadict = functions.get_datadicts() # get default datasets

    datasets_sig = ["wzp6_ee_nunuH_Hbb_ecm240", "wzp6_ee_eeH_Hbb_ecm240", "wzp6_ee_tautauH_Hbb_ecm240", "wzp6_ee_ccH_Hbb_ecm240", "wzp6_ee_bbH_Hbb_ecm240", "wzp6_ee_qqH_Hbb_ecm240", "wzp6_ee_ssH_Hbb_ecm240", "wzp6_ee_mumuH_Hbb_ecm240"]
    datasets_bkg = ["p8_ee_WW_ecm240", "p8_ee_ZZ_ecm240"]# "wzp6_ee_mumu_ecm240", "wzp6_ee_tautau_ecm240", "wzp6_egamma_eZ_Zmumu_ecm240", "wzp6_gammae_eZ_Zmumu_ecm240", "wzp6_gaga_mumu_60_ecm240", "wzp6_gaga_tautau_60_ecm240", "wzp6_ee_nuenueZ_ecm240"]

    datasets_to_run = datasets_sig + datasets_bkg
    result = functions.build_and_run(datadict, datasets_to_run, build_graph, f"output_h_bb.root", args, norm=True, lumi=7200000)
