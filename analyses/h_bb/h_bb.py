import functions
import helpers
import ROOT
import argparse
import logging

import helper_jetclustering
import helper_flavourtagger


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

# setup clustering and flavour tagging
# 2 jets
jet2Cluster = helper_jetclustering.ExclusiveJetClusteringHelper(2, collection="rps_no_leps")
jet2Flavour = helper_flavourtagger.JetFlavourHelper(jet2Cluster.jets, jet2Cluster.constituents)
# 4 jets
jet4Cluster = helper_jetclustering.ExclusiveJetClusteringHelper(4, collection="ReconstructedParticles")
jet4Flavour = helper_flavourtagger.JetFlavourHelper(jet4Cluster.jets, jet4Cluster.constituents)

path = "data/flavourtagger/fccee_flavtagging_edm4hep_wc_v1"
jet2Flavour.load(f"{path}.json", f"{path}.onnx")
jet4Flavour.load(f"{path}.json", f"{path}.onnx")


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
    results.append(df.Histo1D(("cutFlow_mumu", "", *bins_count), "cut0"))
    results.append(df.Histo1D(("cutFlow_ee", "", *bins_count), "cut0"))
    results.append(df.Histo1D(("cutFlow_nunu", "", *bins_count), "cut0"))
    results.append(df.Histo1D(("cutFlow_qq", "", *bins_count), "cut0"))
    
    
    #########
    ### CUT 1: cos theta(miss)
    #########
    df = df.Define("missingEnergy_rp", "FCCAnalyses::missingEnergy(240., ReconstructedParticles)")
    df = df.Define("missingEnergy", "missingEnergy_rp[0].energy")
    df = df.Define("cosTheta_miss", "FCCAnalyses::get_cosTheta_miss(missingEnergy_rp)")
    
    results.append(df.Histo1D(("missingEnergy", "", *missEnergy), "missingEnergy"))
    results.append(df.Histo1D(("cosThetaMiss_nOne", "", *bins_cosThetaMiss), "cosTheta_miss"))
    
    df = df.Filter("cosTheta_miss < 0.96")

    results.append(df.Histo1D(("cutFlow_mumu", "", *bins_count), "cut1"))
    results.append(df.Histo1D(("cutFlow_ee", "", *bins_count), "cut1"))
    results.append(df.Histo1D(("cutFlow_nunu", "", *bins_count), "cut1"))
    results.append(df.Histo1D(("cutFlow_qq", "", *bins_count), "cut1"))
    

    #########
    ### CUT 2: select Z decay product
    #########
    select_mumu = "muons_no == 2 && electrons_no == 0 && missingEnergy < 30"
    select_ee   = "muons_no == 0 && electrons_no == 2 && missingEnergy < 30"
    select_nunu = "muons_no == 0 && electrons_no == 0 && missingEnergy > 95 && missingEnergy < 115"
    select_qq = f"! ( ({select_mumu}) || ({select_ee}) || ({select_nunu}) ) && muons_no == 0 && electrons_no == 0"
    
    df_mumu   = df.Filter(select_mumu)
    df_ee     = df.Filter(select_ee)
    df_nunu   = df.Filter(select_nunu)
    df_quarks = df.Filter(select_qq)
    
    results.append(df_mumu.Histo1D(("cutFlow_mumu", "", *bins_count), "cut2"))
    results.append(df_ee.Histo1D(("cutFlow_ee", "", *bins_count), "cut2"))
    results.append(df_nunu.Histo1D(("cutFlow_nunu", "", *bins_count), "cut2"))
    results.append(df_quarks.Histo1D(("cutFlow_bb", "", *bins_count), "cut2"))
    results.append(df_quarks.Histo1D(("cutFlow_cc", "", *bins_count), "cut2"))
    results.append(df_quarks.Histo1D(("cutFlow_ss", "", *bins_count), "cut2"))
    results.append(df_quarks.Histo1D(("cutFlow_qq", "", *bins_count), "cut2"))


    #########
    ### CUT 3: we want to detect Z->mumu/ee (so we don't have Z->qq interfering with measurement of H->bb)
    ###        Z->qq and Z->nunu do not get the next few cuts
    #########

    # build the Z resonance based on the available leptons. Returns the best lepton pair compatible with the Z mass and recoil at 125 GeV
    # technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system (Z), index and 2 the leptons of the pair
    
    # muons
    df_mumu = df_mumu.Define("hbuilder_result", "FCCAnalyses::resonanceBuilder_mass_recoil(91.2, 125, 0, 240, false)(muons, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    
    df_mumu = df_mumu.Filter("hbuilder_result.size() > 0")
    
    results.append(df_mumu.Histo1D(("cutFlow_mumu", "", *bins_count), "cut3"))

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
    
    results.append(df_ee.Histo1D(("cutFlow_ee", "", *bins_count), "cut3"))

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
    df_mumu = df_mumu.Filter("zmumu_recoil_m > 122")
    results.append(df_mumu.Histo1D(("cutFlow_mumu", "", *bins_count), "cut4"))
    
    results.append(df_ee.Histo1D(("ee_recoil_m_nOne", "", *bins_m), "zee_recoil_m"))
    df_ee = df_ee.Filter("zee_recoil_m > 115 && zee_recoil_m < 135")
    results.append(df_ee.Histo1D(("cutFlow_ee", "", *bins_count), "cut4"))
    
    # graphs for inspecting the WW background
    #results.append(df_mumu.Graph("missingEnergy", "zmumu_recoil_m"))
    #results.append(df_mumu.Graph("zmumu_recoil_m", "missingEnergy"))
    

    #########
    ### CUT 5: momentum
    #########
    results.append(df_mumu.Histo1D(("mumu_p_nOne", "", *bins_p), "zmumu_p"))
    df_mumu = df_mumu.Filter("zmumu_p > 35 && zmumu_p < 60")
    results.append(df_mumu.Histo1D(("cutFlow_mumu", "", *bins_count), "cut5"))

    results.append(df_ee.Histo1D(("ee_p_nOne", "", *bins_p), "zee_p"))
    df_ee = df_ee.Filter("zee_p > 35 && zee_p < 60")
    results.append(df_ee.Histo1D(("cutFlow_ee", "", *bins_count), "cut5"))

    #########
    ### CUT 6: acolinearity
    #########
    df_mumu = df_mumu.Define("acoplanarity", "FCCAnalyses::acoplanarity(zmumu_leps)")
    df_mumu = df_mumu.Define("acolinearity", "FCCAnalyses::acolinearity(zmumu_leps)")
    results.append(df_mumu.Histo1D(("acolinearity_mumu", "", *bins_aco), "acolinearity"))

    df_mumu = df_mumu.Filter("acolinearity > 0.05")
    results.append(df_mumu.Histo1D(("cutFlow_mumu", "", *bins_count), "cut6"))


    df_ee = df_ee.Define("acoplanarity", "FCCAnalyses::acoplanarity(zee_leps)")
    df_ee = df_ee.Define("acolinearity", "FCCAnalyses::acolinearity(zee_leps)")
    results.append(df_ee.Histo1D(("acolinearity_ee", "", *bins_aco), "acolinearity"))
    
    df_ee = df_ee.Filter("acolinearity > 0.05")
    results.append(df_ee.Histo1D(("cutFlow_ee", "", *bins_count), "cut6"))


    #########
    ### CUT 7: cut on Z mass
    #########
    results.append(df_mumu.Histo1D(("zmumu_m_nOne", "", *bins_m), "zmumu_m"))
    df_mumu = df_mumu.Filter("zmumu_m > 85 && zmumu_m < 95")
    results.append(df_mumu.Histo1D(("cutFlow_mumu", "", *bins_count), "cut7"))

    results.append(df_ee.Histo1D(("zee_m_nOne", "", *bins_m), "zee_m"))
    df_ee = df_ee.Filter("zee_m > 85 && zee_m < 95")
    results.append(df_ee.Histo1D(("cutFlow_ee", "", *bins_count), "cut7"))


    # jet analysis for the case of 2 jets (Z -> leps)
    for leps, df in [("muons", df_mumu), ("electrons", df_ee), ("neutrinos", df_nunu)]:
        # define PF candidates collection by removing the leptons
        if leps != "neutrinos":
            df = df.Define("rps_no_leps", f"FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, {leps})")
        else:
            df = df.Alias("rps_no_leps", "ReconstructedParticles")
        
        # clustering and flavour tagging
        df = jet2Cluster.define(df)
        df = jet2Flavour.define_and_inference(df)
        df = df.Define("jet_tlv", "FCCAnalyses::makeLorentzVectors(jet_px, jet_py, jet_pz, jet_e)")
        
        # cut on b jet confidence
        df = df.Filter("recojet_isB[0] > 0.97 && recojet_isB[1] > 0.97")
        if leps != "neutrinos":
            results.append(df.Histo1D((f"cutFlow_{'mumu' if leps == 'muons' else 'ee'}", "", *bins_count), "cut8"))
        else:
            results.append(df.Histo1D(("cutFlow_nunu", "", *bins_count), "cut3"))
        
        # calculate dijet m and p
        df = df.Define("dijet", "jet_tlv[0] + jet_tlv[1]")
        df = df.Define("dijet_m", "dijet.M()")
        df = df.Define("dijet_p", "dijet.P()")
        
        results.append(df.Histo1D((f"h{leps}_m", "", *bins_m), "dijet_m"))
        results.append(df.Histo1D((f"h{leps}_p", "", *bins_m), "dijet_p"))
        
        # for neutrinos, cut on Higgs mass
        if leps == "neutrinos":
            df = df.Filter("dijet_m > 120 && dijet_m < 128")
            results.append(df.Histo1D(("cutFlow_nunu", "", *bins_count), "cut4"))
    
    # Z->qq analyses
    # clustering and flavour tagging
    df_quarks = jet4Cluster.define(df_quarks)
    df_quarks = jet4Flavour.define_and_inference(df_quarks)
    df_quarks = df_quarks.Define("jet_tlv", "FCCAnalyses::makeLorentzVectors(jet_px, jet_py, jet_pz, jet_e)")
    
    df_quarks = df_quarks.Define("bbidx", "FCCAnalyses::getMaxAndSecondMaxIdx(recojet_isB)")
    df_quarks = df_quarks.Define("ssidx", "FCCAnalyses::getMaxAndSecondMaxIdx(recojet_isS)")
    df_quarks = df_quarks.Define("ccidx", "FCCAnalyses::getMaxAndSecondMaxIdx(recojet_isC)")
    df_quarks = df_quarks.Define("qqidx", "FCCAnalyses::getMaxAndSecondMaxIdx(recojet_isQ)")
    
    # sort by tag
    # special case ZH->bbbb
    df_bb = df_quarks.Filter("recojet_isB[0] > 0.97 && recojet_isB[1] > 0.97 && recojet_isB[2] > 0.97 && recojet_isB[3] > 0.97")
    # make sure that there are 2 b jets
    df_quarks = df_quarks.Filter("recojet_isB[bbidx[0]] > 0.97 && recojet_isB[bbidx[1]] > 0.97")
    # filter by other jet types
    df_cc = df_quarks.Filter("recojet_isC[ccidx[0]] > 0.97 && recojet_isC[ccidx[1]] > 0.97")
    df_ss = df_quarks.Filter("recojet_isS[ssidx[0]] > 0.97 && recojet_isS[ssidx[1]] > 0.97")
    df_qq = df_quarks.Filter("recojet_isQ[qqidx[0]] > 0.97 && recojet_isQ[qqidx[1]] > 0.97")
    
    results.append(df_bb.Histo1D(("cutFlow_bb", "", *bins_count), "cut3"))
    results.append(df_cc.Histo1D(("cutFlow_cc", "", *bins_count), "cut3"))
    results.append(df_ss.Histo1D(("cutFlow_ss", "", *bins_count), "cut3"))
    results.append(df_qq.Histo1D(("cutFlow_qq", "", *bins_count), "cut3"))
    
    # Z->cc/ss/qq case
    for q, df in [("cc", df_cc), ("ss", df_ss), ("qq", df_qq)]:
        # compute Z and H masses and momenta
        df = df.Define("z_dijet", f"jet_tlv[{q}idx[0]] + jet_tlv[{q}idx[1]]")
        df = df.Define("h_dijet", f"jet_tlv[bbidx[0]] + jet_tlv[bbidx[1]]")
        
        df = df.Define("z_dijet_m", "z_dijet.M()")
        df = df.Define("z_dijet_p", "z_dijet.P()")
        df = df.Define("h_dijet_m", "h_dijet.M()")
        df = df.Define("h_dijet_p", "h_dijet.P()")
        
        results.append(df.Histo1D((f"z{q}_m", "", *bins_m), "z_dijet_m"))
        results.append(df.Histo1D((f"h{q}_m", "", *bins_m), "h_dijet_m"))
        results.append(df.Histo1D((f"z{q}_p", "", *bins_m), "z_dijet_p"))
        results.append(df.Histo1D((f"h{q}_p", "", *bins_m), "h_dijet_p"))
        
        # filter on Z mass
        df = df.Filter("z_dijet_m > 85 && z_dijet_m < 95")
        results.append(df.Histo1D((f"cutFlow_{q}", "", *bins_count), "cut4"))
        
        # filter on H mass
        df = df.Filter("h_dijet_m > 115 && h_dijet_m < 135")
        results.append(df.Histo1D((f"cutFlow_{q}", "", *bins_count), "cut5"))

    # Z->bb case
    # TODO pair jets based on momentum
    

    return results, weightsum


if __name__ == "__main__":

    datadict = functions.get_datadicts() # get default datasets

    datasets_sig = ["wzp6_ee_nunuH_Hbb_ecm240", "wzp6_ee_eeH_Hbb_ecm240", "wzp6_ee_tautauH_Hbb_ecm240", "wzp6_ee_ccH_Hbb_ecm240", "wzp6_ee_bbH_Hbb_ecm240", "wzp6_ee_qqH_Hbb_ecm240", "wzp6_ee_ssH_Hbb_ecm240", "wzp6_ee_mumuH_Hbb_ecm240"]
    datasets_bkg = ["p8_ee_WW_ecm240", "p8_ee_ZZ_ecm240"]# "wzp6_ee_mumu_ecm240", "wzp6_ee_tautau_ecm240", "wzp6_egamma_eZ_Zmumu_ecm240", "wzp6_gammae_eZ_Zmumu_ecm240", "wzp6_gaga_mumu_60_ecm240", "wzp6_gaga_tautau_60_ecm240", "wzp6_ee_nuenueZ_ecm240"]

    datasets_to_run = datasets_sig + datasets_bkg
    result = functions.build_and_run(datadict, datasets_to_run, build_graph, f"output_h_bb.root", args, norm=True, lumi=7200000)
