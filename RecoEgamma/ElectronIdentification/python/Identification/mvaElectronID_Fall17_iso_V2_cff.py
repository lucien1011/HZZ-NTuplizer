import FWCore.ParameterSet.Config as cms

# Documentation of the MVA
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentificationRun2
# https://rembserj.web.cern.ch/rembserj/notes/Electron_MVA_ID_2017_documentation

#
# In this file we define the locations of the MVA weights, cuts on the MVA values
# for specific working points, and configure those cuts in VID
#

# This MVA implementation class name
mvaFall17ClassName = "ElectronMVAEstimatorRun2Fall17Iso"
# The tag is an extra string attached to the names of the products
# such as ValueMaps that needs to distinguish cases when the same MVA estimator
# class is used with different tuning/weights
mvaTag = "V2"

# The parameters according to which the training bins are split:
ptSplit = 10.      # we have above and below 10 GeV categories
ebSplit = 0.800    # barrel is split into two regions
ebeeSplit = 1.479  # division between barrel and endcap

# There are 6 categories in this MVA. They have to be configured in this strict order
# (cuts and weight files order):
#   0   EB1 (eta<0.8)  pt 5-10 GeV     |   pt < ptSplit && |eta| < ebSplit
#   1   EB2 (eta>=0.8) pt 5-10 GeV     |   pt < ptSplit && |eta| >= ebSplit && |eta| < ebeeSplit
#   2   EE             pt 5-10 GeV     |   pt < ptSplit && |eta| >= ebeeSplit
#   3   EB1 (eta<0.8)  pt 10-inf GeV   |   pt >= ptSplit && |eta| < ebSplit
#   4   EB2 (eta>=0.8) pt 10-inf GeV   |   pt >= ptSplit && |eta| >= ebSplit && |eta| < ebeeSplit
#   5   EE             pt 10-inf GeV   |   pt >= ptSplit && |eta| >= ebeeSplit


mvaFall17WeightFiles_V2 = cms.vstring(
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_iso_V2_EB1_5.weights.xml.gz",
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_iso_V2_EB2_5.weights.xml.gz",
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_iso_V2_EE_5.weights.xml.gz",
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_iso_V2_EB1_10.weights.xml.gz",
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_iso_V2_EB2_10.weights.xml.gz",
    "RecoEgamma/ElectronIdentification/data/Fall17/electronID_mva_Fall17_iso_V2_EE_10.weights.xml.gz"
    )

# Load some common definitions for MVA machinery
from RecoEgamma.ElectronIdentification.Identification.mvaElectronID_tools \
    import (EleMVA_WP,
            configureVIDMVAEleID_V1)

# The locatoins of value maps with the actual MVA values and categories
# for all particles.
# The names for the maps are "<module name>:<MVA class name>Values"
# and "<module name>:<MVA class name>Categories"
mvaProducerModuleLabel = "electronMVAValueMapProducer"
mvaValueMapName        = mvaProducerModuleLabel + ":" + mvaFall17ClassName + mvaTag + "Values"
mvaCategoriesMapName   = mvaProducerModuleLabel + ":" + mvaFall17ClassName + mvaTag + "Categories"

## The working point for this MVA that is expected to have about 90% signal
# WP tuned to give about 90 and 80% signal efficiecny for electrons from Drell-Yan with pT > 25 GeV
# The working point for the low pt categories is just taken over from the high pt
idName90 = "mvaEleID-Fall17-iso-V2-wp90"
MVA_WP90 = EleMVA_WP(
    idName = idName90,
    mvaValueMapName = mvaValueMapName,           # map with MVA values for all particles
    mvaCategoriesMapName = mvaCategoriesMapName, # map with category index for all particles
    cutCategory0_C0 = 0.9296353578292107,  # EB1 low pt
    cutCategory0_C1 = 2.70587899253673,
    cutCategory0_C2 = 0.9350596744699461,
    cutCategory1_C0 = 0.8802860483482772,  # EB2 low pt
    cutCategory1_C1 = 2.0294814189708106,
    cutCategory1_C2 = 0.7164022704821694,
    cutCategory2_C0 = 0.7723145105695042,  # EE low pt
    cutCategory2_C1 = 1.2376055989522428,
    cutCategory2_C2 = 19.494936881558534,
    cutCategory3_C0 = 0.9679627068677376,  # EB1
    cutCategory3_C1 = 9.017734241637768,
    cutCategory3_C2 = 2.2733580727642253,
    cutCategory4_C0 = 0.9386903800951701,  # EB2
    cutCategory4_C1 = 8.978063234072666,
    cutCategory4_C2 = 2.7533412287363532,
    cutCategory5_C0 = 0.8855859008355582,  # EE
    cutCategory5_C1 = 10.968321599721932,
    cutCategory5_C2 = 4.049909321201802
)

idName80 = "mvaEleID-Fall17-iso-V2-wp80"
MVA_WP80 = EleMVA_WP(
    idName = idName80,
    mvaValueMapName = mvaValueMapName,           # map with MVA values for all particles
    mvaCategoriesMapName = mvaCategoriesMapName, # map with category index for all particles
    cutCategory0_C0 = 0.9664176694347087, # EB1 low pt
    cutCategory0_C1 = 2.942813218971673,
    cutCategory0_C2 = 0.27802152770462213,
    cutCategory1_C0 = 0.9466283733433212, # EB2 low pt
    cutCategory1_C1 = 2.020856368381169,
    cutCategory1_C2 = 0.3440849560318207,
    cutCategory2_C0 = 0.9175812962031529, # EE low pt
    cutCategory2_C1 = 1.109971297629491,
    cutCategory2_C2 = 11.76158799411055,
    cutCategory3_C0 = 0.987447686987344,  # EB1
    cutCategory3_C1 = 10.17080888826826,
    cutCategory3_C2 = 0.501698069257147,
    cutCategory4_C0 = 0.9773274263001541, # EB2
    cutCategory4_C1 = 9.112455878045488,
    cutCategory4_C2 = 0.9430928489900158,
    cutCategory5_C0 = 0.9539002157907921, # EE
    cutCategory5_C1 = 8.716459513452563,
    cutCategory5_C2 = 2.728601529879707
)

### WP tuned for HZZ analysis with very high efficiency (about 98%)
# The working points were found by requiring the same signal efficiencies in
# each category as for the Spring 16 HZZ ID
# (see RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Spring16_HZZ_V2_cff.py)
idNamewpLoose = "mvaEleID-Fall17-iso-V2-wpLoose"
MVA_WPLoose = EleMVA_WP(
    idName = idNamewpLoose,
    mvaValueMapName = mvaValueMapName,           # map with MVA values for all particles
    mvaCategoriesMapName = mvaCategoriesMapName, # map with category index for all particles
    cutCategory0 = -0.20915430749445943, # EB1 low pt
    cutCategory1 = -0.4157355054931669,  # EB2 low pt
    cutCategory2 = -0.17859352060344227, # EE low pt
    cutCategory3 = -0.852951639116353,   # EB1
    cutCategory4 = -0.7960197700002585,  # EB2
    cutCategory5 = -0.7450408698516471   # EE
    )

### WP tuned for HZZ analysis to match the Spring16 ID efficiency plus iso cut of < 0.35
# The working points were found by requiring the same signal efficiencies in
# each category as for the Spring 16 HZZ ID
# (see RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Spring16_HZZ_V1_cff.py)
idNamewpHZZ = "mvaEleID-Fall17-iso-V2-wpHZZ"
MVA_WPHZZ = EleMVA_WP(
    idName = idNamewpHZZ,
    mvaValueMapName = mvaValueMapName,           # map with MVA values for all particles
    mvaCategoriesMapName = mvaCategoriesMapName, # map with category index for all particles
    cutCategory0 = 0.5739521065342641,    # EB1 low pt
    cutCategory1 = 0.5504628790992929,   # EB2 low pt
    cutCategory2 = 0.5924627534389098,   # EE low pt
    cutCategory3 = -0.03391387993354392,  # EB1
    cutCategory4 = -0.018451958064666783,  # EB2
    cutCategory5 = -0.38565459150737535   # EE
    )

#
# Configure variable names and the values they are clipped to.
# They have to appear in the same order as in the weights xml file
#

#                Name  |  Lower clip value  | upper clip value
variablesInfo = [
                 ("ele_oldsigmaietaieta"              ,  None, None),
                 ("ele_oldsigmaiphiiphi"              ,  None, None),
                 ("ele_oldcircularity"                ,   -1.,   2.),
                 ("ele_oldr9"                         ,  None,   5.),
                 ("ele_scletawidth"                   ,  None, None),
                 ("ele_sclphiwidth"                   ,  None, None),
                 ("ele_oldhe"                         ,  None, None),
                 ("ele_kfhits"                        ,  None, None),
                 ("ele_kfchi2"                        ,  None,  10.),
                 ("ele_gsfchi2"                       ,  None, 200.),
                 ("ele_fbrem"                         ,   -1., None),
                 ("ele_gsfhits"                       ,  None, None),
                 ("ele_expected_inner_hits"           ,  None, None),
                 ("ele_conversionVertexFitProbability",  None, None),
                 ("ele_ep"                            ,  None,  20.),
                 ("ele_eelepout"                      ,  None,  20.),
                 ("ele_IoEmIop"                       ,  None, None),
                 ("ele_deltaetain"                    , -0.06, 0.06),
                 ("ele_deltaphiin"                    ,  -0.6,  0.6),
                 ("ele_deltaetaseed"                  ,  -0.2,  0.2),
                 ("ele_pfPhotonIso"                   ,  None, None), #
                 ("ele_pfChargedHadIso"               ,  None, None), # PF isolations
                 ("ele_pfNeutralHadIso"               ,  None, None), #
                 ("rho"                               ,  None, None),
                 ("ele_psEoverEraw"                   ,  None, None), # EE only
                ]

varNames, clipLower, clipUpper = [list(l) for l in zip(*variablesInfo)]
for i, x in enumerate(clipLower):
    if x == None:
        clipLower[i] = -float('Inf')
for i, x in enumerate(clipUpper):
    if x == None:
        clipUpper[i] =  float('Inf')

#
# Finally, set up VID configuration for all cuts
#

# Create the PSet that will be fed to the MVA value map producer
mvaEleID_Fall17_iso_V2_producer_config = cms.PSet(
    mvaName            = cms.string(mvaFall17ClassName),
    mvaTag             = cms.string(mvaTag),
    # This MVA uses conversion info, so configure several data items on that
    beamSpot           = cms.InputTag('offlineBeamSpot'),
    conversionsAOD     = cms.InputTag('allConversions'),
    conversionsMiniAOD = cms.InputTag('reducedEgamma:reducedConversions'),
    # Category split parameters
    ptSplit            = cms.double(ptSplit),
    ebSplit            = cms.double(ebSplit),
    ebeeSplit          = cms.double(ebeeSplit),
    # Variable clipping parameters
    varNames           = cms.vstring(*varNames),
    clipLower          = cms.vdouble(*clipLower),
    clipUpper          = cms.vdouble(*clipUpper),
    #
    weightFileNames    = mvaFall17WeightFiles_V2
    )
# Create the VPset's for VID cuts
mvaEleID_Fall17_V2_wpLoose = configureVIDMVAEleID_V1( MVA_WPLoose )
mvaEleID_Fall17_V2_wpHZZ = configureVIDMVAEleID_V1( MVA_WPHZZ )
mvaEleID_Fall17_V2_wp90 = configureVIDMVAEleID_V1( MVA_WP90, cutName="GsfEleMVAExpoScalingCut")
mvaEleID_Fall17_V2_wp80 = configureVIDMVAEleID_V1( MVA_WP80, cutName="GsfEleMVAExpoScalingCut")

mvaEleID_Fall17_V2_wpLoose.isPOGApproved = cms.untracked.bool(True)
mvaEleID_Fall17_V2_wp90.isPOGApproved = cms.untracked.bool(True)
mvaEleID_Fall17_V2_wp80.isPOGApproved = cms.untracked.bool(True)
