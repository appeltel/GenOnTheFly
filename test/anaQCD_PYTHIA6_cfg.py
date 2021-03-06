import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

# ========= Import of Standard Configurations ==============

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")

# ================= Input Variable Parsing ======================

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('standard')

options.output = 'QCDAna.root'
options.maxEvents = 100

options.register('tune',
                 "Z2",
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Pythia tune")

options.register('processType',
                 "NSD",
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Pythia process type with pT_hat range")

options.register('sqrtS',
                 2760.0,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.float,
                 "Center-of-mass energy")

options.register('ptHatLow',
                 120,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Minimum pt-hat")

options.register('ptHatHigh',
                 160,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.int,
                 "Maximum pt-hat")

options.parseArguments()


# ===================== General Settings ===================

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    makeTriggerResults = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.output)
)

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()

# ============= PYTHIA Generator =====================

# ============= Pythia setting  ================================
if options.tune == 'Z2':
  from Configuration.Generator.PythiaUEZ2Settings_cfi import *
if options.tune == 'Z2star':
  from Configuration.Generator.PythiaUEZ2starSettings_cfi import *
if options.tune == 'AMBT2':
  from PythiaUEAMBT2Settings_cfi import *

process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    PNInitialState = cms.untracked.bool(True),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    comEnergy = cms.double(5020.0),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
         pythiaUESettingsBlock,
         processParameters = cms.vstring(
                 'MSEL=1         ! High Pt QCD'
                 ),
         parameterSets = cms.vstring('pythiaUESettings',
                                     'processParameters',
                                     #'pthatLow',
                                     #'pthatHigh')
                                    )
         )
)

from customisePYTHIA6_cfi import *
updatePy6ProcParameters(process.generator,options.processType,options.ptHatLow,options.ptHatHigh,options.sqrtS)

print process.generator.PythiaParameters.processParameters


# ============= Gen Jets ================================

process.ak3GenJets = process.ak5GenJets.clone( rParam = 0.3 )
process.ak4GenJets = process.ak5GenJets.clone( rParam = 0.4 )
process.ak7GenJets = process.ak5GenJets.clone( rParam = 0.7 )


# ========== Flavor Matching ==========================

process.myPartons = cms.EDProducer("PartonSelector",
    withLeptons = cms.bool(False),
    src = cms.InputTag("genParticles")
)

process.flavourByRef = cms.EDProducer("JetPartonMatcher",
    jets = cms.InputTag("ak3GenJets"),
    coneSizeToAssociate = cms.double(0.3),
    partons = cms.InputTag("myPartons")
)


# =============== Analysis =============================

process.qcdAna = cms.EDAnalyzer('QCDAnalyzer',
    genJetSrc = cms.InputTag("ak3GenJets"),
    genParticleSrc = cms.InputTag("genParticles"),
    doFlavor = cms.bool(False),
    doSpecies = cms.bool(False),
    onlyDS = cms.bool(False),
    flavorSrc = cms.InputTag("flavourByRef"),
    flavorId = cms.vint32(21), # gluons
    speciesId = cms.vint32(11), # gamma
    useRapidity = cms.bool(False),
    jetEtaMin = cms.double(-0.5),
    jetEtaMax = cms.double(0.5),
    hEtaMin = cms.double(-1.0),
    hEtaMax = cms.double(1.0),
    jetRadius = cms.double(0.3),
    pthatMin = cms.double(options.ptHatLow),
    pthatMax = cms.double(options.ptHatHigh),
    jetPtBins = cms.vdouble(), 
    hPtBins = cms.vdouble(), 
    qScalePtBins = cms.vdouble(), 
    etaBins = cms.vdouble() 
)

# hadron pT bins from 0 to 200 in increments of 2
# jet pT bins from 0 to 1000 in increments of 10
# qscale pT bins from 0 to 1000 in increments of 10
for x in range(0,101):
  process.qcdAna.hPtBins += [ float(x)*2.0 ]
  process.qcdAna.jetPtBins += [ float(x)*10.0 ]
  process.qcdAna.qScalePtBins += [ float(x)*10.0 ]

#eta bins from -2.4 to 2.4 in increments of 0.4
for x in range(0,13):
  process.qcdAna.etaBins += [ float(x)*0.4-2.4 ]

longLivedBaryonIDs = cms.vint32( 2212, -2212, # p (pbar)
                                 2112, -2112, # n
                                 3122, -3122, # lambda
                                 3222, -3222, # sigma+
                                 3112, -3112, # sigma-
                                 3312, -3312, # xi
                                 3334, -3334  # omega
)

longLivedMesonIDs = cms.vint32( 211, -211, # pion
                                130, 310, 321, -321, # K
)

quarks = cms.vint32( -6,-5,-4,-3,-2,-1,1,2,3,4,5,6 )

process.qcdAna_baryons = process.qcdAna.clone(
    doSpecies = cms.bool(True),
    speciesId = longLivedBaryonIDs
)

process.qcdAna_mesons = process.qcdAna.clone(
    doSpecies = cms.bool(True),
    speciesId = longLivedMesonIDs
)                             

process.qcdAna_qJets = process.qcdAna.clone( 
    doFlavor = cms.bool(True),
    flavorId = quarks
)

process.qcdAna_qJets_baryons = process.qcdAna_baryons.clone(
    doFlavor = cms.bool(True),
    flavorId = quarks
)

process.qcdAna_qJets_mesons = process.qcdAna_mesons.clone(
    doFlavor = cms.bool(True),
    flavorId = quarks
)

process.qcdAna_gJets = process.qcdAna.clone(
    doFlavor = cms.bool(True),
    flavorId = cms.vint32(21) 
)

process.qcdAna_gJets_baryons = process.qcdAna_baryons.clone(
    doFlavor = cms.bool(True),
    flavorId = cms.vint32(21) 
)

process.qcdAna_gJets_mesons = process.qcdAna_mesons.clone(
    doFlavor = cms.bool(True),
    flavorId = cms.vint32(21) 
)


# ================== Paths and Schedule ===================

process.gen_step = cms.Path(process.generator * process.genParticles )
process.genjet_step = cms.Path(process.genJetParticles 
                               * process.ak3GenJets
                          #     * process.ak4GenJets
                          #     * process.ak5GenJets
                          #     * process.ak7GenJets
)
process.flavor_step = cms.Path( process.myPartons * process.flavourByRef )
process.ana_step = cms.Path(
    process.qcdAna *
    process.qcdAna_baryons *
    process.qcdAna_mesons *
    process.qcdAna_qJets *
    process.qcdAna_qJets_baryons *
    process.qcdAna_qJets_mesons *
    process.qcdAna_gJets *
    process.qcdAna_gJets_baryons *
    process.qcdAna_gJets_mesons 
)

process.schedule = cms.Schedule(process.gen_step,
                                process.genjet_step,
                                process.flavor_step,
                                process.ana_step)

