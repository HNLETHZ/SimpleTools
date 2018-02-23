import ROOT
import math
from array import array
from collections import OrderedDict
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi, inConeCollection, bestMatch
from tree_variables import branches # here the ntuple branches are defined
from PhysicsTools.Heppy.physicsobjects.GenParticle import GenParticle
from PhysicsTools.Heppy.physicsobjects.Muon import Muon
from PhysicsTools.Heppy.physicsobjects.Electron import Electron
from PhysicsTools.Heppy.physicsobjects.Photon import Photon
from PhysicsTools.Heppy.physicsobjects.Tau import Tau
from PhysicsTools.Heppy.physicsobjects.PhysicsObject import PhysicsObject
from utils import isAncestor, displacement2D, displacement3D, makeRecoVertex # utility functions
from samples import samples

##########################################################################################
# load custom library to ROOT. This contains the kinematic vertex fitter class
ROOT.gSystem.Load('libSimpleToolsHNLNtuple')
from ROOT import HNLKinematicVertexFitter as VertexFitter

##########################################################################################
# initialise output files to save the flat ntuples
outfile_gen = ROOT.TFile('hnl_gen_tuple.root', 'recreate')
ntuple_gen = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
tofill_gen = OrderedDict(zip(branches, [-99.]*len(branches))) # initialise all branches to unphysical -99       

##########################################################################################
# Get ahold of the events
print '... loading files'
events = Events(samples[:2]) # make sure this corresponds to your file name!
# events = Events(samples) # make sure this corresponds to your file name!
print 'done ...'
maxevents = -1 # max events to process
totevents = events.size() # total number of events in the files

##########################################################################################
# instantiate the handles to the relevant collections.
# Do this *outside* the event loop

handles = OrderedDict()
handles['taus'       ] = ( ('slimmedTaus'                  , '', 'PAT' ), Handle('std::vector<pat::Tau>'                         ) )
handles['muons'      ] = ( ('slimmedMuons'                 , '', 'PAT' ), Handle('std::vector<pat::Muon>'                        ) )
handles['electrons'  ] = ( ('slimmedElectrons'             , '', 'PAT' ), Handle('std::vector<pat::Electron>'                    ) )
handles['photons'    ] = ( ('slimmedPhotons'               , '', 'PAT' ), Handle('std::vector<pat::Photon>'                      ) )
handles['jets'       ] = ( ('slimmedJets'                  , '', 'PAT' ), Handle('std::vector<pat::Jet>'                         ) )
handles['genp_pruned'] = ( ('prunedGenParticles'           , '', 'PAT' ), Handle('std::vector<reco::GenParticle>'                ) )
handles['genp_packed'] = ( ('packedGenParticles'           , '', 'PAT' ), Handle('std::vector<pat::PackedGenParticle>'           ) )
handles['pvs'        ] = ( ('offlineSlimmedPrimaryVertices', '', 'PAT' ), Handle('std::vector<reco::Vertex>'                     ) )
handles['svs'        ] = ( ('slimmedSecondaryVertices'     , '', 'PAT' ), Handle('std::vector<reco::VertexCompositePtrCandidate>') )
handles['dsmuons'    ] = ( ('displacedStandAloneMuons'     , '', 'RECO'), Handle('std::vector<reco::Track>'                      ) )
handles['dgmuons'    ] = ( ('displacedGlobalMuons'         , '', 'RECO'), Handle('std::vector<reco::Track>'                      ) )
handles['beamspot'   ] = ( ('offlineBeamSpot'              , '', 'RECO'), Handle('reco::BeamSpot'                                ) )

##########################################################################################
# stuff I need to instantiate only once
vtxfit = VertexFitter()
# create a std::vector<RecoChargedCandidate> to be passed to the fitter 
tofit = ROOT.std.vector('reco::RecoChargedCandidate')()

##########################################################################################
# start looping on the events
for i, ev in enumerate(events):
    
    ######################################################################################
    # controls on the events being processed
    if maxevents>0 and i>maxevents:
        break
        
    if i%100==0:
        print '===> processing %d / %d event' %(i, totevents)
    
    ######################################################################################
    # load the handles
    for k, v in handles.iteritems():
        ev.getByLabel(v[0], v[1])
        setattr(ev, k, v[1].product())
       
    # loosely filter the reco taus 
    ev.taus = [tau for tau in ev.taus if tau.pt()>18.]

    # loosely filter jets
    ev.jets = [jet for jet in ev.jets if jet.pt()>18. and abs(jet.eta())<2.5]
    
    # all gen particles
    ev.genp = [ip for ip in ev.genp_pruned] + [ip for ip in ev.genp_packed]

    # get the heavy neutrino
    the_hns = [ip for ip in ev.genp_pruned if abs(ip.pdgId())==9900012 and ip.isLastCopy()]
    if len(the_hns)!=1:
        import pdb ; pdb.set_trace()

    the_hn = the_hns[0]

    # all W bosons
    try:
        the_w = [ip for ip in ev.genp_pruned if abs(ip.pdgId())==24 and ip.isLastCopy() and 9900012 in [ip.daughter(jj).pdgId() for jj in range(ip.numberOfDaughters())]][0]
        # prompt lepton
        the_pl = [the_w.daughter(jj) for jj in range(the_w.numberOfDaughters()) if abs(the_w.daughter(jj).pdgId()) in [11,13]][0]
    except:
        print 'no W, must be off-shell. Picking the prompt lepton from somewhere else'
        # prompt lepton
        the_pl = map(GenParticle, [ip for ip in ev.genp_pruned if abs(ip.pdgId()) in [11,13] and ip.isPromptFinalState() and not isAncestor(the_hn, ip)])[0]      

    # get the immediate daughters of the heavy neutrino decay
    the_hn.initialdaus = [the_hn.daughter(jj) for jj in range(the_hn.numberOfDaughters())]
    if len(the_hn.initialdaus) != 3:
        import pdb ; pdb.set_trace()

    the_hn.lep1 = max([ii for ii in the_hn.initialdaus if abs(ii.pdgId()) in [11, 13]], key = lambda x : x.pt())
    the_hn.lep2 = min([ii for ii in the_hn.initialdaus if abs(ii.pdgId()) in [11, 13]], key = lambda x : x.pt())
    the_hn.neu  =     [ii for ii in the_hn.initialdaus if abs(ii.pdgId()) in [12, 14]][0] # there can be only one
        
    # check that the two lepton come from the same vertex 
    if the_hn.lep1.vx() != the_hn.lep2.vx():
        import pdb ; pdb.set_trace()
    
    # identify the secondary vertex
    the_sv = the_hn.lep1.vertex()

    # need to analyse the lepton after they radiated / converted
    for ip in [the_hn.lep1, the_hn.lep2, the_pl]:
        finaldaus  = []
        for ipp in ev.genp_packed:
            mother = ipp.mother(0)
            if mother and isAncestor(ip, mother):
                finaldaus.append(ipp)
        ip.finaldaughters = sorted(finaldaus , key = lambda x : x.pt(), reverse = True)
        ip.hasConvOrRad = (len(ip.finaldaughters)>1)
        if len(ip.finaldaughters)>1:
            try:
                ip.finallep = max([ii for ii in ip.finaldaughters if ii.pdgId()==ip.pdgId()], key = lambda x : x.pt())
            except:
                ip.finallep = max([ii for ii in ip.finaldaughters if abs(ii.pdgId()) in [11, 13]], key = lambda x : x.pt())                
        else:
            ip.finallep = ip

    # 4-momentum of the visible part of the HN
    the_hn.vishn = the_hn.lep1.finallep.p4() + the_hn.lep2.finallep.p4()

    # map our objects to convenient Heppy objects
    ev.electrons = map(Electron     , ev.electrons)
    ev.photons   = map(Photon       , ev.photons  )
    ev.muons     = map(Muon         , ev.muons    )
    ev.taus      = map(Tau          , ev.taus     )
    ev.dsmuons   = map(PhysicsObject, ev.dsmuons  )
    ev.dgmuons   = map(PhysicsObject, ev.dgmuons  )

    # impose the muon PDG ID to the displaced objects, that otherwise carry none
    for mm in ev.dsmuons + ev.dgmuons:
        mm.mass   = lambda : 0.10565837
        mm.pdgId  = lambda : -(mm.charge()*13)

    # also append a TrackRef, this will be needed for the refitting 
    for jj, mm in enumerate(ev.dsmuons):
        mm.track = lambda : ROOT.reco.TrackRef(handles['dsmuons'][1].product(), jj)

    for jj, mm in enumerate(ev.dgmuons):
        mm.track = lambda : ROOT.reco.TrackRef(handles['dgmuons'][1].product(), jj)
    
    # all matchable objects
    # matchable = ev.electrons + ev.photons + ev.muons + ev.taus + ev.dsmuons + ev.dgmuons 
    matchable = ev.electrons + ev.photons + ev.muons + ev.dsmuons + ev.dgmuons # better not to use taus for the time being
            
    # match gen to reco
    for ip in [the_hn.lep1.finallep, the_hn.lep2.finallep, the_pl.finallep]:
        ip.bestmatch     = None
        ip.bestmatchtype = None
        ip.matches = inConeCollection(ip, matchable, 0.2, 0.)
        # to find the best match, give precedence to any matched 
        # particle in the matching cone with the correct PDG ID
        # then to the one which is closest
        ip.matches.sort(key = lambda x : (x.pdgId()==ip.pdgId(), -deltaR(x, ip)), reverse = True )
        if len(ip.matches):
            ip.bestmatch = ip.matches[0]
            # remove already matched particles, avoid multiple matches to the same candidate
            matchable.remove(ip.bestmatch)
            # record which is which
            if ip.bestmatch in ev.electrons: ip.bestmatchtype = 0
            if ip.bestmatch in ev.photons  : ip.bestmatchtype = 1
            if ip.bestmatch in ev.muons    : ip.bestmatchtype = 2
            if ip.bestmatch in ev.taus     : ip.bestmatchtype = 3
            if ip.bestmatch in ev.dsmuons  : ip.bestmatchtype = 4
            if ip.bestmatch in ev.dgmuons  : ip.bestmatchtype = 5
    
    # clear it before doing it again
    ev.recoSv = None

    # let's refit the secondary vertex, IF both leptons match to some reco particle
    if not(the_hn.lep1.finallep.bestmatch is None or the_hn.lep2.finallep.bestmatch is None):
        # clear the vector
        tofit.clear()
        # create a RecoChargedCandidate for each reconstructed lepton and flush it into the vector
        for il in [the_hn.lep1.finallep.bestmatch, the_hn.lep2.finallep.bestmatch]:
            # if the reco particle is a displaced thing, it does not have the p4() method, so let's build it 
            myp4 = ROOT.Math.LorentzVector('<ROOT::Math::PxPyPzE4D<double> >')(il.px(), il.py(), il.pz(), math.sqrt(il.mass()**2 + il.px()**2 + il.py()**2 + il.pz()**2))
            ic = ROOT.reco.RecoChargedCandidate() # instantiate a dummy RecoChargedCandidate
            ic.setCharge(il.charge())             # assign the correct charge
            ic.setP4(myp4)                        # assign the correct p4
            ic.setTrack(il.track())               # set the correct TrackRef
            if ic.track().isNonnull():            # check that the track is valid, there are photons around too!
                tofit.push_back(ic)

        # further sanity check: two *distinct* tracks
        if tofit.size()==2 and tofit[0].track() != tofit[1].track():
            # fit it!
            svtree = vtxfit.Fit(tofit) # actual vertex fitting
            # check that the vertex is good
            if not svtree.get().isEmpty() and svtree.get().isValid():
                svtree.movePointerToTheTop()
                sv = svtree.currentDecayVertex().get()
                ev.recoSv = makeRecoVertex(sv, kinVtxTrkSize=2) # need to do some gymastics
        
#     if ev.eventAuxiliary().event() in [159, 163, 160, 161, 166, 165, 157, 162, 169, 170, 168, 172, 180, 176, 178, 181, 182, 188, 185, 190, 191, 186, 192]:
#         if the_hn.lep1.finallep.bestmatch is None:
#             import pdb ; pdb.set_trace()
    
    ######################################################################################
    # fill the ntuple: each gen tau makes an entry
    for k, v in tofill_gen.iteritems(): tofill_gen[k] = -99. # initialise before filling
    
    tofill_gen['run'         ] = ev.eventAuxiliary().run()
    tofill_gen['lumi'        ] = ev.eventAuxiliary().luminosityBlock()
    tofill_gen['event'       ] = ev.eventAuxiliary().event()
    tofill_gen['nvtx'        ] = len(ev.pvs)

    tofill_gen['in_acc'      ] = abs(the_hn.lep1.finallep.eta())<2.5 and \
                                 abs(the_hn.lep1.finallep.eta())<2.5 and \
                                 abs(the_pl.finallep     .eta())<2.5
                                
    tofill_gen['hn_mass'     ] = the_hn.mass()
    tofill_gen['hn_pt'       ] = the_hn.pt()
    tofill_gen['hn_eta'      ] = the_hn.eta()
    tofill_gen['hn_phi'      ] = the_hn.phi()
    tofill_gen['hn_q'        ] = the_hn.charge()
    tofill_gen['hn_vis_mass' ] = the_hn.vishn.mass()
    tofill_gen['hn_vis_pt'   ] = the_hn.vishn.pt()
    tofill_gen['hn_vis_eta'  ] = the_hn.vishn.eta()
    tofill_gen['hn_vis_phi'  ] = the_hn.vishn.phi()
  
    tofill_gen['hn_2d_disp'  ] = displacement2D(the_hn.lep1, the_hn)
    tofill_gen['hn_3d_disp'  ] = displacement3D(the_hn.lep1, the_hn)
  
    tofill_gen['l0_mass'     ] = the_pl.mass()
    tofill_gen['l0_pt'       ] = the_pl.finallep.pt()
    tofill_gen['l0_eta'      ] = the_pl.finallep.eta()
    tofill_gen['l0_phi'      ] = the_pl.finallep.phi()
    tofill_gen['l0_q'        ] = the_pl.charge()
    tofill_gen['l0_pdgid'    ] = the_pl.pdgId()
    tofill_gen['l0_conv_rad' ] = the_pl.hasConvOrRad
    tofill_gen['l0_pt_loss'  ] = the_pl.finallep.pt() / the_pl.pt()
    
    if the_pl.finallep.bestmatch:
        tofill_gen['l0_reco_mass'    ] = the_pl.finallep.bestmatch.mass()
        tofill_gen['l0_reco_pt'      ] = the_pl.finallep.bestmatch.pt()
        tofill_gen['l0_reco_eta'     ] = the_pl.finallep.bestmatch.eta()
        tofill_gen['l0_reco_phi'     ] = the_pl.finallep.bestmatch.phi()
        tofill_gen['l0_reco_q'       ] = the_pl.finallep.bestmatch.charge()
        tofill_gen['l0_reco_pdgid'   ] = the_pl.finallep.bestmatch.pdgId()
        tofill_gen['l0_reco_type'    ] = the_pl.finallep.bestmatchtype
        tofill_gen['l0_reco_vx'      ] = the_pl.finallep.bestmatch.vx()
        tofill_gen['l0_reco_vy'      ] = the_pl.finallep.bestmatch.vy()
        tofill_gen['l0_reco_vz'      ] = the_pl.finallep.bestmatch.vz()
        
    tofill_gen['l1_mass'     ] = the_hn.lep1.mass()
    tofill_gen['l1_pt'       ] = the_hn.lep1.finallep.pt()
    tofill_gen['l1_eta'      ] = the_hn.lep1.finallep.eta()
    tofill_gen['l1_phi'      ] = the_hn.lep1.finallep.phi()
    tofill_gen['l1_q'        ] = the_hn.lep1.charge()
    tofill_gen['l1_pdgid'    ] = the_hn.lep1.pdgId()
    tofill_gen['l1_conv_rad' ] = the_hn.lep1.hasConvOrRad
    tofill_gen['l1_pt_loss'  ] = the_hn.lep1.finallep.pt() / the_hn.lep1.pt()

    if the_hn.lep1.finallep.bestmatch:
        tofill_gen['l1_reco_mass'    ] = the_hn.lep1.finallep.bestmatch.mass()
        tofill_gen['l1_reco_pt'      ] = the_hn.lep1.finallep.bestmatch.pt()
        tofill_gen['l1_reco_eta'     ] = the_hn.lep1.finallep.bestmatch.eta()
        tofill_gen['l1_reco_phi'     ] = the_hn.lep1.finallep.bestmatch.phi()
        tofill_gen['l1_reco_q'       ] = the_hn.lep1.finallep.bestmatch.charge()
        tofill_gen['l1_reco_pdgid'   ] = the_hn.lep1.finallep.bestmatch.pdgId()
        tofill_gen['l1_reco_type'    ] = the_hn.lep1.finallep.bestmatchtype
        tofill_gen['l1_reco_vx'      ] = the_hn.lep1.finallep.bestmatch.vx()
        tofill_gen['l1_reco_vy'      ] = the_hn.lep1.finallep.bestmatch.vy()
        tofill_gen['l1_reco_vz'      ] = the_hn.lep1.finallep.bestmatch.vz()
  
    tofill_gen['l2_mass'     ] = the_hn.lep2.mass()
    tofill_gen['l2_pt'       ] = the_hn.lep2.finallep.pt()
    tofill_gen['l2_eta'      ] = the_hn.lep2.finallep.eta()
    tofill_gen['l2_phi'      ] = the_hn.lep2.finallep.phi()
    tofill_gen['l2_q'        ] = the_hn.lep2.charge()
    tofill_gen['l2_pdgid'    ] = the_hn.lep2.pdgId()
    tofill_gen['l2_conv_rad' ] = the_hn.lep2.hasConvOrRad
    tofill_gen['l2_pt_loss'  ] = the_hn.lep2.finallep.pt() / the_hn.lep2.pt()

    if the_hn.lep2.finallep.bestmatch:
        tofill_gen['l2_reco_mass'    ] = the_hn.lep2.finallep.bestmatch.mass()
        tofill_gen['l2_reco_pt'      ] = the_hn.lep2.finallep.bestmatch.pt()
        tofill_gen['l2_reco_eta'     ] = the_hn.lep2.finallep.bestmatch.eta()
        tofill_gen['l2_reco_phi'     ] = the_hn.lep2.finallep.bestmatch.phi()
        tofill_gen['l2_reco_q'       ] = the_hn.lep2.finallep.bestmatch.charge()
        tofill_gen['l2_reco_pdgid'   ] = the_hn.lep2.finallep.bestmatch.pdgId()
        tofill_gen['l2_reco_type'    ] = the_hn.lep2.finallep.bestmatchtype
        tofill_gen['l2_reco_vx'      ] = the_hn.lep2.finallep.bestmatch.vx()
        tofill_gen['l2_reco_vy'      ] = the_hn.lep2.finallep.bestmatch.vy()
        tofill_gen['l2_reco_vz'      ] = the_hn.lep2.finallep.bestmatch.vz()
   
    tofill_gen['nu_mass'     ] = the_hn.neu.mass()
    tofill_gen['nu_pt'       ] = the_hn.neu.pt()
    tofill_gen['nu_eta'      ] = the_hn.neu.eta()
    tofill_gen['nu_phi'      ] = the_hn.neu.phi()
    tofill_gen['nu_q'        ] = the_hn.neu.charge()
    tofill_gen['nu_pdgid'    ] = the_hn.neu.pdgId()

    tofill_gen['pv_x'        ] = the_hn.vx()
    tofill_gen['pv_y'        ] = the_hn.vy()
    tofill_gen['pv_z'        ] = the_hn.vz()

    tofill_gen['sv_x'        ] = the_sv.x()
    tofill_gen['sv_y'        ] = the_sv.y()
    tofill_gen['sv_z'        ] = the_sv.z()

    if len(ev.pvs)>0:
        tofill_gen['pv_reco_x'] = ev.pvs[0].x()
        tofill_gen['pv_reco_y'] = ev.pvs[0].y()
        tofill_gen['pv_reco_z'] = ev.pvs[0].z()

    if hasattr(ev, 'recoSv') and ev.recoSv:
        tofill_gen['sv_reco_x'   ] = ev.recoSv.x()
        tofill_gen['sv_reco_y'   ] = ev.recoSv.y()
        tofill_gen['sv_reco_z'   ] = ev.recoSv.z()
        tofill_gen['sv_reco_prob'] = ROOT.TMath.Prob(ev.recoSv.chi2(), int(ev.recoSv.ndof()))

        # for the 2D displacement, take the BS coordinates which are more accurate
        tofill_gen['hn_reco_2d_disp'] = math.sqrt( (ev.recoSv.x() - ev.beamspot.x(ev.pvs[0].z()))**2 +\
                                                   (ev.recoSv.y() - ev.beamspot.y(ev.pvs[0].z()))**2  )
        tofill_gen['hn_reco_3d_disp'] = math.sqrt( (ev.recoSv.x() - ev.beamspot.x(ev.pvs[0].z()))**2 +\
                                                   (ev.recoSv.y() - ev.beamspot.y(ev.pvs[0].z()))**2 +\
                                                   (ev.recoSv.z() - ev.pvs[0].z()               )**2  )

    ntuple_gen.Fill(array('f',tofill_gen.values()))

##########################################################################################
# write the ntuples and close the files
outfile_gen.cd()
ntuple_gen.Write()
outfile_gen.Close()
