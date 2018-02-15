import ROOT
from array import array
from collections import OrderedDict
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from PhysicsTools.Heppy.physicsutils.TauDecayModes import tauDecayModes
from tree_variables import branches # here the ntuple branches are defined
from utils import isAncestor, displacement2D, displacement3D # utility functions
from samples import samples

##########################################################################################
# initialise output files to save the flat ntuples
outfile_gen = ROOT.TFile('hnl_gen_tuple.root', 'recreate')
ntuple_gen = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
tofill_gen = OrderedDict(zip(branches, [-99.]*len(branches))) # initialise all branches to unphysical -99       

##########################################################################################
# Get ahold of the events
print '... loading files'
events = Events(samples[:6]) # make sure this corresponds to your file name!
print 'done ...'
maxevents = -1 # max events to process
totevents = events.size() # total number of events in the files

##########################################################################################
# instantiate the handles to the relevant collections.
# Do this *outside* the event loop

handles = OrderedDict()
handles['taus'       ] = ( ('slimmedTaus'                  , '', 'PAT'), Handle('std::vector<pat::Tau>'                         ) )
handles['muons'      ] = ( ('slimmedMuons'                 , '', 'PAT'), Handle('std::vector<pat::Muon>'                        ) )
handles['electrons'  ] = ( ('slimmedElectrons'             , '', 'PAT'), Handle('std::vector<pat::Electron>'                    ) )
handles['jets'       ] = ( ('slimmedJets'                  , '', 'PAT'), Handle('std::vector<pat::Jet>'                         ) )
handles['genp_pruned'] = ( ('prunedGenParticles'           , '', 'PAT'), Handle('std::vector<reco::GenParticle>'                ) )
handles['genp_packed'] = ( ('packedGenParticles'           , '', 'PAT'), Handle('std::vector<pat::PackedGenParticle>'           ) )
handles['pvs'        ] = ( ('offlineSlimmedPrimaryVertices', '', 'PAT'), Handle('std::vector<reco::Vertex>'                     ) )
handles['svs'        ] = ( ('slimmedSecondaryVertices'     , '', 'PAT'), Handle('std::vector<reco::VertexCompositePtrCandidate>') )

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
    except:
        print 'no W, must be off-shell. Passing on'
        continue
#         import pdb ; pdb.set_trace()

    # prompt lepton
    the_pl = [the_w.daughter(jj) for jj in range(the_w.numberOfDaughters())][0]

    # get the immediate daughters of the heavy neutrino decay
    the_hn.initialdaus = [the_hn.daughter(jj) for jj in range(the_hn.numberOfDaughters())]
    if len(the_hn.initialdaus) != 3:
        import pdb ; pdb.set_trace()

    the_hn.lep1 = max([ii for ii in the_hn.initialdaus if abs(ii.pdgId()) in [11, 13]], key = lambda x : x.pt())
    the_hn.lep2 = min([ii for ii in the_hn.initialdaus if abs(ii.pdgId()) in [11, 13]], key = lambda x : x.pt())
    the_hn.neu  =     [ii for ii in the_hn.initialdaus if abs(ii.pdgId()) in [12, 14]][0] # there can be only one
    
    if the_hn.lep1.vx() != the_hn.lep2.vx():
        import pdb ; pdb.set_trace()

    # need to analyse the lepton after they radiated / converted
    for ip in [the_hn.lep1, the_hn.lep2]:
        finaldaus  = []
        for ipp in ev.genp_packed:
            mother = ipp.mother(0)
            if mother and isAncestor(ip, mother):
                finaldaus.append(ipp)
        ip.finaldaughters = sorted(finaldaus , key = lambda x : x.pt(), reverse = True)
        ip.hasConvOrRad = (len(ip.finaldaughters)>0)
        
#     import pdb ; pdb.set_trace()

    ######################################################################################
    # fill the ntuple: each gen tau makes an entry
    for k, v in tofill_gen.iteritems(): tofill_gen[k] = -99. # initialise before filling
    tofill_gen['run'       ] = ev.eventAuxiliary().run()
    tofill_gen['lumi'      ] = ev.eventAuxiliary().luminosityBlock()
    tofill_gen['event'     ] = ev.eventAuxiliary().event()
    tofill_gen['nvtx'      ] = len(ev.pvs)
    tofill_gen['hn_mass'   ] = the_hn.mass()
    tofill_gen['hn_pt'     ] = the_hn.pt()
    tofill_gen['hn_eta'    ] = the_hn.eta()
    tofill_gen['hn_phi'    ] = the_hn.phi()
    tofill_gen['hn_q'      ] = the_hn.charge()

    tofill_gen['hn_2d_disp'] = displacement2D(the_hn.lep1, the_w)
    tofill_gen['hn_3d_disp'] = displacement3D(the_hn.lep1, the_w)

    tofill_gen['l0_mass'   ] = the_pl.mass()
    tofill_gen['l0_pt'     ] = the_pl.pt()
    tofill_gen['l0_eta'    ] = the_pl.eta()
    tofill_gen['l0_phi'    ] = the_pl.phi()
    tofill_gen['l0_q'      ] = the_pl.charge()
    tofill_gen['l0_pdgid'  ] = the_pl.pdgId()
  
    tofill_gen['l1_mass'   ] = the_hn.lep1.mass()
    tofill_gen['l1_pt'     ] = the_hn.lep1.pt()
    tofill_gen['l1_eta'    ] = the_hn.lep1.eta()
    tofill_gen['l1_phi'    ] = the_hn.lep1.phi()
    tofill_gen['l1_q'      ] = the_hn.lep1.charge()
    tofill_gen['l1_pdgid'  ] = the_hn.lep1.pdgId()
    tofill_gen['l1_cov_rad'] = the_hn.lep1.hasConvOrRad
  
    tofill_gen['l2_mass'   ] = the_hn.lep2.mass()
    tofill_gen['l2_pt'     ] = the_hn.lep2.pt()
    tofill_gen['l2_eta'    ] = the_hn.lep2.eta()
    tofill_gen['l2_phi'    ] = the_hn.lep2.phi()
    tofill_gen['l2_q'      ] = the_hn.lep2.charge()
    tofill_gen['l2_pdgid'  ] = the_hn.lep2.pdgId()
    tofill_gen['l2_cov_rad'] = the_hn.lep2.hasConvOrRad
  
    tofill_gen['nu_mass'   ] = the_hn.neu.mass()
    tofill_gen['nu_pt'     ] = the_hn.neu.pt()
    tofill_gen['nu_eta'    ] = the_hn.neu.eta()
    tofill_gen['nu_phi'    ] = the_hn.neu.phi()
    tofill_gen['nu_q'      ] = the_hn.neu.charge()
    tofill_gen['nu_pdgid'  ] = the_hn.neu.pdgId()

    ntuple_gen.Fill(array('f',tofill_gen.values()))

##########################################################################################
# write the ntuples and close the files
outfile_gen.cd()
ntuple_gen.Write()
outfile_gen.Close()
