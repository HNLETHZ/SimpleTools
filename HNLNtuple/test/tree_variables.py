'''
In this file you can define the branch to be included in the flat ntuple.
You can find some basic quantities here, expand with more specific observables,
such as isolation etc...
'''
branches = [
    'run'            ,
    'lumi'           ,
    'event'          ,
    'nvtx'           ,
       
    'in_acc'         , 
       
    'hn_mass'        ,
    'hn_pt'          ,
    'hn_eta'         ,
    'hn_phi'         ,
    'hn_q'           ,
    'hn_vis_mass'    ,
    'hn_vis_pt'      ,
    'hn_vis_eta'     ,
    'hn_vis_phi'     ,
   
    'hn_2d_disp'     ,
    'hn_3d_disp'     ,
   
    'l0_mass'        ,
    'l0_pt'          ,
    'l0_eta'         ,
    'l0_phi'         ,
    'l0_q'           ,
    'l0_pdgid'       ,
    'l0_conv_rad'    ,
    'l0_pt_loss'     ,
  
    'l0_reco_mass'   ,
    'l0_reco_pt'     ,
    'l0_reco_eta'    ,
    'l0_reco_phi'    ,
    'l0_reco_q'      ,
    'l0_reco_pdgid'  ,
    'l0_reco_type'   ,
    'l0_reco_vx'     ,
    'l0_reco_vy'     ,
    'l0_reco_vz'     ,
  
    'l1_mass'        ,
    'l1_pt'          ,
    'l1_eta'         ,
    'l1_phi'         ,
    'l1_q'           ,
    'l1_pdgid'       ,
    'l1_conv_rad'    ,
    'l1_pt_loss'     ,
  
    'l1_reco_mass'   ,
    'l1_reco_pt'     ,
    'l1_reco_eta'    ,
    'l1_reco_phi'    ,
    'l1_reco_q'      ,
    'l1_reco_pdgid'  ,
    'l1_reco_type'   ,
    'l1_reco_vx'     ,
    'l1_reco_vy'     ,
    'l1_reco_vz'     ,
   
    'l2_mass'        ,
    'l2_pt'          ,
    'l2_eta'         ,
    'l2_phi'         ,
    'l2_q'           ,
    'l2_pdgid'       ,
    'l2_conv_rad'    ,
    'l2_pt_loss'     ,
  
    'l2_reco_mass'   ,
    'l2_reco_pt'     ,
    'l2_reco_eta'    ,
    'l2_reco_phi'    ,
    'l2_reco_q'      ,
    'l2_reco_pdgid'  ,
    'l2_reco_type'   ,
    'l2_reco_vx'     ,
    'l2_reco_vy'     ,
    'l2_reco_vz'     ,
   
    'nu_mass'        ,
    'nu_pt'          ,
    'nu_eta'         ,
    'nu_phi'         ,
    'nu_q'           ,
    'nu_pdgid'       ,

    'pv_x'           ,    
    'pv_y'           ,    
    'pv_z'           ,     
   
    'sv_x'           ,    
    'sv_y'           ,    
    'sv_z'           ,     

    'pv_reco_x'      ,    
    'pv_reco_y'      ,    
    'pv_reco_z'      ,     
  
    'sv_reco_x'      ,    
    'sv_reco_y'      ,    
    'sv_reco_z'      ,     
    'sv_reco_prob'   ,     
    
    'hn_reco_2d_disp',
    'hn_reco_3d_disp',

]
