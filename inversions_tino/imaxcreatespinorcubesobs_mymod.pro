;******************************************************************************
; ~/MyProjects/src/idlstuff/idl_jmb/ImaxCreateSpinorCubesObs.pro
;******************************************************************************
; Diese Routine liest die mit ImaxInterpSolarEvolV8.pro erzeugten Cubes oder
; die mit IMaX4DCubes2Fits.gsl erzeugten Stokes-Einzelbilder von
; V8/4-IMaX-Daten ein und speichert sie als 4D-Cubes in einem Format ab, das
; von Spinor für eine Stokes-Inversion gelesen werden kann.
; Es wird dabei das räumliche gemittelte Stokes I,Q,U und V-Profil berechnet
; und ein gegebener Anteil (factor) des mittleren Profils vom Cube abgezogen.
; So soll der Cube grob streulicht-korregiert werden. (Man kann den Faktor Null
; setzen, wenn man keine Streulichtkorrektur dieser Art wünscht.) Anschliessend
; zeigt eine Inversion des streulicht-korregierten Cubes hoffentlich keine
; Feldabschwächung in den Poren. Das mittlere Profil wird als TXT-Datei gespeichert.
;
; Example:
; idl
; IDL> .r ImaxCreateSpinorCubesObs.pro
; IDL> ImaxCreateSpinorCubesObs,0.12
;******************************************************************************
PRO ImaxCreateSpinorCubesObs_mymod,factor

  CLOSE,/ALL

  ;;******************************************************************************
  ;;* 2009 IMaX data
  ;;******************************************************************************
  ;SrcFile = MAKE_ARRAY(1,/STRING)
  ;FileNo = 0
  ;SrcFile[FileNo] = 'imax_177_027_14_27_05_lev2.fits' & FileNo++
  ;NrFiles = FileNo
  ;
  ;SrcDir  = '/data/sunrise/2009/IMaX/level2_TimeInterpolated/'
  ;DstDir  = '/data/slam/home/riethmue/IMaX09/'
  ;
  ;;ROI_X1  = 479    ; BP1Paper, Fig 1: cropped region in horizontally flipped SrcFile: (x,y,w,h) = (479,459,10,9), cropped region (x1,y1)-(x2,y2) = (479,459)-(488,467)
  ;;ROI_Y1  = 936-1-467
  ;;ROI_X2  = 488
  ;;ROI_Y2  = 936-1-459
  ;ROI_X1  = 362    ; BP1Paper, Fig 1: cropped region in horizontally flipped SrcFile: (x,y,w,h) = (362,426,183,147), (x1,y1)-(x2,y2) = (362,426)-(544,572)
  ;ROI_Y1  = 936-1-572
  ;ROI_X2  = 544
  ;ROI_Y2  = 936-1-426
  ;
  ;;* IMaX mode V5/6
  ;wlref      = 5250.2080D0;-0.020D0
  ;Wavelength = double([-0.08D0,-0.04D0,0.04D0,0.08D0,0.227D0])


  ;;******************************************************************************
  ;;* 2013 IMaX data (all 28 maps of the stable series in the first night)
  ;;******************************************************************************
  ;SrcFile = MAKE_ARRAY(28,/STRING)
  ;FileNo = 0
  ;;SrcFile[FileNo] = 'imax_548_001_23_39_10_lev2.fits' & FileNo++
  ;SrcFile[FileNo] = 'imax_548_002_23_39_46_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_003_23_40_23_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_004_23_40_59_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_005_23_41_36_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_006_23_42_12_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_007_23_42_49_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_008_23_43_25_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_009_23_44_02_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_010_23_44_38_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_011_23_45_15_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_012_23_45_51_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_013_23_46_28_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_014_23_47_04_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_015_23_47_41_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_016_23_48_17_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_017_23_48_54_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_018_23_49_30_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_019_23_50_08_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_020_23_50_45_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_021_23_51_21_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_022_23_51_58_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_023_23_52_34_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_024_23_53_11_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_025_23_53_47_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_026_23_54_24_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_027_23_55_00_lev2.fits' & FileNo++
  ;;SrcFile[FileNo] = 'imax_548_028_23_55_37_lev2.fits' & FileNo++
  ;NrFiles = FileNo
  ;
  ;SrcDir  = '/home/tino/nas/2013/IMaX/level2_2014_Apr_TimeInterpolated/12th/'
  ;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN1/Obs548_Cyc002/'
  ;
  ;ROI_X1  = 0
  ;ROI_Y1  = 0
  ;ROI_X2  = 935
  ;ROI_Y2  = 935
  ;
  ;;* IMaX mode V8/4
  ;wlref      = 5250.2080D0
  ;Wavelength = double([-0.12D0,-0.08D0,-0.04D0,0.00D0,0.04D0,0.08D0,0.12D0,0.227D0])


  ;;******************************************************************************
  ;;* 2013 IMaX data reduced by myself
  ;;******************************************************************************
  ;SrcFile = MAKE_ARRAY(34,/STRING)
  ;FileNo = 0
  ;SrcFile[FileNo] = 'imax_548_001_23_39_10_lev1.fits' & FileNo++
  ;NrFiles = FileNo
  ;
  ;SrcDir  = '/home/tino/nas_copy/2013/IMaX/level1.5_2015_Aug_TimeInterpolated/12th/'
  ;DstDir  = '/home/tino/nas_copy/2013/IMaX/level1.5_2015_Aug_TimeInterpolated_SpinorInput/12th/'
  ;;DstDir  = '/home/tino/gwdg/spinor_Michiel/IMaX13/level1_2015_Aug_TimeInterpolated/12th/'
  ;
  ;ROI_X1  = 0
  ;ROI_Y1  = 0
  ;ROI_X2  = 1023
  ;ROI_Y2  = 1023
  ;
  ;;* IMaX mode V8/4
  ;wlref      = 5250.2080D0
  ;Wavelength = double([-0.12D0,-0.08D0,-0.04D0,0.00D0,0.04D0,0.08D0,0.12D0,0.227D0])


  ;;******************************************************************************
  ;;* 2013 IMaX data reduced by Julian (reduced FOV of the first map for test purposes)
  ;;******************************************************************************
  ;SrcFile = MAKE_ARRAY(28,/STRING)
  ;FileNo = 0
  ;SrcFile[FileNo] = 'imax_548_002_23_39_46_lev1.fits' & FileNo++
  ;NrFiles = FileNo
  ;
  ;;SrcDir  = '/data/sunrise/2013/IMaX/level2_2014_Apr_TimeInterpolated/12th/'
  ;;DstDir  = '/data/slam/home/riethmue/IMaX13/'
  ;
  ;;SrcDir  = '/home/tino/nas/2013/IMaX/level2_2014_Apr_TimeInterpolated/12th/'
  ;;DstDir  = '/home/tino/gwdg/spinor_Michiel/IMaX13/level2_2014_Apr_TimeInterpolated/12th/'
  ;
  ;SrcDir  = '/home/tino/nas/2013/IMaX/level1_2014_Apr_TimeInterpolated/12th/'
  ;;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN1_0.35MeanIQUVProfileSubtracted/Obs548_Cyc002/'
  ;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN1_0.35MeanIProfileSubtracted/Obs548_Cyc002/'
  ;
  ;ROI_X1  = 463    ; small test region with only the two small pores clearly showing the field weakening effect
  ;ROI_Y1  = 292
  ;ROI_X2  = ROI_X1+140-1
  ;ROI_Y2  = ROI_Y1+80-1
  ;
  ;;* IMaX mode V8/4
  ;wlref      = 5250.2080D0
  ;Wavelength = double([-0.12D0,-0.08D0,-0.04D0,0.00D0,0.04D0,0.08D0,0.12D0,0.227D0])


  ;;******************************************************************************
  ;;* 2013 IMaX data reduced by Julian (reduced FOV of the first map for test purposes)
  ;;******************************************************************************
  ;SrcFile = MAKE_ARRAY(28,/STRING)
  ;FileNo = 0
  ;SrcFile[FileNo] = 'imax_548_002_23_39_46_lev2.fits' & FileNo++
  ;NrFiles = FileNo
  ;
  ;SrcDir  = '/home/tino/nas/2013/IMaX/level2_2014_Apr_TimeInterpolated/12th/'
  ;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN1_0.35MeanIProfileSubtracted_FullFov/Obs548_Cyc002/'
  ;
  ;ROI_X1  = 62
  ;ROI_Y1  = 62
  ;ROI_X2  = 873
  ;ROI_Y2  = 873
  ;
  ;;* IMaX mode V8/4
  ;wlref      = 5250.2080D0
  ;Wavelength = double([-0.12D0,-0.08D0,-0.04D0,0.00D0,0.04D0,0.08D0,0.12D0,0.227D0])


  ;;******************************************************************************
  ;;* 2013 IMaX data reduced by Julian in April 2014
  ;;******************************************************************************
  ;SrcFile = MAKE_ARRAY(28,/STRING)
  ;FileNo = 0
  ;SrcFile[FileNo] = 'imax_548_002_23_39_46_lev1.fits' & FileNo++
  ;NrFiles = FileNo
  ;
  ;SrcDir  = '/home/tino/nas/2013/IMaX/level1_2014_Apr_TimeInterpolated/12th/'
  ;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN1_0.35MeanIProfileSubtracted_FullFov/Obs548_Cyc002/lev1/'
  ;;SrcDir  = '/data/sunrise/2013/IMaX/level1_2014_Apr_TimeInterpolated/12th/'
  ;;DstDir  = '/data/slam/home/riethmue/IMaX13/RUN1_0.35MeanIProfileSubtracted_FullFov/Obs548_Cyc002/lev1/'
  ;
  ;ROI_X1  = 0
  ;ROI_Y1  = 0
  ;ROI_X2  = 935
  ;ROI_Y2  = 935
  ;
  ;;* IMaX mode V8/4
  ;wlref      = 5250.2080D0
  ;Wavelength = double([-0.12D0,-0.08D0,-0.04D0,0.00D0,0.04D0,0.08D0,0.12D0,0.227D0])


  ;;******************************************************************************
  ;;* 2013 IMaX data reduced by Julian with an improved cross-talk removal to get rid of the field weakening in pores
  ;;******************************************************************************
  ;SrcFile = MAKE_ARRAY(28,/STRING)
  ;FileNo = 0
  ;SrcFile[FileNo] = 'imax_548_028_23_55_37_lev1.fits' & FileNo++
  ;NrFiles = FileNo
  ;
  ;SrcDir  = '/home/tino/nas/2013/IMaX/Julian/ImprovedCrossTalk/'
  ;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN1_ImprovedCrossTalk_0.30MeanIProfileSubtracted/Obs548_Cyc028/'
  ;
  ;ROI_X1  = 0
  ;ROI_Y1  = 0
  ;ROI_X2  = 935
  ;ROI_Y2  = 935
  ;
  ;;* IMaX mode V8/4
  ;wlref      = 5250.2080D0
  ;Wavelength = double([-0.12D0,-0.08D0,-0.04D0,0.00D0,0.04D0,0.08D0,0.12D0,0.227D0])


  ;;******************************************************************************
  ;;* 2013 IMaX data reduced by Julian in April 2014
  ;;******************************************************************************
  ;SrcFile = MAKE_ARRAY(28,/STRING)
  ;FileNo = 0
  ;SrcFile[FileNo] = 'imax_548_002_23_39_46_lev1.fits' & FileNo++
  ;NrFiles = FileNo
  ;
  ;SrcDir  = '/data/sunrise/2013/IMaX/level1_2014_Apr_TimeInterpolated/12th/'
  ;DstDir  = '/data/slam/home/riethmue/IMaX13/RUN1_2d_EfPsf_0.30MeanIProfileSubtracted/Obs548_Cyc002/'
  ;
  ;ROI_X1  = 0
  ;ROI_Y1  = 0
  ;ROI_X2  = 935
  ;ROI_Y2  = 935
  ;
  ;;* IMaX mode V8/4
  ;wlref      = 5250.2080D0
  ;Wavelength = double([-0.12D0,-0.08D0,-0.04D0,0.00D0,0.04D0,0.08D0,0.12D0,0.227D0])


  ;;******************************************************************************
  ;;* 2013 IMaX data reduced by Julian in April 2014 and stray-light PSF ~/nas_copy/2013/IMaX/Tino/StrayLightAnalysis.Fibre/ManuallyTunedStrayLightPsf_936.fits removed from Stokes I only
  ;;******************************************************************************
  ;SrcFile = MAKE_ARRAY(28,/STRING)
  ;FileNo = 0
  ;SrcFile[FileNo] = '548_002/' & FileNo++
  ;NrFiles = FileNo
  ;
  ;SrcDir  = '/home/tino/nas/2013/IMaX/level1_2014_Apr_StraylightRemoved_v3/12th/reduc/'
  ;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN1_ManuallyTunedStraylightRemoved_v3/Obs548_Cyc002/'
  ;
  ;w       = 936
  ;h       = 936
  ;ROI_X1  = 0
  ;ROI_Y1  = 0
  ;ROI_X2  = 935
  ;ROI_Y2  = 935
  ;
  ;;* IMaX mode V8/4
  ;wlref      = 5250.2080D0
  ;Wavelength = double([-0.12D0,-0.08D0,-0.04D0,0.00D0,0.04D0,0.08D0,0.12D0,0.227D0])


  ;;******************************************************************************
  ;;* 2013 IMaX data reduced by Julian in April 2014
  ;;******************************************************************************
  ;SrcFile = MAKE_ARRAY(28,/STRING)
  ;FileNo = 0
  ;SrcFile[FileNo] = '548_002/' & FileNo++
  ;NrFiles = FileNo
  ;
  ;SrcDir  = '/home/tino/nas_copy/2013/IMaX/level2_2014_Apr_TimeInterpolated/12th/reduc/'
  ;DstDir  = '/home/tino/SUNRISE_Data/IMaX/MASI.20LinesV46/091584/lev2/'
  ;
  ;w       = 936
  ;h       = 936
  ;ROI_X1  = 62
  ;ROI_Y1  = 62
  ;ROI_X2  = 873
  ;ROI_Y2  = 873
  ;
  ;;* IMaX mode V8/4
  ;wlref      = 5250.2080D0
  ;Wavelength = double([-0.12D0,-0.08D0,-0.04D0,0.00D0,0.04D0,0.08D0,0.12D0,0.227D0])


  ;;******************************************************************************
  ;;* 2013 IMaX data reduced by Julian in Jan. 2016, two small pores only
  ;;******************************************************************************
  ;SrcFile = MAKE_ARRAY(28,/STRING)
  ;FileNo = 0
  ;SrcFile[FileNo] = '548_002/' & FileNo++
  ;NrFiles = FileNo
  ;
  ;SrcDir  = '/home/tino/nas_copy/2013/IMaX/level1.2_2016_Jan/12th/reduc/'
  ;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN1/Obs548_Cyc002/lev1.2_mojo_GlobSL0.332_140x105/'
  ;
  ;w       = 936
  ;h       = 936
  ;;ROI_X1  = 463    ; small test region with only the two small pores clearly showing the field weakening effect
  ;;ROI_Y1  = 292
  ;;ROI_X2  = ROI_X1+140-1
  ;;ROI_Y2  = ROI_Y1+ 80-1
  ;ROI_X1  = 463    ; small test region with only the two small pores clearly showing the field weakening effect
  ;ROI_Y1  = 279
  ;ROI_X2  = ROI_X1+140-1
  ;ROI_Y2  = ROI_Y1+105-1
  ;
  ;;* IMaX mode V8/4
  ;wlref      = 5250.2080D0
  ;Wavelength = double([-0.12D0,-0.08D0,-0.04D0,0.00D0,0.04D0,0.08D0,0.12D0,0.227D0])


  ;;******************************************************************************
  ;;* 2013 IMaX data reduced by Julian in Jan. 2016, full FOV
  ;;******************************************************************************
  ;SrcFile = MAKE_ARRAY(28,/STRING)
  ;FileNo = 0
  ;SrcFile[FileNo] = '548_028/' & FileNo++
  ;NrFiles = FileNo
  ;
  ;;SrcDir  = '/home/tino/nas_copy/2013/IMaX/level2.5_2016_Jan/12th/reduc/'
  ;;SrcDir  = '/home/tino/nas_copy/2013/IMaX/level1.2_2016_Jan/12th/reduc/'
  ;;SrcDir  = '/home/tino/nas_copy/2013/IMaX/level2.2_2016_Jan/12th/reduc/'
  ;;SrcDir  = '/data/sunrise/2013/IMaX/level2.2_2016_Jan/12th/reduc/'
  ;;SrcDir  = '/data/sunrise/2013/IMaX/level1.2_2016_Jan/12th/reduc/'
  ;;SrcDir  = '/home/tino/nas/2013/IMaX/level2.1_2016_Jan/12th/reduc/'
  ;SrcDir  = '/home/tino/nas/2013/IMaX/level1.1_2016_Jan/12th/reduc/'
  ;
  ;;DstDir  = '/home/tino/gwdg/spinor_Michiel/IMaX13/RUN1/Obs548_Cyc002/lev2.5/'
  ;;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN1/Obs548_Cyc002/lev2.5_mojo_GlobSL0.000_812x812/'
  ;;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN1/Obs548_Cyc002/lev1.2_atlas_GlobSL0.138_936x936/'
  ;;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN3b/Obs548_Cyc002/lev1.2_atlas_GlobSL0.126_936x936/'
  ;;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN1/Obs548_Cyc002/lev2.2_atlas_GlobSL0.350_936x936/'
  ;;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN1/Obs548_Cyc002/lev2.2_atlas_GlobSL0.200_936x936/'
  ;;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN1/Obs548_Cyc002/lev1.2_atlas_GlobSL0.032_936x936/'
  ;;DstDir  = '/home/tino/gwdg/spinor_Michiel/IMaX13/RUN1/Obs548_Cyc002/lev2.2_atlas_GlobSL0.250_936x936/'
  ;;DstDir  = '/home/tino/gwdg/spinor_Michiel/IMaX13/RUN1/Obs548_Cyc002/lev2.2_atlas_GlobSL0.300_936x936/'
  ;;DstDir  = '/home/tino/gwdg/spinor_Michiel/IMaX13/RUN1/Obs548_Cyc002/lev2.2_atlas_GlobSL0.230_936x936/'
  ;;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN1/Obs548_Cyc002/lev1.2_atlas_GlobSL0.165_936x936/'
  ;;DstDir  = '/data/slam/home/riethmue/IMaX13/RUN1/Obs548_Cyc001/lev2.2_atlas_GlobSL0.250_936x936/'
  ;;DstDir  = '/data/slam/home/riethmue/IMaX13/RUN1/Obs548_Cyc002/lev1.2_GlobSL0.350_936x936/'
  ;;DstDir  = '/home/tino/gwdg/spinor_Michiel/IMaX13/RUN5_Andrei_lev2/Obs548_Cyc001/lev2.1_GlobSL0.250_936x936/'
  ;DstDir  = '/home/tino/gwdg/spinor_Michiel/IMaX13/RUN5_Andrei_lev1/Obs548_Cyc028/lev1.1_GlobSL0.250_936x936/'
  ;
  ;w       = 936
  ;h       = 936
  ;;ROI_X1  = 62
  ;;ROI_Y1  = 62
  ;;ROI_X2  = 873
  ;;ROI_Y2  = 873
  ;ROI_X1  = 0
  ;ROI_Y1  = 0
  ;ROI_X2  = 935
  ;ROI_Y2  = 935
  ;
  ;;* IMaX mode V8/4
  ;wlref      = 5250.2080D0
  ;Wavelength = double([-0.12D0,-0.08D0,-0.04D0,0.00D0,0.04D0,0.08D0,0.12D0,0.227D0])

  ;;******************************************************************************
  ;;* 2013 IMaX data reduced by Julian in Apr. 2014, MASI FOV
  ;;******************************************************************************
  ;SrcFile = MAKE_ARRAY(28,/STRING)
  ;FileNo = 0
  ;SrcFile[FileNo] = '548_002/' & FileNo++
  ;NrFiles = FileNo
  ;
  ;;SrcDir  = '/home/tino/nas_copy/2013/IMaX/Apr_2014_Reduction/level2_2014_Apr/12th/reduc/'
  ;;DstDir  = '/home/tino/SUNRISE_Data/IMaX/MASI.20LinesV46/088000/lev1/'
  ;
  ;;SrcDir  = '/home/tino/nas_copy/2013/IMaX/Apr_2014_Reduction/level2_2014_Apr/12th/reduc/'
  ;;DstDir  = '/home/tino/SUNRISE_Data/IMaX/MASI.20LinesV46/091584/lev2.BP3x_0.25GlSl_1Pol_DistWeighted/'
  ;
  ;;SrcDir  = '/home/tino/nas_copy/2013/IMaX/Apr_2014_Reduction/level1_2014_Apr/12th/reduc/'
  ;;DstDir  = '/home/tino/SUNRISE_Data/IMaX/MASI.20LinesV46/091584/lev1.BP3x_0.25GlSl_1Pol_DistWeighted/'
  ;
  ;;SrcDir  = '/home/tino/nas_copy/2013/IMaX/Apr_2014_Reduction/level1_2014_Apr/12th/reduc/'
  ;;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN1/Obs548_Cyc002/level1_2014_Apr_GlobSL0.25_812x812/'
  ;
  ;SrcDir  = '/home/tino/nas_copy/2013/IMaX/Apr_2014_Reduction/level2_2014_Apr/12th/reduc/'
  ;DstDir  = '/home/tino/pulpodata/home/riethmue/IMaX13/RUN1/Obs548_Cyc002/level2_2014_Apr_GlobSL0.00_812x812/'
  ;
  ;w       = 936
  ;h       = 936
  ;ROI_X1  = 62
  ;ROI_Y1  = 62
  ;ROI_X2  = 873
  ;ROI_Y2  = 873
  ;
  ;
  ;;* IMaX mode V8/4
  ;wlref      = 5250.2080D0
  ;Wavelength = double([-0.12D0,-0.08D0,-0.04D0,0.00D0,0.04D0,0.08D0,0.12D0,0.227D0])

  ;;******************************************************************************
  ;;* MASI 2nd iter. data 812x812
  ;;******************************************************************************
  ;SrcFile = MAKE_ARRAY(28,/STRING)
  ;FileNo = 0
  ;SrcFile[FileNo] = '000_0915/' & FileNo++
  ;NrFiles = FileNo
  ;
  ;SrcDir  = '/home/tino/MURaM/poren_unipolar400G_576x384x576/091584/lev2/SPINOR_Synth_ImaxProfile.pos_936_BP3x_0.25GlobalStrayLight/'
  ;DstDir  = '/home/tino/MURaM/poren_unipolar400G_576x384x576/091584/lev2/SPINOR_Synth_ImaxProfile.pos_936_BP3x_0.25GlobalStrayLight/'
  ;
  ;w       = 936
  ;h       = 936
  ;ROI_X1  = 62
  ;ROI_Y1  = 62
  ;ROI_X2  = 873
  ;ROI_Y2  = 873
  ;
  ;       
  ;;* IMaX mode V46
  ;wlref      = 5250.2080D0
  ;Wavelength = double([-0.15D0,-0.14D0,-0.13D0,-0.12D0,-0.11D0,-0.10D0,-0.09D0,-0.08D0,-0.07D0,-0.06D0,-0.05D0,-0.04D0,-0.03D0,-0.02D0,-0.01D0,0.00D0,0.01D0,0.02D0,0.03D0,0.04D0,0.05D0,0.06D0,0.07D0,0.08D0,0.09D0,0.10D0,0.11D0,0.12D0,0.13D0,0.14D0,0.15D0,0.16D0,0.17D0,0.18D0,0.19D0,0.20D0,0.21D0,0.22D0,0.23D0,0.24D0,0.25D0,0.26D0,0.27D0,0.28D0,0.29D0,0.30D0])

  ;******************************************************************************
  ;* 2009 IMaX data, full FOV
  ;******************************************************************************
  SrcFile = MAKE_ARRAY(2,/STRING)
  FileNo = 0
  SrcFile[FileNo] = 'mean_rem_output_21.fits' & FileNo++
  SrcFile[FileNo] = 'mean_rem_output_22.fits' & FileNo++
  NrFiles = FileNo

  SrcDir  = '/scratch/prabhu/backup_workstation/sunrise_holly/work_IMaX_pole_scripts/inversions_tino/'
  DstDir  = '/scratch/prabhu/backup_workstation/sunrise_holly/work_IMaX_pole_scripts/inversions_tino/interp/'
 
  w       = 936
  h       = 936
  ROI_X1  = 490;0
  ROI_Y1  = 490;0
  ROI_X2  = 500;935
  ROI_Y2  = 500;935

  ;* IMaX mode V5/6
  wlref      = 5250.2080D0
  Wavelength = double([-0.08D0,-0.04D0,0.04D0,0.08D0,0.227D0])
  NWL = 5






  ;;******************************************************************************
  ;;* M A I N
  ;;******************************************************************************

  for FileNo=0,NrFiles-1 do begin
    ;create sub-dir
    spawn,'mkdir -p '+DstDir

    ; read input
    if STRMID(SrcFile[FileNo],0,1,/REVERSE_OFFSET) EQ '/' then begin ; single Stokes images created by IMaX4DCubes2Fits.gsl
      ObsCyc = STRMID(SrcFile[FileNo],0,STRLEN(SrcFile[FileNo])-1)
      NrStokes = 4
      if (size(Wavelength))[1] EQ 46 then begin ; V46
        NWL = 46
        ContIdx = 38
        im = fltarr(w,h,NWL,NrStokes)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m150_I.fits.gz' & print,'Reading file ',FileName & im[*,*, 0,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m150_Q.fits.gz' & print,'Reading file ',FileName & im[*,*, 0,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m150_U.fits.gz' & print,'Reading file ',FileName & im[*,*, 0,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m150_V.fits.gz' & print,'Reading file ',FileName & im[*,*, 0,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m140_I.fits.gz' & print,'Reading file ',FileName & im[*,*, 1,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m140_Q.fits.gz' & print,'Reading file ',FileName & im[*,*, 1,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m140_U.fits.gz' & print,'Reading file ',FileName & im[*,*, 1,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m140_V.fits.gz' & print,'Reading file ',FileName & im[*,*, 1,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m130_I.fits.gz' & print,'Reading file ',FileName & im[*,*, 2,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m130_Q.fits.gz' & print,'Reading file ',FileName & im[*,*, 2,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m130_U.fits.gz' & print,'Reading file ',FileName & im[*,*, 2,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m130_V.fits.gz' & print,'Reading file ',FileName & im[*,*, 2,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m120_I.fits.gz' & print,'Reading file ',FileName & im[*,*, 3,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m120_Q.fits.gz' & print,'Reading file ',FileName & im[*,*, 3,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m120_U.fits.gz' & print,'Reading file ',FileName & im[*,*, 3,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m120_V.fits.gz' & print,'Reading file ',FileName & im[*,*, 3,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m110_I.fits.gz' & print,'Reading file ',FileName & im[*,*, 4,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m110_Q.fits.gz' & print,'Reading file ',FileName & im[*,*, 4,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m110_U.fits.gz' & print,'Reading file ',FileName & im[*,*, 4,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m110_V.fits.gz' & print,'Reading file ',FileName & im[*,*, 4,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m100_I.fits.gz' & print,'Reading file ',FileName & im[*,*, 5,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m100_Q.fits.gz' & print,'Reading file ',FileName & im[*,*, 5,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m100_U.fits.gz' & print,'Reading file ',FileName & im[*,*, 5,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m100_V.fits.gz' & print,'Reading file ',FileName & im[*,*, 5,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m090_I.fits.gz' & print,'Reading file ',FileName & im[*,*, 6,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m090_Q.fits.gz' & print,'Reading file ',FileName & im[*,*, 6,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m090_U.fits.gz' & print,'Reading file ',FileName & im[*,*, 6,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m090_V.fits.gz' & print,'Reading file ',FileName & im[*,*, 6,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m080_I.fits.gz' & print,'Reading file ',FileName & im[*,*, 7,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m080_Q.fits.gz' & print,'Reading file ',FileName & im[*,*, 7,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m080_U.fits.gz' & print,'Reading file ',FileName & im[*,*, 7,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m080_V.fits.gz' & print,'Reading file ',FileName & im[*,*, 7,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m070_I.fits.gz' & print,'Reading file ',FileName & im[*,*, 8,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m070_Q.fits.gz' & print,'Reading file ',FileName & im[*,*, 8,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m070_U.fits.gz' & print,'Reading file ',FileName & im[*,*, 8,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m070_V.fits.gz' & print,'Reading file ',FileName & im[*,*, 8,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m060_I.fits.gz' & print,'Reading file ',FileName & im[*,*, 9,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m060_Q.fits.gz' & print,'Reading file ',FileName & im[*,*, 9,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m060_U.fits.gz' & print,'Reading file ',FileName & im[*,*, 9,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m060_V.fits.gz' & print,'Reading file ',FileName & im[*,*, 9,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m050_I.fits.gz' & print,'Reading file ',FileName & im[*,*,10,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m050_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,10,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m050_U.fits.gz' & print,'Reading file ',FileName & im[*,*,10,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m050_V.fits.gz' & print,'Reading file ',FileName & im[*,*,10,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m040_I.fits.gz' & print,'Reading file ',FileName & im[*,*,11,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m040_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,11,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m040_U.fits.gz' & print,'Reading file ',FileName & im[*,*,11,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m040_V.fits.gz' & print,'Reading file ',FileName & im[*,*,11,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m030_I.fits.gz' & print,'Reading file ',FileName & im[*,*,12,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m030_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,12,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m030_U.fits.gz' & print,'Reading file ',FileName & im[*,*,12,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m030_V.fits.gz' & print,'Reading file ',FileName & im[*,*,12,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m020_I.fits.gz' & print,'Reading file ',FileName & im[*,*,13,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m020_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,13,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m020_U.fits.gz' & print,'Reading file ',FileName & im[*,*,13,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m020_V.fits.gz' & print,'Reading file ',FileName & im[*,*,13,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m010_I.fits.gz' & print,'Reading file ',FileName & im[*,*,14,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m010_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,14,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m010_U.fits.gz' & print,'Reading file ',FileName & im[*,*,14,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m010_V.fits.gz' & print,'Reading file ',FileName & im[*,*,14,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p000_I.fits.gz' & print,'Reading file ',FileName & im[*,*,15,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p000_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,15,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p000_U.fits.gz' & print,'Reading file ',FileName & im[*,*,15,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p000_V.fits.gz' & print,'Reading file ',FileName & im[*,*,15,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p010_I.fits.gz' & print,'Reading file ',FileName & im[*,*,16,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p010_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,16,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p010_U.fits.gz' & print,'Reading file ',FileName & im[*,*,16,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p010_V.fits.gz' & print,'Reading file ',FileName & im[*,*,16,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p020_I.fits.gz' & print,'Reading file ',FileName & im[*,*,17,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p020_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,17,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p020_U.fits.gz' & print,'Reading file ',FileName & im[*,*,17,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p020_V.fits.gz' & print,'Reading file ',FileName & im[*,*,17,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p030_I.fits.gz' & print,'Reading file ',FileName & im[*,*,18,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p030_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,18,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p030_U.fits.gz' & print,'Reading file ',FileName & im[*,*,18,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p030_V.fits.gz' & print,'Reading file ',FileName & im[*,*,18,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p040_I.fits.gz' & print,'Reading file ',FileName & im[*,*,19,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p040_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,19,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p040_U.fits.gz' & print,'Reading file ',FileName & im[*,*,19,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p040_V.fits.gz' & print,'Reading file ',FileName & im[*,*,19,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p050_I.fits.gz' & print,'Reading file ',FileName & im[*,*,20,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p050_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,20,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p050_U.fits.gz' & print,'Reading file ',FileName & im[*,*,20,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p050_V.fits.gz' & print,'Reading file ',FileName & im[*,*,20,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p060_I.fits.gz' & print,'Reading file ',FileName & im[*,*,21,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p060_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,21,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p060_U.fits.gz' & print,'Reading file ',FileName & im[*,*,21,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p060_V.fits.gz' & print,'Reading file ',FileName & im[*,*,21,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p070_I.fits.gz' & print,'Reading file ',FileName & im[*,*,22,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p070_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,22,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p070_U.fits.gz' & print,'Reading file ',FileName & im[*,*,22,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p070_V.fits.gz' & print,'Reading file ',FileName & im[*,*,22,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p080_I.fits.gz' & print,'Reading file ',FileName & im[*,*,23,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p080_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,23,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p080_U.fits.gz' & print,'Reading file ',FileName & im[*,*,23,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p080_V.fits.gz' & print,'Reading file ',FileName & im[*,*,23,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p090_I.fits.gz' & print,'Reading file ',FileName & im[*,*,24,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p090_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,24,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p090_U.fits.gz' & print,'Reading file ',FileName & im[*,*,24,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p090_V.fits.gz' & print,'Reading file ',FileName & im[*,*,24,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p100_I.fits.gz' & print,'Reading file ',FileName & im[*,*,25,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p100_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,25,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p100_U.fits.gz' & print,'Reading file ',FileName & im[*,*,25,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p100_V.fits.gz' & print,'Reading file ',FileName & im[*,*,25,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p110_I.fits.gz' & print,'Reading file ',FileName & im[*,*,26,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p110_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,26,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p110_U.fits.gz' & print,'Reading file ',FileName & im[*,*,26,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p110_V.fits.gz' & print,'Reading file ',FileName & im[*,*,26,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p120_I.fits.gz' & print,'Reading file ',FileName & im[*,*,27,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p120_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,27,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p120_U.fits.gz' & print,'Reading file ',FileName & im[*,*,27,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p120_V.fits.gz' & print,'Reading file ',FileName & im[*,*,27,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p130_I.fits.gz' & print,'Reading file ',FileName & im[*,*,28,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p130_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,28,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p130_U.fits.gz' & print,'Reading file ',FileName & im[*,*,28,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p130_V.fits.gz' & print,'Reading file ',FileName & im[*,*,28,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p140_I.fits.gz' & print,'Reading file ',FileName & im[*,*,29,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p140_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,29,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p140_U.fits.gz' & print,'Reading file ',FileName & im[*,*,29,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p140_V.fits.gz' & print,'Reading file ',FileName & im[*,*,29,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p150_I.fits.gz' & print,'Reading file ',FileName & im[*,*,30,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p150_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,30,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p150_U.fits.gz' & print,'Reading file ',FileName & im[*,*,30,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p150_V.fits.gz' & print,'Reading file ',FileName & im[*,*,30,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p160_I.fits.gz' & print,'Reading file ',FileName & im[*,*,31,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p160_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,31,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p160_U.fits.gz' & print,'Reading file ',FileName & im[*,*,31,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p160_V.fits.gz' & print,'Reading file ',FileName & im[*,*,31,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p170_I.fits.gz' & print,'Reading file ',FileName & im[*,*,32,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p170_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,32,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p170_U.fits.gz' & print,'Reading file ',FileName & im[*,*,32,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p170_V.fits.gz' & print,'Reading file ',FileName & im[*,*,32,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p180_I.fits.gz' & print,'Reading file ',FileName & im[*,*,33,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p180_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,33,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p180_U.fits.gz' & print,'Reading file ',FileName & im[*,*,33,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p180_V.fits.gz' & print,'Reading file ',FileName & im[*,*,33,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p190_I.fits.gz' & print,'Reading file ',FileName & im[*,*,34,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p190_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,34,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p190_U.fits.gz' & print,'Reading file ',FileName & im[*,*,34,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p190_V.fits.gz' & print,'Reading file ',FileName & im[*,*,34,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p200_I.fits.gz' & print,'Reading file ',FileName & im[*,*,35,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p200_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,35,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p200_U.fits.gz' & print,'Reading file ',FileName & im[*,*,35,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p200_V.fits.gz' & print,'Reading file ',FileName & im[*,*,35,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p210_I.fits.gz' & print,'Reading file ',FileName & im[*,*,36,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p210_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,36,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p210_U.fits.gz' & print,'Reading file ',FileName & im[*,*,36,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p210_V.fits.gz' & print,'Reading file ',FileName & im[*,*,36,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p220_I.fits.gz' & print,'Reading file ',FileName & im[*,*,37,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p220_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,37,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p220_U.fits.gz' & print,'Reading file ',FileName & im[*,*,37,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p220_V.fits.gz' & print,'Reading file ',FileName & im[*,*,37,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p230_I.fits.gz' & print,'Reading file ',FileName & im[*,*,38,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p230_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,38,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p230_U.fits.gz' & print,'Reading file ',FileName & im[*,*,38,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p230_V.fits.gz' & print,'Reading file ',FileName & im[*,*,38,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p240_I.fits.gz' & print,'Reading file ',FileName & im[*,*,39,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p240_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,39,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p240_U.fits.gz' & print,'Reading file ',FileName & im[*,*,39,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p240_V.fits.gz' & print,'Reading file ',FileName & im[*,*,39,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p250_I.fits.gz' & print,'Reading file ',FileName & im[*,*,40,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p250_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,40,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p250_U.fits.gz' & print,'Reading file ',FileName & im[*,*,40,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p250_V.fits.gz' & print,'Reading file ',FileName & im[*,*,40,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p260_I.fits.gz' & print,'Reading file ',FileName & im[*,*,41,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p260_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,41,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p260_U.fits.gz' & print,'Reading file ',FileName & im[*,*,41,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p260_V.fits.gz' & print,'Reading file ',FileName & im[*,*,41,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p270_I.fits.gz' & print,'Reading file ',FileName & im[*,*,42,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p270_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,42,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p270_U.fits.gz' & print,'Reading file ',FileName & im[*,*,42,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p270_V.fits.gz' & print,'Reading file ',FileName & im[*,*,42,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p280_I.fits.gz' & print,'Reading file ',FileName & im[*,*,43,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p280_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,43,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p280_U.fits.gz' & print,'Reading file ',FileName & im[*,*,43,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p280_V.fits.gz' & print,'Reading file ',FileName & im[*,*,43,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p290_I.fits.gz' & print,'Reading file ',FileName & im[*,*,44,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p290_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,44,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p290_U.fits.gz' & print,'Reading file ',FileName & im[*,*,44,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p290_V.fits.gz' & print,'Reading file ',FileName & im[*,*,44,3] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p300_I.fits.gz' & print,'Reading file ',FileName & im[*,*,45,0] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p300_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,45,1] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p300_U.fits.gz' & print,'Reading file ',FileName & im[*,*,45,2] = readfits(FileName)
        FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p300_V.fits.gz' & print,'Reading file ',FileName & im[*,*,45,3] = readfits(FileName)
      endif else begin ; V8/4 or V5/6
        if (size(Wavelength))[1] EQ 8 then begin ; V8/4
          NWL = 8
          ContIdx = NWL-1
          im = fltarr(w,h,NWL,NrStokes)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m120_I.fits.gz' & print,'Reading file ',FileName & im[*,*,0,0] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m120_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,0,1] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m120_U.fits.gz' & print,'Reading file ',FileName & im[*,*,0,2] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m120_V.fits.gz' & print,'Reading file ',FileName & im[*,*,0,3] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m080_I.fits.gz' & print,'Reading file ',FileName & im[*,*,1,0] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m080_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,1,1] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m080_U.fits.gz' & print,'Reading file ',FileName & im[*,*,1,2] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m080_V.fits.gz' & print,'Reading file ',FileName & im[*,*,1,3] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m040_I.fits.gz' & print,'Reading file ',FileName & im[*,*,2,0] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m040_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,2,1] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m040_U.fits.gz' & print,'Reading file ',FileName & im[*,*,2,2] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m040_V.fits.gz' & print,'Reading file ',FileName & im[*,*,2,3] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p000_I.fits.gz' & print,'Reading file ',FileName & im[*,*,3,0] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p000_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,3,1] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p000_U.fits.gz' & print,'Reading file ',FileName & im[*,*,3,2] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p000_V.fits.gz' & print,'Reading file ',FileName & im[*,*,3,3] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p040_I.fits.gz' & print,'Reading file ',FileName & im[*,*,4,0] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p040_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,4,1] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p040_U.fits.gz' & print,'Reading file ',FileName & im[*,*,4,2] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p040_V.fits.gz' & print,'Reading file ',FileName & im[*,*,4,3] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p080_I.fits.gz' & print,'Reading file ',FileName & im[*,*,5,0] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p080_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,5,1] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p080_U.fits.gz' & print,'Reading file ',FileName & im[*,*,5,2] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p080_V.fits.gz' & print,'Reading file ',FileName & im[*,*,5,3] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p120_I.fits.gz' & print,'Reading file ',FileName & im[*,*,6,0] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p120_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,6,1] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p120_U.fits.gz' & print,'Reading file ',FileName & im[*,*,6,2] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p120_V.fits.gz' & print,'Reading file ',FileName & im[*,*,6,3] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p227_I.fits.gz' & print,'Reading file ',FileName & im[*,*,7,0] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p227_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,7,1] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p227_U.fits.gz' & print,'Reading file ',FileName & im[*,*,7,2] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p227_V.fits.gz' & print,'Reading file ',FileName & im[*,*,7,3] = readfits(FileName)
        endif else begin ; V5/6
          NWL = 5
          ContIdx = NWL-1
          im = fltarr(w,h,NWL,NrStokes)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m080_I.fits.gz' & print,'Reading file ',FileName & im[*,*,0,0] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m080_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,0,1] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m080_U.fits.gz' & print,'Reading file ',FileName & im[*,*,0,2] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m080_V.fits.gz' & print,'Reading file ',FileName & im[*,*,0,3] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m040_I.fits.gz' & print,'Reading file ',FileName & im[*,*,1,0] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m040_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,1,1] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m040_U.fits.gz' & print,'Reading file ',FileName & im[*,*,1,2] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_m040_V.fits.gz' & print,'Reading file ',FileName & im[*,*,1,3] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p040_I.fits.gz' & print,'Reading file ',FileName & im[*,*,2,0] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p040_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,2,1] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p040_U.fits.gz' & print,'Reading file ',FileName & im[*,*,2,2] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p040_V.fits.gz' & print,'Reading file ',FileName & im[*,*,2,3] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p080_I.fits.gz' & print,'Reading file ',FileName & im[*,*,3,0] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p080_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,3,1] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p080_U.fits.gz' & print,'Reading file ',FileName & im[*,*,3,2] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p080_V.fits.gz' & print,'Reading file ',FileName & im[*,*,3,3] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p227_I.fits.gz' & print,'Reading file ',FileName & im[*,*,4,0] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p227_Q.fits.gz' & print,'Reading file ',FileName & im[*,*,4,1] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p227_U.fits.gz' & print,'Reading file ',FileName & im[*,*,4,2] = readfits(FileName)
          FileName = SrcDir+SrcFile[FileNo]+'reduc_rnr_'+ObsCyc+'_p227_V.fits.gz' & print,'Reading file ',FileName & im[*,*,4,3] = readfits(FileName)
        endelse
      endelse
      MeanProfileFileName = DstDir+'MeanProfile.txt'
      DstFileName=DstDir+'imax.fits'
      IContFileName=DstDir+'ICont.fits'
    endif else begin                                                 ; 4D cubes created by ImaxInterpSolarEvolV8.pro
      ContIdx = NWL-1
      BaseName = STRMID(SrcFile[FileNo],0,18)                       ; 'imax_177_033' aus 'imax_177_033_14_30_25_lev2.fits' ausschneiden
      print,'Reading file ',FileNo+1,' of ',NrFiles,': ',SrcDir+SrcFile[FileNo]
      im = readfits(SrcDir+SrcFile[FileNo])
      MeanProfileFileName = DstDir+'MeanProfile'+BaseName+'.txt'
      DstFileName=DstDir+BaseName+'_spincube.fits'
      IContFileName=DstDir+'ICont_'+BaseName+'.fits'
    endelse

    dim      = size(im)
    width    = dim[1];
    height   = dim[2];
    NWL      = dim[3] ; number of wavelength points
    NrStokes = dim[4] ; number of Stokes parameters
    if NWL NE (size(Wavelength))[1] then begin
      print,'Error: NWL != (size(Wavelength))[1] (',NWL,' != ',(size(Wavelength))[1],')!'
      stop
    endif

    ; calculate the spatial mean profile for each Stokes parameter in the full FOV
    MeanProfile = total(total(im,1),1)/(width*height)
    ;  plot,MeanProfile[*,0]


    ;open txt file for writing
    print,'Writing the mean profile'
    openw,lun,MeanProfileFileName,/GET_LUN

    ;write header
    printf,lun,'Gen.NrChans                             5'
    printf,lun,'Chan.Calibration                        1                 1                 1                 1                 1'
    printf,lun,'Chan.ChannelDataType                    ''float64''         ''float64''         ''float64''         ''float64''         ''float64'''
    printf,lun,'Chan.ChannelName                        ''Wavelength''      ''Stokes I''        ''Stokes Q''        ''Stokes U''        ''Stokes V'''
    printf,lun,'Chan.ChannelNumber                      0                 1                 2                 3                 4'
    printf,lun,'Chan.ChannelOffset                      0                 0                 0                 0                 0'
    printf,lun,'Chan.ChannelSampleFreq                  1                 1                 1                 1                 1'
    printf,lun,'Chan.PhysQuantity                       ''Wavelength''      ''Stokes I''        ''Stokes Q''        ''Stokes U''        ''Stokes V'''
    printf,lun,'Chan.Unit.Ampere                        0                 0                 0                 0                 0'
    printf,lun,'Chan.Unit.Candela                       0                 0                 0                 0                 0'
    printf,lun,'Chan.Unit.Kelvin                        0                 0                 0                 0                 0'
    printf,lun,'Chan.Unit.Kilogram                      0                 0                 0                 0                 0'
    printf,lun,'Chan.Unit.LogRef                        1                 1                 1                 1                 1'
    printf,lun,'Chan.Unit.LogScale                      20                20                20                20                20'
    printf,lun,'Chan.Unit.Meter                         1                 0                 0                 0                 0'
    printf,lun,'Chan.Unit.Name                          ''\xc5''            ''''                ''''                ''''                '''''
    printf,lun,'Chan.Unit.Offset                        0                 0                 0                 0                 0'
    printf,lun,'Chan.Unit.Radian                        0                 0                 0                 0                 0'
    printf,lun,'Chan.Unit.Scale                         1e+10             1                 1                 1                 1'
    printf,lun,'Chan.Unit.Second                        0                 0                 0                 0                 0'

    for l=0,NWL-1 do begin ;loop over all wavelength points
      printf,lun,FORMAT='(E20.10, A, E20.10, A, E20.10, A, E20.10, A, E20.10)',Wavelength[l]/1.0e10,' ',MeanProfile[l,0],' ',MeanProfile[l,1],' ',MeanProfile[l,2],' ',MeanProfile[l,3]
    endfor

    ;close txt file
    free_lun,lun

    print,'Reformating and cropping data'
    cube = fltarr(NWL,NrStokes,ROI_X2-ROI_X1+1,ROI_Y2-ROI_Y1+1)
    for l=0,NWL-1 do begin ;loop over all wavelength points
      for s=0,NrStokes-1 do begin ;loop over all Stokes parameters
        cube[l,s,*,*] = im[ROI_X1:ROI_X2,ROI_Y1:ROI_Y2,l,s]
      end
    end

    ;   ; do the straylight correction for all 4 Stokes parameters
    ;   print,'Correcting for straylight'
    ;   for l=0,NWL-1 do begin ;loop over all wavelength points
    ;      for s=0,NrStokes-1 do begin ;loop over all Stokes parameters
    ;         cube[l,s,*,*] -= (MeanProfile[l,s]*factor)
    ;      endfor
    ;   endfor

    ; do the straylight correction only for Stokes I because in 2013 the IMaX internal stray light was depolarized
    print,'Correcting for straylight'
    for l=0,NWL-1 do begin ;loop over all wavelength points
      cube[l,0,*,*] -= (MeanProfile[l,0]*factor)
    endfor

    ; determine the Stokes I continuum image
    icont = reform(cube[ContIdx,0,*,*])

    ; write output file
    wavel = double(wlref) + Wavelength ; apply WL offset

    sxaddpar,h2,'NSR',1
    sxaddpar,h2,'WLI1',1,FORMAT='(I)'

    print,'Writing data to file ',DstFileName
    ftswr,cube,DstFileName,h2,/create
    ftswr,wavel,DstFileName,h2
    ftswr,icont,DstFileName,h2

    writefits,IContFileName,icont
  end ; for

end ; procedure
