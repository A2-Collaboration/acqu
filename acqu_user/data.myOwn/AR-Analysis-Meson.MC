#------------------------------------------------------------------------------
#                 Configuration file for TA2SaschaPhysics
#------------------------------------------------------------------------------
#
Reconstruct:
#
# CaLib configuration
#
# General:      dbhost                  dbname   dbuser      dbpass        calibration
Use-CaLib:      a2calib.online.a2.kph   CaLib_calib    calib       basel10       ETAP_Aug_12
#Use-CaLib:      a2calib.online.a2.kph   2004_01-2012_08    calib       basel10       2012_08_EtaP_Ott
#
# miscellaneous calibrations    # target position  (0=off/1=on)
Use-CaLib-Misc:                 0
# tagger calib:     time  (0=off/1=on)
Use-CaLib-TAGG:     1
# CB calib:         energy  time    timewalk    quadr. energy  (0=off/1=on)
Use-CaLib-CB:       1       1       0           1
# TAPS calib:       energy  time    quadr. energy   LED thr.  (0=off/1=on)
Use-CaLib-TAPS:     0       1       0               0
# PID calib:        phi     droop   energy  time  (0=off/1=on)
Use-CaLib-PID:      1       0       0       1
# Veto calib:       energy  time  (0=off/1=on)
Use-CaLib-Veto:     0       1
#
#Use-CaLib-CBEnergyPerRun:       data/CBEnergyCorrectionPerRun.txt
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Target information: position (x,y,z) and mass in MeV
Target: 0. 0. 0. 938.272
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Maximum number of detected particles to be accepted from different apparati.
#This defines sizes of internal arrays; so don't select these numbers too
#small...
#             Tagger  Ball  TAPS
Permutations: 128     32    16
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Standard AcquRoot option
ParticleID: 176
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Enable Trigger decoding for multiplicities M1 to M4
#This decodes either hardware information from ADCs 0 and 1, or
#evaluates software trigger implementation, if acitvated.
#        M1  M2  M3  M4
Trigger:  0   1   1   1
#
#Software implementation/simulation of L1 and L2 trigger, options are:
#1.       CB energy sum threshold in MeV (L1 trigger)
#2.       CB energy sum threshold deviation in MeV (L1 trigger)
#3.       Number of inner TAPS rings not contributing to multiplicity,
#         0 accepts all,
#         1 disables the 1st inner ring (6x4 PbWO4s)
#         2 disables also the 2nd inner ring (12x4 PbWO4s)
#         3 and higher completely disables TAPS contribution
#4. to 7. Prescaling factors for M1-4 events, i.e. skipped events
#         between two accepted events (L2 trigger)
#         0     accepts all events of this multiplicity,
#         65535 supresses this multiplicity as the hardware does
#4.       M1 prescale
#5.       M2 prescale
#6.       M3 prescale
#7.       M4 prescale
#8.       TAPS crystals per block (64 or 73)
#               1.       2.      3.     4.      5.     6.     7.     8.
TriggerParam:   500.0    1.0     2      65535   0      0      0	     73
#
#Text files describing the single crystal behaviour fo L1 and L2 triggers:
#ThresBaF2: points to a text file containing the LED discriminator thresholds
#           for each BaF2 crystal (needed for L2 trigger)
#ThresNaI:  points to a text file containing the discriminator thresholds
#           (high threshold) for each NaI(Tl) crystal (needed for L2 trigger)
#ThresNaI:  points to a text file containing the relative energy calibration
#           for each NaI(Tl) crystal (needed for L1 trigger)
*ThresBaF2: data/TAPS_BaF2_Thres.txt
*ThresNaI:  data/CB_NaI_Thres.txt
*GainsNaI:  data/CB_NaI_Gain.txt
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Times for the prompt and random windows used in the analysis
#
Prompt: -7. 6.
Random: -40. -10. 10. 40.
#
# Save Prompt-Random-Ratio as histogram for later subtraction
# Maybe it's better to place this keyword after Initialise:
#  but normally this should be catched within TA2DataManager
#
Display: 1D promptRandomRatio 1000 0. 1.
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# Print some additional information for debugging purposes
#Debug:
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# Display speed information of the analysis
SpeedInfo:
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# File to read TCutG from used for energy and momentum balance cut
#           path to file                                   cut name
BalanceCut: /home/wagners/acqu/acqu_user/data/myCuts.root  balanceCut
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Options to control single meson analysis (cut windows, etc.):
#- Invariant mass (2 photons) range to be accepted as meson (pi0/eta)
#- Missing mass (from beam, proton target and reconstructed meson) range
#  to be accepted as recoil nucleon
#- Beam energy threshold for meson production reaction. Do not use the
#  exact value but stay a bit below, because of tagger energy resolution
#- Particle ID code for reconstructed meson (111=pi0, 221=eta)
#Example pi0->2gamma analysis:
*Misc: InvariantMass 115.0 150.0
*Misc: MissingMass 900.0 980.0
*Misc: Threshold 140.0
*Misc: MesonID 111
#Example eta->2gamma analysis:
*Misc: InvariantMass 507.0 587.0
*Misc: MissingMass 900.0 980.0
*Misc: Threshold 705.0
*Misc: MesonID 221
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Enable output of TA2Particles in ROOT tree file
*Misc: OutputFile A2Events.root
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Tagger prompt/random handling: Text file containing time window definitions
#for prompt and random hits. These times are valid for the EventTime histogram
#(see below) which gives the time difference between tagged photon time and
#detected time in CB/TAPS
*Misc: TimeWindows data/Meson_Cut_Time.txt
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Standard AcquRoot option
# Place keyword AFTER all the config and misc stuff to execute SetConfig after
# all other config data was processed
#
Initialise:
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Trigger spectra: CB energy sum (L1), CB/TAPS combined multiplicity (L2),
#                 energy in each of 45 CB blocks (16 crystals each)
Display: 1D ESum 750 0 1500
Display: 1D Mult 20 0 20
*Display: 1D EDisc 200 0 200 0 44
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Event time as calculated from timr difference between tagged photon and
#detected photons in CB/TAPS
Display: 1D EventTime 600 -150 150
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Invariant mass of 2 photons
*Display: 1D MesonMinv   800 0 800
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Missing mass from beam, proton target and reconstructed meson
#W1 and W2 indicate different (prompt/random) time windows
*Display: 1D MesonMmissW1  1000 0 2000
*Display: 1D MesonMmissW2  1000 0 2000
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Some properties (polar angle, cluster size, central detector, veto detector,
#veto energy) of photons and protons detected in CB and TAPS. Demonstrates
#some possibilities of TA2Particle
*Display: 1D PhotonThetaCB 180 0 180 0 1
*Display: 1D ProtonThetaCB 180 0 180
*Display: 1D PhotonThetaTAPS 180 0 180 0 1
*Display: 1D ProtonThetaTAPS 180 0 180
*Display: 1D PhotonSizeCB 30 0 30 0 1
*Display: 1D ProtonSizeCB 30 0 30
*Display: 1D PhotonSizeTAPS 30 0 30 0 1
*Display: 1D ProtonSizeTAPS 30 0 30
*Display: 1D PhotonCentralCB 720 0 720 0 1
*Display: 1D ProtonCentralCB 720 0 720
*Display: 1D PhotonCentralTAPS 384 0 384 0 1
*Display: 1D ProtonCentralTAPS 384 0 384
*Display: 1D ProtonVetoCB 24 0 24
*Display: 1D ProtonVetoTAPS 384 0 384
*Display: 1D ProtonDeltaECB 100 0 10
*Display: 1D ProtonDeltaETAPS 100 0 10
#
*Display: 2D PhotonEnergyCB 400 0 800  0 PhotonSizeCB 30 0 30  0
*Display: 2D PhotonEnergyCB 400 0 800  1 PhotonSizeCB 30 0 30  1
*Display: 2D ProtonEnergyCB 400 0 800 -1 ProtonSizeCB 30 0 30 -1
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Some distributions for reconstructed g p -> p meson reactions in CM frame
#(proton & meson energies, azimuthal angles, polar angles, total CM energy W)
#W1 and W2 indicate different (prompt/random) time windows
*Display: 1D CMEProtonW1  2000 0 2000
*Display: 1D CMEProtonW2  2000 0 2000
*Display: 1D CMThetaProtonW1  180 0 180
*Display: 1D CMThetaProtonW2  180 0 180
*Display: 1D CMPhiProtonW1  360 -180 180
*Display: 1D CMPhiProtonW2  360 -180 180
*Display: 1D CMEMesonW1  800 0 800
*Display: 1D CMEMesonW2  800 0 800
*Display: 1D CMThetaMesonW1  180 0 180
*Display: 1D CMThetaMesonW2  180 0 180
*Display: 1D CMPhiMesonW1  360 -180 180
*Display: 1D CMPhiMesonW2  360 -180 180
*Display: 1D CMETotalW1  1000 1000 2000
*Display: 1D CMETotalW2  1000 1000 2000
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Some distributions for reconstructed g p -> p meson reactions in Lab frame
#(proton & meson energies, azimuthal angles, polar angles, beam energy)
#W1 and W2 indicate different (prompt/random) time windows
*Display: 1D LabEProtonW1  2000 0 2000
*Display: 1D LabEProtonW2  2000 0 2000
*Display: 1D LabThetaProtonW1  180 0 180
*Display: 1D LabThetaProtonW2  180 0 180
*Display: 1D LabPhiProtonW1  360 -180 180
*Display: 1D LabPhiProtonW2  360 -180 180
*Display: 1D LabEMesonW1  800 0 800
*Display: 1D LabEMesonW2  800 0 800
*Display: 1D LabThetaMesonW1  180 0 180
*Display: 1D LabThetaMesonW2  180 0 180
*Display: 1D LabPhiMesonW1  360 -180 180
*Display: 1D LabPhiMesonW2  360 -180 180
*Display: 1D LabEBeamW1  1500 0 1500
*Display: 1D LabEBeamW2  1500 0 1500
*Display: 1D LabChBeamW1  352 0 352
*Display: 1D LabChBeamW2  352 0 352
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Two-dimensional plot of meson CM polar angle vs Lab beam energy
*Display: 2D CMThetaMesonW1 90 0 180 -1 LabEBeamW1 150 0 1500 -1
*Display: 2D CMThetaMesonW2 90 0 180 -1 LabEBeamW2 150 0 1500 -1
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Spectra for true information from Pluto
Display: 1D trueNPart 20 0 20
Display: 1D trueIDPart 20 0 20
Display: 1D trueNProt 10 0 10
Display: 1D trueProtMass 10 930 940
Display: 1D trueProtEnergy 300 0 1200
*Display: 1D trueProtTheta 90 0 180
Display: 1D trueProtTheta 100 0 60
Display: 1D trueProtPhi 180 -180 180
Display: 1D trueNGamma 10 0 10
Display: 1D trueGammaMass 10 -0.5 0.5
Display: 1D trueGammaEnergy 700 0 1400
Display: 1D trueGammaTheta 90 0 180
Display: 1D trueGammaPhi 180 -180 180
Display: 1D trueNPosi 10 0 10
Display: 1D truePosiMass 20 0 1
Display: 1D truePosiEnergy 500 0 1000
Display: 1D truePosiTheta 90 0 180
Display: 1D truePosiPhi 180 -180 180
Display: 1D trueNElec 10 0 10
Display: 1D trueElecMass 20 0 1
Display: 1D trueElecEnergy 500 0 1000
Display: 1D trueElecTheta 90 0 180
Display: 1D trueElecPhi 180 -180 180
Display: 1D trueNMuPls 10 0 10
*Display: 1D trueMuPlsMass 20 105 106
Display: 1D trueMuPlsMass 400 0 200
Display: 1D trueMuPlsEnergy 500 0 1000
Display: 1D trueMuPlsTheta 90 0 180
Display: 1D trueMuPlsPhi 180 -180 180
Display: 1D trueNMuMns 10 0 10
*Display: 1D trueMuMnsMass 20 105 106
Display: 1D trueMuMnsMass 400 0 200
Display: 1D trueMuMnsEnergy 500 0 1000
Display: 1D trueMuMnsTheta 90 0 180
Display: 1D trueMuMnsPhi 180 -180 180
Display: 1D trueNPiPls 10 0 10
*Display: 1D truePiPlsMass 20 139 140
Display: 1D truePiPlsMass 500 0 250
Display: 1D truePiPlsEnergy 500 0 1000
Display: 1D truePiPlsTheta 90 0 180
Display: 1D truePiPlsPhi 180 -180 180
Display: 1D trueNPiMns 10 0 10
*Display: 1D truePiMnsMass 20 139 140
Display: 1D truePiMnsMass 500 0 250
Display: 1D truePiMnsEnergy 500 0 1000
Display: 1D truePiMnsTheta 90 0 180
Display: 1D truePiMnsPhi 180 -180 180

Display: 1D trueInvM 2400 0 1200

Display: 1D invM_allPart 1600 0 1600
Display: 1D invM_2charged 1600 0 1600
Display: 1D timeProton 500 -250 250
Display: 1D timePhoton 500 -250 250
Display: 1D timeFS 500 -250 250
Display: 1D timeLeptons 500 -250 250
Display: 1D timeTagger 500 -250 250
Display: 1D dTimeFS 500 -250 250
Display: 1D coplanarity 360 0 360
Display: 1D protEnergyReconstr 800 0 800

Display: 1D invMass2g 1200 0 1200
Display: 1D invMass2g1p 1200 0 1200
Display: 1D invMass6g 1200 0 1200
Display: 1D invMass6g1p 1200 0 1200
Display: 1D invMass2CB 1200 0 1200
Display: 1D invMass2CB1TAPS 1200 0 1200
Display: 1D invMass6CB 1200 0 1200
Display: 1D invMass6CB1TAPS 1200 0 1200

# Prompt
Display: 1D invM_protE_prompt 1600 0 1600
Display: 1D invM_copl_prompt 1600 0 1600
Display: 1D invM_balance_prompt 1600 0 1600
Display: 1D invM_dAlpha_prompt 1600 0 1600
Display: 1D invM_missM_prompt 1600 0 1600
Display: 1D invM_copl_balance_prompt 1600 0 1600
Display: 1D invM_balance_missM_prompt 1600 0 1600
Display: 1D invM_protE_copl_prompt 1600 0 1600
Display: 1D invM_copl_missM_prompt 1600 0 1600
Display: 1D invM_balance_dAlpha_prompt 1600 0 1600
Display: 1D invM_copl_dAlpha_prompt 1600 0 1600
Display: 1D invM_protE_copl_balance_prompt 1600 0 1600
Display: 1D invM_copl_balance_dAlpha_prompt 1600 0 1600
Display: 1D invM_allCuts_prompt 1600 0 1600
Display: 1D missM_protE_prompt 1600 200 1800
Display: 1D missM_copl_prompt 1600 200 1800
Display: 1D missM_balance_prompt 1600 200 1800
Display: 1D missM_dAlpha_prompt 1600 200 1800
Display: 1D missM_invM_prompt 1600 200 1800
Display: 1D missM_copl_balance_prompt 1600 200 1800
Display: 1D missM_protE_copl_prompt 1600 200 1800
Display: 1D missM_balance_dAlpha_prompt 1600 200 1800
Display: 1D missM_copl_dAlpha_prompt 1600 200 1800
Display: 1D missM_protE_copl_balance_prompt 1600 200 1800
Display: 1D missM_copl_balance_dAlpha_prompt 1600 200 1800
Display: 1D missM_allCuts_prompt 1600 200 1800
Display: 1D balancePx_prompt 4000 -2000 2000
Display: 1D balancePy_prompt 4000 -2000 2000
Display: 1D balancePz_prompt 4000 -2000 2000
Display: 1D balanceE_prompt 4000 -2000 2000
Display: 1D missM_prompt 1000 500 1500
Display: 1D protonEnergyExpect_prompt 800 0 800
Display: 1D protDAlphaTAPSCl_prompt 180 0 90
# Random
Display: 1D invM_protE_random 1600 0 1600
Display: 1D invM_copl_random 1600 0 1600
Display: 1D invM_balance_random 1600 0 1600
Display: 1D invM_dAlpha_random 1600 0 1600
Display: 1D invM_missM_random 1600 0 1600
Display: 1D invM_copl_balance_random 1600 0 1600
Display: 1D invM_balance_missM_random 1600 0 1600
Display: 1D invM_protE_copl_random 1600 0 1600
Display: 1D invM_copl_missM_random 1600 0 1600
Display: 1D invM_balance_dAlpha_random 1600 0 1600
Display: 1D invM_copl_dAlpha_random 1600 0 1600
Display: 1D invM_protE_copl_balance_random 1600 0 1600
Display: 1D invM_copl_balance_dAlpha_random 1600 0 1600
Display: 1D invM_allCuts_random 1600 0 1600
Display: 1D missM_protE_random 1600 200 1800
Display: 1D missM_copl_random 1600 200 1800
Display: 1D missM_balance_random 1600 200 1800
Display: 1D missM_dAlpha_random 1600 200 1800
Display: 1D missM_invM_random 1600 200 1800
Display: 1D missM_copl_balance_random 1600 200 1800
Display: 1D missM_protE_copl_random 1600 200 1800
Display: 1D missM_balance_dAlpha_random 1600 200 1800
Display: 1D missM_copl_dAlpha_random 1600 200 1800
Display: 1D missM_protE_copl_balance_random 1600 200 1800
Display: 1D missM_copl_balance_dAlpha_random 1600 200 1800
Display: 1D missM_allCuts_random 1600 200 1800
Display: 1D balancePx_random 4000 -2000 2000
Display: 1D balancePy_random 4000 -2000 2000
Display: 1D balancePz_random 4000 -2000 2000
Display: 1D balanceE_random 4000 -2000 2000
Display: 1D missM_random 1000 500 1500
Display: 1D protonEnergyExpect_random 800 0 800
Display: 1D protDAlphaTAPSCl_random 180 0 90
