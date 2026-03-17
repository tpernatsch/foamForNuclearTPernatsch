
% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2024 15:50:46' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Fuel Pin' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1:  4])  = 'main' ;
WORKING_DIRECTORY         (idx, [1:103])  = '/home/thomas-guilbaud/Software/foamForNuclear/pythonapi/tutorials/tests/test_fuelPin_monteCarlo/Serpent' ;
HOSTNAME                  (idx, [1:  7])  = 'TS-P350' ;
CPU_TYPE                  (idx, [1: 46])  = '11th Gen Intel(R) Core(TM) i9-11900K @ 3.50GHz' ;
CPU_MHZ                   (idx, 1)        = 98.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Dec  2 16:37:57 2025' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Dec  2 16:40:08 2025' ;

% Run parameters:

POP                       (idx, 1)        = 10000 ;
CYCLES                    (idx, 1)        = 1000 ;
SKIP                      (idx, 1)        = 100 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1764689877254 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 15 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  15]) = [  1.01577E+00  1.02975E+00  1.00301E+00  1.01183E+00  9.99208E-01  1.00158E+00  9.83312E-01  9.89737E-01  9.90049E-01  9.74149E-01  9.93896E-01  9.94224E-01  1.00857E+00  9.92815E-01  1.01211E+00  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 71])  = '/home/thomas-guilbaud/Software/Serpent/xsdata/endfb7/sss_endfb7u.xsdata' ;
DECAY_DATA_FILE_PATH      (idx, [1:  3])  = 'N/A' ;
SFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
NFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  1.04694E-02 0.00070  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.89531E-01 7.5E-06  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  7.51604E-01 3.8E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  7.52090E-01 3.8E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.76776E+00 0.00014  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  3.19299E+01 0.00019  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  3.19299E+01 0.00019  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  1.05251E+01 0.00025  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  3.65280E-01 0.00077  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 1000 ;
SIMULATED_HISTORIES       (idx, 1)        = 10001055 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  1.00011E+04 0.00051 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  1.00011E+04 0.00051 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  2.67639E+01 ;
RUNNING_TIME              (idx, 1)        =  2.19012E+00 ;
INIT_TIME                 (idx, [1:  2])  = [  1.19833E-02  1.19833E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.33332E-04  2.33332E-04 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  2.17790E+00  2.17790E+00  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  2.18980E+00  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 12.22029 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  1.31870E+01 0.00622 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  8.79040E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 63913.86 ;
ALLOC_MEMSIZE             (idx, 1)        = 408.02;
MEMSIZE                   (idx, 1)        = 264.81;
XS_MEMSIZE                (idx, 1)        = 145.32;
MAT_MEMSIZE               (idx, 1)        = 14.72;
RES_MEMSIZE               (idx, 1)        = 4.27;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 100.50;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 143.20;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 7 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 128116 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.00000E+37 ;
URES_EMAX                 (idx, 1)        = -1.00000E+37 ;
URES_AVAIL                (idx, 1)        = 7 ;
URES_USED                 (idx, 1)        = 0 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 13 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 13 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 0 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 347 ;
TOT_TRANSMU_REA           (idx, 1)        = 0 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 0 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 1 ;
TMS_MODE                  (idx, 1)        = 0 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  0.00000E+00 ;
TOT_DECAY_HEAT            (idx, 1)        =  0.00000E+00 ;
TOT_SF_RATE               (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  0.00000E+00 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  0.00000E+00 ;
INHALATION_TOXICITY       (idx, 1)        =  0.00000E+00 ;
INGESTION_TOXICITY        (idx, 1)        =  0.00000E+00 ;
ACTINIDE_INH_TOX          (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ING_TOX          (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  0.00000E+00 ;
SR90_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
TE132_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
I131_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
I132_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
CS134_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
CS137_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  0.00000E+00 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  0.00000E+00 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  0.00000E+00 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  0.00000E+00 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  9.98582E-05 0.00027  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  5.39002E-01 0.00076 ];
U235_FISS                 (idx, [1:   4]) = [  4.35748E-01 0.00041  9.38022E-01 0.00012 ];
U238_FISS                 (idx, [1:   4]) = [  2.87640E-02 0.00187  6.19124E-02 0.00176 ];
U235_CAPT                 (idx, [1:   4]) = [  9.27791E-02 0.00102  1.73244E-01 0.00095 ];
U238_CAPT                 (idx, [1:   4]) = [  2.83368E-01 0.00061  5.29109E-01 0.00042 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 10001055 1.00000E+07 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 1.51611E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 10001055 1.00152E+07 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 5355383 5.36308E+06 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 4645672 4.65208E+06 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 10001055 1.00152E+07 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -1.21072E-07 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.50775E-11 0.00014 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  1.14288E+00 0.00014 ];
TOT_FISSRATE              (idx, [1:   2]) = [  4.64575E-01 0.00014 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  5.35425E-01 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.98582E-01 0.00027 ];
TOT_FLUX                  (idx, [1:   2]) = [  4.51379E+01 0.00022 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  3.19307E+01 0.00017 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.76384E+00 0.00028 ];
SIX_FF_F                  (idx, [1:   2]) = [  7.95396E-01 0.00020 ];
SIX_FF_P                  (idx, [1:   2]) = [  6.61915E-01 0.00023 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  1.23249E+00 0.00023 ];
SIX_FF_LF                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_LT                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.14448E+00 0.00035 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.14448E+00 0.00035 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.46004E+00 8.7E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02564E+02 9.1E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.14440E+00 0.00036  1.13662E+00 0.00035  7.85261E-03 0.00558 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.14461E+00 0.00014 ];
COL_KEFF                  (idx, [1:   2]) = [  1.14459E+00 0.00033 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.14461E+00 0.00014 ];
ABS_KINF                  (idx, [1:   2]) = [  1.14461E+00 0.00014 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.76226E+01 0.00013 ];
IMP_ALF                   (idx, [1:   2]) = [  1.76240E+01 5.4E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  4.45383E-07 0.00224 ];
IMP_EALF                  (idx, [1:   2]) = [  4.43859E-07 0.00096 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  2.08096E-01 0.00196 ];
IMP_AFGE                  (idx, [1:   2]) = [  2.07406E-01 0.00081 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.24107E-03 0.00389  1.74110E-04 0.02285  9.89280E-04 0.00910  9.79162E-04 0.00900  2.85504E-03 0.00568  9.36400E-04 0.00939  3.07070E-04 0.01759 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  8.07935E-01 0.00904  1.07795E-02 0.01261  3.16553E-02 0.00015  1.10172E-01 0.00019  3.20588E-01 0.00016  1.34554E+00 0.00012  8.60233E+00 0.00596 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  6.88925E-03 0.00548  1.94780E-04 0.03335  1.08813E-03 0.01373  1.08468E-03 0.01342  3.14323E-03 0.00807  1.03296E-03 0.01466  3.45472E-04 0.02525 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  8.17511E-01 0.01316  1.24908E-02 1.9E-06  3.16565E-02 0.00021  1.10193E-01 0.00028  3.20681E-01 0.00023  1.34554E+00 0.00016  8.89360E+00 0.00145 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  2.46166E-05 0.00072  2.46044E-05 0.00073  2.63206E-05 0.00745 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  2.81677E-05 0.00063  2.81537E-05 0.00063  3.01157E-05 0.00742 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.86045E-03 0.00563  1.92990E-04 0.03374  1.08306E-03 0.01382  1.08385E-03 0.01359  3.12644E-03 0.00834  1.02969E-03 0.01442  3.44413E-04 0.02658 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  8.14833E-01 0.01397  1.24908E-02 2.4E-06  3.16659E-02 0.00024  1.10153E-01 0.00030  3.20560E-01 0.00025  1.34579E+00 0.00018  8.88137E+00 0.00167 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  2.46358E-05 0.00157  2.46255E-05 0.00158  2.63151E-05 0.01735 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  2.81894E-05 0.00153  2.81776E-05 0.00153  3.01110E-05 0.01733 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  6.92380E-03 0.01680  1.85253E-04 0.10796  1.16040E-03 0.04184  1.10877E-03 0.04196  3.13077E-03 0.02501  1.00716E-03 0.04436  3.31450E-04 0.07858 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  7.52734E-01 0.03852  1.24908E-02 5.5E-06  3.16718E-02 0.00053  1.10150E-01 0.00070  3.20789E-01 0.00069  1.34491E+00 0.00044  8.93125E+00 0.00386 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  6.98345E-03 0.01646  1.90092E-04 0.10745  1.17569E-03 0.04111  1.12502E-03 0.04106  3.13617E-03 0.02413  1.01618E-03 0.04297  3.40288E-04 0.07664 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  7.59419E-01 0.03775  1.24908E-02 5.5E-06  3.16721E-02 0.00053  1.10162E-01 0.00071  3.20768E-01 0.00068  1.34499E+00 0.00043  8.93495E+00 0.00387 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -2.81848E+02 0.01690 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  2.46332E-05 0.00047 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  2.81865E-05 0.00030 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  6.91349E-03 0.00326 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -2.80749E+02 0.00333 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.05784E-07 0.00033 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  2.90081E-06 0.00029  2.90087E-06 0.00029  2.89182E-06 0.00352 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  3.12270E-05 0.00039  3.12287E-05 0.00039  3.09638E-05 0.00450 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  6.62311E-01 0.00023  6.61556E-01 0.00023  7.94700E-01 0.00604 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.02411E+01 0.00905 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  3.19299E+01 0.00019  3.51211E+01 0.00026 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '0' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  6.79050E+04 0.00223  2.74411E+05 0.00096  5.72159E+05 0.00058  6.25496E+05 0.00058  5.82568E+05 0.00040  6.34398E+05 0.00048  4.32400E+05 0.00042  3.84657E+05 0.00054  2.94273E+05 0.00053  2.40412E+05 0.00057  2.07493E+05 0.00058  1.87361E+05 0.00061  1.72764E+05 0.00064  1.64248E+05 0.00065  1.59955E+05 0.00060  1.38226E+05 0.00066  1.36673E+05 0.00074  1.35384E+05 0.00063  1.33107E+05 0.00062  2.60105E+05 0.00048  2.50799E+05 0.00055  1.81460E+05 0.00063  1.17515E+05 0.00064  1.35372E+05 0.00070  1.27610E+05 0.00068  1.16452E+05 0.00076  1.89855E+05 0.00054  4.34676E+04 0.00115  5.45192E+04 0.00114  4.94822E+04 0.00118  2.87548E+04 0.00160  5.00248E+04 0.00090  3.40209E+04 0.00133  2.89510E+04 0.00131  5.51939E+03 0.00302  5.49269E+03 0.00242  5.63893E+03 0.00267  5.77512E+03 0.00251  5.72821E+03 0.00267  5.64893E+03 0.00288  5.85103E+03 0.00275  5.47305E+03 0.00226  1.03516E+04 0.00183  1.65086E+04 0.00142  2.10752E+04 0.00159  5.55513E+04 0.00090  5.81809E+04 0.00101  6.27323E+04 0.00102  4.16391E+04 0.00107  3.05667E+04 0.00082  2.32994E+04 0.00123  2.70326E+04 0.00103  5.00414E+04 0.00095  6.58681E+04 0.00075  1.22840E+05 0.00068  1.82777E+05 0.00057  2.64666E+05 0.00056  1.65624E+05 0.00066  1.16947E+05 0.00060  8.32445E+04 0.00064  7.40467E+04 0.00068  7.24501E+04 0.00071  6.01007E+04 0.00066  4.03248E+04 0.00081  3.69708E+04 0.00087  3.26306E+04 0.00112  2.73903E+04 0.00099  2.13541E+04 0.00097  1.41022E+04 0.00128  4.88222E+03 0.00193 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.14459E+00 0.00030 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  3.67523E+01 0.00026  8.38662E+00 0.00022 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.57420E-01 6.4E-05  1.36472E+00 0.00010 ];
INF_CAPT                  (idx, [1:   4]) = [  6.93541E-03 0.00038  3.34539E-02 0.00010 ];
INF_ABS                   (idx, [1:   4]) = [  9.21921E-03 0.00029  7.88450E-02 0.00016 ];
INF_FISS                  (idx, [1:   4]) = [  2.28380E-03 0.00025  4.53911E-02 0.00022 ];
INF_NSF                   (idx, [1:   4]) = [  5.86000E-03 0.00025  1.10604E-01 0.00022 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.56590E+00 3.8E-05  2.43670E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.03898E+02 3.9E-06  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  5.96618E-08 0.00025  2.46083E-06 9.2E-05 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.48193E-01 6.8E-05  1.28591E+00 0.00011 ];
INF_SCATT1                (idx, [1:   4]) = [  2.47927E-01 0.00010  3.35126E-01 0.00025 ];
INF_SCATT2                (idx, [1:   4]) = [  9.75018E-02 0.00020  8.18830E-02 0.00074 ];
INF_SCATT3                (idx, [1:   4]) = [  7.52487E-03 0.00227  2.48471E-02 0.00155 ];
INF_SCATT4                (idx, [1:   4]) = [ -1.03134E-02 0.00152 -6.47429E-03 0.00690 ];
INF_SCATT5                (idx, [1:   4]) = [  2.58664E-04 0.04202  5.13894E-03 0.00818 ];
INF_SCATT6                (idx, [1:   4]) = [  5.16060E-03 0.00218 -1.33239E-02 0.00312 ];
INF_SCATT7                (idx, [1:   4]) = [  7.38704E-04 0.01341 -2.42958E-05 1.00000 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.48234E-01 6.8E-05  1.28591E+00 0.00011 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.47927E-01 0.00010  3.35126E-01 0.00025 ];
INF_SCATTP2               (idx, [1:   4]) = [  9.75021E-02 0.00020  8.18830E-02 0.00074 ];
INF_SCATTP3               (idx, [1:   4]) = [  7.52485E-03 0.00227  2.48471E-02 0.00155 ];
INF_SCATTP4               (idx, [1:   4]) = [ -1.03135E-02 0.00152 -6.47429E-03 0.00690 ];
INF_SCATTP5               (idx, [1:   4]) = [  2.58758E-04 0.04195  5.13894E-03 0.00818 ];
INF_SCATTP6               (idx, [1:   4]) = [  5.16053E-03 0.00218 -1.33239E-02 0.00312 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.38648E-04 0.01341 -2.42958E-05 1.00000 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.28155E-01 0.00021  9.09044E-01 0.00014 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.46100E+00 0.00021  3.66686E-01 0.00014 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  9.17801E-03 0.00029  7.88450E-02 0.00016 ];
INF_REMXS                 (idx, [1:   4]) = [  2.75058E-02 0.00015  8.01014E-02 0.00034 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  0.00000E+00 ;
PM148_LAMBDA              (idx, 1)        =  0.00000E+00 ;
PM148M_LAMBDA             (idx, 1)        =  0.00000E+00 ;
PM149_LAMBDA              (idx, 1)        =  0.00000E+00 ;
I135_LAMBDA               (idx, 1)        =  0.00000E+00 ;
XE135_LAMBDA              (idx, 1)        =  0.00000E+00 ;
XE135M_LAMBDA             (idx, 1)        =  0.00000E+00 ;
I135_BR                   (idx, 1)        =  0.00000E+00 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.29914E-01 6.5E-05  1.82786E-02 0.00022  1.28801E-03 0.00378  1.28462E+00 0.00011 ];
INF_S1                    (idx, [1:   8]) = [  2.42603E-01 0.00010  5.32366E-03 0.00062  6.07080E-04 0.00535  3.34519E-01 0.00025 ];
INF_S2                    (idx, [1:   8]) = [  9.90917E-02 0.00019 -1.58993E-03 0.00157  3.23735E-04 0.00753  8.15593E-02 0.00074 ];
INF_S3                    (idx, [1:   8]) = [  9.40773E-03 0.00179 -1.88286E-03 0.00144  1.17271E-04 0.01453  2.47298E-02 0.00156 ];
INF_S4                    (idx, [1:   8]) = [ -9.70039E-03 0.00161 -6.13036E-04 0.00373  9.65419E-07 1.00000 -6.47525E-03 0.00688 ];
INF_S5                    (idx, [1:   8]) = [  2.29515E-04 0.04743  2.91486E-05 0.06326 -4.82435E-05 0.02695  5.18718E-03 0.00815 ];
INF_S6                    (idx, [1:   8]) = [  5.30683E-03 0.00209 -1.46237E-04 0.01353 -5.78311E-05 0.02035 -1.32661E-02 0.00315 ];
INF_S7                    (idx, [1:   8]) = [  9.18489E-04 0.01063 -1.79785E-04 0.01008 -5.29896E-05 0.01896  2.86938E-05 1.00000 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.29955E-01 6.5E-05  1.82786E-02 0.00022  1.28801E-03 0.00378  1.28462E+00 0.00011 ];
INF_SP1                   (idx, [1:   8]) = [  2.42603E-01 0.00010  5.32366E-03 0.00062  6.07080E-04 0.00535  3.34519E-01 0.00025 ];
INF_SP2                   (idx, [1:   8]) = [  9.90921E-02 0.00019 -1.58993E-03 0.00157  3.23735E-04 0.00753  8.15593E-02 0.00074 ];
INF_SP3                   (idx, [1:   8]) = [  9.40771E-03 0.00179 -1.88286E-03 0.00144  1.17271E-04 0.01453  2.47298E-02 0.00156 ];
INF_SP4                   (idx, [1:   8]) = [ -9.70047E-03 0.00161 -6.13036E-04 0.00373  9.65419E-07 1.00000 -6.47525E-03 0.00688 ];
INF_SP5                   (idx, [1:   8]) = [  2.29610E-04 0.04735  2.91486E-05 0.06326 -4.82435E-05 0.02695  5.18718E-03 0.00815 ];
INF_SP6                   (idx, [1:   8]) = [  5.30677E-03 0.00209 -1.46237E-04 0.01353 -5.78311E-05 0.02035 -1.32661E-02 0.00315 ];
INF_SP7                   (idx, [1:   8]) = [  9.18433E-04 0.01064 -1.79785E-04 0.01008 -5.29896E-05 0.01896  2.86938E-05 1.00000 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  2.41510E-01 0.00038  8.27022E-01 0.00140 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  2.41571E-01 0.00061  8.29876E-01 0.00246 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  2.41456E-01 0.00054  8.30401E-01 0.00270 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  2.41510E-01 0.00070  8.21370E-01 0.00233 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  1.38022E+00 0.00038  4.03091E-01 0.00140 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  1.37988E+00 0.00060  4.01785E-01 0.00245 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  1.38053E+00 0.00054  4.01554E-01 0.00267 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  1.38024E+00 0.00070  4.05934E-01 0.00234 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  6.88925E-03 0.00548  1.94780E-04 0.03335  1.08813E-03 0.01373  1.08468E-03 0.01342  3.14323E-03 0.00807  1.03296E-03 0.01466  3.45472E-04 0.02525 ];
LAMBDA                    (idx, [1:  14]) = [  8.17511E-01 0.01316  1.24908E-02 1.9E-06  3.16565E-02 0.00021  1.10193E-01 0.00028  3.20681E-01 0.00023  1.34554E+00 0.00016  8.89360E+00 0.00145 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2024 15:50:46' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Fuel Pin' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1:  4])  = 'main' ;
WORKING_DIRECTORY         (idx, [1:103])  = '/home/thomas-guilbaud/Software/foamForNuclear/pythonapi/tutorials/tests/test_fuelPin_monteCarlo/Serpent' ;
HOSTNAME                  (idx, [1:  7])  = 'TS-P350' ;
CPU_TYPE                  (idx, [1: 46])  = '11th Gen Intel(R) Core(TM) i9-11900K @ 3.50GHz' ;
CPU_MHZ                   (idx, 1)        = 98.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Dec  2 16:37:57 2025' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Dec  2 16:40:08 2025' ;

% Run parameters:

POP                       (idx, 1)        = 10000 ;
CYCLES                    (idx, 1)        = 1000 ;
SKIP                      (idx, 1)        = 100 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1764689877254 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 15 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  15]) = [  1.01577E+00  1.02975E+00  1.00301E+00  1.01183E+00  9.99208E-01  1.00158E+00  9.83312E-01  9.89737E-01  9.90049E-01  9.74149E-01  9.93896E-01  9.94224E-01  1.00857E+00  9.92815E-01  1.01211E+00  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 71])  = '/home/thomas-guilbaud/Software/Serpent/xsdata/endfb7/sss_endfb7u.xsdata' ;
DECAY_DATA_FILE_PATH      (idx, [1:  3])  = 'N/A' ;
SFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
NFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  1.04694E-02 0.00070  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.89531E-01 7.5E-06  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  7.51604E-01 3.8E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  7.52090E-01 3.8E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.76776E+00 0.00014  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  3.19299E+01 0.00019  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  3.19299E+01 0.00019  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  1.05251E+01 0.00025  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  3.65280E-01 0.00077  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 1000 ;
SIMULATED_HISTORIES       (idx, 1)        = 10001055 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  1.00011E+04 0.00051 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  1.00011E+04 0.00051 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  2.67642E+01 ;
RUNNING_TIME              (idx, 1)        =  2.19015E+00 ;
INIT_TIME                 (idx, [1:  2])  = [  1.19833E-02  1.19833E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.33332E-04  2.33332E-04 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  2.17790E+00  2.17790E+00  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  2.18980E+00  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 12.22024 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  1.31870E+01 0.00622 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  8.79027E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 63913.86 ;
ALLOC_MEMSIZE             (idx, 1)        = 408.02;
MEMSIZE                   (idx, 1)        = 264.81;
XS_MEMSIZE                (idx, 1)        = 145.32;
MAT_MEMSIZE               (idx, 1)        = 14.72;
RES_MEMSIZE               (idx, 1)        = 4.27;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 100.50;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 143.20;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 7 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 128116 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.00000E+37 ;
URES_EMAX                 (idx, 1)        = -1.00000E+37 ;
URES_AVAIL                (idx, 1)        = 7 ;
URES_USED                 (idx, 1)        = 0 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 13 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 13 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 0 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 347 ;
TOT_TRANSMU_REA           (idx, 1)        = 0 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 0 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 1 ;
TMS_MODE                  (idx, 1)        = 0 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  0.00000E+00 ;
TOT_DECAY_HEAT            (idx, 1)        =  0.00000E+00 ;
TOT_SF_RATE               (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  0.00000E+00 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  0.00000E+00 ;
INHALATION_TOXICITY       (idx, 1)        =  0.00000E+00 ;
INGESTION_TOXICITY        (idx, 1)        =  0.00000E+00 ;
ACTINIDE_INH_TOX          (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ING_TOX          (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  0.00000E+00 ;
SR90_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
TE132_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
I131_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
I132_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
CS134_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
CS137_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  0.00000E+00 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  0.00000E+00 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  0.00000E+00 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  0.00000E+00 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  9.98582E-05 0.00027  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  5.39002E-01 0.00076 ];
U235_FISS                 (idx, [1:   4]) = [  4.35748E-01 0.00041  9.38022E-01 0.00012 ];
U238_FISS                 (idx, [1:   4]) = [  2.87640E-02 0.00187  6.19124E-02 0.00176 ];
U235_CAPT                 (idx, [1:   4]) = [  9.27791E-02 0.00102  1.73244E-01 0.00095 ];
U238_CAPT                 (idx, [1:   4]) = [  2.83368E-01 0.00061  5.29109E-01 0.00042 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 10001055 1.00000E+07 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 1.51611E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 10001055 1.00152E+07 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 5355383 5.36308E+06 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 4645672 4.65208E+06 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 10001055 1.00152E+07 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -1.21072E-07 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.50775E-11 0.00014 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  1.14288E+00 0.00014 ];
TOT_FISSRATE              (idx, [1:   2]) = [  4.64575E-01 0.00014 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  5.35425E-01 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.98582E-01 0.00027 ];
TOT_FLUX                  (idx, [1:   2]) = [  4.51379E+01 0.00022 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  3.19307E+01 0.00017 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.76384E+00 0.00028 ];
SIX_FF_F                  (idx, [1:   2]) = [  7.95396E-01 0.00020 ];
SIX_FF_P                  (idx, [1:   2]) = [  6.61915E-01 0.00023 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  1.23249E+00 0.00023 ];
SIX_FF_LF                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_LT                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.14448E+00 0.00035 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.14448E+00 0.00035 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.46004E+00 8.7E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02564E+02 9.1E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.14440E+00 0.00036  1.13662E+00 0.00035  7.85261E-03 0.00558 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.14461E+00 0.00014 ];
COL_KEFF                  (idx, [1:   2]) = [  1.14459E+00 0.00033 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.14461E+00 0.00014 ];
ABS_KINF                  (idx, [1:   2]) = [  1.14461E+00 0.00014 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.76226E+01 0.00013 ];
IMP_ALF                   (idx, [1:   2]) = [  1.76240E+01 5.4E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  4.45383E-07 0.00224 ];
IMP_EALF                  (idx, [1:   2]) = [  4.43859E-07 0.00096 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  2.08096E-01 0.00196 ];
IMP_AFGE                  (idx, [1:   2]) = [  2.07406E-01 0.00081 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.24107E-03 0.00389  1.74110E-04 0.02285  9.89280E-04 0.00910  9.79162E-04 0.00900  2.85504E-03 0.00568  9.36400E-04 0.00939  3.07070E-04 0.01759 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  8.07935E-01 0.00904  1.07795E-02 0.01261  3.16553E-02 0.00015  1.10172E-01 0.00019  3.20588E-01 0.00016  1.34554E+00 0.00012  8.60233E+00 0.00596 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  6.88925E-03 0.00548  1.94780E-04 0.03335  1.08813E-03 0.01373  1.08468E-03 0.01342  3.14323E-03 0.00807  1.03296E-03 0.01466  3.45472E-04 0.02525 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  8.17511E-01 0.01316  1.24908E-02 1.9E-06  3.16565E-02 0.00021  1.10193E-01 0.00028  3.20681E-01 0.00023  1.34554E+00 0.00016  8.89360E+00 0.00145 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  2.46166E-05 0.00072  2.46044E-05 0.00073  2.63206E-05 0.00745 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  2.81677E-05 0.00063  2.81537E-05 0.00063  3.01157E-05 0.00742 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.86045E-03 0.00563  1.92990E-04 0.03374  1.08306E-03 0.01382  1.08385E-03 0.01359  3.12644E-03 0.00834  1.02969E-03 0.01442  3.44413E-04 0.02658 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  8.14833E-01 0.01397  1.24908E-02 2.4E-06  3.16659E-02 0.00024  1.10153E-01 0.00030  3.20560E-01 0.00025  1.34579E+00 0.00018  8.88137E+00 0.00167 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  2.46358E-05 0.00157  2.46255E-05 0.00158  2.63151E-05 0.01735 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  2.81894E-05 0.00153  2.81776E-05 0.00153  3.01110E-05 0.01733 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  6.92380E-03 0.01680  1.85253E-04 0.10796  1.16040E-03 0.04184  1.10877E-03 0.04196  3.13077E-03 0.02501  1.00716E-03 0.04436  3.31450E-04 0.07858 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  7.52734E-01 0.03852  1.24908E-02 5.5E-06  3.16718E-02 0.00053  1.10150E-01 0.00070  3.20789E-01 0.00069  1.34491E+00 0.00044  8.93125E+00 0.00386 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  6.98345E-03 0.01646  1.90092E-04 0.10745  1.17569E-03 0.04111  1.12502E-03 0.04106  3.13617E-03 0.02413  1.01618E-03 0.04297  3.40288E-04 0.07664 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  7.59419E-01 0.03775  1.24908E-02 5.5E-06  3.16721E-02 0.00053  1.10162E-01 0.00071  3.20768E-01 0.00068  1.34499E+00 0.00043  8.93495E+00 0.00387 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -2.81848E+02 0.01690 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  2.46332E-05 0.00047 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  2.81865E-05 0.00030 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  6.91349E-03 0.00326 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -2.80749E+02 0.00333 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.05784E-07 0.00033 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  2.90081E-06 0.00029  2.90087E-06 0.00029  2.89182E-06 0.00352 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  3.12270E-05 0.00039  3.12287E-05 0.00039  3.09638E-05 0.00450 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  6.62311E-01 0.00023  6.61556E-01 0.00023  7.94700E-01 0.00604 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.02411E+01 0.00905 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  3.19299E+01 0.00019  3.51211E+01 0.00026 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  5])  = 'uFuel' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  2.14727E+04 0.00285  8.72833E+04 0.00146  1.82189E+05 0.00096  1.98928E+05 0.00087  1.85496E+05 0.00063  2.02025E+05 0.00081  1.36351E+05 0.00068  1.20961E+05 0.00094  9.13985E+04 0.00103  7.40332E+04 0.00098  6.34439E+04 0.00105  5.71727E+04 0.00119  5.23840E+04 0.00115  4.98514E+04 0.00117  4.83308E+04 0.00120  4.17304E+04 0.00123  4.10716E+04 0.00147  4.06794E+04 0.00134  3.98788E+04 0.00114  7.80088E+04 0.00094  7.45624E+04 0.00104  5.36401E+04 0.00126  3.47575E+04 0.00133  3.91946E+04 0.00135  3.65909E+04 0.00124  3.49726E+04 0.00147  5.38829E+04 0.00096  1.30809E+04 0.00222  1.64858E+04 0.00235  1.49791E+04 0.00259  8.68310E+03 0.00262  1.51921E+04 0.00205  1.03132E+04 0.00252  8.67745E+03 0.00299  1.63922E+03 0.00591  1.64685E+03 0.00538  1.68532E+03 0.00511  1.72166E+03 0.00605  1.71861E+03 0.00589  1.70153E+03 0.00714  1.75808E+03 0.00508  1.63815E+03 0.00535  3.11568E+03 0.00362  4.95138E+03 0.00305  6.32916E+03 0.00365  1.66892E+04 0.00201  1.74159E+04 0.00207  1.86914E+04 0.00195  1.23280E+04 0.00190  9.00309E+03 0.00210  6.85962E+03 0.00231  7.90454E+03 0.00273  1.46286E+04 0.00197  1.93131E+04 0.00173  3.60264E+04 0.00117  5.34520E+04 0.00093  7.69136E+04 0.00088  4.76547E+04 0.00102  3.33961E+04 0.00101  2.36319E+04 0.00125  2.08767E+04 0.00114  2.02534E+04 0.00098  1.66794E+04 0.00128  1.10747E+04 0.00151  1.00261E+04 0.00155  8.75865E+03 0.00186  7.21805E+03 0.00214  5.45784E+03 0.00195  3.40750E+03 0.00229  1.03758E+03 0.00361 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.35704E+00 0.00032 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  1.13454E+01 0.00032  2.40662E+00 0.00023 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  4.18075E-01 7.5E-05  6.15272E-01 5.7E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  2.05995E-02 0.00044  6.03705E-02 0.00014 ];
INF_ABS                   (idx, [1:   4]) = [  2.79976E-02 0.00033  2.18550E-01 0.00015 ];
INF_FISS                  (idx, [1:   4]) = [  7.39816E-03 0.00024  1.58179E-01 0.00016 ];
INF_NSF                   (idx, [1:   4]) = [  1.89829E-02 0.00024  3.85436E-01 0.00016 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.56590E+00 3.8E-05  2.43670E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.03898E+02 3.9E-06  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  5.75370E-08 0.00047  2.41325E-06 0.00014 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  3.90050E-01 7.8E-05  3.96817E-01 0.00018 ];
INF_SCATT1                (idx, [1:   4]) = [  4.35269E-02 0.00098  7.95062E-03 0.00940 ];
INF_SCATT2                (idx, [1:   4]) = [  2.33899E-02 0.00117  3.73745E-04 0.18314 ];
INF_SCATT3                (idx, [1:   4]) = [  1.25120E-02 0.00155  1.41750E-04 0.32481 ];
INF_SCATT4                (idx, [1:   4]) = [  7.96893E-03 0.00197  1.02107E-06 1.00000 ];
INF_SCATT5                (idx, [1:   4]) = [  4.04251E-03 0.00396  2.96159E-05 1.00000 ];
INF_SCATT6                (idx, [1:   4]) = [  2.02661E-03 0.00746  9.59105E-06 1.00000 ];
INF_SCATT7                (idx, [1:   4]) = [  8.54247E-04 0.01669  4.73578E-05 0.73156 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  3.90177E-01 7.7E-05  3.96817E-01 0.00018 ];
INF_SCATTP1               (idx, [1:   4]) = [  4.35277E-02 0.00098  7.95062E-03 0.00940 ];
INF_SCATTP2               (idx, [1:   4]) = [  2.33908E-02 0.00117  3.73745E-04 0.18314 ];
INF_SCATTP3               (idx, [1:   4]) = [  1.25119E-02 0.00154  1.41750E-04 0.32481 ];
INF_SCATTP4               (idx, [1:   4]) = [  7.96863E-03 0.00198  1.02107E-06 1.00000 ];
INF_SCATTP5               (idx, [1:   4]) = [  4.04275E-03 0.00395  2.96159E-05 1.00000 ];
INF_SCATTP6               (idx, [1:   4]) = [  2.02645E-03 0.00749  9.59105E-06 1.00000 ];
INF_SCATTP7               (idx, [1:   4]) = [  8.54086E-04 0.01677  4.73578E-05 0.73156 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  3.03532E-01 0.00026  5.85945E-01 0.00014 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.09818E+00 0.00026  5.68882E-01 0.00014 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  2.78703E-02 0.00034  2.18550E-01 0.00015 ];
INF_REMXS                 (idx, [1:   4]) = [  2.88581E-02 0.00054  2.19305E-01 0.00035 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  0.00000E+00 ;
PM148_LAMBDA              (idx, 1)        =  0.00000E+00 ;
PM148M_LAMBDA             (idx, 1)        =  0.00000E+00 ;
PM149_LAMBDA              (idx, 1)        =  0.00000E+00 ;
I135_LAMBDA               (idx, 1)        =  0.00000E+00 ;
XE135_LAMBDA              (idx, 1)        =  0.00000E+00 ;
XE135M_LAMBDA             (idx, 1)        =  0.00000E+00 ;
I135_BR                   (idx, 1)        =  0.00000E+00 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  3.89217E-01 7.8E-05  8.32821E-04 0.00314  8.50514E-04 0.00827  3.95966E-01 0.00018 ];
INF_S1                    (idx, [1:   8]) = [  4.37358E-02 0.00097 -2.08829E-04 0.00712 -3.73026E-05 0.07212  7.98793E-03 0.00939 ];
INF_S2                    (idx, [1:   8]) = [  2.34092E-02 0.00116 -1.93041E-05 0.05898 -3.92378E-05 0.06393  4.12982E-04 0.16621 ];
INF_S3                    (idx, [1:   8]) = [  1.25148E-02 0.00153 -2.87197E-06 0.32261 -1.38130E-05 0.13488  1.55563E-04 0.29470 ];
INF_S4                    (idx, [1:   8]) = [  7.97049E-03 0.00197 -1.56246E-06 0.56805 -9.11040E-06 0.18644  1.01315E-05 1.00000 ];
INF_S5                    (idx, [1:   8]) = [  4.04334E-03 0.00394 -8.29964E-07 0.94616 -8.13502E-06 0.19351  3.77510E-05 1.00000 ];
INF_S6                    (idx, [1:   8]) = [  2.02737E-03 0.00743 -7.62424E-07 0.72127 -1.49877E-06 1.00000  1.10898E-05 1.00000 ];
INF_S7                    (idx, [1:   8]) = [  8.54732E-04 0.01672 -4.85303E-07 1.00000 -2.64321E-07 1.00000  4.76221E-05 0.72869 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.89344E-01 7.7E-05  8.32821E-04 0.00314  8.50514E-04 0.00827  3.95966E-01 0.00018 ];
INF_SP1                   (idx, [1:   8]) = [  4.37366E-02 0.00097 -2.08829E-04 0.00712 -3.73026E-05 0.07212  7.98793E-03 0.00939 ];
INF_SP2                   (idx, [1:   8]) = [  2.34101E-02 0.00116 -1.93041E-05 0.05898 -3.92378E-05 0.06393  4.12982E-04 0.16621 ];
INF_SP3                   (idx, [1:   8]) = [  1.25148E-02 0.00153 -2.87197E-06 0.32261 -1.38130E-05 0.13488  1.55563E-04 0.29470 ];
INF_SP4                   (idx, [1:   8]) = [  7.97020E-03 0.00198 -1.56246E-06 0.56805 -9.11040E-06 0.18644  1.01315E-05 1.00000 ];
INF_SP5                   (idx, [1:   8]) = [  4.04358E-03 0.00393 -8.29964E-07 0.94616 -8.13502E-06 0.19351  3.77510E-05 1.00000 ];
INF_SP6                   (idx, [1:   8]) = [  2.02721E-03 0.00746 -7.62424E-07 0.72127 -1.49877E-06 1.00000  1.10898E-05 1.00000 ];
INF_SP7                   (idx, [1:   8]) = [  8.54571E-04 0.01680 -4.85303E-07 1.00000 -2.64321E-07 1.00000  4.76221E-05 0.72869 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  4.35378E-01 0.00103  1.55650E+00 0.00302 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  4.34424E-01 0.00135  1.61985E+00 0.00546 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  4.35363E-01 0.00171  1.62715E+00 0.00540 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  4.36431E-01 0.00149  1.44198E+00 0.00566 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  7.65659E-01 0.00103  2.14249E-01 0.00296 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  7.67369E-01 0.00135  2.06072E-01 0.00531 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  7.65754E-01 0.00171  2.05150E-01 0.00541 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  7.63854E-01 0.00149  2.31526E-01 0.00564 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  6.88925E-03 0.00548  1.94780E-04 0.03335  1.08813E-03 0.01373  1.08468E-03 0.01342  3.14323E-03 0.00807  1.03296E-03 0.01466  3.45472E-04 0.02525 ];
LAMBDA                    (idx, [1:  14]) = [  8.17511E-01 0.01316  1.24908E-02 1.9E-06  3.16565E-02 0.00021  1.10193E-01 0.00028  3.20681E-01 0.00023  1.34554E+00 0.00016  8.89360E+00 0.00145 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2024 15:50:46' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Fuel Pin' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1:  4])  = 'main' ;
WORKING_DIRECTORY         (idx, [1:103])  = '/home/thomas-guilbaud/Software/foamForNuclear/pythonapi/tutorials/tests/test_fuelPin_monteCarlo/Serpent' ;
HOSTNAME                  (idx, [1:  7])  = 'TS-P350' ;
CPU_TYPE                  (idx, [1: 46])  = '11th Gen Intel(R) Core(TM) i9-11900K @ 3.50GHz' ;
CPU_MHZ                   (idx, 1)        = 98.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Dec  2 16:37:57 2025' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Dec  2 16:40:08 2025' ;

% Run parameters:

POP                       (idx, 1)        = 10000 ;
CYCLES                    (idx, 1)        = 1000 ;
SKIP                      (idx, 1)        = 100 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1764689877254 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 15 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  15]) = [  1.01577E+00  1.02975E+00  1.00301E+00  1.01183E+00  9.99208E-01  1.00158E+00  9.83312E-01  9.89737E-01  9.90049E-01  9.74149E-01  9.93896E-01  9.94224E-01  1.00857E+00  9.92815E-01  1.01211E+00  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 71])  = '/home/thomas-guilbaud/Software/Serpent/xsdata/endfb7/sss_endfb7u.xsdata' ;
DECAY_DATA_FILE_PATH      (idx, [1:  3])  = 'N/A' ;
SFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
NFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  1.04694E-02 0.00070  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.89531E-01 7.5E-06  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  7.51604E-01 3.8E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  7.52090E-01 3.8E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.76776E+00 0.00014  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  3.19299E+01 0.00019  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  3.19299E+01 0.00019  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  1.05251E+01 0.00025  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  3.65280E-01 0.00077  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 1000 ;
SIMULATED_HISTORIES       (idx, 1)        = 10001055 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  1.00011E+04 0.00051 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  1.00011E+04 0.00051 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  2.67645E+01 ;
RUNNING_TIME              (idx, 1)        =  2.19018E+00 ;
INIT_TIME                 (idx, [1:  2])  = [  1.19833E-02  1.19833E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.33332E-04  2.33332E-04 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  2.17790E+00  2.17790E+00  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  2.18980E+00  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 12.22019 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  1.31870E+01 0.00622 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  8.79013E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 63913.86 ;
ALLOC_MEMSIZE             (idx, 1)        = 408.02;
MEMSIZE                   (idx, 1)        = 264.81;
XS_MEMSIZE                (idx, 1)        = 145.32;
MAT_MEMSIZE               (idx, 1)        = 14.72;
RES_MEMSIZE               (idx, 1)        = 4.27;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 100.50;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 143.20;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 7 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 128116 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.00000E+37 ;
URES_EMAX                 (idx, 1)        = -1.00000E+37 ;
URES_AVAIL                (idx, 1)        = 7 ;
URES_USED                 (idx, 1)        = 0 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 13 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 13 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 0 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 347 ;
TOT_TRANSMU_REA           (idx, 1)        = 0 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 0 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 1 ;
TMS_MODE                  (idx, 1)        = 0 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  0.00000E+00 ;
TOT_DECAY_HEAT            (idx, 1)        =  0.00000E+00 ;
TOT_SF_RATE               (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  0.00000E+00 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  0.00000E+00 ;
INHALATION_TOXICITY       (idx, 1)        =  0.00000E+00 ;
INGESTION_TOXICITY        (idx, 1)        =  0.00000E+00 ;
ACTINIDE_INH_TOX          (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ING_TOX          (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  0.00000E+00 ;
SR90_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
TE132_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
I131_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
I132_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
CS134_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
CS137_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  0.00000E+00 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  0.00000E+00 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  0.00000E+00 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  0.00000E+00 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  9.98582E-05 0.00027  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  5.39002E-01 0.00076 ];
U235_FISS                 (idx, [1:   4]) = [  4.35748E-01 0.00041  9.38022E-01 0.00012 ];
U238_FISS                 (idx, [1:   4]) = [  2.87640E-02 0.00187  6.19124E-02 0.00176 ];
U235_CAPT                 (idx, [1:   4]) = [  9.27791E-02 0.00102  1.73244E-01 0.00095 ];
U238_CAPT                 (idx, [1:   4]) = [  2.83368E-01 0.00061  5.29109E-01 0.00042 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 10001055 1.00000E+07 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 1.51611E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 10001055 1.00152E+07 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 5355383 5.36308E+06 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 4645672 4.65208E+06 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 10001055 1.00152E+07 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -1.21072E-07 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.50775E-11 0.00014 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  1.14288E+00 0.00014 ];
TOT_FISSRATE              (idx, [1:   2]) = [  4.64575E-01 0.00014 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  5.35425E-01 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.98582E-01 0.00027 ];
TOT_FLUX                  (idx, [1:   2]) = [  4.51379E+01 0.00022 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  3.19307E+01 0.00017 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.76384E+00 0.00028 ];
SIX_FF_F                  (idx, [1:   2]) = [  7.95396E-01 0.00020 ];
SIX_FF_P                  (idx, [1:   2]) = [  6.61915E-01 0.00023 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  1.23249E+00 0.00023 ];
SIX_FF_LF                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_LT                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.14448E+00 0.00035 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.14448E+00 0.00035 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.46004E+00 8.7E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02564E+02 9.1E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.14440E+00 0.00036  1.13662E+00 0.00035  7.85261E-03 0.00558 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.14461E+00 0.00014 ];
COL_KEFF                  (idx, [1:   2]) = [  1.14459E+00 0.00033 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.14461E+00 0.00014 ];
ABS_KINF                  (idx, [1:   2]) = [  1.14461E+00 0.00014 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.76226E+01 0.00013 ];
IMP_ALF                   (idx, [1:   2]) = [  1.76240E+01 5.4E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  4.45383E-07 0.00224 ];
IMP_EALF                  (idx, [1:   2]) = [  4.43859E-07 0.00096 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  2.08096E-01 0.00196 ];
IMP_AFGE                  (idx, [1:   2]) = [  2.07406E-01 0.00081 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.24107E-03 0.00389  1.74110E-04 0.02285  9.89280E-04 0.00910  9.79162E-04 0.00900  2.85504E-03 0.00568  9.36400E-04 0.00939  3.07070E-04 0.01759 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  8.07935E-01 0.00904  1.07795E-02 0.01261  3.16553E-02 0.00015  1.10172E-01 0.00019  3.20588E-01 0.00016  1.34554E+00 0.00012  8.60233E+00 0.00596 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  6.88925E-03 0.00548  1.94780E-04 0.03335  1.08813E-03 0.01373  1.08468E-03 0.01342  3.14323E-03 0.00807  1.03296E-03 0.01466  3.45472E-04 0.02525 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  8.17511E-01 0.01316  1.24908E-02 1.9E-06  3.16565E-02 0.00021  1.10193E-01 0.00028  3.20681E-01 0.00023  1.34554E+00 0.00016  8.89360E+00 0.00145 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  2.46166E-05 0.00072  2.46044E-05 0.00073  2.63206E-05 0.00745 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  2.81677E-05 0.00063  2.81537E-05 0.00063  3.01157E-05 0.00742 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.86045E-03 0.00563  1.92990E-04 0.03374  1.08306E-03 0.01382  1.08385E-03 0.01359  3.12644E-03 0.00834  1.02969E-03 0.01442  3.44413E-04 0.02658 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  8.14833E-01 0.01397  1.24908E-02 2.4E-06  3.16659E-02 0.00024  1.10153E-01 0.00030  3.20560E-01 0.00025  1.34579E+00 0.00018  8.88137E+00 0.00167 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  2.46358E-05 0.00157  2.46255E-05 0.00158  2.63151E-05 0.01735 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  2.81894E-05 0.00153  2.81776E-05 0.00153  3.01110E-05 0.01733 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  6.92380E-03 0.01680  1.85253E-04 0.10796  1.16040E-03 0.04184  1.10877E-03 0.04196  3.13077E-03 0.02501  1.00716E-03 0.04436  3.31450E-04 0.07858 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  7.52734E-01 0.03852  1.24908E-02 5.5E-06  3.16718E-02 0.00053  1.10150E-01 0.00070  3.20789E-01 0.00069  1.34491E+00 0.00044  8.93125E+00 0.00386 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  6.98345E-03 0.01646  1.90092E-04 0.10745  1.17569E-03 0.04111  1.12502E-03 0.04106  3.13617E-03 0.02413  1.01618E-03 0.04297  3.40288E-04 0.07664 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  7.59419E-01 0.03775  1.24908E-02 5.5E-06  3.16721E-02 0.00053  1.10162E-01 0.00071  3.20768E-01 0.00068  1.34499E+00 0.00043  8.93495E+00 0.00387 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -2.81848E+02 0.01690 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  2.46332E-05 0.00047 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  2.81865E-05 0.00030 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  6.91349E-03 0.00326 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -2.80749E+02 0.00333 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.05784E-07 0.00033 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  2.90081E-06 0.00029  2.90087E-06 0.00029  2.89182E-06 0.00352 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  3.12270E-05 0.00039  3.12287E-05 0.00039  3.09638E-05 0.00450 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  6.62311E-01 0.00023  6.61556E-01 0.00023  7.94700E-01 0.00604 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.02411E+01 0.00905 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  3.19299E+01 0.00019  3.51211E+01 0.00026 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  5])  = 'uClad' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  7.40343E+03 0.00460  3.00417E+04 0.00215  6.27132E+04 0.00156  6.88024E+04 0.00128  6.42254E+04 0.00120  6.99596E+04 0.00102  4.74491E+04 0.00118  4.22538E+04 0.00143  3.22868E+04 0.00135  2.63074E+04 0.00147  2.26876E+04 0.00171  2.05492E+04 0.00185  1.89410E+04 0.00217  1.79962E+04 0.00254  1.74748E+04 0.00186  1.50825E+04 0.00233  1.49223E+04 0.00200  1.48497E+04 0.00220  1.46349E+04 0.00257  2.84618E+04 0.00165  2.73974E+04 0.00190  1.98522E+04 0.00223  1.28078E+04 0.00256  1.48501E+04 0.00242  1.39332E+04 0.00241  1.27131E+04 0.00252  2.07252E+04 0.00201  4.73950E+03 0.00443  5.97153E+03 0.00333  5.41731E+03 0.00374  3.13532E+03 0.00563  5.47158E+03 0.00450  3.73399E+03 0.00413  3.16984E+03 0.00488  6.13243E+02 0.01109  6.02130E+02 0.00916  6.17952E+02 0.00974  6.29682E+02 0.00928  6.27671E+02 0.01101  6.17047E+02 0.01026  6.43693E+02 0.01027  6.02701E+02 0.01071  1.11802E+03 0.00807  1.79301E+03 0.00663  2.30977E+03 0.00572  6.07693E+03 0.00335  6.34663E+03 0.00333  6.86540E+03 0.00328  4.50978E+03 0.00408  3.31724E+03 0.00403  2.52364E+03 0.00507  2.95842E+03 0.00500  5.45314E+03 0.00401  7.12846E+03 0.00294  1.33685E+04 0.00240  1.98429E+04 0.00192  2.87493E+04 0.00149  1.79707E+04 0.00198  1.26549E+04 0.00228  8.97094E+03 0.00286  7.96278E+03 0.00388  7.81805E+03 0.00325  6.45486E+03 0.00370  4.31930E+03 0.00378  3.94462E+03 0.00407  3.47265E+03 0.00519  2.88641E+03 0.00438  2.25872E+03 0.00494  1.46582E+03 0.00672  4.93097E+02 0.01052 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.03034E+00 0.00044  9.07389E-01 0.00069 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  3.13507E-01 0.00012  2.98986E-01 9.5E-06 ];
INF_CAPT                  (idx, [1:   4]) = [  1.53408E-03 0.00156  4.25042E-03 0.00036 ];
INF_ABS                   (idx, [1:   4]) = [  1.53408E-03 0.00156  4.25042E-03 0.00036 ];
INF_FISS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NSF                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NUBAR                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  5.94901E-08 0.00083  2.44957E-06 0.00036 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  3.11973E-01 0.00012  2.94764E-01 6.3E-05 ];
INF_SCATT1                (idx, [1:   4]) = [  4.18951E-02 0.00126  2.22679E-03 0.04480 ];
INF_SCATT2                (idx, [1:   4]) = [  1.96347E-02 0.00208  1.34220E-04 0.52951 ];
INF_SCATT3                (idx, [1:   4]) = [  5.37268E-03 0.00651  1.23047E-04 0.51401 ];
INF_SCATT4                (idx, [1:   4]) = [  4.10477E-03 0.00782 -3.12111E-05 1.00000 ];
INF_SCATT5                (idx, [1:   4]) = [  1.89176E-03 0.01397 -2.12138E-05 1.00000 ];
INF_SCATT6                (idx, [1:   4]) = [  9.11833E-04 0.03108 -9.01575E-05 0.55240 ];
INF_SCATT7                (idx, [1:   4]) = [  1.96548E-04 0.12221  3.47618E-05 1.00000 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  3.11990E-01 0.00012  2.94764E-01 6.3E-05 ];
INF_SCATTP1               (idx, [1:   4]) = [  4.18956E-02 0.00126  2.22679E-03 0.04480 ];
INF_SCATTP2               (idx, [1:   4]) = [  1.96353E-02 0.00208  1.34220E-04 0.52951 ];
INF_SCATTP3               (idx, [1:   4]) = [  5.37267E-03 0.00650  1.23047E-04 0.51401 ];
INF_SCATTP4               (idx, [1:   4]) = [  4.10482E-03 0.00782 -3.12111E-05 1.00000 ];
INF_SCATTP5               (idx, [1:   4]) = [  1.89195E-03 0.01397 -2.12138E-05 1.00000 ];
INF_SCATTP6               (idx, [1:   4]) = [  9.11651E-04 0.03114 -9.01575E-05 0.55240 ];
INF_SCATTP7               (idx, [1:   4]) = [  1.96492E-04 0.12220  3.47618E-05 1.00000 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.31136E-01 0.00035  2.96670E-01 0.00034 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.44216E+00 0.00035  1.12359E+00 0.00034 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  1.51687E-03 0.00162  4.25042E-03 0.00036 ];
INF_REMXS                 (idx, [1:   4]) = [  1.86835E-03 0.00442  4.85935E-03 0.00432 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  0.00000E+00 ;
PM148_LAMBDA              (idx, 1)        =  0.00000E+00 ;
PM148M_LAMBDA             (idx, 1)        =  0.00000E+00 ;
PM149_LAMBDA              (idx, 1)        =  0.00000E+00 ;
I135_LAMBDA               (idx, 1)        =  0.00000E+00 ;
XE135_LAMBDA              (idx, 1)        =  0.00000E+00 ;
XE135M_LAMBDA             (idx, 1)        =  0.00000E+00 ;
I135_BR                   (idx, 1)        =  0.00000E+00 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  3.11638E-01 0.00012  3.34792E-04 0.00893  6.37432E-04 0.01188  2.94127E-01 7.0E-05 ];
INF_S1                    (idx, [1:   8]) = [  4.19749E-02 0.00125 -7.97988E-05 0.02084 -8.79685E-05 0.05498  2.31476E-03 0.04307 ];
INF_S2                    (idx, [1:   8]) = [  1.96407E-02 0.00209 -6.00040E-06 0.20733 -2.56330E-05 0.14050  1.59853E-04 0.44722 ];
INF_S3                    (idx, [1:   8]) = [  5.37660E-03 0.00650 -3.91563E-06 0.29029 -7.02833E-06 0.45475  1.30075E-04 0.48534 ];
INF_S4                    (idx, [1:   8]) = [  4.10524E-03 0.00780 -4.74573E-07 1.00000 -1.05507E-06 1.00000 -3.01560E-05 1.00000 ];
INF_S5                    (idx, [1:   8]) = [  1.89216E-03 0.01393 -4.00338E-07 1.00000 -8.17302E-06 0.28678 -1.30408E-05 1.00000 ];
INF_S6                    (idx, [1:   8]) = [  9.11556E-04 0.03116  2.76932E-07 1.00000  9.67531E-07 1.00000 -9.11250E-05 0.54987 ];
INF_S7                    (idx, [1:   8]) = [  1.96560E-04 0.12192 -1.16783E-08 1.00000 -1.76509E-06 1.00000  3.65269E-05 1.00000 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.11656E-01 0.00012  3.34792E-04 0.00893  6.37432E-04 0.01188  2.94127E-01 7.0E-05 ];
INF_SP1                   (idx, [1:   8]) = [  4.19754E-02 0.00125 -7.97988E-05 0.02084 -8.79685E-05 0.05498  2.31476E-03 0.04307 ];
INF_SP2                   (idx, [1:   8]) = [  1.96413E-02 0.00209 -6.00040E-06 0.20733 -2.56330E-05 0.14050  1.59853E-04 0.44722 ];
INF_SP3                   (idx, [1:   8]) = [  5.37658E-03 0.00649 -3.91563E-06 0.29029 -7.02833E-06 0.45475  1.30075E-04 0.48534 ];
INF_SP4                   (idx, [1:   8]) = [  4.10529E-03 0.00780 -4.74573E-07 1.00000 -1.05507E-06 1.00000 -3.01560E-05 1.00000 ];
INF_SP5                   (idx, [1:   8]) = [  1.89235E-03 0.01394 -4.00338E-07 1.00000 -8.17302E-06 0.28678 -1.30408E-05 1.00000 ];
INF_SP6                   (idx, [1:   8]) = [  9.11374E-04 0.03122  2.76932E-07 1.00000  9.67531E-07 1.00000 -9.11250E-05 0.54987 ];
INF_SP7                   (idx, [1:   8]) = [  1.96504E-04 0.12190 -1.16783E-08 1.00000 -1.76509E-06 1.00000  3.65269E-05 1.00000 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  8.30889E-01 0.00360  4.82062E+00 0.00948 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  8.30103E-01 0.00472  4.88474E+00 0.01718 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  8.29143E-01 0.00401  5.13422E+00 0.01647 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  8.34532E-01 0.00519  4.59803E+00 0.01643 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  4.01430E-01 0.00358  6.94490E-02 0.00939 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  4.01996E-01 0.00473  6.91589E-02 0.01605 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  4.02340E-01 0.00403  6.57390E-02 0.01554 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  3.99953E-01 0.00519  7.34492E-02 0.01628 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
LAMBDA                    (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Jun  3 2024 15:50:46' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Fuel Pin' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1:  4])  = 'main' ;
WORKING_DIRECTORY         (idx, [1:103])  = '/home/thomas-guilbaud/Software/foamForNuclear/pythonapi/tutorials/tests/test_fuelPin_monteCarlo/Serpent' ;
HOSTNAME                  (idx, [1:  7])  = 'TS-P350' ;
CPU_TYPE                  (idx, [1: 46])  = '11th Gen Intel(R) Core(TM) i9-11900K @ 3.50GHz' ;
CPU_MHZ                   (idx, 1)        = 98.0 ;
START_DATE                (idx, [1: 24])  = 'Tue Dec  2 16:37:57 2025' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Tue Dec  2 16:40:08 2025' ;

% Run parameters:

POP                       (idx, 1)        = 10000 ;
CYCLES                    (idx, 1)        = 1000 ;
SKIP                      (idx, 1)        = 100 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1764689877254 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 15 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  15]) = [  1.01577E+00  1.02975E+00  1.00301E+00  1.01183E+00  9.99208E-01  1.00158E+00  9.83312E-01  9.89737E-01  9.90049E-01  9.74149E-01  9.93896E-01  9.94224E-01  1.00857E+00  9.92815E-01  1.01211E+00  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 71])  = '/home/thomas-guilbaud/Software/Serpent/xsdata/endfb7/sss_endfb7u.xsdata' ;
DECAY_DATA_FILE_PATH      (idx, [1:  3])  = 'N/A' ;
SFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
NFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  1.04694E-02 0.00070  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.89531E-01 7.5E-06  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  7.51604E-01 3.8E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  7.52090E-01 3.8E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.76776E+00 0.00014  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  3.19299E+01 0.00019  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  3.19299E+01 0.00019  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  1.05251E+01 0.00025  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  3.65280E-01 0.00077  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 1000 ;
SIMULATED_HISTORIES       (idx, 1)        = 10001055 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  1.00011E+04 0.00051 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  1.00011E+04 0.00051 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  2.67648E+01 ;
RUNNING_TIME              (idx, 1)        =  2.19022E+00 ;
INIT_TIME                 (idx, [1:  2])  = [  1.19833E-02  1.19833E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.33332E-04  2.33332E-04 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  2.17790E+00  2.17790E+00  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  2.18980E+00  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 12.22014 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  1.31870E+01 0.00622 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  8.79000E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 63913.86 ;
ALLOC_MEMSIZE             (idx, 1)        = 408.02;
MEMSIZE                   (idx, 1)        = 264.81;
XS_MEMSIZE                (idx, 1)        = 145.32;
MAT_MEMSIZE               (idx, 1)        = 14.72;
RES_MEMSIZE               (idx, 1)        = 4.27;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 100.50;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 143.20;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 7 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 128116 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.00000E+37 ;
URES_EMAX                 (idx, 1)        = -1.00000E+37 ;
URES_AVAIL                (idx, 1)        = 7 ;
URES_USED                 (idx, 1)        = 0 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 13 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 13 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 0 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 347 ;
TOT_TRANSMU_REA           (idx, 1)        = 0 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 0 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 1 ;
TMS_MODE                  (idx, 1)        = 0 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  0.00000E+00 ;
TOT_DECAY_HEAT            (idx, 1)        =  0.00000E+00 ;
TOT_SF_RATE               (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  0.00000E+00 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  0.00000E+00 ;
INHALATION_TOXICITY       (idx, 1)        =  0.00000E+00 ;
INGESTION_TOXICITY        (idx, 1)        =  0.00000E+00 ;
ACTINIDE_INH_TOX          (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ING_TOX          (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  0.00000E+00 ;
SR90_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
TE132_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
I131_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
I132_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
CS134_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
CS137_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  0.00000E+00 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  0.00000E+00 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  0.00000E+00 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  0.00000E+00 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  9.98582E-05 0.00027  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  5.39002E-01 0.00076 ];
U235_FISS                 (idx, [1:   4]) = [  4.35748E-01 0.00041  9.38022E-01 0.00012 ];
U238_FISS                 (idx, [1:   4]) = [  2.87640E-02 0.00187  6.19124E-02 0.00176 ];
U235_CAPT                 (idx, [1:   4]) = [  9.27791E-02 0.00102  1.73244E-01 0.00095 ];
U238_CAPT                 (idx, [1:   4]) = [  2.83368E-01 0.00061  5.29109E-01 0.00042 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 10001055 1.00000E+07 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 1.51611E+04 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 10001055 1.00152E+07 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 5355383 5.36308E+06 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 4645672 4.65208E+06 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 10001055 1.00152E+07 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -1.21072E-07 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.50775E-11 0.00014 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  1.14288E+00 0.00014 ];
TOT_FISSRATE              (idx, [1:   2]) = [  4.64575E-01 0.00014 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  5.35425E-01 0.00012 ];
TOT_ABSRATE               (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.98582E-01 0.00027 ];
TOT_FLUX                  (idx, [1:   2]) = [  4.51379E+01 0.00022 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  3.19307E+01 0.00017 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.76384E+00 0.00028 ];
SIX_FF_F                  (idx, [1:   2]) = [  7.95396E-01 0.00020 ];
SIX_FF_P                  (idx, [1:   2]) = [  6.61915E-01 0.00023 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  1.23249E+00 0.00023 ];
SIX_FF_LF                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_LT                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.14448E+00 0.00035 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.14448E+00 0.00035 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.46004E+00 8.7E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02564E+02 9.1E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.14440E+00 0.00036  1.13662E+00 0.00035  7.85261E-03 0.00558 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.14461E+00 0.00014 ];
COL_KEFF                  (idx, [1:   2]) = [  1.14459E+00 0.00033 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.14461E+00 0.00014 ];
ABS_KINF                  (idx, [1:   2]) = [  1.14461E+00 0.00014 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.76226E+01 0.00013 ];
IMP_ALF                   (idx, [1:   2]) = [  1.76240E+01 5.4E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  4.45383E-07 0.00224 ];
IMP_EALF                  (idx, [1:   2]) = [  4.43859E-07 0.00096 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  2.08096E-01 0.00196 ];
IMP_AFGE                  (idx, [1:   2]) = [  2.07406E-01 0.00081 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 6 ;
FWD_ANA_BETA_ZERO         (idx, [1:  14]) = [  6.24107E-03 0.00389  1.74110E-04 0.02285  9.89280E-04 0.00910  9.79162E-04 0.00900  2.85504E-03 0.00568  9.36400E-04 0.00939  3.07070E-04 0.01759 ];
FWD_ANA_LAMBDA            (idx, [1:  14]) = [  8.07935E-01 0.00904  1.07795E-02 0.01261  3.16553E-02 0.00015  1.10172E-01 0.00019  3.20588E-01 0.00016  1.34554E+00 0.00012  8.60233E+00 0.00596 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  14]) = [  6.88925E-03 0.00548  1.94780E-04 0.03335  1.08813E-03 0.01373  1.08468E-03 0.01342  3.14323E-03 0.00807  1.03296E-03 0.01466  3.45472E-04 0.02525 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  14]) = [  8.17511E-01 0.01316  1.24908E-02 1.9E-06  3.16565E-02 0.00021  1.10193E-01 0.00028  3.20681E-01 0.00023  1.34554E+00 0.00016  8.89360E+00 0.00145 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  2.46166E-05 0.00072  2.46044E-05 0.00073  2.63206E-05 0.00745 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  2.81677E-05 0.00063  2.81537E-05 0.00063  3.01157E-05 0.00742 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  14]) = [  6.86045E-03 0.00563  1.92990E-04 0.03374  1.08306E-03 0.01382  1.08385E-03 0.01359  3.12644E-03 0.00834  1.02969E-03 0.01442  3.44413E-04 0.02658 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  14]) = [  8.14833E-01 0.01397  1.24908E-02 2.4E-06  3.16659E-02 0.00024  1.10153E-01 0.00030  3.20560E-01 0.00025  1.34579E+00 0.00018  8.88137E+00 0.00167 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  2.46358E-05 0.00157  2.46255E-05 0.00158  2.63151E-05 0.01735 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  2.81894E-05 0.00153  2.81776E-05 0.00153  3.01110E-05 0.01733 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  14]) = [  6.92380E-03 0.01680  1.85253E-04 0.10796  1.16040E-03 0.04184  1.10877E-03 0.04196  3.13077E-03 0.02501  1.00716E-03 0.04436  3.31450E-04 0.07858 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  14]) = [  7.52734E-01 0.03852  1.24908E-02 5.5E-06  3.16718E-02 0.00053  1.10150E-01 0.00070  3.20789E-01 0.00069  1.34491E+00 0.00044  8.93125E+00 0.00386 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  14]) = [  6.98345E-03 0.01646  1.90092E-04 0.10745  1.17569E-03 0.04111  1.12502E-03 0.04106  3.13617E-03 0.02413  1.01618E-03 0.04297  3.40288E-04 0.07664 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  14]) = [  7.59419E-01 0.03775  1.24908E-02 5.5E-06  3.16721E-02 0.00053  1.10162E-01 0.00071  3.20768E-01 0.00068  1.34499E+00 0.00043  8.93495E+00 0.00387 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -2.81848E+02 0.01690 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  2.46332E-05 0.00047 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  2.81865E-05 0.00030 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  6.91349E-03 0.00326 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -2.80749E+02 0.00333 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  5.05784E-07 0.00033 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  2.90081E-06 0.00029  2.90087E-06 0.00029  2.89182E-06 0.00352 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  3.12270E-05 0.00039  3.12287E-05 0.00039  3.09638E-05 0.00450 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  6.62311E-01 0.00023  6.61556E-01 0.00023  7.94700E-01 0.00604 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.02411E+01 0.00905 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  3.19299E+01 0.00019  3.51211E+01 0.00026 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  6])  = 'uWater' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  3.90289E+04 0.00259  1.57086E+05 0.00119  3.27257E+05 0.00052  3.57766E+05 0.00060  3.32847E+05 0.00045  3.62413E+05 0.00043  2.48600E+05 0.00042  2.21442E+05 0.00045  1.70588E+05 0.00046  1.40071E+05 0.00053  1.21361E+05 0.00053  1.09639E+05 0.00049  1.01439E+05 0.00048  9.64001E+04 0.00055  9.41491E+04 0.00048  8.14133E+04 0.00051  8.06791E+04 0.00057  7.98547E+04 0.00052  7.85937E+04 0.00053  1.53634E+05 0.00038  1.48839E+05 0.00037  1.07968E+05 0.00046  6.99498E+04 0.00044  8.13275E+04 0.00052  7.70861E+04 0.00056  6.87665E+04 0.00055  1.15247E+05 0.00045  2.56472E+04 0.00075  3.20618E+04 0.00075  2.90858E+04 0.00066  1.69364E+04 0.00130  2.93612E+04 0.00068  1.99737E+04 0.00125  1.71037E+04 0.00118  3.26692E+03 0.00272  3.24371E+03 0.00255  3.33566E+03 0.00213  3.42378E+03 0.00227  3.38193E+03 0.00233  3.33035E+03 0.00217  3.44926E+03 0.00265  3.23219E+03 0.00208  6.11792E+03 0.00168  9.76420E+03 0.00132  1.24363E+04 0.00128  3.27852E+04 0.00084  3.44184E+04 0.00090  3.71754E+04 0.00075  2.48013E+04 0.00101  1.82464E+04 0.00090  1.39161E+04 0.00124  1.61696E+04 0.00092  2.99596E+04 0.00087  3.94265E+04 0.00070  7.34452E+04 0.00060  1.09482E+05 0.00058  1.59003E+05 0.00052  9.99984E+04 0.00056  7.08965E+04 0.00062  5.06417E+04 0.00063  4.52073E+04 0.00066  4.43787E+04 0.00066  3.69664E+04 0.00079  2.49308E+04 0.00063  2.30001E+04 0.00073  2.03993E+04 0.00085  1.72858E+04 0.00096  1.36375E+04 0.00075  9.22886E+03 0.00091  3.35154E+03 0.00122 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  2.13765E+01 0.00024  5.07261E+00 0.00028 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  6.77363E-01 7.7E-05  1.91092E+00 4.6E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  7.01738E-04 0.00023  2.59077E-02 7.6E-05 ];
INF_ABS                   (idx, [1:   4]) = [  7.01738E-04 0.00023  2.59077E-02 7.6E-05 ];
INF_FISS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NSF                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NUBAR                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  6.08219E-08 0.00021  2.48543E-06 7.6E-05 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  6.76662E-01 7.8E-05  1.88502E+00 4.9E-05 ];
INF_SCATT1                (idx, [1:   4]) = [  3.95255E-01 9.0E-05  5.49899E-01 0.00021 ];
INF_SCATT2                (idx, [1:   4]) = [  1.51517E-01 0.00017  1.35177E-01 0.00066 ];
INF_SCATT3                (idx, [1:   4]) = [  5.28379E-03 0.00492  4.09910E-02 0.00157 ];
INF_SCATT4                (idx, [1:   4]) = [ -2.27350E-02 0.00105 -1.06990E-02 0.00624 ];
INF_SCATT5                (idx, [1:   4]) = [ -2.05747E-03 0.00779  8.48604E-03 0.00744 ];
INF_SCATT6                (idx, [1:   4]) = [  7.62500E-03 0.00205 -2.20170E-02 0.00303 ];
INF_SCATT7                (idx, [1:   4]) = [  7.79575E-04 0.01962 -6.88368E-05 0.73749 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  6.76662E-01 7.8E-05  1.88502E+00 4.9E-05 ];
INF_SCATTP1               (idx, [1:   4]) = [  3.95255E-01 9.0E-05  5.49899E-01 0.00021 ];
INF_SCATTP2               (idx, [1:   4]) = [  1.51517E-01 0.00017  1.35177E-01 0.00066 ];
INF_SCATTP3               (idx, [1:   4]) = [  5.28379E-03 0.00492  4.09910E-02 0.00157 ];
INF_SCATTP4               (idx, [1:   4]) = [ -2.27350E-02 0.00105 -1.06990E-02 0.00624 ];
INF_SCATTP5               (idx, [1:   4]) = [ -2.05747E-03 0.00779  8.48604E-03 0.00744 ];
INF_SCATTP6               (idx, [1:   4]) = [  7.62500E-03 0.00205 -2.20170E-02 0.00303 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.79575E-04 0.01962 -6.88368E-05 0.73749 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  1.75954E-01 0.00026  1.16469E+00 0.00014 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.89444E+00 0.00026  2.86199E-01 0.00014 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  7.01738E-04 0.00023  2.59077E-02 7.6E-05 ];
INF_REMXS                 (idx, [1:   4]) = [  3.16217E-02 0.00022  2.75176E-02 0.00079 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  0.00000E+00 ;
PM148_LAMBDA              (idx, 1)        =  0.00000E+00 ;
PM148M_LAMBDA             (idx, 1)        =  0.00000E+00 ;
PM149_LAMBDA              (idx, 1)        =  0.00000E+00 ;
I135_LAMBDA               (idx, 1)        =  0.00000E+00 ;
XE135_LAMBDA              (idx, 1)        =  0.00000E+00 ;
XE135M_LAMBDA             (idx, 1)        =  0.00000E+00 ;
I135_BR                   (idx, 1)        =  0.00000E+00 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  6.45742E-01 7.5E-05  3.09208E-02 0.00022  1.61194E-03 0.00378  1.88341E+00 4.9E-05 ];
INF_S1                    (idx, [1:   8]) = [  3.85976E-01 9.2E-05  9.27873E-03 0.00058  1.03713E-03 0.00465  5.48862E-01 0.00021 ];
INF_S2                    (idx, [1:   8]) = [  1.54239E-01 0.00017 -2.72216E-03 0.00160  5.58442E-04 0.00649  1.34619E-01 0.00066 ];
INF_S3                    (idx, [1:   8]) = [  8.51869E-03 0.00298 -3.23489E-03 0.00144  2.01700E-04 0.01297  4.07893E-02 0.00157 ];
INF_S4                    (idx, [1:   8]) = [ -2.16820E-02 0.00109 -1.05306E-03 0.00371  6.11259E-06 0.34027 -1.07051E-02 0.00622 ];
INF_S5                    (idx, [1:   8]) = [ -2.10810E-03 0.00771  5.06299E-05 0.06298 -7.44423E-05 0.02695  8.56048E-03 0.00744 ];
INF_S6                    (idx, [1:   8]) = [  7.87608E-03 0.00192 -2.51071E-04 0.01328 -9.50764E-05 0.01912 -2.19219E-02 0.00306 ];
INF_S7                    (idx, [1:   8]) = [  1.08842E-03 0.01353 -3.08841E-04 0.01013 -8.71664E-05 0.01830  1.83295E-05 1.00000 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  6.45742E-01 7.5E-05  3.09208E-02 0.00022  1.61194E-03 0.00378  1.88341E+00 4.9E-05 ];
INF_SP1                   (idx, [1:   8]) = [  3.85976E-01 9.2E-05  9.27873E-03 0.00058  1.03713E-03 0.00465  5.48862E-01 0.00021 ];
INF_SP2                   (idx, [1:   8]) = [  1.54239E-01 0.00017 -2.72216E-03 0.00160  5.58442E-04 0.00649  1.34619E-01 0.00066 ];
INF_SP3                   (idx, [1:   8]) = [  8.51869E-03 0.00298 -3.23489E-03 0.00144  2.01700E-04 0.01297  4.07893E-02 0.00157 ];
INF_SP4                   (idx, [1:   8]) = [ -2.16820E-02 0.00109 -1.05306E-03 0.00371  6.11259E-06 0.34027 -1.07051E-02 0.00622 ];
INF_SP5                   (idx, [1:   8]) = [ -2.10810E-03 0.00771  5.06299E-05 0.06298 -7.44423E-05 0.02695  8.56048E-03 0.00744 ];
INF_SP6                   (idx, [1:   8]) = [  7.87608E-03 0.00192 -2.51071E-04 0.01328 -9.50764E-05 0.01912 -2.19219E-02 0.00306 ];
INF_SP7                   (idx, [1:   8]) = [  1.08842E-03 0.01353 -3.08841E-04 0.01013 -8.71664E-05 0.01830  1.83295E-05 1.00000 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  1.76282E-01 0.00044  6.03542E-01 0.00163 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  1.76433E-01 0.00071  6.01773E-01 0.00285 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  1.76254E-01 0.00064  6.01066E-01 0.00297 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  1.76168E-01 0.00080  6.08296E-01 0.00257 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  1.89093E+00 0.00044  5.52367E-01 0.00163 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  1.88933E+00 0.00071  5.54139E-01 0.00284 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  1.89125E+00 0.00063  5.54806E-01 0.00293 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  1.89220E+00 0.00080  5.48156E-01 0.00257 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
LAMBDA                    (idx, [1:  14]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

