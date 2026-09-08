# Multiscale Model of Thrombolysis and Clot Lysis

MATLAB implementation of a mechanistic multiscale model for thrombolysis, coupling systemic thrombolytic pharmacology with local transport, fibrinolytic reactions, clot degradation, and blood flow through an occluded cerebral artery.

The model was developed to investigate how administration of tissue plasminogen activator (tPA) leads to changes in circulating fibrinolytic proteins, transport of thrombolytic agents to an intravascular clot, biochemical fibrin degradation, and eventual restoration of blood flow.

## Overview

Thrombolysis is a multiscale process involving events ranging from systemic drug administration to molecular reactions within an occluding blood clot.

This repository contains two connected modelling components:

```text
tPA administration
       │
       ▼
Systemic pharmacokinetic / biochemical model
       │
       │ circulating concentrations
       ▼
Arterial inlet concentrations
       │
       ▼
1D convection–diffusion–reaction model
       │
       ▼
Transport into the clot
       │
       ▼
tPA–PLG–PLS fibrinolytic reactions
       │
       ▼
Fibrin degradation
       │
       ▼
↑ porosity / ↑ permeability
       │
       ▼
↓ clot resistance
       │
       ▼
Restoration of blood flow
```

The framework therefore links **systemic thrombolytic treatment** with **local clot lysis and haemodynamics**.

## Repository structure

```text
clotlysis/
│
└── clotLysis/
    │
    ├── TwoCompartCoupledLysisModel/
    │   ├── compartmentalModel_NewMain.m
    │   ├── compartmentalODE_NewFunc.m
    │   ├── dosageRegimen_func.m
    │   ├── infusionRate_func_general.m
    │   ├── elimination_para.m
    │   ├── kinetic_para.m
    │   ├── initConc_values.m
    │   ├── compare_model_expData_main.m
    │   ├── paraEst_main.m
    │   ├── paraEst_expData_func.m
    │   ├── calculateError_btw_simNexp.m
    │   └── experimental data (*.mat)
    │
    └── integratedModel_coupledWithCompart/
        ├── integratedModel_main.m
        ├── integratedModel_new_main.m
        ├── integratedModel_new_main1.m
        ├── integratedModel_new_main2.m
        ├── integratedModel_new_main3.m
        ├── integratedModel_new_main4.m
        ├── solveArteries_func.m
        ├── convectionDiffusionReaction_ode_func.m
        ├── convectionDiffusionReaction_mat_func.m
        ├── reactionRate_cal_func.m
        ├── fibrinolytic_para.m
        ├── calculateClotProperties_func.m
        ├── permeability_func.m
        ├── DarcyLaw_func.m
        ├── clotProperty_func.m
        ├── matrix_cal_func.m
        └── plot3D.m
```

## 1. Systemic thrombolysis model

The `TwoCompartCoupledLysisModel` directory contains a two-compartment pharmacokinetic and biochemical model describing systemic changes following administration of tPA.

The model includes:

- tPA administration and infusion regimens
- central and peripheral tPA compartments
- plasminogen (PLG)
- plasmin (PLS)
- fibrinogen (FBG)
- α2-antiplasmin (AP)
- α2-macroglobulin (MG)
- plasminogen activator inhibitor (PAI)
- elimination and inter-compartmental transport
- systemic fibrinolytic reactions

A typical clinical alteplase regimen can be defined through:

```matlab
tPAdose = 0.9*80;            % total tPA dose [mg]

N_stage = 2;

t_infusion = [1,60]*60;      % infusion duration [s]

perc_infusion = [0.1,0.9];   % fraction delivered at each stage
```

The systemic equations are solved using MATLAB stiff ODE solvers such as `ode15s` and `ode23s`.

### Experimental comparison

The repository includes digitised literature datasets used for comparison with the systemic model:

```text
Collen1986.mat
Noe1987.mat
Tanswell1989.mat
Tebbe1989.mat
Tienfenbrunn1986.mat
Verstraete1985.mat
```

`compare_model_expData_main.m` can be used to compare simulated circulating concentrations with experimental data reported in the literature.

The comparison includes quantities such as:

- tPA
- plasminogen
- fibrinogen
- α2-antiplasmin

Parameter-estimation utilities are also provided.

## 2. Local clot-lysis model

The `integratedModel_coupledWithCompart` directory contains a spatially resolved one-dimensional model of thrombolysis in a cerebral artery.

The model represents the arterial geometry along its axial direction and describes:

- blood flow
- thrombolytic protein transport
- diffusion/dispersion
- biochemical fibrinolysis
- fibrin-bound reactions
- clot degradation
- changes in clot porosity
- changes in clot permeability
- changes in hydraulic resistance
- restoration of arterial flow

The model includes representative dimensions for the:

- internal carotid artery (ICA)
- anterior cerebral artery (ACA)
- middle cerebral artery (MCA)

and allows an occluding clot to be positioned within the MCA.

## Protein transport

Transport of fibrinolytic species is represented using convection–diffusion–reaction equations of the general form

```text
concentration change
        =
convection
        +
diffusion
        +
biochemical reaction
```

The free-phase species represented in the extended model include:

| Species | Description |
|---|---|
| tPA | tissue plasminogen activator |
| PLG | plasminogen |
| PLS | plasmin |
| FBG | fibrinogen |
| AP | α2-antiplasmin |
| MG | α2-macroglobulin |
| PAI | plasminogen activator inhibitor |
| PLS–AP | plasmin–antiplasmin complex |

Fibrin-associated species are represented separately from free plasma species.

## Fibrinolytic reaction model

Within the clot, tPA, plasminogen and plasmin interact with fibrin binding sites.

The model includes:

```text
free tPA
   ⇅
fibrin-bound tPA

free PLG
   ⇅
fibrin-bound PLG

free PLS
   ⇅
fibrin-bound PLS
```

Bound tPA activates bound plasminogen to form plasmin.

```text
tPA + PLG
     │
     ▼
   plasmin
     │
     ▼
fibrin cleavage
```

Plasmin-mediated cleavage progressively decreases the concentration of intact fibrin binding sites.

The fibrinolysis kinetics include:

- adsorption/desorption of tPA
- adsorption/desorption of plasminogen
- adsorption/desorption of plasmin
- plasminogen activation
- plasmin-mediated fibrin degradation
- release of plasmin from lysed binding sites

## Plasma fibrinolysis

The model can additionally include reactions occurring in the plasma phase.

These include:

```text
tPA + PLG → PLS + tPA

PLS + FBG → fibrin degradation products

PLS + AP ⇌ PLS–AP → inactive products

PLS + MG → inactive products

tPA + PAI → inactive products
```

This allows systemic depletion and inhibition of fibrinolytic proteins to be represented alongside clot-specific reactions.

## Dynamic clot properties

Fibrin degradation changes the structural properties of the clot.

The local extent of lysis is calculated from the remaining fibrin binding-site concentration:

```text
Extent of lysis
       │
       ▼
decrease in fibrin volume fraction
       │
       ▼
increase in porosity
       │
       ▼
increase in permeability
       │
       ▼
decrease in hydraulic resistance
```

A clot segment is treated as effectively lysed when its extent of lysis exceeds a prescribed threshold.

For example:

```matlab
EL_crit = 0.95;
```

corresponds to a 95% local lysis threshold.

## Permeability and blood flow

Clot permeability is calculated as a function of fibrin fibre radius and fibrin volume fraction.

The permeability changes dynamically as fibrin is degraded.

Flow through the remaining clot is then estimated using Darcy's law:

```text
          A
Q = ------------- (dp/dx)
       R μ
```

where:

- `Q` is volumetric flow rate
- `A` is arterial cross-sectional area
- `R` is clot hydraulic resistance
- `μ` is blood viscosity
- `dp/dx` is the pressure gradient

As fibrinolysis progresses,

```text
fibrin ↓
   ↓
clot permeability ↑
   ↓
clot resistance ↓
   ↓
blood flow ↑
```

providing feedback between biochemical clot degradation and haemodynamics.

## Example simulation

A developed example of the one-dimensional model is provided in:

```matlab
integratedModel_new_main4.m
```

The example defines parameters including:

```matlab
L_MCA   = 15e-3;    % MCA length [m]

L_clot  = 5e-3;     % clot length [m]

R_f0    = 300;      % fibrin fibre radius [nm]

epsilon_0 = 0.95;   % initial clot porosity

dx = 0.01e-3;       % spatial discretisation [m]

dt = 0.05;          % time step [s]

tf = 60*180;        % maximum simulation time [s]
```

These values are research-model parameters and should be checked or modified for the intended application.

## Running the model

### Requirements

The repository is written entirely in MATLAB.

A recent MATLAB installation is recommended.

### Clone the repository

```bash
git clone https://github.com/keviccc/clotlysis.git
```

Open MATLAB and navigate to:

```text
clotlysis/clotLysis/
```

### Run the systemic model

Navigate to:

```text
TwoCompartCoupledLysisModel/
```

and run:

```matlab
compartmentalModel_NewMain
```

To compare model predictions against the included literature datasets:

```matlab
compare_model_expData_main
```

### Run the 1D clot-lysis model

Navigate to:

```text
integratedModel_coupledWithCompart/
```

and run an integrated-model script, for example:

```matlab
integratedModel_new_main4
```

Several `integratedModel_new_main*.m` files are retained because they correspond to different stages and configurations of model development.

## Model outputs

Depending on the selected script, the model can generate spatial and temporal predictions of:

### Circulating species

- tPA
- plasminogen
- plasmin
- fibrinogen
- α2-antiplasmin
- α2-macroglobulin
- PAI

### Local thrombolysis

- free tPA concentration
- free plasminogen concentration
- free plasmin concentration
- fibrin-bound tPA
- fibrin-bound plasminogen
- fibrin-bound plasmin
- remaining fibrin binding sites
- extent of clot lysis

### Clot properties

- fibrin volume fraction
- porosity
- permeability
- hydraulic resistance
- remaining clot region

### Haemodynamics

- flow through the occluded artery
- change in flow during clot degradation
- predicted time to clot dissolution / recanalisation

## Numerical framework

The complete model combines several numerical approaches:

```text
Two-compartment ODE model
           │
           ▼
systemic protein concentrations
           │
           ▼
1D finite-difference transport
           │
           ▼
convection–diffusion–reaction
           │
           ▼
fibrinolysis kinetics
           │
           ▼
dynamic clot properties
           │
           ▼
Darcy flow
```

This framework enables thrombolytic treatment to be investigated across multiple spatial and temporal scales.

## Research use

The code in this repository was developed as research code for mechanistic investigation of thrombolysis.

It is intended primarily for:

- computational studies of fibrinolysis
- thrombolytic drug-delivery modelling
- investigation of clot permeability and recanalisation
- multiscale modelling of thrombolytic treatment
- development and testing of new thrombolytic treatment strategies

The repository contains several historical model variants. Users should therefore check parameter definitions, units, boundary conditions, and numerical settings before applying the model to new scenarios.

## Related topics

`thrombolysis` · `fibrinolysis` · `stroke` · `tPA` · `pharmacokinetics` · `PKPD` · `drug-delivery` · `reaction-transport` · `Darcy-flow` · `biomedical-modelling` · `multiscale-modeling` · `MATLAB`

## License

No software license is currently specified.


## Contact

For questions regarding the model or code, please open an issue in this repository.# clotlysis
