# **POLYRETINA Retinal Implant: Investigating the Impact of Field of View on Cognitive Load  and Performance**
 
## **1. Introduction**
This repository contains an EEG analysis pipeline developed for investigating the impact of varying fields of view on the effectiveness of the POLYRETINA retinal implant. The pipeline is based on the BeMoBIL framework and is designed to handle preprocessing, analysis, and visualization of EEG data. The POLYRETINA project aims to assess how different visual angles influence cognitive load and user performance, with the ultimate goal of enhancing autonomy for individuals with acquired blindness.

## **2. Installation**

Follow these steps to set up the project environment:

### Clone the Repository
To start, clone the repository to your local machine:

```console
git clone https://github.com/paulmoreau40/EEG-POLYRETINA.git && cd EEG-POLYRETINA
```

### Create a Virtual Environment
To create a virtual environment, run the following command:

```console
python -m venv polyretina
```

### Activate the Virtual Environment
Activate the virtual environment using one of the following commands based on your operating system:

**Windows:**
```console
.\polyretina\Scripts\activate
```
**Linux / Mac:**
```console
source polyretina/bin/activate
```

### Install Dependencies
Next, install the necessary dependencies by running:

```console
pip install -r requirements.txt
```



## **3. Project Structure**
Important codes:
- Behavioural data: **EEG_behavioural_analysis.ipynb** (time taken to do the task or "reaction time", and success rate) for each Field of View.
- EEG data:
    - **ConfigEEGPOL.m** (study configuration and parameters)
    - **main_preprocessing_POL.m** (XDF files import, data preparation, data preprocessing, ICA algorithm, Dipole Fitting, ICs autolabelling, ICs manual selection)
    - **keptComponents_POL.m** (list of Independent Components manually selected)
    - **main_single_participant_analysis_4Polyretina.m** (main analysis and statistics with data visualisation)
 
    - Data preparation and preprocessing:
        - **xdf2set.m**
        - **export_events_EEGPOL.m** (import and adapt the structure of events to create a table with the events name, type, time, latency, blocks and trials)
        - **events_check_EEGPOL.m** (check the integrity of the EEG data for each event) 

## **4. Authors, Credits, and Acknowledgments**
- Authors: Paul Moreau, Antonin Duret, Sandrine Hinrichs
- Research Labs:
    - [Aging in Vision and Action, Sorbonne Université / Institut de la Vison / INSERM / CNRS](https://www.institut-vision.org/en/research/aging-vision-and-action#:~:text=Our%20team%20analyzes%20the%20aging,research%20and%20innovative%20technology%20transfer.)
    - [Ophthalmic and neural technologies Laboratory](https://ghezzi-lab.org/)
    - [Ecole Polytechnique Fédérale de Lausanne, EPFL](https://www.epfl.ch/fr/)
- Credits: The EEG analysis pipeline is based on the BeMoBIL framework, with adaptations for the POLYRETINA field of view study.

## **5. References**
Bibliography: 
- [BeMoBIL: Brain Mobile Brain/Body Imaging Lab](https://www.tu.berlin/en/bpn/research/berlin-mobile-brain-body-imaging-lab)
- [POLYRETINA restores light responses in vivo in blind Göttingen minipigs (Vagni et al., 2022)](https://www.nature.com/articles/s41467-022-31180-z)


