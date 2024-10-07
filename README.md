# **POLYRETINA Retinal Implant: Investigating the Impact of the Field of View on Cognitive Load  and Performance**
 
## **1. Introduction**
This repository contains an EEG analysis pipeline developed for investigating the impact of varying fields of view on the effectiveness of the POLYRETINA retinal implant. The pipeline is based on the BeMoBIL framework and is designed to handle preprocessing, analysis, and visualisation of EEG data. The POLYRETINA project aims to assess how different visual angles influence cognitive load and user performance, with the ultimate goal of enhancing autonomy for individuals with acquired blindness.

The project is conducted at the Institut de la Vision (Paris, FRANCE), within the Aging in Vision and Action Laboratory (AVA Lab).

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
python -m venv polyretina-env
```

### Activate the Virtual Environment
Activate the virtual environment using one of the following commands based on your operating system:

**Windows:**
```console
.\polyretina-env\Scripts\activate
```
**Linux / Mac:**
```console
source polyretina-env/bin/activate
```

### Install Dependencies
Next, install the necessary dependencies by running:

```console
pip install -r requirements.txt
```

### Additional Requirements

**EEGLAB 2024.0**:  
   - Download EEGLAB version 2024.0 from the official website: https://sccn.ucsd.edu/eeglab/download.php.
   - Unzip the downloaded archive.
   - Place the **EEGLAB** folder in the **toolboxes/** directory of the project:
     ```
     /toolboxes/eeglab2024.0/
     ```



## **3. Project Structure**

**Summary of the folders structure**
```
├── /data/
│   └── /analysis/             # Data for EEG analysis (raw, preprocessed, processed, analysis)
│   └── /behavioural/          # Data for behavioural analysis (trial time, success rate)
│
├── /polyretina/               # Functions and codes
│   ├── /behavioural/          # Contains Jupyter notebooks for behavioral analysis
│   │   └── analysis_notebook.ipynb
│   ├── /eeg/                  # Contains EEG-related functions and MATLAB scripts
│       ├── main_preprocessing.m
│       └── main_single_participant_analysis.m
│       └── ...
│
├── /figures/                 # Figures produced during the processing steps
│   ├── /AutoBadSampsCleaning/
│   ├── /ChannelsImport/
│   └── /PREP_distributions/
│
├── /toolboxes/                # Contains external toolboxes and dependencies
│   ├── /bemobil_pipeline0.2/  # Functions from the BeMoBIL pipeline
│   ├── /eeglab2024.0/         # EEGLAB toolbox for EEG analysis
│   └── /ParforProgMon/        # Parfor Progress Monitor toolbox
│
├── requirements.txt           # Python dependencies file
```

**Important codes for the analysis**
1. **Behavioural data**

    **`EEG_behavioural_analysis.ipynb`**: Analyses time taken to complete the task and success rate for each Field of View (20°/45°).

2. **EEG data**
   
    Mains:
    - **`main_preprocessing.m`**: Imports XDF files, prepares data and preprocesses it, performs ICA and Dipole Fitting, and labels Independent Components (ICs) automatically and manually.
    - **`main_single_participant_analysis.m`**: Conducts main analysis and statistics with data visualisation.

    Configuration:
    - **`ConfigEEGPOL.m`**: Contains study parameters (threshold filters, channels ELC filenames, types of data etc.). Can be modified according to the needs of the study.
    - **`getMainFoldersNames.m`**:  Manages paths to fit the directory structure (data, figures, etc.). Can be personalised if necessary.

    Data preparation and preprocessing:
    - **`xdf2set.m`**: Converts raw XDF files to EEGLAB .set format, processing EEG, motion capture and eye tracking data, events, and channels.
    - **`export_events_EEGPOL.m`**: Imports the events and creates a table with event names, types, latencies, blocks, and trials.
    - **`events_check_EEGPOL.m`**: Checks the integrity of the EEG data for each event.
    - **`keptComponents_POL.m`**: Lists Independent Components that were manually selected (filled by the user).


## **4. Authors, Credits, and Acknowledgments**
- Authors and collaborators: Paul Moreau, Antonin Duret, Sandrine Hinrichs, Denis Sheynikhovich, Angelo Arleo, Diego Ghezzi
- Research Labs:
    - [Aging in Vision and Action, Sorbonne Université / Institut de la Vison / INSERM / CNRS](https://www.institut-vision.org/en/research/aging-vision-and-action#:~:text=Our%20team%20analyzes%20the%20aging,research%20and%20innovative%20technology%20transfer.)
    - [Ophthalmic and neural technologies Laboratory](https://ghezzi-lab.org/)
    - [Ecole Polytechnique Fédérale de Lausanne, EPFL](https://www.epfl.ch/fr/)
- Credits: The EEG analysis pipeline is based on the [BeMoBIL framework](https://www.tu.berlin/en/bpn/research/berlin-mobile-brain-body-imaging-lab), with adaptations for the POLYRETINA field of view study.

## **5. References**
Bibliography: 
- [The BeMoBIL Pipeline for automated analyses of multimodal mobile brain and body imaging data (Klug et al., 2022)](https://www.biorxiv.org/content/10.1101/2022.09.29.510051v2.abstract)
- [POLYRETINA restores light responses in vivo in blind Göttingen minipigs (Vagni et al., 2022)](https://www.nature.com/articles/s41467-022-31180-z)


