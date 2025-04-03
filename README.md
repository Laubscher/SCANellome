# SCANellome

SCANellome is a software application for Linux and macOS that allows users to perform anellome analysis of biological sequencing data in the FASTQ file format. The software is designed to be user-friendly and provides a simple interface for analyzing anellome.

# Installation

To install SCANellome on linux, simply download the executable file from https://laubscher.github.io/Anelloviruses/SCANellome .  
In the *Permisions* tab of the downloaded file check the box : *allow executing file as program*.

Generate axecutable from source code with pyinstaller: 
> pyinstaller --onefile --hidden-import=PIL._tkinter_finder SCANellome-2.x.x.py

 With files Dicodb.py and img.py in the same folder.

# Usage

To use SCANellome, first create or select a project. This will allow you to organize your data and analysis results in a convenient and structured manner. Once you have created or selected a project, you can import your FASTQ files from your computer.
Once you have configured the analysis parameters, click the "Run" button to begin the analysis.

# Author

SCANellome was created by Florian Laubscher. If you have any questions or feedback, you can contact the author at florian.laubscher@hcuge.ch

If you use SCANellome in your work, please cite:
>Laubscher, F., Kaiser, L., & Cordey, S. (2023). SCANellome: Analysis of the Genomic Diversity of Human and Non-Human Primate Anelloviruses from Metagenomics Data.
><br/>*Viruses*, **15**(7), 1575. [doi: 10.3390/v15071575][doi]

# License

SCANellome is distributed under the GNU GENERAL PUBLIC LICENSE Version 3. See the LICENSE file for more information.

# Disclaimer

SCANellome is provided "as is" without warranty of any kind.

[doi]: https://doi.org/10.3390/v15071575



