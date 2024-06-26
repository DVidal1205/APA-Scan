# APA-Scan
APA-Scan detects the potential Alternative Polyadenylation (APA) events in the downstream 3'-UTR of a gene among two different biological conditions.

# CHANGES IN VERSION 2.0

## Paralellization

1. APA-Scan now supports parallel processing to speed up the quantification process. Users can specify the number of **cores** to use for parallel processing in the configuration.ini file. The default value is NULL, which will run the program in sequential mode. Users can also specify a MAX value to use all available cores on the machine.

## Change in Samtools

Debian distributions have issues with easily configuring execution permissions for the samtools binary. To address this issue, I have installed the samtools binary to my system in the /usr/local/bin directory. To accomodate this change, I had to modify samtools_directory in the methods.py file to /usr/local/bin/samtools. This change will allow the program to run on Debian distributions without any issues. (Note: the version of samtools downloaded on my system is identical to the one in the repository. I am looking into a more permanent solution for this issue to add in the documentation.)

# Installation
APA-Scan is a python tool which can be downloaded directly from github. Python (version 3.0 or higher) is required to be installed on users machine to run APA-Scan. It can work on Windows, Linux and Mac platform.

Users need to update the configuration.ini file and run the specific file APA-Scan.py to demonstrate the APA events.

$python3 APA-Scan.py

To generate plots for any selected 3'-UTR region, users need to run makePlots.py and provide input in a specific pattern: chromName:geneName:start-end

$python3 Make-plots.py

# User manual
User manual for APA-Scan is available on github https://github.com/compbiolabucf/APA-Scan/blob/master/APA-Scan%20User%20Manual.pdf

# Citation
Please use the following information to cite...

# Contact the Author
Naima Ahmed Fahmi: fnaima@knights.ucf.edu
Wei Zhang: wzhang.cs@ucf.edu
Jeongsik Yong: jyong@umn.edu
