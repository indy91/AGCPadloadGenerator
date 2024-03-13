# AGCPadloadGenerator

Program to generate erasable memory loads (called padloads) for the Apollo Guidance Computer.

Instructions for CMC and LGC:

First select a mission that is the closest to the desired custom mission to be flown. Then select the AGC version (rope). In case a custom rope is used then select the data set for the planetary inertial orientation subroutine (PIOS). The data sets can be found in the PIOSDataSets.txt file. The data sets are typically valid for a year ranging from July 1st to July 1st. An exception to this is the year 1970, for which the file has three versions (NBY1970_V1, NBY1970_V2, NBY1970_V3). NBY1970_V1 is the dataset flown on Apollo 11 in Comanche 55 and Luminary 99. These numbers were in error and got corrected for Apollo 12 and also used on Apollo 13 (NBY1970_V2). Custom ropes have been created using Artemis and Luminary 1E. These ropes have some parameters in fixed memory instead of erasable memory. For this the third data set (NBY1970_V3) for the same year got added.

Rope Constants:

New data for more launch years can be created with this tool as well. Use the rope constants button, enter the desired year and a RopeConstants.txt file will be written. This file contains everything that needs to be changed in the AGC source code to make the rope valid for the selected year. This data includes star vectors and Earth and Moon PIOS fixed constants. The LGC also has the solar and lunar epehemerides in fixed memory. These can not be calculated yet. So only Artemis 72 could be updated so far, no Luminary version.

The text file also contains the PIOS data set which can be added to the PIOSDataSets.txt file to support padloads for the modified rope.

More instructions to follow...