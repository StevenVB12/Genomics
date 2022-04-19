#!/bin/bash
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=72:00:00
#SBATCH --job-name=sweepE
#SBATCH --error=sweepE
#SBATCH --output=sweepE
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=1

module load python2
#module load jdk/11.0.7

ID=$((SLURM_ARRAY_TASK_ID -1))

scafs=(Herato0101 Herato0201 Herato0202 Herato0203 Herato0204 Herato0205 Herato0206 Herato0207 Herato0208 Herato0209 Herato0210 Herato0211 Herato0212 Herato0213 Herato0214 Herato0215 Herato0301 Herato0302 Herato0303 Herato0304 Herato0305 Herato0306 Herato0307 Herato0308 Herato0309 Herato0310 Herato0401 Herato0402 Herato0403 Herato0404 Herato0405 Herato0406 Herato0407 Herato0408 Herato0409 Herato0410 Herato0411 Herato0412 Herato0413 Herato0414 Herato0415 Herato0416 Herato0417 Herato0418 Herato0419 Herato0501 Herato0502 Herato0503 Herato0504 Herato0505 Herato0506 Herato0507 Herato0508 Herato0509 Herato0510 Herato0511 Herato0601 Herato0602 Herato0603 Herato0604 Herato0605 Herato0606 Herato0607 Herato0608 Herato0609 Herato0701 Herato0801 Herato0802 Herato0803 Herato0804 Herato0805 Herato0806 Herato0807 Herato0808 Herato0809 Herato0810 Herato0811 Herato0812 Herato0813 Herato0814 Herato0815 Herato0816 Herato0817 Herato0818 Herato0819 Herato0820 Herato0821 Herato0901 Herato0902 Herato0903 Herato0904 Herato1001 Herato1002 Herato1003 Herato1004 Herato1005 Herato1006 Herato1007 Herato1101 Herato1102 Herato1103 Herato1104 Herato1105 Herato1106 Herato1107 Herato1108 Herato1109 Herato1110 Herato1111 Herato1112 Herato1113 Herato1114 Herato1115 Herato1116 Herato1201 Herato1202 Herato1301 Herato1401 Herato1402 Herato1403 Herato1404 Herato1405 Herato1406 Herato1407 Herato1408 Herato1409 Herato1410 Herato1411 Herato1501 Herato1502 Herato1503 Herato1504 Herato1505 Herato1506 Herato1507 Herato1508 Herato1509 Herato1510 Herato1511 Herato1512 Herato1513 Herato1514 Herato1515 Herato1516 Herato1517 Herato1518 Herato1519 Herato1520 Herato1521 Herato1522 Herato1523 Herato1524 Herato1601 Herato1602 Herato1603 Herato1604 Herato1605 Herato1701 Herato1702 Herato1703 Herato1704 Herato1705 Herato1706 Herato1707 Herato1708 Herato1709 Herato1710 Herato1711 Herato1712 Herato1713 Herato1714 Herato1715 Herato1716 Herato1717 Herato1718 Herato1719 Herato1801 Herato1802 Herato1803 Herato1804 Herato1805 Herato1806 Herato1807 Herato1901 Herato1902 Herato1903 Herato1904 Herato1905 Herato1906 Herato1907 Herato1908 Herato1909 Herato1910 Herato2001 Herato2101)


/work/rpapa/share/programs/SF2/SweepFinder2 -sg 2000 /work/rpapa/share/sweepfinder/sweep_in_erato/demophoon.demophoon_$(echo "${scafs[ID]}").sweepfinder.input /work/rpapa/share/sweepfinder/sweep_out_erato/$(echo "${scafs[ID]}").sf2.out 