# PlotBio V1.0

A tool to view biologic informations.

This tool, shows biologic informations, such as, GC content, width and dinucleotides, using Fasta files.
The user can provide up to 10 Fasta files.  

# AUTHOR/SUPPORT
Vitor Gregorio - vitorgregorio@alunos.utfpr.edu.br;</br>
* Alexandre Rossi Paschoal - paschoal@utfpr.edu.br;</br>

# HARDWARE/SOFTWARE REQUIREMENTS
Operating system(s): Linux/Unix or similar</br>
Programming language: Python3</br>
Library: Biopython, Matplotlib</br>
License: GNU GPL</br>

# EXAMPLE

In this example, there are six fasta files in the same folder. The user need to specify the number of files (-N).


python3 PlotBio.py -N 6 Arabidopsis_thaliana_GreeNC.fasta Arabidopsis_thaliana_CANTATA.fasta Oryza_sativa_Japonica_Group_GreeNC.fasta Oryza_sativa_Japonica_Group_CANTATA.fasta Zea_mays_GreeNC.fasta Zea_mays_CANTATA.fasta

The fasta files were get from CANTATAdb (<http://cantata.amu.edu.pl/>) and GreeNC (<http://greenc.sciencedesigners.com/wiki/Main_Page>)
