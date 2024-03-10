#!/usr/bin/python3

#Modules
import os
import sys
import subprocess

#######################  OPTIONS  #######################
max_number_of_sequences=1000
availthreads = subprocess.check_output("nproc", shell=True).decode("utf-8").rstrip()
#availthreads = 64
rmkeywords = ("associated", "predicted", "isoform")


#######################  FUNCTIONS  #######################
#for calling bash commands for each individual sequence
def indivbash(bashline, outfileloci, outfileformat, dirorheader=1, goutfileformat =""):
  os.system("mkdir " + outfileloci)
  for header in seqdict.keys():
    with open("temp.fasta", "w") as infile:
      infile.write(seqdict.get(header))
    if (dirorheader == 1):
      try:
        headerfileformat = header.replace(">", "").replace(" ", "_").replace("[", "").replace("]", "").replace(",", "").replace(":", "")
        os.system(bashline + outfileloci + headerfileformat + outfileformat)
      except:
        print("Error with " + header)
    if (dirorheader == 2):
      try:
        headerfileformat = header.replace(">", "").replace(" ", "_").replace("[", "").replace("]", "").replace(",", "").replace(":", "")
        os.system(bashline + outfileloci + headerfileformat + outfileformat + " -goutfile " + outfileloci + headerfileformat + goutfileformat)
      except:
        print("Error with " + header)
    else:
      try:
        os.system(bashline + outfileloci + outfileformat)
      except:
        print("Error with " + header)


#######################  Folders  #######################
os.system("mkdir results")
os.system("mkdir temp")






##########################################################
#######################  Pipeline  #######################
##########################################################



#######################  1. Gather user input for taxonomic group and protein family  #######################
print("-------------------------------------------")
print("-------------------------------------------")

taxo = input("Enter Taxonomic Group:\n")
pfam = input("Enter Protein Family:\n")
print("-------------------------------------------")
while True:
  pfampartialflag = input("Remove partial sequences? (NOT PARTIAL) (y/n)\n")
  if (pfampartialflag == "y" or pfampartialflag == "n"):
    break
  else:
    print("Not valid response. Please retry")
  
print("-------------------------------------------")

while True:
  print("Remove sequences with keywords in rmkeywords? \n(" + str(" ".join(rmkeywords)) + ")")
  rmkeywordflag = input("Filter against keywords? (y/n)\n")
  if (rmkeywordflag == "y" or rmkeywordflag == "n"):
    break
  else:
    print("Not valid response. Please retry")


print("-------------------------------------------")
print("-------------------------------------------")



#######################  2. Gather desired protein sequences  #######################
#2a. Query for taxonID
print("Gathering taxonID for " + taxo + "\n...")
esearchTaxoquery = "esearch -db taxonomy -spell -query \"" + taxo + "\" | efetch -format uid"
try:
  esearchTaxoUID = subprocess.check_output(esearchTaxoquery, shell=True).decode("utf-8").rstrip() #QUERY TAXONID
except:
  print("Invalid taxonomic group (Please avoid plurals)\nExiting Program...")
  exit()

print("Done" + "\nTAXONID: " + str(esearchTaxoUID) )

if (esearchTaxoUID == ""): #error check
  print("Invalid taxonomic group (Please avoid plurals)\nExiting program...")
  exit()
  
print("-------------------------------------------")

#2b. Query for protein sequences filtered with taxonID with :exp to get taxo subtree groups
print("Gathering protein sequences of " + pfam + " in " + taxo + " txid:" +esearchTaxoUID + "\n...")

if pfampartialflag == "n":
  esearchProtquery = "esearch -db protein -spell -query \"txid" + esearchTaxoUID + "[Organism:exp]" + " AND " + pfam + "[PROT]" + "\" | efetch -format fasta > seq.fasta"
elif pfampartialflag == "y":
  esearchProtquery = "esearch -db protein -spell -query \"txid" + esearchTaxoUID + "[Organism:exp]" + " AND " + pfam + "[PROT]" + " NOT PARTIAL\" | efetch -format fasta > seq.fasta"

if (pfam[-1].lower() == "s"): #plural catch
  print("PLURAL DETECTED IN PROTEIN FAMILY QUERY\nQuery protein family may be in plural form which may affect query sequences for a 0 hit result. Suggested to remove plural and write singular protein family. (eg. ABC transporter rather than ABC transporters)")
  usercont = input("Do you wish to continue? (y/n)")
  if (usercont == "y"):
    print("Continuing...")
  elif (usercont == "n"):
    print("Exiting...")
    exit()
  else:
    print("Not valid response (y/n)\nExiting program")
    exit()

try: #error catch
  esearchProtfasta = os.system(esearchProtquery)
except:
  print("Invalid protein family (Please avoid plurals)\n Exiting program")
  exit()

print("Done")

#2c. Give user more information about queried sequences (no. of seq)
with open("seq.fasta") as infile:
  infileread = infile.read()
  seqcount = infileread.count(">")

print("Number of sequences gathered: " + str(seqcount))

if (seqcount == 0): #error check for empty or >maxseqcount seqdata
  print("Invalid protein family (Please avoid plurals)\nExiting program...")
  exit()
elif (seqcount > max_number_of_sequences):
  print("Exceeds threshold max number of sequences, " + str(max_number_of_sequences))
  maxnumbercheck = input("Do you wish to continue with " + seqcount + " sequences? (y/n)")
  if (maxnumbercheck ==  "y"):
    print("Continuing...")
  elif(maxnumbercheck == "n"):
    print("Exiting...")
    exit()
  else:
    print("Not valid response (y/n)\nExiting Program")
    exit()

print("-------------------------------------------")



#######################  3. Sequence data prep  #######################
print("Sequence data preparation\n...")
#3a. Separate sequences in seq.fasta
with open("seq.fasta") as infile:
  allseq = infile.read().rstrip().split(">")
  allseq.pop(0)
  allseq1 = [">" + seqelement for seqelement in allseq]

#3bi. Create list of headers
seqheaders = []
for myseq in allseq1:
  headerend = myseq.find("\n")
  header = myseq[:headerend]
  seqheaders.append(header)
  
#3bii. Create dictionary of header:fasta
seqdict = {}
i=0
for header in seqheaders:
  seqdict[header] = allseq1[i]
  i=i+1

#3d. Remove headers with keywords
open("finalseq.fasta", "w").close()

if rmkeywordflag == "y":
  print("Removing keyword sequences")
  for rmkeyword in rmkeywords:
    for header in seqheaders:
      if (header.lower().find(rmkeyword) > -1):
        del seqdict[header]
elif rmkeywordflag == "n":
  print("Keeping keyword sequences")


for header in seqdict.keys():
  with open("finalseq.fasta", "a") as outfile:
    outfile.write(seqdict[header])
    
with open("finalseq.fasta") as infile:
  infileread = infile.read()
  seqcountfinal = infileread.count(">")

print(str(seqcount - seqcountfinal) + " sequences removed")

#3d. Get list of species
seqspecies = []
for header in seqdict.keys():
  startpos = header.find("[")
  endpos = header.find("]")
  seqspecies.append(header[startpos+1:endpos])

print("Dataset contains sequences from " + str(len(set(seqspecies))) + " unique species")
print("Final dataset contains " + str(seqcountfinal) + " sequences")

print("Done")
print("-------------------------------------------")


#######################  SEQUENCE ANALYSIS  #######################
##4. CLUSTALO for alignment, EMBOSS for plotcon for sequence conservation plot
#4a. ClustalO
print("Aligning sequences via ClustalO with: " + availthreads + " threads\n...")
os.system("clustalo -i finalseq.fasta -o ./results/aligned.fasta --force --threads=" + str(availthreads))
print("Done")
print("-------------------------------------------")

#4b.  plotcon
print("Plotting convservation of sequence alignment\n...")
os.system("plotcon -sequences ./results/aligned.fasta -winsize=4 -graph png -sprotein1 -gdirectory results")
print("Done")
print("-------------------------------------------")

##5 Run patmatmotifs on each seq in allseq1 (Scan for PROSITE motifs in sequences patmatmotifs)
print("Scanning sequences for motifs\n...")
indivbash("patmatmotifs -auto -full -raccshow2 -rstrandshow2 -rusashow2 -rdesshow2 -rscoreshow2 -sequence temp.fasta -outfile ", "./results/motifs/", ".patmatmotifs", 1)
print("Done - Results in ./results/motifs/")
print("-------------------------------------------")

##6. EMBOSS Analysis 1 - sigcleave - signal sequence cleavage site
print("Scanning sequences for signal peptide cleavage sites\n...")
indivbash("sigcleave -rdesshow2 -rscoreshow2 -rusashow2 -auto -sequence temp.fasta -outfile " , "./results/sigcleave/", ".sigcleave", 1)
print("Done - Results in ./results/sigcleave/")
print("-------------------------------------------")

##7. EMBOSS Analysis 2 - charge - protein charge plot
print("Gathering protein charge\n...")
indivbash("charge -auto -seqall temp.fasta -outfile ", "./results/charge/", ".charge", 1)
print("Done - Results in ./results/charge/")
print("-------------------------------------------")

##8. EMBOSS Analysis 3 - freak - resuidue frequency plot
print("Plotting residue frequency\n...")
indivbash("freak -auto -graph svg -seqall temp.fasta -odirectory ", "./results/freak/", "", 0)
print("Done - Results in ./results/freak/")
print("-------------------------------------------")

##9. EMBOSS Analysis 4 - helixturnhelix - helix turn helix motif searching for nucleic acid binding motifs
print("Scanning for nucleic acid binding site motifs\n...")
indivbash("helixturnhelix -sprotein1 -warning FALSE -rdesshow2 -auto -sequence temp.fasta -outfile ", "./results/helixturnhelix/", ".helixturnhelix", 1)
print("Done - Results in ./results/helixturnhelix/")
print("-------------------------------------------")

##10. EMBOSS Analysis 5 - hmoment - calculate and plot hydrophobic moment 
print("Plotting hydrophobic moments\n...")
indivbash("hmoment -auto -seqall temp.fasta -outfile ", "./results/hmoment/", ".hmoment", 1)
print("Done - Results in ./results/hmoment/")
print("-------------------------------------------")

##11. EMBOSS Analysis 6 - iep - isoelectric point
print("Calculating isoelectric point\n...")
indivbash("iep -auto -sequence temp.fasta -outfile ", "./results/iep/", ".iep", 1)
print("Done - Results in ./results/iep/")
print("-------------------------------------------")

##12. EMBOSS Analysis 7 - pepcoil - coiled coil regions
print("Predicting coiled coil regions\n...")
indivbash("pepcoil -auto -rdesshow2 -sequence temp.fasta -outfile ", "./results/pepcoil/", ".pepcoil", 1)
print("Done - Results in ./results/pepcoil/")
print("-------------------------------------------")

##13. EMBOSS Analysis 9 - pepstats - Molecular weight, Number of residues ,Average residue weight, Charge ,Isoelectric point,For each type of amino acid: number, molar percent, DayhoffStat,For each physico-chemical class of amino acid: number, molar percent,Probability of protein expression in E. coli inclusion bodies,Molar extinction coefficient (A280),Extinction coefficient at 1 mg/ml (A280)
print("Gathering pepstats data\n...")
indivbash("pepstats -auto -sequence temp.fasta -outfile ", "./results/pepstats/", ".pepstats", 1)
print("Done - Results in ./results/pepstats/")
print("-------------------------------------------")
