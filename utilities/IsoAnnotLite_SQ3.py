#Script to generate a GFF3 file from SQANTI3 output and using a tappAS GFF3 as reference.

import argparse
import math
from pickle import FALSE
import sys
import os
import bisect

#Global Variables
verbose = False
USE_GFF3 = False
USE_NAME = False
USE_STDOUT = False
ALL_AS_NOVELS = False
INTRONIC = True
STATS = False
SAVE_PROB_TRANSCRIPTS = False

version = "2.7.3"
CLASS_COLUMN_USED = [0,1,2,3,5,6,7,30,32,33]
CLASS_COLUMN_NAME = ["isoform", "chrom", "strand", "length", "structural_category", "associated_gene", "associated_transcript", "ORF_length","CDS_start", "CDS_end"]

LST_TRANSCRIPTFEATURES_NOTIN_CDS = ["uORF", "miRNA_Binding", "PAS", "3UTRmotif", "5UTRmotif"]
LST_TRANSCRIPTSOURCES_INTRONIC = ["PAR-clip"]

#Functions
def createGTFFromSqanti(file_exons, file_trans, file_junct, dc_gene_description, filename):
    global verbose
    res = open(filename,"w+")
    source = "tappAS"
    feature = ""
    start = ""
    end = ""
    aux = "."
    strand = ""
    desc = ""

    dc_coding = {}
    dc_gene = {}
    dc_SQstrand = {}
    f = open(file_trans)

    #check header
    global CLASS_COLUMN_USED
    global CLASS_COLUMN_NAME

    header = next(f)
    fields = header.split("\t")
    #index = 0
    #for column in CLASS_COLUMN_NAME: #check all the columns we used
    #    if column not in fields[CLASS_COLUMN_USED[index]]: #if now in the correct possition...
    #        print("File classification does not have the correct structure. The column \"" + column + "\" is not in the possition " + str(CLASS_COLUMN_USED[index]) + " in the classification file. We have found the column \"" + str(fields[CLASS_COLUMN_USED[index]]) + "\".")
    #        sys.exit()
    #    else:
    #        index = index + 1

    for column in CLASS_COLUMN_NAME: #check all the columns we used
        if column not in fields: #we do not find the column
            print("We can not find the column \"" + column + "\" is the SQANTI file.")
            sys.exit()
        else:
            #update position
            CLASS_COLUMN_USED[CLASS_COLUMN_NAME.index(column)] = fields.index(column)

    #positions for each column
    sq_isoform = CLASS_COLUMN_USED[CLASS_COLUMN_NAME.index("isoform")]
    sq_chrom = CLASS_COLUMN_USED[CLASS_COLUMN_NAME.index("chrom")]
    sq_strand = CLASS_COLUMN_USED[CLASS_COLUMN_NAME.index("strand")]
    sq_length = CLASS_COLUMN_USED[CLASS_COLUMN_NAME.index("length")]
    sq_structural_category = CLASS_COLUMN_USED[CLASS_COLUMN_NAME.index("structural_category")]
    sq_associated_gene = CLASS_COLUMN_USED[CLASS_COLUMN_NAME.index("associated_gene")]
    sq_associated_transcript = CLASS_COLUMN_USED[CLASS_COLUMN_NAME.index("associated_transcript")]
    sq_ORF_length = CLASS_COLUMN_USED[CLASS_COLUMN_NAME.index("ORF_length")]
    sq_CDS_start = CLASS_COLUMN_USED[CLASS_COLUMN_NAME.index("CDS_start")]
    sq_CDS_end = CLASS_COLUMN_USED[CLASS_COLUMN_NAME.index("CDS_end")]

    #add transcript, gene and CDS
    for line in f:
        fields = line.split("\t")

        #trans
        transcript = fields[sq_isoform]
        #source        
        feature = "transcript"
        start = "1"
        end = fields[sq_length]
        #aux
        strand = fields[sq_strand]

        dc_SQstrand.update({str(transcript) : strand}) #saving strand

        desc = "ID="+fields[sq_associated_transcript]+"; primary_class="+fields[sq_structural_category]+"\n"
        res.write("\t".join([transcript, source, feature, start, end, aux, strand, aux, desc]))
        #gene
        transcript = fields[sq_isoform]
        #source        
        feature = "gene"
        start = "1"
        end = fields[sq_length]
        #aux
        strand = fields[sq_strand]
        if dc_gene_description.get(fields[sq_associated_gene]):
            desc = dc_gene_description.get(fields[sq_associated_gene]) #gene name
        else:
            desc = "ID="+fields[sq_associated_gene]+"; Name="+fields[sq_associated_gene]+"; Desc=" + fields[sq_associated_gene] + "\n"
        res.write("\t".join([transcript, source, feature, start, end, aux, strand, aux, desc]))
        #CDS
        transcript = fields[sq_isoform]
        #source        
        feature = "CDS"
        start = fields[sq_CDS_start] #30
        end = fields[sq_CDS_end] #31
        #aux
        strand = fields[sq_strand]
        desc = "ID=Protein_"+transcript+"; Name=Protein_"+transcript+"; Desc=Protein_"+transcript+"\n"
        if start != "NA":
            res.write("\t".join([transcript, source, feature, start, end, aux, strand, aux, desc]))
            res.write("\t".join([transcript, source, "protein", "1", str(int(math.ceil((int(end)-int(start)-1)/3))), aux, strand, aux, desc]))
        else:
            res.write("\t".join([transcript, source, feature, ".", ".", aux, strand, aux, desc]))
        #genomic    
        desc = "Chr="+fields[sq_chrom]+"\n"
        
        #Gene
        gene = fields[sq_associated_gene]
        category = fields[sq_structural_category]
        transAssociated = fields[sq_associated_transcript]

        if transAssociated.startswith("ENS"):
            transAssociated = transAssociated.split(".") #ENSMUS213123.1 -> #ENSMUS213123
            transAssociated = transAssociated[0]

        if(not dc_gene.get(transcript)):
            dc_gene.update({str(transcript) : [gene, category, transAssociated]})
        else:
            dc_gene.update({str(transcript) : dc_gene.get(transcript) + [gene, category, transAssociated]})
        
        #Coding Dictionary
        CDSstart = fields[sq_CDS_start] #30
        CDSend = fields[sq_CDS_end] #31
        orf = fields[sq_ORF_length] #28
        
        if(not dc_coding.get(transcript)):
            dc_coding.update({str(transcript) : [CDSstart, CDSend, orf]})
        else:
            dc_coding.update({str(transcript) : dc_coding.get(transcript) + [CDSstart, CDSend, orf]})

        res.write("\t".join([transcript, source, "genomic", "1", "1", aux, strand, aux, desc]))

        #Write TranscriptAttributes
        sourceAux = "TranscriptAttributes"
        lengthTranscript = fields[sq_length]
        if not CDSstart == "NA":
            #3'UTR
            feature = "3UTR_Length"
            start = int(CDSend) + 1
            end = lengthTranscript
            desc = "ID=3UTR_Length; Name=3UTR_Length; Desc=3UTR_Length\n"
            res.write("\t".join([transcript, sourceAux, feature, str(start), str(end), aux, strand, aux, desc]))
            #5'UTR
            feature = "5UTR_Length"
            start = 1
            end = int(fields[sq_CDS_start]) - 1 + 1 #30
            desc = "ID=5UTR_Length; Name=5UTR_Length; Desc=5UTR_Length\n"
            res.write("\t".join([transcript, sourceAux, feature, str(start), str(end), aux, strand, aux, desc]))
            #CDS
            feature = "CDS"
            start = CDSstart
            end = CDSend
            desc = "ID=CDS; Name=CDS; Desc=CDS\n"
            res.write("\t".join([transcript, sourceAux, feature, str(start), str(end), aux, strand, aux, desc]))
            #polyA
            feature = "polyA_Site"
            start = lengthTranscript
            end = lengthTranscript
            desc = "ID=polyA_Site; Name=polyA_Site; Desc=polyA_Site\n"
            res.write("\t".join([transcript, sourceAux, feature, str(start), str(end), aux, strand, aux, desc]))

    f.close()

    f = open(file_exons)
    dc_exons = {}
    dc_geneID2geneName = {}
    dc_geneName2geneID = {}
    #add exons
    for line in f:
        fields = line.split("\t")
        if len(fields) == 9:
            transcript = fields[8].split('"')[1].strip()
            #source
            feature = fields[sq_strand]
            if(feature == "transcript"): #just want exons
                split_length = len(fields[8].split('"'))
                gID = fields[8].split('"')[3].strip()
                if split_length>5: #it can be possible some transcripts do not have Gene Name
                    gName = fields[8].split('"')[5].strip()
                else:
                    gName = ""

                # if verbose:
                #     print("Reading gene ID and gene Name from junction file: " + gID + " and " + gName)

                #GeneNameID Dictionary
                if(not dc_geneID2geneName.get(gID)):
                    dc_geneID2geneName.update({str(gID) : gName}) #geneName can be NA
                else:
                    dc_geneID2geneName.update({str(gID) : gName})

                if(not dc_geneName2geneID.get(gName)):
                    dc_geneName2geneID.update({str(gName) : gID})
                else:
                    dc_geneName2geneID.update({str(gName) : gID})

                continue

            start = int(fields[sq_length])
            end = int(fields[4])
            #aux
            strand = fields[sq_associated_gene]
            #desc = fields[8]
            desc = "Chr=" + str(fields[sq_isoform]) + "\n"

            #Exons Dictionary
            if(not dc_exons.get(transcript)):
                dc_exons.update({str(transcript) : [[start,end]]})
            else:
                dc_exons.update({str(transcript) : dc_exons.get(transcript) + [[start,end]]})
            
            res.write("\t".join([transcript, source, feature, str(start), str(end), aux, strand, aux, desc]))
        else:
            print("File corrected doesn't have the correct number of columns (9).")
    f.close()

    #add junctions
    f = open(file_junct)
    #header
    header = next(f)
    for line in f:
        fields = line.split("\t")
        #Junctions file can have a dvierse number of columns, not only 19 but 0-14 are allways the same
        transcript = fields[sq_isoform]
        #source        
        feature = "splice_junction"
        start = fields[4]
        end = fields[sq_structural_category]
        #aux
        strand = fields[sq_strand]
        desc = "ID="+fields[sq_length]+"_"+fields[14]+"; Chr="+fields[sq_chrom]+"\n"

        res.write("\t".join([transcript, source, feature, start, end, aux, strand, aux, desc]))
    f.close()
    res.close()

    return dc_exons, dc_coding, dc_gene, dc_SQstrand, dc_geneID2geneName, dc_geneName2geneID

def createDCgeneTrans(dc_SQtransGene):
    global verbose
    #create dc_SQgeneTrans
    dc_SQgeneTrans = {}
    for trans, values in dc_SQtransGene.items():
        gene = values[0] #gene
        if(not dc_SQgeneTrans.get(gene)):
            dc_SQgeneTrans.update({str(gene) : [trans]})
        else:
            dc_SQgeneTrans.update({str(gene) : dc_SQgeneTrans.get(gene) + [trans]})

    return dc_SQgeneTrans

def readGFF(gff3):
    global verbose
    f = open(gff3)
    #create dictionary for each transcript and dictionary for exons
    dc_GFF3 = {}
    dc_GFF3_raw_annot = {}
    dc_GFF3exonsTrans = {}
    dc_GFF3transExons = {}
    dc_GFF3coding = {}
    dc_GFF3strand = {}
    dc_GFF3geneTrans = {}
    dc_GFF3transGene = {}
    dc_GFF3geneNameTrans = {}
    dc_gene_description = {}
    for line in f:
        fields = line.split("\t")
        if len(fields) == 9:
            #feature (transcript, gene, exons...)
            transcript = fields[0]
            feature = fields[2]
            start = fields[3]
            end = fields[4]
            strand = fields[6]

            if not strand == ".":
                dc_GFF3strand.update({str(transcript) : strand}) #saving strand

            text = fields[8].split(" ")
            if not text[-1].endswith("\n"):
                line = line + "\n"

            if feature == "gene":
                g = text[0][3:-1] #delete "ID=" and final ";" == GENE
                #save description
                if not dc_gene_description.get(str(g)):
                    dc_gene_description.update({str(g) : fields[8]})
                else:                
                    dc_gene_description.update({str(g) : fields[8]})

                #keep the Gene ID:
                if not dc_GFF3geneTrans.get(str(g)):
                    dc_GFF3geneTrans.update({str(g) : [transcript]})
                else:                
                    dc_GFF3geneTrans.update({str(g) : dc_GFF3geneTrans.get(str(g)) + [transcript]})

                #trans - trans-Gene
                if not dc_GFF3transGene.get(str(transcript)):
                    dc_GFF3transGene.update({str(transcript) : str(g)})
                else:                
                    dc_GFF3transGene.update({str(transcript) : str(g)})

                #keep the Gene NAME:
                g_name = text[1][5:-1] #delete "NAME=" and final ";" == GENE

                # if verbose:
                #     print(g_name)

                if not dc_GFF3geneNameTrans.get(str(g_name)):
                    dc_GFF3geneNameTrans.update({str(g_name) : [transcript]})
                else:                
                    dc_GFF3geneNameTrans.update({str(g_name) : dc_GFF3geneNameTrans.get(str(g_name)) + [transcript]})

            if feature == "exon":
                if not dc_GFF3transExons.get(str(transcript)):
                    dc_GFF3transExons.update({str(transcript) : [[int(start), int(end)]]})
                else:                
                    dc_GFF3transExons.update({str(transcript) : dc_GFF3transExons.get(str(transcript)) + [[int(start), int(end)]]})

                if not dc_GFF3exonsTrans.get(int(start)):
                    dc_GFF3exonsTrans.update({int(start) : [transcript]})
                else:
                    dc_GFF3exonsTrans.update({int(start) : dc_GFF3exonsTrans.get(int(start)) + [transcript]})
            elif feature == "CDS":
                if not dc_GFF3coding.get(str(transcript)):
                    dc_GFF3coding.update({str(transcript) : [int(start), int(end)]}) #not int bc can be NA
                else:                
                    dc_GFF3coding.update({str(transcript) : dc_GFF3coding.get(str(transcript)) + [int(start), int(end)]})

            elif feature in ["splice_junction","transcript","gene","protein","genomic"] :
                if not dc_GFF3_raw_annot.get(transcript):
                    dc_GFF3_raw_annot.update({str(transcript) : [line]})
                else:
                    dc_GFF3_raw_annot.update({str(transcript) : dc_GFF3_raw_annot.get(transcript) + [line]})

            else:
                if not dc_GFF3.get(transcript):
                    dc_GFF3.update({str(transcript) : [[start, end, line]]})
                else:
                    dc_GFF3.update({str(transcript) : dc_GFF3.get(transcript) + [[start, end, line]]})
        else:
            print("File GFF3 doesn't have the correct number of columns (9).")

    sorted(dc_GFF3exonsTrans.keys())
    return dc_GFF3, dc_GFF3exonsTrans, dc_GFF3transExons, dc_GFF3coding, dc_GFF3strand, dc_GFF3geneTrans, dc_GFF3transGene, dc_GFF3geneNameTrans, dc_GFF3_raw_annot, dc_gene_description

def unique(list1): 
    global verbose
    # intilize a null list 
    unique_list = [] 
      
    # traverse for all elements 
    for x in list1: 
        # check if exists in unique_list or not 
        if x not in unique_list: 
            unique_list.append(x) 

def isGenomicPosition(start, end):
    global verbose
    res = False
    if int(start)>=100000 or int(end)>=100000:
        res = True
    
    return res

def transformTransFeaturesToGenomic(dc_GFF3, dc_GFF3transExons, dc_GFF3coding, dc_GFF3strand):
    global verbose
    newdc_GFF3 = {}
    bnegative = False
    transcriptsAnnotated = 0

    for trans in dc_GFF3transExons.keys():
        #Print Length Transformation
        transcriptsAnnotated = transcriptsAnnotated + 1
        perct = transcriptsAnnotated/len(dc_GFF3transExons)*100
        print("\t" + "%.2f" % perct + " % of features transformed...", end='\r')

        if not dc_GFF3.get(trans):
            continue

        annot = dc_GFF3.get(trans)

        for values in annot:
            bProt = False
            line = values[2]
            fields = line.split("\t")
            text = fields[8].split(" ")
            strand = dc_GFF3strand.get(trans)

            start = 0
            end = 0
            startG = 0
            endG = 0

            fields = values[2].split("\t")

            #Transcript calculate normal - include CDS
            if text[-1].endswith("T\n") and not fields[3] == ".":
                start = int(fields[3])
                end = int(fields[4])
            elif text[-1].endswith("P\n") and not fields[3] == ".":
                start = int(fields[3])
                end = int(fields[4])
                bProt = True
                #amoniacids to nucleotics
                start = start * 3
                end = end * 3
            else:
                if not newdc_GFF3.get(trans):
                    newdc_GFF3.update({str(trans) : [values]})
                    continue
                else:
                    newdc_GFF3.update({str(trans) : newdc_GFF3.get(trans) + [values]})
                    continue

	        # PAR-clip | Genomic as Transcript
            alreadyGenomic = False
            alreadyGenomic = isGenomicPosition(start, end)

            if alreadyGenomic:
                if not newdc_GFF3.get(trans):
                    newdc_GFF3.update({str(trans) : [values]})
                    continue
                else:
                    newdc_GFF3.update({str(trans) : newdc_GFF3.get(trans) + [values]})
                    continue

            if bProt: #prot
                totalDiff = end - start + 3 #(total diff can be 0 for features in the same location)
            else:
                totalDiff = end - start + 1 #(total diff can be 0 for features in the same location)

            if not bProt:
                allExons = dc_GFF3transExons.get(trans)
            else:
                allExons = dc_GFF3coding.get(trans)
                if not allExons:
                    continue

            if strand == "+":
                allExons = sorted(allExons)
            else:
                allExons = sorted(allExons, reverse = True)

            bstart = False
            bend = False

            for exon in allExons:
                if totalDiff < 0:
                    bnegative = True
                    break

                #START already found
                if bstart:
                    if strand == "+":
                        if exon[0]+totalDiff-1 <= exon[1]: #pos ends here
                            endG = exon[0]+totalDiff-1
                            bend = True
                        else: #pos ends in other exon and we add the final exon
                            totalDiff = totalDiff - (exon[1]-exon[0]+1)
                    else:
                        if exon[1]-totalDiff+1 >= exon[0]: #pos ends here
                            endG = exon[1]-totalDiff+1
                            bend = True
                        else: #pos ends in other exon and we add the final exon
                            totalDiff = totalDiff - (exon[1]-exon[0]+1)

                #Search for START
                if exon[1]-exon[0]+1 >= start and not bstart: #pos starts here
                    if strand == "+":
                        startG = exon[0]+int(start)-1
                        bstart = True
                        if startG+totalDiff-1 <= exon[1]: #pos ends here
                            endG = startG+totalDiff-1
                            bend = True
                        else: #pos ends in other exon and we add the final exon
                            totalDiff = totalDiff - (exon[1]-startG+1)
                    else:
                        startG = exon[1]-int(start)+1
                        bstart = True
                        if startG-totalDiff+1 >= exon[0]: #pos ends here
                            endG = startG-totalDiff+1
                            bend = True
                        else: #pos ends in other exon and we add the final exon
                            totalDiff = totalDiff - (startG-exon[0]+1)
                else:
                    #not in first exon, update the start and end pos substrating exon length
                    start = start - (exon[1]-exon[0]+1)
                    end = end - (exon[1]-exon[0]+1)
                    #dist = (exon[1]-exon[0]+1)
                if bend:
                    if strand == "-": #reorder
                        aux = startG
                        startG = endG
                        endG = aux
                    if not bProt:
                        newline = fields[0] + "\t" + fields[1] + "\t" + fields[2] + "\t" + str(startG) + "\t" + str(endG) + "\t" + fields[5] + "\t" + fields[6] + "\t" + fields[7] + "\t" + fields[8]
                        if not newdc_GFF3.get(trans):
                            newdc_GFF3.update({str(trans) : [[startG, endG, newline]]})
                            break
                        else:
                            newdc_GFF3.update({str(trans) : newdc_GFF3.get(trans) + [[startG, endG, newline]]})
                            break
                    else:
                        if not newdc_GFF3.get(trans):
                            newdc_GFF3.update({str(trans) : [[startG, endG, values[2]]]})
                            break
                        else:
                            newdc_GFF3.update({str(trans) : newdc_GFF3.get(trans) + [[startG, endG, values[2]]]})
                            break
            
            if bnegative:
                break

    return newdc_GFF3

def transformTransFeaturesToLocale(dc_GFF3, dc_SQexons):
    global verbose
    dc_newGFF3 = {}
    for trans in dc_GFF3.keys():
        annot = dc_GFF3.get(trans)
        line = annot[0]
        line = line.split("\t")
        strand = line[6]

        exons = dc_SQexons.get(trans)
        if strand == "+":
            exons = sorted(exons)
        else:
            exons = sorted(exons, reverse = True)

        start = 0
        end = 0
        for line in annot:
            fields = line.split("\t")
            text = fields[8].split(" ")
            if fields[1] == "tappAS":
                continue

            elif text[-1].endswith("T\n"):
                isGenomic = False
                if(not fields[3] == "." and not fields[4] == "."):
                    isGenomic = isGenomicPosition(fields[3], fields[4])

                if not isGenomic: #not need the transformation to a local position
                    if not dc_newGFF3.get(trans):
                        dc_newGFF3.update({str(trans) : [line]})
                        continue
                    else:
                        dc_newGFF3.update({str(trans) : dc_newGFF3.get(trans) + [line]})
                        continue

                if strand == "+":
                    startG = fields[3]
                    endG = fields[4]
                else:
                    startG = fields[4] #end
                    endG = fields[3] #start
                bstart = False
                bend = False
                distance = 0 #other exons
                for ex in exons:
                    if not startG=="." or not endG==".":
                        #SEARCH FOR START
                        if int(ex[0])<=int(startG) and int(startG)<=int(ex[1]) and not bstart: #start
                            if strand == "+":
                                start = (int(startG)-int(ex[0]) + 1) + distance
                                bstart = True
                                if int(ex[0])<=int(endG) and int(endG)<=int(ex[1]):
                                    end = start + (int(endG)-int(startG) + 1) - 1
                                    bend = True
                                    break
                                else:
                                    distance = int(ex[1]) - int(startG) + 1
                                    continue
                            else: #negative strand
                                start = (int(ex[1])-int(startG) + 1) + distance
                                bstart = True
                                if int(ex[0])<=int(endG) and int(endG)<=int(ex[1]):
                                    end = start + (int(startG)-int(endG) + 1) - 1
                                    bend = True
                                    break
                                else:
                                    distance = int(startG)-int(ex[0]) + 1
                                    continue
                                
                        elif not bstart:
                            distance = distance + (int(ex[1])-int(ex[0]) + 1)

                        #SEARCH FOR END
                        if bstart:
                            if int(ex[0])<=int(endG) and int(endG)<=int(ex[1]):
                                if strand == "+":
                                    end = start + distance + (int(endG)-int(ex[0]) + 1) -1 #no afecta el -1
                                    bend = True
                                    break
                                else:
                                    end = start + distance + (int(ex[1]-int(endG)) + 1) -1 #no afecta el -1
                                    bend = True
                                    break
                            else:
                                distance = distance + (int(ex[1]) - int(ex[0]) + 1)
                    else:
                        start = startG
                        end = endG
                        bend = True
                        break
                if bend: #to be sure in full-spliced match cases
                    newline = fields[0] + "\t" + fields[1] + "\t" + fields[2] + "\t" + str(start) + "\t" + str(end) + "\t" + fields[5] + "\t" + fields[6] + "\t" + fields[7] + "\t" + fields[8]
                    if not dc_newGFF3.get(trans):
                        dc_newGFF3.update({str(trans) : [newline]})
                    else:
                        dc_newGFF3.update({str(trans) : dc_newGFF3.get(trans) + [newline]})
            else:
                if not dc_newGFF3.get(trans):
                    dc_newGFF3.update({str(trans) : [line]})
                else:
                    dc_newGFF3.update({str(trans) : dc_newGFF3.get(trans) + [line]})
    return dc_newGFF3

def transformProtFeaturesToLocale(dc_GFF3, dc_SQexons, dc_SQcoding):
    transInterest = "NM_001081080.1"
    featureOfInterest = "DOMAIN"
    global verbose

    dc_newGFF3 = {}
    for trans in dc_GFF3.keys():
        annot = dc_GFF3.get(trans)
        line = annot[0]
        line = line.split("\t")
        strand = line[6]

        exons = dc_SQexons.get(trans)
        if strand == "+":
            exons = sorted(exons)
        else:
            exons = sorted(exons, reverse = True)
        
        annot = dc_GFF3.get(trans)

        start = 0
        end = 0
        if not dc_SQcoding.get(trans):
            continue

        startcoding = dc_SQcoding.get(trans)[0]
        startcoding = startcoding[0]
        if startcoding=="NA":
            continue
            
        for line in annot:
            fields = line.split("\t")
            text = fields[8].split(" ")

            if fields[1] == "tappAS":
                continue
            if text[-1].endswith("P\n"):
                startG = fields[3]
                endG = fields[4]

                if trans == transInterest and fields[2]==featureOfInterest and verbose:
                    print(trans)
                    print(fields[2])
                    print("****originalLine****")
                    print(line)

                bstart = False
                CDSstart = False
                distance = 0 #other exons
                for ex in exons:
                    if not startG=="." or not endG==".":
                        if not CDSstart: #CDS start
                            if strand == "+":
                                if int(ex[0])<=int(startcoding) and int(startcoding)<=int(ex[1]) and not bstart: #start
                                    start = int(startG)-int(ex[0]) + 1 + distance #CDSstart
                                    CDSstart = True
                                else:
                                    distance = distance + int(ex[1]) - int(startG) + 1
                            else: #egative strand
                                if int(ex[0])<=int(startcoding) and int(startcoding)<=int(ex[1]) and not bstart: #start
                                    start = int(ex[1])-int(startG) + 1 + distance #CDSstart
                                    CDSstart = True
                                else:
                                    distance = distance + int(startG) - int(ex[0]) + 1

                        if int(ex[0])<=int(startG) and int(startG)<=int(ex[1]) and CDSstart and not bstart: #start
                            if strand == "+":
                                start = int(startG)-int(start) + 1 + distance #diff between genomic pos and CDSstart
                                bstart = True
                                if int(endG)>=int(ex[0]) and int(endG)<=int(ex[1]):
                                    end = start + int(endG)-int(startG) + 1
                                    break
                                else:
                                    distance = distance + int(ex[1]) - int(startG) + 1
                            else: #negative strand
                                start = int(startG)-int(start) + 1 + distance #diff between genomic pos and CDSstart
                                bstart = True
                                if int(endG)>=int(ex[0]) and int(endG)<=int(ex[1]):
                                    end = start + int(startG)-int(endG) + 1
                                    break
                                else:
                                    distance = distance + int(startG) - int(ex[0])+ 1
                        else:
                            distance = int(ex[1])-int(ex[0]) + 1 
                        if bstart and CDSstart:
                            if strand == "+":
                                if int(endG)>=int(ex[0]) and int(endG)<=int(ex[1]):
                                    end = int(endG)-int(ex[0]) + 1 + distance
                                    break
                                else:
                                    distance = distance + (int(ex[1]) - int(ex[0]) + 1)
                            else: #negative strand
                                if int(endG)>=int(ex[0]) and int(endG)<=int(ex[1]):
                                    end = int(ex[1]) - int(endG) + 1 + distance
                                    break
                                else:
                                    distance = distance + (int(ex[1]) - int(ex[0]) + 1)
                    else:
                        start = startG
                        end = endG
                newline = fields[0] + "\t" + fields[1] + "\t" + fields[2] + "\t" + str(start) + "\t" + str(end) + "\t" + fields[5] + "\t" + fields[6] + "\t" + fields[7] + "\t" + fields[8]
                if not dc_newGFF3.get(trans):
                    dc_newGFF3.update({str(trans) : [newline]})
                else:
                    dc_newGFF3.update({str(trans) : dc_newGFF3.get(trans) + [newline]})
            else:
                if not dc_newGFF3.get(trans):
                    dc_newGFF3.update({str(trans) : [line]})
                else:
                    dc_newGFF3.update({str(trans) : dc_newGFF3.get(trans) + [line]})
            
            if trans == transInterest and fields[2]==featureOfInterest and verbose:
                print("****newLine****")
                print(new_line)

    return dc_newGFF3

def transformCDStoGenomic(dc_SQcoding, dc_SQexons, dc_SQstrand):
    #toTest
    transInterest = "PB.14090.2" #sometimes we can found problems with CDS positions in classification file
    
    global verbose
    newdc_coding = {}
    bnegative = False

    for trans in dc_SQcoding.keys(): #for each coding trans...
        newCDS = []
        aux = []
        CDS = dc_SQcoding.get(trans)
        
        if CDS[0] == "NA": #if NA we have to added as well... But we can not convert the position to genomic... but we do not want to lose the transcripts. Then go for next transcript.
            if not newdc_coding.get(str(trans)):
                newdc_coding.update({str(trans) : [CDS]})
            else:                
                newdc_coding.update({str(trans) : newdc_coding.get(str(trans)) + [CDS]})
            continue

        totalDiff = int(CDS[1]) - int(CDS[0]) # number of aminoacids between CDS region

        #Exons manipulation
        allExons = dc_SQexons.get(trans)
        if not allExons:
            continue

        #Sort it depending of strand
        if dc_SQstrand.get(trans) == "+":
            allExons = sorted(allExons)
        else:
            allExons = sorted(allExons, reverse = True)

        bstart = False
        bend = False
        start = 0
        end = 0
        otherExonsDistance = 0

        if(trans==transInterest and verbose):
            print(trans)
            print(allExons)
            print(totalDiff)
            print(CDS)

        for exon in allExons: #for each exon in all the transcript...
            if totalDiff < 0: #[3]
                print("The difference can't be negative.")
                bnegative = True
                break

            #END
            if bstart:
                if dc_SQstrand.get(trans) == "+":
                    if exon[0]+totalDiff-1 <= exon[1]: #CDS ends here
                        end = exon[0]+totalDiff-1
                        aux = [[exon[0],end]]
                        newCDS = newCDS + aux
                        bend = True
                    else: #CDS ends in other exon and we add the final exon
                        aux = [[exon[0], exon[1]]]
                        newCDS = newCDS + aux
                        totalDiff = totalDiff - (exon[1]-exon[0]+1)
                else:
                    if exon[1]-totalDiff+1 >= exon[0]: #CDS ends here
                        end = exon[1]-totalDiff+1
                        aux = [[end,exon[1]]]
                        newCDS = newCDS + aux
                        bend = True
                    else: #CDS ends in other exon and we add the final exon
                        aux = [[exon[0], exon[1]]]
                        newCDS = newCDS + aux
                        totalDiff = totalDiff - (exon[1]-exon[0]+1)

            #START
            if exon[1]-exon[0]+1 >= (int(CDS[0]) - otherExonsDistance) and not bstart: #CDS starts here : We are looking for the first position on the exons, CDS is in local position and we have to look for the exon with the genome position number CDS[0]
                if dc_SQstrand.get(trans) == "+":
                    start = exon[0]+(int(CDS[0]) - otherExonsDistance) #genomic position to start isexonIni + (CDS - 1)
                    bstart = True
                    if start+totalDiff-1 <= exon[1]: #CDS ends here
                        end = start+totalDiff-1
                        aux = [[start,end]]
                        newCDS = newCDS + aux
                        bend = True
                    else: #CDS ends in other exon and we add the final exon
                        aux = [[start, exon[1]]]
                        newCDS = newCDS + aux
                        totalDiff = totalDiff - (exon[1]-start+1) #we have to update the length of CDS we already use it - we do not have to add one, bc then later we can add without adding extra numbers
                else:
                    start = exon[1]-(int(CDS[0]) - otherExonsDistance) #end genomic position minus start isexonIni + (CDS - 1)
                    bstart = True
                    if start-totalDiff+1 >= exon[0]: #CDS ends here
                        end = start-totalDiff+1
                        aux = [[start,end]]
                        newCDS = newCDS + aux
                        bend = True
                    else: #CDS ends in other exon and we add the final exon
                        aux = [[exon[0], start]]
                        newCDS = newCDS + aux
                        totalDiff = totalDiff - (start-exon[0]+1) #we have to update the length of CDS we already use it - we do not have to add one, bc then later we can add without adding extra numbers

            elif not bstart: #CDS start in other exon...
                otherExonsDistance = otherExonsDistance + (exon[1]-exon[0]+1) # we have to substract each exon we look for bc is too short to start here

            if bend:
                if not newdc_coding.get(str(trans)):
                    newdc_coding.update({str(trans) : newCDS})
                else:                
                    newdc_coding.update({str(trans) : newdc_coding.get(str(trans)) + newCDS})

                if(trans==transInterest and verbose):
                    print(newCDS)

                break

        if bnegative:
            break

    return newdc_coding

def checkSameCDS(dc_SQcoding, dc_GFF3coding, transSQ, transGFF3, strand):
    global verbose
    coding = True
    semicoding = True
    total_semi = 0
    total_annot = 0

    if dc_SQcoding.get(transSQ) and dc_GFF3coding.get(transGFF3): #both coding
        if not dc_SQcoding.get(transSQ)[0][0]=="NA":
            #We already have intervale exons range
            #   If all of them match, is coding
            #   If all of them match except sub-exons (begining or end) is semi-coding
            allExonsGFF3 = dc_GFF3coding.get(transGFF3)
            
            if strand == "+":
                allExonsGFF3 = sorted(allExonsGFF3)
            else:
                allExonsGFF3 = sorted(allExonsGFF3, reverse = True)

            for ex in allExonsGFF3:
                allExonsSQ = dc_SQcoding.get(transSQ)
                if strand == "+":
                    allExonsSQ = sorted(allExonsSQ)
                else:
                    allExonsSQ = sorted(allExonsSQ, reverse = True)

                if ex in allExonsSQ:
                    total_annot = total_annot + 1
                    continue
                else:
                    coding = False
                    semicoding = False #Check if we found semicoding
                    for exSQ in allExonsSQ:
                        if ex[0] <= exSQ[0] and exSQ[1] <= ex[1]: #Region inside
                            total_semi = total_semi + 1
                            semicoding = True
                            break
                        elif exSQ[0] <= ex[0] and exSQ[1] <= ex[1]: #or region bigger by left
                            total_semi = total_semi + 1
                            semicoding = True
                            break
                        elif ex[0] <= exSQ[0] and ex[1] <= exSQ[1]: #or region bigger by right
                            total_semi = total_semi + 1
                            semicoding = True
                            break
                        elif exSQ[0] <= ex[0] and ex[1] <= exSQ[1]: #or region bigger by both sides
                            total_semi = total_semi + 1
                            semicoding = True
                            break

        if total_annot == len(dc_GFF3coding.get(transGFF3)) and not dc_GFF3coding.get(transGFF3):
            coding = True
        elif total_annot > 0 or total_semi > 0:
            semicoding = True
        else:
            coding = False
            semicoding = False
    else: #if any not coding false
        coding = False
        semicoding = False
    return coding, semicoding

def checkFeatureInCDS(dc_SQcoding, dc_GFF3coding, transSQ, transGFF3, start, end, strand):
    #toTest
    transInterest = "PB.14090.2" #sometimes we can found problems with CDS positions in classification file
    startInterest = 166180380
    
    global verbose
    bstart = False
    if not dc_SQcoding.get(transSQ)[0]=="NA" and dc_SQcoding.get(transSQ) and dc_GFF3coding.get(transGFF3):
        #We already have intervale exons range
        #   If all of them match, is coding
        #   If all of them match except sub-exons (begining or end) is semi-coding
        allExonsGFF3 = dc_GFF3coding.get(transGFF3)
        allExonsSQ = dc_SQcoding.get(transSQ)

        if strand == "+":
            allExonsGFF3 = sorted(allExonsGFF3)
            allExonsSQ = sorted(allExonsSQ)
        else:
            allExonsGFF3 = sorted(allExonsGFF3, reverse = True)
            allExonsSQ = sorted(allExonsSQ, reverse = True)
        
        if strand == "-":
            aux = end
            end = start
            start = aux

        # if transSQ == transInterest and verbose:
        #     print("Trans: " + transSQ)
        #     print("All features positions")
        #     print(str(start), str(end))

        if transSQ == transInterest and start == startInterest and verbose:
            print("CheckFeatureCDS")
            print("Exons GFF3")
            print(allExonsGFF3)
            print("Exons SQ")
            print(allExonsSQ)
            print(str(start), str(end))

        for ex in allExonsGFF3:
            #########
            #  END  #
            #########
            if bstart:
                if strand=="+":
                    if ex[0] <= end and end <= ex[1]:#end in exon
                        if ex in allExonsSQ: #if exon exist
                            return True
                        else:
                            for exSQ in allExonsSQ: #we just need end subexon
                                if exSQ[0] == ex[0] and end <= exSQ[1] and end >= exSQ[0]: #and feature in range
                                    if transSQ == transInterest and start == startInterest and verbose:
                                        print("Found")
                                    return True
                            return False #does not find the feture in same exon
                    elif end > ex[1]: #we have to check the next exon in list if this exon is equal to one in SQlist
                        if not ex in allExonsSQ:
                            return False #end in another exons and we don't have that intermediate in SQ
                        else:
                            if transSQ == transInterest and start == startInterest and verbose:
                                print("Después de encontrar el start, seguimos buscando el exon correcto")
                            continue
                else:
                    if ex[0] <= end and end <= ex[1]:#end in exon
                        if ex in allExonsSQ: #if exon exist
                            return True
                        else:
                            for exSQ in allExonsSQ: #we just need end subexon
                                if exSQ[1] == ex[1] and end <= exSQ[1] and end >= exSQ[0]: #and feature in range
                                    if transSQ == transInterest and start == startInterest and verbose:
                                        print("Encontrado")
                                    return True
                            return False #does not find the feture in same exon
                    elif end < ex[0]: #we have to check the next exon in list if this exon is equal to one in SQlist
                        if not ex in allExonsSQ:
                            return False #end in another exons and we don't have that intermediate in SQ
                        else:
                            if transSQ == transInterest and start == startInterest and verbose:
                                print("Después de encontrar el start, seguimos buscando el exon correcto")
                            continue

            #########
            # START #
            #########
            if ex[0]<= start and start <= ex[1] and not bstart: #start in exon
                if strand=="+":
                    if transSQ == transInterest and start == startInterest and verbose:
                        print("Start found:")
                        print(start, ex)
                        print("Where is the end? " + str(end))
                    if ex[0] <= end and end <= ex[1]: #end in exon
                        if ex in allExonsSQ: #if exon exist
                            return True
                        else:
                            for exSQ in allExonsSQ: #we just need start and end in subexon on SQANTI
                                if exSQ[0]<= start and start <= exSQ[1] and exSQ[0]<= end and end <= exSQ[1]: #start and end in one exonSQ
                                    return True
                            return False #does not find the feture in same exon

                    else: #we need an exSQ that ends in same position to continue
                        for exSQ in allExonsSQ:
                            if exSQ[0]<= start and start <= exSQ[1] and ex[1]==exSQ[1]: #it can start before but end in same place
                                if transSQ == transInterest and start == startInterest and verbose:
                                    print("Continue looking in next exon...")
                                bstart = True
                else: #negative strand
                    if transSQ == transInterest and start == startInterest and verbose:
                        print("Start found:")
                        print(start, ex)
                        print("Where is the end? " + str(end))
                    if ex[0] <= end and end <= ex[1]: #end in exon
                        if ex in allExonsSQ: #if exon exist
                            return True
                        else:
                            for exSQ in allExonsSQ: #we just need start and end in subexon on SQANTI
                                if exSQ[0]<= start and start <= exSQ[1] and exSQ[0]<= end and end <= exSQ[1]: #start and end in one exonSQ
                                    return True
                            return False #does not find the feture in same exon

                    else: #we need an exSQ that ends in same position to continue
                        for exSQ in allExonsSQ:
                            if exSQ[0]<= start and start <= exSQ[1] and ex[0]==exSQ[0]: #it can start before but end in same place
                                if transSQ == transInterest and start == startInterest and verbose:
                                    print("Continue looking in next exon...")
                                bstart = True
    return False

def transformGenomicToLocale(dc_SQcoding, transSQ, start, end, strand, PROT):
    transInterest = "NM_001081080.1"
    global verbose
    bstart = False
    if not dc_SQcoding.get(transSQ)[0]=="NA" and dc_SQcoding.get(transSQ):
        # We already have intervale exons range
        #   If all of them match, is coding
        #   If all of them match except sub-exons (begining or end) is semi-coding
        allExonsSQ = dc_SQcoding.get(transSQ)

        if strand == "+":
            allExonsSQ = sorted(allExonsSQ)
        else:
            allExonsSQ = sorted(allExonsSQ, reverse = True)
        
        if strand == "-":
            aux = end
            end = start
            start = aux

        localStart = 0
        localEnd = 0
        acumulated = 0

        if transSQ==transInterest and verbose:
            print(transSQ)
            print(PROT)
            print(start)
            print(end)
            print(allExonsSQ)

        for ex in allExonsSQ:
            if(ex[0]=="NA"):
                break #not coding
            #########
            #  END  #
            #########
            if bstart:
                if strand=="+":
                    if ex[0] <= end and end <= ex[1]:#end in exon
                        #localStart
                        localEnd = end-ex[0]+1 + acumulated + localStart
                        break
                        #return(localStart, localEnd)

                    elif end > ex[1]: #we have to check the next exon in list if this exon is equal to one in SQlist
                        bstart = True
                        acumulated = acumulated + ex[1]-ex[0]+1
                else:
                    if ex[0] <= end and end <= ex[1]: #end in exon
                        #localStart
                        localEnd = ex[1]-end+1 + acumulated + localStart
                        break
                        #return(localStart, localEnd)

                    else: #we need an exSQ that ends in same position to continue
                        bstart = True
                        acumulated = acumulated + ex[1]-ex[0]+1

            #########
            # START #
            #########
            if ex[0]<= start and start <= ex[1] and not bstart: #start in exon
                if strand=="+":
                    if ex[0] <= end and end <= ex[1]: #end in exon
                        localStart = start-ex[0]+1+acumulated
                        localEnd = end-start+1 + localStart
                        break
                        #return(int(localStart/3), int(localEnd/3))

                    else: #get acumulative sum and look for the next one
                        bstart = True
                        localStart = start-ex[0]+1+acumulated
                        acumulated = ex[1]-start+1 #now for the end
                else: #negative strand
                    if ex[0] <= end and end <= ex[1]: #end in exon
                        localStart = ex[1]-start+1+acumulated
                        localEnd = start-end + localStart
                        break
                        #return(localStart, localEnd)

                    else: #we need an exSQ that ends in same position to continue
                        bstart = True
                        localStart = ex[1]-start+1+acumulated
                        acumulated = start-ex[0]+1 #now for the end
            elif not bstart: #not found in 1st exn, get acumulated
                    acumulated = acumulated + ex[1]-ex[0]+1

        if transSQ==transInterest and verbose:
            print(transSQ)
            print(PROT)
            print(localStart)
            print(localEnd)

        if(PROT):
            return(int(localStart/3), int(localEnd/3)) #to transform into aminoacids
        else:
            return(int(localStart), int(localEnd))

def checkFeatureInTranscript(dc_SQexons, dc_GFF3transExons, transSQ, transGFF3, start, end, strand, dc_SQcoding, dc_GFF3coding):
    #toTest
    transInterest = "ENST00000486627" #sometimes we can found problems with CDS positions in classification file
    startInterest = 0
    global verbose
    bstart = False
    bnotMiddleExon = False
    res = False
    #check if the feature is inside both CDS regions or outside for both CDS regions #si Codificante entonces estará, si no no
    codingSQ = False
    codingGFF3 = False

    if(transSQ==transInterest and verbose):
        print("****START*****")
        print(transSQ, transInterest)
        print(start, end)

    #NON-CODING CASE
    if not dc_SQcoding.get(transSQ) and not dc_GFF3coding.get(transGFF3):
        codingSQ = False
        codingGFF3 = False
    #CODING CASE
    if dc_SQcoding.get(transSQ) and dc_GFF3coding.get(transGFF3):
        CDSSQ =  dc_SQcoding.get(transSQ) #list of CDS exons
        CDSGFF3 = dc_GFF3coding.get(transGFF3) #list of CDS exons
        
        CDSSQ_start = 100000000000
        CDSSQ_end = 0
        for ex in CDSSQ:
            aux_s = min(ex)
            aux_m = max(ex)
            if aux_s == "NA" or aux_m == "NA":
                continue 
            if aux_s < CDSSQ_start:
                 CDSSQ_start = aux_s
            if aux_m > CDSSQ_end:
                CDSSQ_end = aux_m
        
        CDSGFF3_start = 100000000000
        CDSGFF3_end = 0
        for ex in CDSGFF3:
            aux_s = min(ex)
            aux_m = max(ex)
            if aux_s < CDSGFF3_start:
                 CDSGFF3_start = aux_s #min
            if aux_m > CDSGFF3_end:
                CDSGFF3_end = aux_m #max

        if CDSSQ_start <= start <= CDSSQ_end and CDSSQ_start <= end <= CDSSQ_end:
            codingSQ = True
        if CDSGFF3_start <= start <= CDSGFF3_end and CDSGFF3_start <= end <= CDSGFF3_end:
            codingGFF3 = True

    #SEMI CODING CASE A:
    if not dc_SQcoding.get(transSQ) and dc_GFF3coding.get(transGFF3):
        CDSGFF3 = dc_GFF3coding.get(transGFF3) #list of CDS exons
        CDSGFF3_start = 100000000000
        CDSGFF3_end = 0
        for ex in CDSGFF3:
            aux_s = min(ex)
            aux_m = max(ex)
            if aux_s == "NA" or aux_m == "NA":
                continue 
            if aux_s < CDSGFF3_start:
                 CDSGFF3_start = aux_s #min
            if aux_m > CDSGFF3_end:
                CDSGFF3_end = aux_m #max

        if CDSGFF3_start <= start <= CDSGFF3_end and CDSGFF3_start <= end <= CDSGFF3_end:
            codingGFF3 = True

    #SEMI CODING CASE B:
    if dc_SQcoding.get(transSQ) and not dc_GFF3coding.get(transGFF3):
        CDSSQ =  dc_SQcoding.get(transSQ) #list of CDS exons
        CDSSQ_start = 100000000000
        CDSSQ_end = 0
        for ex in CDSSQ:
            aux_s = min(ex)
            aux_m = max(ex)
            if aux_s == "NA" or aux_m == "NA":
                continue 
            if aux_s < CDSSQ_start:
                 CDSSQ_start = aux_s
            if aux_m > CDSSQ_end:
                CDSSQ_end = aux_m

        if CDSSQ_start <= start <= CDSSQ_end and CDSSQ_start <= end <= CDSSQ_end:
            codingSQ = True

    if(transSQ==transInterest and verbose):
        print("CODING")
        print(codingSQ, codingGFF3)

    if not codingSQ == codingGFF3:
        if(transSQ==transInterest and verbose):
            print("Entro en coding ==??")
        return res

    if dc_SQexons.get(transSQ) and dc_GFF3transExons.get(transGFF3):
        #We already have intervale exons range
        #   If all of them match, is coding
        #   If all of them match except sub-exons (begining or end) is semi-coding
        allExonsGFF3 = dc_GFF3transExons.get(transGFF3)
        allExonsSQ = dc_SQexons.get(transSQ)
        if strand == "+":
            allExonsGFF3 = sorted(allExonsGFF3)
            allExonsSQ = sorted(allExonsSQ)
        else:
            allExonsGFF3 = sorted(allExonsGFF3, reverse = True)
            allExonsSQ = sorted(allExonsSQ, reverse = True)
            a = start
            start = end
            end = a

        #we already have genomic position, just check is inside SQ trasncript exons
        for ex in allExonsGFF3:
            ##1
            if(transSQ==transInterest and start==startInterest and verbose):
                print(start, end)
                print(ex)
                print(ex[0] <= start and start <= ex[1])
                print("****")
            if ex[0] <= start and start <= ex[1] and not bstart: #Annot in exon
                if(transSQ==transInterest and start==startInterest and verbose):
                    print("strand" + strand)
                if strand == "+":
                    exSQ = [0,0]
                    for ex_aux_sq in allExonsSQ: #Look for start in SQ
                        if(transSQ==transInterest and start==startInterest and verbose):
                            print(ex_aux_sq)
                        if ex_aux_sq[0] <= start and start <= ex_aux_sq[1]:
                            exSQ = ex_aux_sq
                            break

                    if(transSQ==transInterest and start==startInterest and verbose):
                        print(exSQ)

                    if exSQ[0] <= end and end <= exSQ[1]: #also end it's here
                        res = True
                        break
                    elif ex[1] == exSQ[1]: #end in another exon, we need same ending
                        if(transSQ==transInterest and start==startInterest and verbose):
                            print("In another exon")
                        bstart = True
                    else: #different termination
                        break
                else:
                    exSQ = [0,0]
                    for ex_aux_sq in allExonsSQ: #Look for start in SQ
                        if ex_aux_sq[0] <= start and start <= ex_aux_sq[1]:
                            exSQ = ex_aux_sq
                            break
                    if(transSQ==transInterest and start==startInterest and verbose):
                        print("strand negativo")
                        print("ini" + str(start) + " fin: " + str(end))
                        print(exSQ, ex)
                        print(ex[0] <= start and start <= ex[1])
                        print("_______")
                    if exSQ[0] <= end and end <= exSQ[1]: #also end it's here
                        res = True
                        break
                    elif ex[0] == exSQ[0]: #end in another exon, we need same ending
                        if(transSQ==transInterest and start==startInterest and verbose):
                            print("In another exon")
                        bstart = True
                    else: #different termination
                        break
            ##2 final exon
            elif bstart and ex[0] <= end and end <= ex[1]: #End Annot in exon
                if strand == "+":
                    exSQ = [0,0]
                    for ex_aux_sq in allExonsSQ: #Look for start in SQ
                        if ex_aux_sq[0] <= end and end <= ex_aux_sq[1]:
                            exSQ = ex_aux_sq
                            break
                    if(transSQ==transInterest and start==startInterest and verbose):
                        print("End is here")
                        print(ex)
                        print(exSQ)
                        print("end: " + str(end))
                        
                    if exSQ[0] <= end and end <= exSQ[1] and exSQ[0]==ex[0]: #also end it's here
                        res = True
                        break
                    else: #SQexon shorter
                        break

                else:
                    exSQ = [0,0]
                    for ex_aux_sq in allExonsSQ: #Look for start in SQ
                        if ex_aux_sq[0] <= end and end <= ex_aux_sq[1]:
                            exSQ = ex_aux_sq
                            break
                    if(transSQ==transInterest and start==startInterest and verbose):
                        print("End is here")
                        print(ex)
                        print(exSQ)
                        print("end: " + str(end))
                        
                    if exSQ[0] <= end and end <= exSQ[1] and exSQ[1]==ex[1]: #also end it's here
                        res = True
                        break
                    else: #SQexon shorter
                        break

            elif(bstart): #we jump a complete exon (that is also present in SQtrans)
                found = False
                for ex_aux_sq in allExonsSQ:
                    if(ex==ex_aux_sq):
                        found = True
                if(found):
                    if(transSQ==transInterest and start==startInterest and verbose):
                        print("saltamos un exon...")
                    continue
                else:
                    break
            else:
                continue #if start is not found, we keep looking
            
            if bnotMiddleExon:
                break

    if(transSQ==transInterest and start==startInterest and verbose):
        print("*********")
        print(start, end)
        print(allExonsGFF3)
        print(res)
    return res

def getTranscriptomePerGene(dc_SQgeneTrans, dc_GFF3_Genomic):
    global verbose
    dc_GFF3Gene_Genomic = {}
    for gene in dc_SQgeneTrans.keys():
        #get trans
        lst_trans = dc_SQgeneTrans.get(gene)
        for trans in lst_trans:
            if dc_GFF3_Genomic.get(trans):
                lst_lines = dc_GFF3_Genomic.get(trans) # trans:[[s,e,info], [s,e,info]] - we want only the lines, no list of lines
                if(not dc_GFF3Gene_Genomic.get(gene)):
                    dc_GFF3Gene_Genomic.update({str(gene) : lst_lines})
                else:
                    dc_GFF3Gene_Genomic.update({str(gene) : dc_GFF3Gene_Genomic.get(gene) + lst_lines})
    
    return dc_GFF3Gene_Genomic

def range_subset(range1, range2):
    """Whether range1 is a subset of range2."""
    if not range1:
        return True  # empty range is subset of anything
    if not range2:
        return False  # non-empty range can't be subset of empty range
    if len(range1) > 1 and range1.step % range2.step:
        return False  # must have a single value or integer multiple step
    return range1.start in range2 and range1[-1] in range2

def addTranscriptWithAnnot(feature, transSQ, AnnotTranscriptCounts_dc, TotalTranscripts = False):
    #check if trans already in transAnnot - count transcript once, if it is a new trans, then added
    hs_check = AnnotTranscriptCounts_dc.get(feature)[0] #annotated - hashmaps
    hs_checkTotal = AnnotTranscriptCounts_dc.get(feature)[1] #total - hashmaps
    
    if not TotalTranscripts:
        # print("Add to the left: " + transSQ + "con id" + feature)
        # print(hs_check)
        if hs_check=={}:
            hs_check = {}
            hs_check.update({transSQ : ""})
            AnnotTranscriptCounts_dc.update({feature : [hs_check, hs_checkTotal]})
        elif not hs_check.get(transSQ):
            hs_check.update({transSQ : ""})
            AnnotTranscriptCounts_dc.update({feature : [hs_check, hs_checkTotal]})
        # print("*****")
        # print(hs_check)
    else:
        # print("Add to the right: " + transSQ + "con id" + feature)
        # print(hs_checkTotal)
        if hs_checkTotal=={}:
            hs_checkTotal = {}
            hs_checkTotal.update({transSQ : ""})
            AnnotTranscriptCounts_dc.update({feature : [hs_check, hs_checkTotal]})
        elif not hs_checkTotal.get(transSQ):
            hs_checkTotal.update({transSQ : ""})
            AnnotTranscriptCounts_dc.update({feature : [hs_check, hs_checkTotal]})
    #     print("*****")
    #     print(hs_checkTotal)
    # print("_________")
    return AnnotTranscriptCounts_dc

#same for Features 
def addFeatureWithAnnot(feature, key, AnnotFeatureCounts_dc, TotalFeatures = False):
    #check if trans already in transAnnot - count transcript once, if it is a new trans, then added
    hs_check = AnnotFeatureCounts_dc.get(feature)[0] #annotated - hashmaps
    hs_checkTotal = AnnotFeatureCounts_dc.get(feature)[1] #total - hashmaps
    
    if not TotalFeatures:
        if hs_check=={}:
            hs_check = {}
            hs_check.update({key : ""})
            AnnotFeatureCounts_dc.update({feature : [hs_check, hs_checkTotal]})
        elif not hs_check.get(key):
            hs_check.update({key : ""})
            AnnotFeatureCounts_dc.update({feature : [hs_check, hs_checkTotal]})
    else:
        if hs_checkTotal=={}:
            hs_checkTotal = {}
            hs_checkTotal.update({key : ""})
            AnnotFeatureCounts_dc.update({feature : [hs_check, hs_checkTotal]})
        elif not hs_checkTotal.get(key):
            hs_checkTotal.update({key : ""})
            AnnotFeatureCounts_dc.update({feature : [hs_check, hs_checkTotal]})
    return AnnotFeatureCounts_dc

def mappingFeatures(dc_SQexons, dc_SQcoding, dc_SQtransGene, dc_SQgeneTrans, dc_SQstrand, dc_GFF3exonsTrans, dc_GFF3transExons, dc_GFF3, dc_GFF3gene, dc_GFF3coding, dc_GFF3geneTrans, dc_geneID2geneName, dc_geneName2geneID, dc_GFF3transGene, filename, filename_prints):
    #toTest
    transInterest = "ENST00000486627"
    featureOfInterest = "miRNA_Binding"
    global verbose
    global USE_STDOUT
    global LST_TRANSCRIPTFEATURES_NOTIN_CDS
    global LST_TRANSCRIPTSOURCES_INTRONIC
    global ALL_AS_NOVELS
    global INTRONIC
    global STATS
    global SAVE_PROB_TRANSCRIPTS

    if SAVE_PROB_TRANSCRIPTS:
        file_trans_not_annot_by_PF = open("file_trans_not_annot_by_PF.txt","w+")
        file_novel_not_annot_by_PF = open("file_novel_not_annot_by_by_PF.txt","w+")
        file_reference_gene_not_annot = open("file_reference_gene_not_annot.txt","w+")
        file_reference_transcript_not_annot = open("file_reference_transcript_not_annot.txt","w+")
        file_trans_not_gene_ID = open("file_transcript_wo_gene_ID.txt","w+")

    f = open(filename,"a+")
    print("\n")

    realNewTrans = False

    novelTranscripts = 0
    novelTranscriptsAnnotated = 0
    novelTranscriptsAnnotated_wo_features = 0
    novelTranscriptsNotAnnotated = 0
    novelTranscriptsNotGeneMatch = 0
    
    transcriptsAnnotated = 0
    transcriptsAnnotated_wo_features = 0
    transcriptsNotAnnotated = 0
    transcriptsNotAnnotatedNotGeneInformation = 0
    transcriptsNotGeneID = 0

    featuresAnnotated = 0

    anotationsChecked = 0
    transcriptsChecked = 0
    perct = 0

    #countDomain = 0

    transAnnotTotal_dc = {}
    protAnnotTotal_dc = {}
    geneAnnotTotal_dc = {}

    transAnnotTranscriptCounts_dc = {}
    protAnnotTranscriptCounts_dc = {}
    geneAnnotTranscriptCounts_dc = {}

    annotLines = [] #for novel transcripts
    
    for transSQ in dc_SQexons.keys():
        transcriptsChecked = transcriptsChecked + 1
        realNewTrans = False
        if ALL_AS_NOVELS:
            newTrans = True
        else:
            newTrans = False

        anyAnnot = False
        anyNovelAnnotation = False
        isNovel = False
        referenceTransNoAnnot = True
        referenceGeneNoAnnot = True

        #Be carefully - not all tranSQ must be in SQtransGene [we do not have gene and transcript related information]
        if not dc_SQtransGene.get(str(transSQ)):
            if SAVE_PROB_TRANSCRIPTS:
                file_trans_not_gene_ID.write(transSQ + "\n")
            transcriptsNotGeneID = transcriptsNotGeneID + 1
            continue
        
        perct = transcriptsChecked/len(dc_SQexons)*100
        m0 = "\n\n\t" + "%.2f" % perct + " % of transcripts annotated."
        print("\t" + "%.2f" % perct + " % of GTF transcripts checked...", end="\r")

        #######################
        #IF FULL-SPLICED-MATCH#
        #######################
        infoGenomic = dc_SQtransGene.get(transSQ)
        transGFF3 = infoGenomic[2] #transAssociated
        geneSQ = infoGenomic[0] #gene

        ###########################
        #IF NOT FULL-SPLICED-MATCH#
        ###########################
        val = ""

        if dc_GFF3.get(transSQ):
            transGFF3 = transSQ
            val = dc_GFF3.get(transGFF3)
        elif dc_GFF3.get(transGFF3): #dc_GFF3 key=gene
            val = dc_GFF3.get(transGFF3)
        elif(dc_GFF3transExons.get(transGFF3)): #if transcript have exons, is not novel, it means the reference have no annotations but can be annotated as a novel
            newTrans = True
        else: #Novel Transcript won't be annoted
            newTrans = True

        if(dc_GFF3transGene.get(transSQ) or dc_GFF3transGene.get(transGFF3)):
            realNewTrans = False #trans or reference trans exits in GFF3 (is not a novel transcript)
        elif(dc_GFF3geneTrans.get(geneSQ)):
            realNewTrans = True #gene exits in GFF3 (with or without annot)
        else: #not transcript and not gene exists in GFF3, novel gene?
            if SAVE_PROB_TRANSCRIPTS:
                file_reference_gene_not_annot.write(transSQ + "\n")
            novelTranscripts = novelTranscripts + 1
            novelTranscriptsNotGeneMatch = novelTranscriptsNotGeneMatch + 1
            continue

        if(transSQ==transInterest and verbose):
            print("First print:")
            print(val)
            print(newTrans)
            print(referenceTransNoAnnot)

        if(newTrans==False):
            strand = dc_SQstrand.get(transSQ)
            #Check if we had same CDS to add Protein information
            coding, semicoding = checkSameCDS(dc_SQcoding, dc_GFF3coding, transSQ, transGFF3, strand)

            for values in dc_GFF3.get(transGFF3):
                fields = values[2].split("\t")
                text = fields[8].split(" ")
                if fields[1] == "tappAS":
                    continue
                anotationsChecked = anotationsChecked + 1

                #################
                #GENE ANNOTATION#
                #################
                if (text[-1].endswith("G\n") or text[-1].endswith("N\n")): #gene - always put if the gene is the same
                    
                    referenceTransNoAnnot = False #exists at least one feature

                    #update total prot annotation
                    key = fields[2]+fields[3]+str(values[0])+str(values[1])+transSQ #miRNAmi-hsa-3p345234PB.2.2

                    if(not geneAnnotTotal_dc.get(fields[2])):
                        countFeatureAnnot_dc = {}
                        countFeatureTotalAnnot_dc = {}
                        countFeatureTotalAnnot_dc[key] = ""
                        geneAnnotTotal_dc.update({fields[2] : [countFeatureAnnot_dc, countFeatureTotalAnnot_dc]})

                        countTranscriptsAnnot_dc = {}
                        countTranscriptsTotalAnnot_dc = {}
                        countTranscriptsTotalAnnot_dc[transSQ] = ""
                        geneAnnotTranscriptCounts_dc.update({fields[2] : [countTranscriptsAnnot_dc, countTranscriptsTotalAnnot_dc]})
                    else:
                        geneAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, geneAnnotTotal_dc, True)
                        geneAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, geneAnnotTranscriptCounts_dc, True)

                    #check if trans already in transAnnot - count transcript once, if it is a new trans, then added
                    geneAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, geneAnnotTranscriptCounts_dc, False)
                    geneAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, geneAnnotTotal_dc, False) #sum 1

                    featuresAnnotated = featuresAnnotated + 1
                    anyAnnot = True
                    
                    index = values[2].find("\t")
                    line_text = values[2].split("\t")
                    if values[2].endswith("\n"):
                        new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:]))])) #write line
                    else:
                        new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:])), "\n"])) #write line
                    annotLines.append(new_line)

                ####################
                #PROTEIN ANNOTATION#
                ####################
                elif (text[-1].endswith("P\n")): #protein

                    referenceTransNoAnnot = False #exists at least one feature

                    if(transGFF3==transInterest and fields[2] == featureOfInterest and verbose):
                        print("\n\nTrans:" + transSQ + " con feature: " + featureOfInterest)
                        print("start: " + str(values[0]) + " end: " + str(values[1]))
                        print(strand)

                    #update total prot annotation
                    key = fields[2]+fields[3]+str(values[0])+str(values[1])+transSQ #miRNAmi-hsa-3p345234PB.2.2
                    if(not protAnnotTotal_dc.get(fields[2])):
                        countFeatureAnnot_dc = {}
                        countFeatureTotalAnnot_dc = {}
                        countFeatureTotalAnnot_dc[key] = ""
                        protAnnotTotal_dc.update({fields[2] : [countFeatureAnnot_dc, countFeatureTotalAnnot_dc]})

                        countTranscriptsAnnot_dc = {}
                        countTranscriptsTotalAnnot_dc = {}
                        countTranscriptsTotalAnnot_dc[transSQ] = ""
                        protAnnotTranscriptCounts_dc.update({fields[2] : [countTranscriptsAnnot_dc, countTranscriptsTotalAnnot_dc]})
                    else:
                        protAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, protAnnotTotal_dc, True)
                        protAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, protAnnotTranscriptCounts_dc, True)

                    if coding:
                        index = values[2].find("\t")
                        line_text = values[2].split("\t")

                        #check if trans already in transAnnot - count transcript once, if it is a new trans, then added
                        protAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, protAnnotTranscriptCounts_dc, False)
                        protAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, protAnnotTotal_dc, False) #sum 1
                        
                        featuresAnnotated = featuresAnnotated + 1
                        anyAnnot = True
                        if values[2].endswith("\n"):
                            new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:]))])) #write line
                        else:
                            new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:])), "\n"])) #write line
                        annotLines.append(new_line)

                    elif semicoding and not values[0]=="." and not values[1] == ".":
                        bannot = False
                        #function: match annot to its our CDSexons and match to CDSexonsSQ
                        bannot = checkFeatureInCDS(dc_SQcoding, dc_GFF3coding, transSQ, transGFF3, int(values[0]), int(values[1]), strand)

                        if bannot:
                            index = values[2].find("\t")
                            line_text = values[2].split("\t")

                            #check if trans already in transAnnot - count transcript once, if it is a new trans, then added
                            protAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, protAnnotTranscriptCounts_dc, False)
                            protAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, protAnnotTotal_dc, False) #sum 1
                        
                            featuresAnnotated = featuresAnnotated + 1
                            anyAnnot = True
                            if values[2].endswith("\n"):
                                new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:]))])) #write line
                            else:
                                new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:])), "\n"])) #write line
                            annotLines.append(new_line)
                        else: #not bannot
                            if(fields[2] == featureOfInterest) and verbose:
                                print(featureOfInterest + " (interest):" + transSQ)

                    elif semicoding and values[0]=="." and values[1] == ".":
                        index = values[2].find("\t")
                        line_text = values[2].split("\t")

                        #check if trans already in transAnnot - count transcript once, if it is a new trans, then added
                        protAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, protAnnotTranscriptCounts_dc, False)
                        protAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, protAnnotTotal_dc, False) #sum 1
                        
                        featuresAnnotated = featuresAnnotated + 1
                        anyAnnot = True
                        if values[2].endswith("\n"):
                            new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:]))])) #write line
                        else:
                            new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:])), "\n"])) #write line
                        annotLines.append(new_line)
                    # else: #not annotated
                    #     if(fields[2] == "DOMAIN"):
                    #         countDomain = countDomain + 1
                    #         print("DOMAIN (interest):" + transSQ)
                    #         print(countDomain)

                #######################
                #TRANSCRIPT ANNOTATION#
                #######################
                elif not values[0]=="." and not values[1] == "." and text[-1].endswith("T\n"):

                    referenceTransNoAnnot = False #exists at least one feature

                    #update total trans annotation
                    key = fields[2]+fields[3]+str(values[0])+str(values[1])+transSQ #miRNAmi-hsa-3p345234PB.2.2
                    if(not transAnnotTotal_dc.get(fields[2])):
                        countFeatureAnnot_dc = {}
                        countFeatureTotalAnnot_dc = {}
                        countFeatureTotalAnnot_dc[key] = ""
                        transAnnotTotal_dc.update({fields[2] : [countFeatureAnnot_dc, countFeatureTotalAnnot_dc]})

                        countTranscriptsAnnot_dc = {}
                        countTranscriptsTotalAnnot_dc = {}
                        countTranscriptsTotalAnnot_dc[transSQ] = ""
                        transAnnotTranscriptCounts_dc.update({fields[2] : [countTranscriptsAnnot_dc, countTranscriptsTotalAnnot_dc]})
                    else:
                        transAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, transAnnotTotal_dc, True)
                        transAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, transAnnotTranscriptCounts_dc, True)

                    bannot = False
                    if(fields[1] in LST_TRANSCRIPTSOURCES_INTRONIC and INTRONIC):
                        exons = dc_SQexons.get(transSQ)
                        start = exons[0][0]
                        end = exons[-1][1]
                        if start <= int(values[0]) <= end and start <= int(values[1]) <= end:
                            bannot = True #check if feature inside ini-end gene position of the transcript #do not check for strand because don't care
                        else:
                            bannot = False
                    else:
                        bannot = checkFeatureInTranscript(dc_SQexons, dc_GFF3transExons, transSQ, transGFF3, int(values[0]), int(values[1]), strand, dc_SQcoding, dc_GFF3coding)

                    if transSQ == transInterest and verbose:
                        print("Transcript Feature:")
                        print(values)
                        print(dc_SQexons)
                        print(bannot)

                    if bannot:
                        index = values[2].find("\t")
                        line_text = values[2].split("\t")

                        #check if trans already in transAnnot - count transcript once, if it is a new trans, then added
                        transAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, transAnnotTranscriptCounts_dc, False)
                        transAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, transAnnotTotal_dc, False)

                        featuresAnnotated = featuresAnnotated + 1
                        anyAnnot = True
                        #sometime (PAR-CLIP) we have not strand for transcript information, in that case, we need to add it
                        if line_text[6]!=".": #normal
                            if values[2].endswith("\n"):
                                new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:]))])) #write line
                            else:
                                new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:])), "\n"])) #write line
                            annotLines.append(new_line)
                        else: #no strand
                            if values[2].endswith("\n"):
                                new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:6])), strand ,'\t'.join(map(str, line_text[7:]))])) #write line
                                #print(new_line)
                            else:
                                new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:6])), strand ,'\t'.join(map(str, line_text[7:])), "\n"])) #write line
                                #print(new_line)
                            annotLines.append(new_line)
                    else:
                        #no annotate trans
                        # if(fields[2] == "3UTRmotif"):
                        #     print("3UTRmotif (interest):" + transSQ + " using " + transGFF3)
                        # if(fields[2] == "miRNA"):
                        #     print("miRNA (interest):" + transSQ + " using " + transGFF3)
                        if(fields[2] == featureOfInterest) and verbose:
                            print(featureOfInterest + " (interest):" + transSQ + " using " + transGFF3)
                
                elif values[0]=="." and values[1] == "." and text[-1].endswith("T\n"): #NMD
                    #update total trans annotation
                    key = fields[2]+fields[3]+str(values[0])+str(values[1])+transSQ #miRNAmi-hsa-3p345234PB.2.2
                    if(not transAnnotTotal_dc.get(fields[2])):
                        countFeatureAnnot_dc = {}
                        countFeatureTotalAnnot_dc = {}
                        countFeatureTotalAnnot_dc[key] = ""
                        transAnnotTotal_dc.update({fields[2] : [countFeatureAnnot_dc, countFeatureTotalAnnot_dc]})
                        
                        countTranscriptsAnnot_dc = {}
                        countTranscriptsTotalAnnot_dc = {}
                        countTranscriptsTotalAnnot_dc[transSQ] = ""
                        transAnnotTranscriptCounts_dc.update({fields[2] : [countTranscriptsAnnot_dc, countTranscriptsTotalAnnot_dc]})
                    else:
                        transAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, transAnnotTotal_dc, True)
                        transAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, transAnnotTranscriptCounts_dc, True)


                    #check if trans already in transAnnot - count transcript once, if it is a new trans, then added
                    transAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, transAnnotTranscriptCounts_dc, False)
                    transAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, transAnnotTotal_dc, False)

                    featuresAnnotated = featuresAnnotated + 1
                    anyAnnot = True
                    
                    index = values[2].find("\t")

                    if values[2].endswith("\n"):
                        new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:]))])) #write line
                    else:
                        new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:])), "\n"])) #write line
                        annotLines.append(new_line)

                else:
                    if verbose:
                        print("We are not looking for this annotation")
                        print(transSQ + " en feature " + fields[2])
                        print(fields)

        else: #new Transcript and we have to run possible matches for all transcripts of the gene
            if(realNewTrans):
                novelTranscripts = novelTranscripts + 1 #real novelTranscript, can be here also by force novel transferred annotation

            gene = infoGenomic[0] #geneAssociated
            alternative_gene = ""
            anyNovelAnnotation = False
            isNovel = True

            l_gene = 0
            l_alternative_gene = 0

            # if verbose:
            #     print("Gene associated to map features is: " + gene)
            if len(gene) != 0:
                if dc_geneID2geneName.get(gene):
                    alternative_gene = dc_geneID2geneName.get(gene)
                elif dc_geneName2geneID.get(gene):
                    alternative_gene = dc_geneName2geneID.get(gene)
                elif alternative_gene == "" and dc_GFF3transGene.get(transGFF3): #if still NULL
                    alternative_gene = dc_GFF3transGene.get(transGFF3) #uses associated transcript in the GFF3 and get it's gene

                if verbose:
                    print("1º: " + transSQ)
                    print("Gene: " + gene)
                    print("Gene alt: " + alternative_gene)
                    print(dc_GFF3geneTrans.get(gene))
                    print(dc_GFF3geneTrans.get(alternative_gene))

                if dc_GFF3geneTrans.get(gene):
                    l_gene = len(dc_GFF3geneTrans.get(gene))
                if dc_GFF3geneTrans.get(alternative_gene):
                    l_alternative_gene = len(dc_GFF3geneTrans.get(alternative_gene))

                if l_gene < l_alternative_gene:
                    if verbose:
                        print("se cambia el gen!!\n")
                        print("El transcript de referencia original es: " + transSQ + " y el de asociación " + transGFF3)
                        print("gene: " + gene + " alternative_gene: " + alternative_gene + "\n")
                    gene = alternative_gene
                    if verbose:
                        print(dc_GFF3geneTrans.get(gene))

            #it is possible some transcripts do not have gene associated but its transcript it is annotated
            #or that associated gene is not present in the dictionaries - in that case look for the real gene can be slow if a lot of transcripts
            #we have to look for the gene id in another dictionary trans -> raw_gene 
            elif l_alternative_gene==0 and l_gene==0:
                if dc_GFF3transGene.get(transGFF3):
                    gene = dc_GFF3transGene.get(transGFF3)
                    if verbose:
                        print("2º")
                        print("Gene has not been defined, using the classification information. New gene is: " + gene)
                # for gene_test, lst_transcripts in dc_GFF3geneTrans.items():
                #     if transGFF3 in lst_transcripts:
                #         gene = gene_test
                #         break

            if verbose:
                print(transSQ)
                print(gene)
                print(dc_GFF3geneTrans.get(gene))
                print("\n")

            if dc_GFF3geneTrans.get(gene):
                
                for transGFF3 in dc_GFF3geneTrans.get(gene): #Each trans in GFF3 inside gene

                    if transSQ == transInterest and verbose:
                        print("Estamos revisando el transcrito del GFF3... " + transGFF3)

                    if dc_GFF3.get(transGFF3):
                        val = dc_GFF3.get(transGFF3)
                    else:
                        continue #if we do not have values, next transcript
                    
                    line = val[0][2].split("\t")
                    strand = dc_SQstrand.get(transSQ)
                    
                    for values in dc_GFF3.get(transGFF3):
                        fields = values[2].split("\t")
                        text = fields[8].split(" ")
                        if fields[1] == "tappAS":
                            continue
                        anotationsChecked = anotationsChecked + 1

                        #################
                        #GENE ANNOTATION#
                        #################
                        if (text[-1].endswith("G\n") or text[-1].endswith("N\n")): #gene - always put if the gene is the same

                            referenceGeneNoAnnot = False #exists at least one feature in one transcript with annot

                            #update total prot annotation
                            key = fields[2]+fields[3]+str(values[0])+str(values[1])+transSQ #miRNAmi-hsa-3p345234PB.2.2
                            if(not geneAnnotTotal_dc.get(fields[2])):
                                countFeatureAnnot_dc = {}
                                countFeatureTotalAnnot_dc = {}
                                countFeatureTotalAnnot_dc[key] = ""
                                geneAnnotTotal_dc.update({fields[2] : [countFeatureAnnot_dc, countFeatureTotalAnnot_dc]})
                        
                                countTranscriptsAnnot_dc = {}
                                countTranscriptsTotalAnnot_dc = {}
                                countTranscriptsTotalAnnot_dc[transSQ] = ""
                                geneAnnotTranscriptCounts_dc.update({fields[2] : [countTranscriptsAnnot_dc, countTranscriptsTotalAnnot_dc]})
                            else:
                                geneAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, geneAnnotTotal_dc, True)
                                geneAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, geneAnnotTranscriptCounts_dc, True)

                            #check if trans already in transAnnot - count transcript once, if it is a new trans, then added
                            geneAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, geneAnnotTranscriptCounts_dc, False)
                            geneAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, geneAnnotTotal_dc, False)

                            featuresAnnotated = featuresAnnotated + 1
                            anyNovelAnnotation  = True

                            line_text = values[2].split("\t")
                            if values[2].endswith("\n"):
                                new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:3])), values[0], values[1], '\t'.join(map(str, line_text[5:]))]))
                            else:
                                new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:3])), values[0], values[1], '\t'.join(map(str, line_text[5:])), "\n"]))
                            annotLines.append(new_line)

                        ####################
                        #PROTEIN ANNOTATION#
                        ####################
                        elif (text[-1].endswith("P\n")): #protein

                            referenceGeneNoAnnot = False #exists at least one feature in one transcript with annot

                            #Check if we had same CDS to add Protein information
                            coding, semicoding = checkSameCDS(dc_SQcoding, dc_GFF3coding, transSQ, transGFF3, strand)

                            if(coding==False and semicoding==False):
                                continue

                            #update total prot annotation
                            key = fields[2]+fields[3]+str(values[0])+str(values[1])+transSQ #miRNAmi-hsa-3p345234PB.2.2
                            if(not protAnnotTotal_dc.get(fields[2])):
                                countFeatureAnnot_dc = {}
                                countFeatureTotalAnnot_dc = {}
                                countFeatureTotalAnnot_dc[key] = ""
                                protAnnotTotal_dc.update({fields[2] : [countFeatureAnnot_dc, countFeatureTotalAnnot_dc]})
                        
                                countTranscriptsAnnot_dc = {}
                                countTranscriptsTotalAnnot_dc = {}
                                countTranscriptsTotalAnnot_dc[transSQ] = ""
                                protAnnotTranscriptCounts_dc.update({fields[2] : [countTranscriptsAnnot_dc, countTranscriptsTotalAnnot_dc]})
                            else:
                                protAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, protAnnotTotal_dc, True)
                                protAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, protAnnotTranscriptCounts_dc, True)

                            if values[0]!="." and values[1] != ".":
                                localePosition = transformGenomicToLocale(dc_SQcoding, transSQ, int(values[0]), int(values[1]), strand, True) #last argument indicate if prot or not
                            else:
                                localePosition = [".", "."]

                            if coding:
                                index = values[2].find("\t")

                                #check if trans already in transAnnot - count transcript once, if it is a new trans, then added
                                protAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, protAnnotTranscriptCounts_dc, False)
                                protAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, protAnnotTotal_dc, False) #sum 1

                                featuresAnnotated = featuresAnnotated + 1
                                anyNovelAnnotation = True
                                line_text = values[2].split("\t")
                                if values[2].endswith("\n"):
                                    # if verbose:
                                    #     print(transSQ)
                                    #     print(dc_SQcoding.get(transSQ))
                                    #     print(int(values[0]), int(values[1]))
                                    #     print(localePosition)
                                    #     print(coding, semicoding)
                                    #     print("______")

                                    new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:3])), localePosition[0], localePosition[1], '\t'.join(map(str, line_text[5:]))]))
                                else:
                                    new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:3])), localePosition[0], localePosition[1], '\t'.join(map(str, line_text[5:])), "\n"]))
                                annotLines.append(new_line)

                            elif semicoding and not values[0]=="." and not values[1] == ".":
                                bannot = False
                                #funcion match annot to its our CDSexons and match to CDSexonsSQ
                                bannot = checkFeatureInCDS(dc_SQcoding, dc_GFF3coding, transSQ, transGFF3, int(values[0]), int(values[1]), strand)

                                # if verbose:
                                #     print(transSQ)
                                #     print(dc_SQcoding.get(transSQ))
                                #     print(localePosition)
                                #     print(coding, semicoding)
                                #     print("______")

                                if bannot:
                                    index = values[2].find("\t")

                                    #check if trans already in transAnnot - count transcript once, if it is a new trans, then added
                                    protAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, protAnnotTranscriptCounts_dc, False)
                                    protAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, protAnnotTotal_dc, False) #sum 1

                                    featuresAnnotated = featuresAnnotated + 1
                                    anyNovelAnnotation = True
                                    line_text = values[2].split("\t")
                                    if values[2].endswith("\n"):
                                        new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:3])), localePosition[0], localePosition[1], '\t'.join(map(str, line_text[5:]))]))
                                    else:
                                        new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:3])), localePosition[0], localePosition[1], '\t'.join(map(str, line_text[5:])), "\n"]))
                                    annotLines.append(new_line)

                            elif semicoding and values[0]=="." and values[1] == ".":
                                index = values[2].find("\t")

                                #check if trans already in transAnnot - count transcript once, if it is a new trans, then added
                                protAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, protAnnotTranscriptCounts_dc, False)
                                protAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, protAnnotTotal_dc, False) #sum 1

                                featuresAnnotated = featuresAnnotated + 1
                                anyNovelAnnotation = True
                                line_text = values[2].split("\t")
                                if values[2].endswith("\n"):
                                    new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:3])), localePosition[0], localePosition[1], '\t'.join(map(str, line_text[5:]))]))
                                else:
                                    new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:3])), localePosition[0], localePosition[1], '\t'.join(map(str, line_text[5:])), "\n"]))
                                annotLines.append(new_line)

                        #######################
                        #TRANSCRIPT ANNOTATION#
                        #######################

                        elif not values[0]=="." and not values[1] == "." and text[-1].endswith("T\n"):

                            referenceGeneNoAnnot = False #exists at least one feature in one transcript with annot

                            #update total trans annotation
                            key = fields[2]+fields[3]+str(values[0])+str(values[1])+transSQ #miRNAmi-hsa-3p345234PB.2.2
                            if(not transAnnotTotal_dc.get(fields[2])):
                                countFeatureAnnot_dc = {}
                                countFeatureTotalAnnot_dc = {}
                                countFeatureTotalAnnot_dc[key] = ""
                                transAnnotTotal_dc.update({fields[2] : [countFeatureAnnot_dc, countFeatureTotalAnnot_dc]})
                        
                                countTranscriptsAnnot_dc = {}
                                countTranscriptsTotalAnnot_dc = {}
                                countTranscriptsTotalAnnot_dc[transSQ] = ""
                                transAnnotTranscriptCounts_dc.update({fields[2] : [countTranscriptsAnnot_dc, countTranscriptsTotalAnnot_dc]})
                            else:
                                transAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, transAnnotTotal_dc, True)
                                transAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, transAnnotTranscriptCounts_dc, True)

                            bannot = False
                            if(fields[1] in LST_TRANSCRIPTSOURCES_INTRONIC and INTRONIC):
                                exons = dc_SQexons.get(transSQ)
                                start = exons[0][0]
                                end = exons[-1][1]
                                if start <= int(values[0]) <= end and start <= int(values[1]) <= end:
                                    bannot = True #check if feature inside ini-end gene position of the transcript #do not check for strand because don't care
                                else:
                                    bannot = False
                            else:
                                bannot = checkFeatureInTranscript(dc_SQexons, dc_GFF3transExons, transSQ, transGFF3, int(values[0]), int(values[1]), strand, dc_SQcoding, dc_GFF3coding)
    
                            val = values[2].split("\t")
                            shouldAnnotate=True

                            if transSQ == transInterest and transGFF3=="NM_001081080.1"  and verbose:
                                print("Annotation for studied transcript:")
                                print(val)

                            if val[2] in LST_TRANSCRIPTFEATURES_NOTIN_CDS: #if annotation is not inside CDS, check we have that positions outside our CDS
                                if(dc_SQcoding.get(transSQ)): 
                                    cdsExons = dc_SQcoding.get(transSQ)
                                    for ex in cdsExons:
                                        if ex[0]=="NA": #if coding with values, check it 
                                            continue
                                        elif(int(ex[0]) <= int(values[0]) <= int(ex[1]) or int(ex[0]) <= int(values[1]) <= int(ex[1])): #if the annotation is inside, we cannot annoted it
                                            shouldAnnotate = False
                                            break

                            if transSQ == transInterest and transGFF3=="NM_001081080.1" and verbose:
                                print("If it is a 3UTR feature, should we try to annot?" + str(shouldAnnotate))
                                print("annotated: " + str(bannot))

                            if shouldAnnotate:
                                if bannot:
                                    index = values[2].find("\t")

                                    #check if trans already in transAnnot - count transcripts just one time, if it is new them add it
                                    transAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, transAnnotTranscriptCounts_dc, False)
                                    transAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, transAnnotTotal_dc, False)

                                    featuresAnnotated = featuresAnnotated + 1
                                    anyNovelAnnotation = True
                                    line_text = values[2].split("\t")
                                    if values[2].endswith("\n"):
                                        new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:3])), values[0], values[1], '\t'.join(map(str, line_text[5:]))]))
                                    else:
                                        new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:3])), values[0], values[1], '\t'.join(map(str, line_text[5:])), "\n"]))
                                    annotLines.append(new_line)

                            if transSQ == transInterest and transGFF3=="NM_001081080.1" and verbose:
                                print(annotLines[-1])

                        elif values[0]=="." and values[1] == "." and text[-1].endswith("T\n"): #NMD
                            #update total trans annotation
                            key = fields[2]+fields[3]+str(values[0])+str(values[1])+transSQ #miRNAmi-hsa-3p345234PB.2.2
                            if(not transAnnotTotal_dc.get(fields[2])):
                                countFeatureAnnot_dc = {}
                                countFeatureTotalAnnot_dc = {}
                                countFeatureTotalAnnot_dc[key] = ""
                                transAnnotTotal_dc.update({fields[2] : [countFeatureAnnot_dc, countFeatureTotalAnnot_dc]})
                        
                                countTranscriptsAnnot_dc = {}
                                countTranscriptsTotalAnnot_dc = {}
                                countTranscriptsTotalAnnot_dc[transSQ] = ""
                                transAnnotTranscriptCounts_dc.update({fields[2] : [countTranscriptsAnnot_dc, countTranscriptsTotalAnnot_dc]})
                            else:
                                transAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, transAnnotTotal_dc, True)
                                transAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, transAnnotTranscriptCounts_dc, True)

                            #check if trans already in transAnnot - count transcripts just one time, if it is new them add it
                            transAnnotTranscriptCounts_dc = addTranscriptWithAnnot(fields[2], transSQ, transAnnotTranscriptCounts_dc, False)
                            transAnnotTotal_dc = addFeatureWithAnnot(fields[2], key, transAnnotTotal_dc, False)

                            featuresAnnotated = featuresAnnotated + 1
                            anyNovelAnnotation = True
                            line_text = values[2].split("\t")
                            if values[2].endswith("\n"):
                                new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:3])), values[0], values[1], '\t'.join(map(str, line_text[5:]))]))
                            else:
                                new_line = '\t'.join(map(str, [transSQ, '\t'.join(map(str, line_text[1:3])), values[0], values[1], '\t'.join(map(str, line_text[5:])), "\n"]))
                            annotLines.append(new_line)

                        else:
                            if transSQ == transInterest and verbose:
                                print("We are not looking for this annotation")
                                print(transSQ + " en feature " + fields[2])
                                print(fields)

            # else: #gene not associated in novelTrascript
            #     if(verbose):
            #         print("Gene not associated with no-match in novel Transcript -- Gene " + gene + " con el trans" + transSQ)

        if(anyAnnot): #REFERENCE ANNOT
            transcriptsAnnotated = transcriptsAnnotated + 1
        elif(not realNewTrans):
            if(referenceTransNoAnnot):
                if SAVE_PROB_TRANSCRIPTS:
                    file_reference_transcript_not_annot.write(transSQ + "\n")
                transcriptsAnnotated_wo_features = transcriptsAnnotated_wo_features + 1 #transcript has not info to transfered, so mark as checked it
            else:
                if SAVE_PROB_TRANSCRIPTS:
                    file_trans_not_annot_by_PF.write(transSQ + "\n")
                transcriptsNotAnnotated = transcriptsNotAnnotated + 1 #by positional transference

        elif(realNewTrans):
            if(anyNovelAnnotation): #NOVEL ANNOT
                #if(referenceGeneNoAnnot): #if novel annotated with another gene transcript
                novelTranscriptsAnnotated = novelTranscriptsAnnotated + 1
            elif(referenceGeneNoAnnot):
                if SAVE_PROB_TRANSCRIPTS:
                    file_reference_gene_not_annot.write(transSQ + "\n")
                novelTranscriptsAnnotated_wo_features = novelTranscriptsAnnotated_wo_features + 1 #transcript has not info to transfered, so mark as checked it
            else: #features exits in gene
                if SAVE_PROB_TRANSCRIPTS:
                    file_novel_not_annot_by_PF.write(transSQ + "\n")
                novelTranscriptsNotAnnotated = novelTranscriptsNotAnnotated + 1

    # new dictionary to remove duplicates
    realAnnotLines = {}
    annotInterest = ""
    for annot in annotLines:
        annot_aux = annot.split("\t")
        extra_info = ""
        if annot_aux[1] == "GeneOntology":
            extra_info = str(annot_aux[8].split(";")[0])
        text_annot = str(annot_aux[0:7]) + extra_info #check also attribute field
        # if annot_aux[0] == transInterest and verbose:
        #     if annot_aux[2] == "3UTRmotif":
        #         annotInterest = text_annot
        #         print("linea de interés guardada ;)")
        #     print(annot_aux)
        #     print(text_annot)
        #     print(not realAnnotLines.get(text_annot))
        if(not realAnnotLines.get(text_annot)):
            realAnnotLines.update({str(text_annot) : annot})

    # get prot annotation
    # we have to be aware in proteins because can appear domains that include smaller domains
    checkProteinAnnot = {}
    for annot in realAnnotLines.keys():
        annot_aux = realAnnotLines.get(annot)
        annot_aux = annot_aux.split("\t")

        if(annot_aux[-1].endswith("P\n")):
            prot_annot = str(annot_aux[0:3] + annot_aux[5:9])
            if(not annot_aux[3:4][0] == "." and not annot_aux[4:5][0] == "."):
                prot_pos_ini = int(annot_aux[3:4][0])
                prot_pos_fin = int(annot_aux[4:5][0])
            else: #just care about positions
                continue
            if(prot_pos_ini!=prot_pos_fin):
                if(not checkProteinAnnot.get(prot_annot)):
                    checkProteinAnnot.update({str(prot_annot) : [[prot_pos_ini, prot_pos_fin]]})
                else:
                    checkProteinAnnot.update({str(prot_annot) : checkProteinAnnot.get(prot_annot) + [[prot_pos_ini, prot_pos_fin]]})

    lst_to_delete = []
    for annot in checkProteinAnnot.keys():
        #get multiples values
        annot_to_save = annot.split("'")
        annot_to_save = [annot_to_save[1],annot_to_save[3],annot_to_save[5]]
        annot_aux = checkProteinAnnot.get(annot)
        if(len(annot_aux)>1):
            same_interval = False
            for val in annot_aux: #to study
                for aux_val in annot_aux:
                    if(val==aux_val):
                        next
                    else:
                        if(range_subset(range(val[0], val[1]), range(aux_val[0], aux_val[1]))):
                            same_interval = True
                            lst_to_delete.append(str(annot_to_save + [str(val[0])] + [str(val[1])]))
    
    if(len(lst_to_delete)>0):
        for rm in lst_to_delete:
            if(realAnnotLines.get(rm)):
                del realAnnotLines[rm] #remove if is the short range

    annotLines = []
    for annot in realAnnotLines.keys():
        annotLines.append(realAnnotLines.get(annot))

    for line in annotLines:
        f.write(line)

    f.close()

    if SAVE_PROB_TRANSCRIPTS:
        file_novel_not_annot_by_PF.close()
        file_reference_gene_not_annot.close()
        file_reference_transcript_not_annot.close()

    print("\n\n")

    #statistics
    # transcriptsAnnotated
    # transcriptsAnnotated_wo_features
    # transcriptsNotAnnotated

    # novelTranscriptsAnnotated
    # novelTranscriptsAnnotated_wo_features
    # novelTranscriptsNotAnnotated
    # transcriptsNotGeneID
    
    #TRANS_REFERENCE
    totalTranscripts = len(dc_SQexons)

    perctAnnotedMatchID = round((transcriptsAnnotated/totalTranscripts*100),2)
    annotedMatchID = transcriptsAnnotated

    perctAnnotedMatchID_wo_features = round((transcriptsAnnotated_wo_features/totalTranscripts*100),2)
    annotedMatchID_wo_features = transcriptsAnnotated_wo_features

    percTranscriptsNotAnnotated = round((transcriptsNotAnnotated/totalTranscripts*100),2)

    #NOVEL_TRANS
    percNovelTranscripts = round((novelTranscripts/totalTranscripts*100),2)

    percnovelTranscriptsAnnotated = round((novelTranscriptsAnnotated/totalTranscripts*100),2)
    percnovelTranscripts_wo_features = round((novelTranscriptsAnnotated_wo_features/totalTranscripts*100),2)
    percNovelTranscriptsNotAnnotated = round((novelTranscriptsNotAnnotated/totalTranscripts*100),2)
    percnovelTranscripts_no_gene_name = round((transcriptsNotGeneID/totalTranscripts*100),2)
    percnovelTranscripts_no_gene_match = round((novelTranscriptsNotGeneMatch/totalTranscripts*100),2)
    
    
    percTranscriptsAnnotated = round((transcriptsAnnotated/totalTranscripts*100),2)
    percTranscriptsNotAnnotatedNotGeneInformation = round((transcriptsNotAnnotatedNotGeneInformation/totalTranscripts*100),2)

    if(verbose):
        print(novelTranscripts, novelTranscriptsAnnotated, novelTranscriptsNotAnnotated, transcriptsAnnotated, transcriptsNotAnnotated, transcriptsNotAnnotatedNotGeneInformation)
    
    transcriptSection = "\tTRANSCRIPT-LEVEL SUMMARY\n"
    
    info_transAnnot = "\t\t·A total of " + "%.2f" % perctAnnotedMatchID + " % of transcripts in the GTF target annotation file (" + str(annotedMatchID) + " of " + str(totalTranscripts) + ") were annotated because they matched a transcript ID in the GFF3 source annotation file (SQANTI FSM and ISM)."
    info_transAnnot_wo_features = "\t\t·A total of " + "%.2f" % perctAnnotedMatchID_wo_features + " % of transcripts in the GTF target annotation file (" + str(transcriptsAnnotated_wo_features) + " of " + str(totalTranscripts) + ") were not annotated because reference has no annotations."
    info_transAnnot_notAnnot = "\t\t·A total of " + "%.2f" % percTranscriptsNotAnnotated + " % of transcripts in the GTF target annotation file (" + str(transcriptsNotAnnotated) + " of " + str(totalTranscripts) + ") were not annotated by positional feature transference."
    
    info_novelTrans = "\n\t\t·A total of " + "%.2f" % percNovelTranscripts + " % of transcripts in the GTF target annotation file (" + str(novelTranscripts) + " of " + str(totalTranscripts) + ") do not match any of the transcript IDs in the GFF3 file (SQANTI novel transcripts)."
    
    info_novelTransAnnot = "\t\t\t·A total of " + "%.2f" % percnovelTranscriptsAnnotated + " % novel transcripts (" + str(novelTranscriptsAnnotated) + " of " + str(totalTranscripts) + ") were annotated by positional feature transference."
    info_novelTransAnnot_wo_features = "\t\t\t·A total of " + "%.2f" % percnovelTranscripts_wo_features + " % of novel transcripts in the GTF target annotation file (" + str(novelTranscriptsAnnotated_wo_features) + " of " + str(totalTranscripts) + ") were not annotated because reference gene has no annotations in any transcript."
    info_novelTransAnnot_notAnnot = "\t\t\t·A total of " + "%.2f" % percNovelTranscriptsNotAnnotated + " % of novel transcripts (" + str(novelTranscriptsNotAnnotated) + " of " + str(totalTranscripts) + ") were not annotated by any positional feature."
    info_novelTransAnnot_noGeneName = "\t\t\t·A total of " + "%.2f" % percnovelTranscripts_no_gene_name + " % of novel transcripts in the GTF target annotation file (" + str(transcriptsNotGeneID) + " of " + str(totalTranscripts) + ") were not annotated because no gene name was returned by SQANTI output."
    info_novelTransAnnot_noGeneMatch = "\t\t\t·A total of " + "%.2f" % percnovelTranscripts_no_gene_match + " % of novel transcripts in the GTF target annotation file (" + str(novelTranscriptsNotGeneMatch) + " of " + str(totalTranscripts) + ") were not annotated because gene ID was not found in GFF3 source annotation file."

    globalTranscripts = "\n\t\tGlobal statistics:\n"

    perctTotalAnnotated = round(((transcriptsAnnotated+novelTranscriptsAnnotated)/totalTranscripts*100),2)
    total_transcripts_annotated = "\t\t·A total of " + "%.2f" % perctTotalAnnotated + " % of transcripts in the GTF target annotation file (" + str(int(transcriptsAnnotated+novelTranscriptsAnnotated))  + " of " + str(totalTranscripts) + ") were annotated (at least one feature transferred)."
    perctTotalNotAnnotated_no_features = round(((transcriptsAnnotated_wo_features+novelTranscriptsAnnotated_wo_features)/totalTranscripts*100),2)
    total_transcripts_wo_features = "\t\t·Did not annotate " + "%.2f" % perctTotalNotAnnotated_no_features + " % of transcripts (" + str(transcriptsAnnotated_wo_features + novelTranscriptsAnnotated_wo_features) + " of " + str(totalTranscripts) +  ") because reference matches have not annotation to transferred."
    perctTotalNotAnnotated_no_pf = round(((transcriptsNotAnnotated+novelTranscriptsNotAnnotated)/totalTranscripts*100),2)
    total_transcripts_no_pf = "\t\t·Did not annotate " + "%.2f" % perctTotalNotAnnotated_no_pf + " % of transcripts (" + str(transcriptsNotAnnotated + novelTranscriptsNotAnnotated) + " of " + str(totalTranscripts) + ") because no annotations were retrieved by positional feature transference with an existing TRANSCRIPT or GENE reference."
    total_transcripts_no_gene_name = "\t\t·Did not annotate " + "%.2f" % percnovelTranscripts_no_gene_name + " % of transcripts (" + str(transcriptsNotGeneID) + " of " + str(totalTranscripts) + ") because no gene name was returned by SQANTI output."
    total_transcripts_no_gene_match = "\t\t·Did not annotate " + "%.2f" % percnovelTranscripts_no_gene_match + " % of transcripts (" + str(novelTranscriptsNotGeneMatch) + " of " + str(totalTranscripts) + ") because gene ID was not found in GFF3 source annotation file."

    featureTransferenceSection = "\n\n\tFEATURE-TRANSFERENCE SUMMARY"
    transfTrans = ""
    for feat in transAnnotTranscriptCounts_dc.keys():
        #transfTrans = transfTrans +  "\t\t·" + str(feat)+": annotated " + "%.2f" % (len(transAnnotTranscriptCounts_dc.get(feat)[0]) / len(transAnnotTranscriptCounts_dc.get(feat)[1])*100) + " % of transcripts (" + str(len(transAnnotTranscriptCounts_dc.get(feat)[0])) + " of " + str(len(transAnnotTranscriptCounts_dc.get(feat)[1])) + ").\n"
        transfTrans = transfTrans +  "\t\t·" + str(feat)+": annotated " + "%.2f" % (len(transAnnotTranscriptCounts_dc.get(feat)[0].keys()) / totalTranscripts*100) + " % of transcripts (" + str(len(transAnnotTranscriptCounts_dc.get(feat)[0].keys())) + " of " + str(totalTranscripts) + ") with at least one " + feat + " feature.\n"

    transfProt = ""
    for feat in protAnnotTranscriptCounts_dc.keys():
        #transfProt = transfProt +  "\t\t·" + str(feat)+": annotated " + "%.2f" % (len(protAnnotTranscriptCounts_dc.get(feat)[0]) / len(protAnnotTranscriptCounts_dc.get(feat)[1])*100) + " % of transcripts (" + str(len(protAnnotTranscriptCounts_dc.get(feat)[0])) + " of " + str(len(transAnnotTranscriptCounts_dc.get(feat)[1])) + ").\n"
        transfProt = transfProt +  "\t\t·" + str(feat)+": annotated " + "%.2f" % (len(protAnnotTranscriptCounts_dc.get(feat)[0].keys()) / totalTranscripts*100) + " % of transcripts (" + str(len(protAnnotTranscriptCounts_dc.get(feat)[0].keys())) + " of " + str(totalTranscripts) + ") with at least one " + feat + " feature.\n"

    transfGene = ""
    for feat in geneAnnotTranscriptCounts_dc.keys():
        if(feat in ["C", "P", "F"]):
            #transfGene = transfGene +  "\t\t·" + str(feat)+" (GeneOntology): annotated " + "%.2f" % (len(geneAnnotTranscriptCounts_dc.get(feat)[0]) / len(geneAnnotTranscriptCounts_dc.get(feat)[1])*100) + " % of transcripts (" + str(len(geneAnnotTranscriptCounts_dc.get(feat)[0])) + " of " + str(len(transAnnotTranscriptCounts_dc.get(feat)[1])) + ").\n"
            transfGene = transfGene +  "\t\t·" + str(feat)+ " (GeneOntology): annotated " + "%.2f" % (len(geneAnnotTranscriptCounts_dc.get(feat)[0]) / totalTranscripts*100) + " % of transcripts (" + str(len(geneAnnotTranscriptCounts_dc.get(feat)[0])) + " of " + str(totalTranscripts) + ") with at least one " + feat + " feature.\n"
        else:
            #transfGene = transfGene +  "\t\t·" + str(feat)+": annotated " + "%.2f" % (len(geneAnnotTranscriptCounts_dc.get(feat)[0]) / len(geneAnnotTranscriptCounts_dc.get(feat)[1])*100) + " % of transcripts (" + str(len(geneAnnotTranscriptCounts_dc.get(feat)[0])) + " of " + str(len(transAnnotTranscriptCounts_dc.get(feat)[1])) + ").\n"
            transfGene = transfGene +  "\t\t·" + str(feat)+": annotated " + "%.2f" % (len(geneAnnotTranscriptCounts_dc.get(feat)[0].keys()) / totalTranscripts*100) + " % of transcripts (" + str(len(geneAnnotTranscriptCounts_dc.get(feat)[0].keys())) + " of " + str(totalTranscripts) + ") with at least one " + feat + " feature.\n"


    featureLevelSection = "\n\tFEATURE-LEVEL SUMMARY"
    mTrans = ""

    realFeaturesAnnoted = 0
    realTotalFeaturesChecked = 0
    for feat in transAnnotTotal_dc.keys():
        mTrans = mTrans +  "\t\t·" + str(feat)+": transferred " + "%.2f" % (len(transAnnotTotal_dc.get(feat)[0].keys())/len(transAnnotTotal_dc.get(feat)[1].keys())*100) + " % of features (" + str(len(transAnnotTotal_dc.get(feat)[0].keys())) + " of " + str(len(transAnnotTotal_dc.get(feat)[1].keys())) + ").\n"
        realFeaturesAnnoted = realFeaturesAnnoted + len(transAnnotTotal_dc.get(feat)[0].keys())
        realTotalFeaturesChecked = realTotalFeaturesChecked + len(transAnnotTotal_dc.get(feat)[1].keys())

    mProt = ""
    for feat in protAnnotTotal_dc.keys():
        mProt = mProt +  "\t\t·" + str(feat) + ": transferred " + "%.2f" % (len(protAnnotTotal_dc.get(feat)[0].keys())/len(protAnnotTotal_dc.get(feat)[1].keys())*100) + " % of features (" + str(len(protAnnotTotal_dc.get(feat)[0].keys())) + " of " + str(len(protAnnotTotal_dc.get(feat)[1].keys())) + ").\n"
        realFeaturesAnnoted = realFeaturesAnnoted + len(protAnnotTotal_dc.get(feat)[0].keys())
        realTotalFeaturesChecked = realTotalFeaturesChecked + len(protAnnotTotal_dc.get(feat)[1].keys())

    mGene = ""
    for feat in geneAnnotTotal_dc.keys():
        if(feat in ["C", "P", "F"]):
            mGene = mGene +  "\t\t·" + str(feat) + " (GeneOntology): transferred " + "%.2f" % (len(geneAnnotTotal_dc.get(feat)[0].keys())/len(geneAnnotTotal_dc.get(feat)[1].keys())*100) + " % of features (" + str(len(geneAnnotTotal_dc.get(feat)[0].keys())) + " of " + str(len(geneAnnotTotal_dc.get(feat)[1].keys())) + ").\n"
            realFeaturesAnnoted = realFeaturesAnnoted + len(geneAnnotTotal_dc.get(feat)[0].keys())
            realTotalFeaturesChecked = realTotalFeaturesChecked + len(geneAnnotTotal_dc.get(feat)[1].keys())
        else:
            mGene = mGene +  "\t\t·" + str(feat) + ": transferred " + "%.2f" % (len(geneAnnotTotal_dc.get(feat)[0].keys())/len(geneAnnotTotal_dc.get(feat)[1].keys())*100) + " % of features (" + str(len(geneAnnotTotal_dc.get(feat)[0].keys())) + " of " + str(len(geneAnnotTotal_dc.get(feat)[1].keys())) + ").\n"
            realFeaturesAnnoted = realFeaturesAnnoted + len(geneAnnotTotal_dc.get(feat)[0].keys())
            realTotalFeaturesChecked = realTotalFeaturesChecked + len(geneAnnotTotal_dc.get(feat)[1].keys())

    anyNotAnnotated = False
    if((transcriptsNotAnnotated + novelTranscriptsNotAnnotated) > 0):
        anyNotAnnotated = True
    
    anyNovel = False
    if(novelTranscripts!=0):
        anyNovel = True

    anyNovelNotAnnot = False
    if(percNovelTranscriptsNotAnnotated!=0):
        anyNovelNotAnnot = True

    anyTransWithoutFeatures = False
    if(percTranscriptsNotAnnotated!=0):
        anyTransWithoutFeatures = True

    print(transcriptSection)
    print(info_transAnnot)
    print(info_transAnnot_wo_features)
    print(info_transAnnot_notAnnot)
    if(anyNovel):
        print(info_novelTrans)
        print(info_novelTransAnnot)
        print(info_novelTransAnnot_wo_features)
        print(info_novelTransAnnot_notAnnot)
        print(info_novelTransAnnot_noGeneName)
        print(info_novelTransAnnot_noGeneMatch)

    print(globalTranscripts)
    print(total_transcripts_annotated)
    if(anyNotAnnotated):
        print(total_transcripts_wo_features)
        print(total_transcripts_no_pf)
        print(total_transcripts_no_gene_match)
    
    if(percTranscriptsAnnotated!=0):
        #Transference section
        print(featureTransferenceSection)
        m7 = "\n\t\tTranscript Annotation:\n"
        print(m7)
        print(transfTrans)

        m8 = "\t\tProtein Annotation:\n"
        print(m8)
        print(transfProt)

        m9 = "\t\tGene Annotation:\n"
        print(m9)
        print(transfGene)

        if STATS:
            print(featureLevelSection)
            m7 = "\n\t\tTranscript Annotation:\n"
            print(m7)
            print(mTrans)

            m8 = "\t\tProtein Annotation:\n"
            print(m8)
            print(mProt)

            m9 = "\t\tGene Annotation:\n"
            print(m9)
            print(mGene)

            # if(anotationsChecked==0):
            #     perct_annot = 0
            # else:
            #     perct_annot = featuresAnnotated/anotationsChecked*100

            if(realTotalFeaturesChecked==0):
                perct_annotReal = 0
            else:
                perct_annotReal = round(realFeaturesAnnoted/realTotalFeaturesChecked*100,2)

            globalFeatures = "\t\tGlobal statistics:"
            m10 = "\t\t·Annotated a total of " + "%.2f" % perct_annotReal + " % (" + str(realFeaturesAnnoted) + " of " + str(realTotalFeaturesChecked) + ") features from the reference GFF3 file for the study transcripts."
            print(globalFeatures)
            print(m10 + "\n\n")

        else:
            print("\n")
    else:
        print("\n")

    if USE_STDOUT:
        f = open(filename_prints, "w")
        f.write(transcriptSection)
        f.write("\n" + info_transAnnot + "\n")
        f.write(info_transAnnot_wo_features + "\n")
        f.write(info_transAnnot_notAnnot + "\n")
        f.write(info_novelTrans + "\n")
        f.write(info_novelTransAnnot + "\n")
        f.write(info_novelTransAnnot_notAnnot + "\n")
        f.write(info_novelTransAnnot_noGeneName + "\n")
        f.write(info_novelTransAnnot_noGeneMatch + "\n")
        f.write(globalTranscripts)
        f.write(total_transcripts_annotated + "\n")
        if(anyNotAnnotated):
            f.write(total_transcripts_wo_features + "\n")
            f.write(total_transcripts_no_pf + "\n")
            f.write(total_transcripts_no_gene_name + "\n")
            f.write(total_transcripts_no_gene_match + "\n")

        if(percTranscriptsAnnotated!=0):
            #transference
            f.write (featureTransferenceSection)
            f.write("\n" + m7 + "\n")
            f.write(transfTrans + "\n")
            f.write(m8 + "\n")
            f.write(transfProt + "\n")
            f.write(m9 + "\n")
            f.write(transfGene + "\n")
            if STATS:
                f.write (featureLevelSection)
                f.write("\n" + m7 + "\n")
                f.write(mTrans + "\n")
                f.write(m8 + "\n")
                f.write(mProt + "\n")
                f.write(m9 + "\n")
                f.write(mGene + "\n")
                f.write(globalFeatures)
                f.write("\n" + m10)
    
#UPDATE GFF3 - new columns information
def addPosType(res, line, posType):
    global verbose
    if line.endswith(";"):
        res.write(line + " PosType=" + posType + "\n")
    else:
        res.write(line[:-1] + "; PosType=" + posType + "\n")

def updateGTF(filename, filenameMod):
    global verbose
    # open new file
    res = open(filenameMod, "w")
    # open annotation file and process all data
    with open(filename, 'r') as f:
        # process all entries - no header line in file
        for line in f:
            if len(line) == 0:
                break
            else:
                if line and line[0] != "#":
                    fields = line.split("\t")
                    if len(fields) == 9:
                        
                        text = fields[8].split(" ")
                        if text[-1].startswith("PosType"):
                            res.write(line)

                        elif fields[1] == "tappAS":
                            if fields[2] == 'transcript':
                                addPosType(res, line, "T")
                            elif fields[2] == 'gene':
                                addPosType(res, line, "T")
                            elif fields[2] == 'CDS':
                                addPosType(res, line, "T")   
                            elif fields[2] == 'genomic':
                                addPosType(res, line, "G")
                            elif fields[2] == 'exon':
                                addPosType(res, line, "G")
                            elif fields[2] == 'splice_junction':
                                addPosType(res, line, "G")
                            elif fields[2] == 'protein':
                                addPosType(res, line, "P")
                            else:
                                print(line)
                                break

                        elif fields[1] == "COILS":   
                            if fields[2] == 'COILED':
                                addPosType(res, line, "P")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using P type to annotate.")
                                addPosType(res, line, "P")
                                #break

                        elif fields[1] == "GeneOntology":
                            if fields[2] in ('C', 'cellular_component'):
                                addPosType(res, line, "N")
                            elif fields[2] in ('F', 'molecular_function'):
                                addPosType(res, line, "N")
                            elif fields[2] in ('P', 'biological_process'):
                                addPosType(res, line, "N")
                            elif fields[2] in ('eco'):
                                addPosType(res, line, "N") # Fran tomato annot
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using N type to annotate.")
                                addPosType(res, line, "N")
                                ##break

                        elif fields[1] == "MOBIDB_LITE": 
                            if fields[2] == 'DISORDER' :
                                addPosType(res, line, "P")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using P type to annotate.")
                                addPosType(res, line, "P")
                                #break

                        elif fields[1] == "NMD": 
                            if fields[2] == 'NMD':
                                addPosType(res, line, "T")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using T type to annotate.")
                                addPosType(res, line, "T")
                                #break

                        elif fields[1] in ("PAR-CLIP", "PAR-clip"):
                            if fields[2] in ('RNA_binding', 'RNA_Binding_Protein', 'RBP_Binding') or fields[2].startswith('RNA_binding_'):
                                addPosType(res, line, "T")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using T type to annotate.")
                                addPosType(res, line, "T")
                                #break

                        elif fields[1] == "PFAM":    
                            if fields[2] == 'DOMAIN':
                                addPosType(res, line, "P")
                            elif fields[2] in ("CLAN","clan"):
                                addPosType(res, line, "N")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using N type to annotate.")
                                addPosType(res, line, "N")
                                #break

                        elif fields[1] == "Provean": 
                            if fields[2] == 'FunctionalImpact':
                                addPosType(res, line, "N")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using N type to annotate.")
                                addPosType(res, line, "N")
                                #break

                        elif fields[1] in ("REACTOME","Reactome"):
                            if fields[2] in ('PATHWAY','pathway', 'Pathway'):
                                addPosType(res, line, "N")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using N type to annotate.")
                                addPosType(res, line, "N")
                                #break

                        elif fields[1] == "RepeatMasker":
                            if fields[2] == 'repeat':
                                addPosType(res, line, "T")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using T type to annotate.")
                                addPosType(res, line, "T")
                                #break

                        elif fields[1] == "SIGNALP_EUK": 
                            if fields[2] == 'SIGNAL':
                                addPosType(res, line, "P")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using P type to annotate.")
                                addPosType(res, line, "P")
                                #break

                        elif fields[1] == "TMHMM":   
                            if fields[2] == 'TRANSMEM':
                                addPosType(res, line, "P")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using P type to annotate.")
                                addPosType(res, line, "P")
                                #break

                        elif fields[1] == "TranscriptAttributes":
                            addPosType(res, line, "T")

                        elif fields[1] == "UTRsite": 
                            if fields[2] == 'uORF':
                                addPosType(res, line, "T")
                            elif fields[2] == '5UTRmotif':
                                addPosType(res, line, "T")
                            elif fields[2] == 'PAS':
                                addPosType(res, line, "T")
                            elif fields[2] == '3UTRmotif':
                                addPosType(res, line, "T")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using T type to annotate.")
                                addPosType(res, line, "T")
                                #break

                        elif fields[1] in ("UniProtKB/Swiss-Prot_Phosphosite", "Swissprot_Phosphosite"):
                            if fields[2] == 'ACT_SITE':
                                addPosType(res, line, "P")
                            elif fields[2] == 'BINDING':
                                addPosType(res, line, "P")
                            elif fields[2] == 'PTM':
                                addPosType(res, line, "P")
                            elif fields[2] == 'MOTIF':
                                addPosType(res, line, "P")
                            elif fields[2] == 'COILED':
                                addPosType(res, line, "P")
                            elif fields[2] == 'TRANSMEM':
                                addPosType(res, line, "P")
                            elif fields[2] == 'COMPBIAS':
                                addPosType(res, line, "P")
                            elif fields[2] == 'INTRAMEM':
                                addPosType(res, line, "P")
                            elif fields[2] == 'NON_STD':
                                addPosType(res, line, "P")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using P type to annotate.")
                                addPosType(res, line, "P")
                                #break

                        elif fields[1] in ("cNLS_mapper", "NLS_mapper"): 
                            if fields[2] == 'MOTIF': 
                                addPosType(res, line, "P")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using P type to annotate.")
                                addPosType(res, line, "P")
                                #break

                        elif fields[1] in ("miRWalk", "mirWalk"): 
                            if fields[2] in('miRNA', 'miRNA_Binding'):
                                addPosType(res, line, "T")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using T type to annotate.")
                                addPosType(res, line, "T")
                                #break

                        elif fields[1] == "scanForMotifs":   
                            if fields[2] == 'PAS':
                                addPosType(res, line, "T")
                            elif fields[2] in ('3UTRmotif', "3'UTRmotif"):
                                addPosType(res, line, "T")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using T type to annotate.")
                                addPosType(res, line, "T")
                                #break

                        elif fields[1] == "MetaCyc":
                            if fields[2] == 'pathway':
                                addPosType(res, line, "N")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using N type to annotate.")
                                addPosType(res, line, "N")
                                #break

                        elif fields[1] == "KEGG":
                            if fields[2] in ('pathway','Pathway'):
                                addPosType(res, line, "N")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using N type to annotate.")
                                addPosType(res, line, "N")
                                #break

                        elif fields[1] == "SUPERFAMILY":
                            if fields[2] == 'DOMAIN':
                                addPosType(res, line, "P")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using P type to annotate.")
                                addPosType(res, line, "P")
                                #break

                        elif fields[1] == "SMART":
                            if fields[2] == 'DOMAIN':
                                addPosType(res, line, "P")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using P type to annotate.")
                                addPosType(res, line, "P")
                                #break

                        elif fields[1] == "TIGRFAM":
                            if fields[2] == 'DOMAIN':
                                addPosType(res, line, "P")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using P type to annotate.")
                                addPosType(res, line, "P")
                                #break

                        elif fields[1] == "psRNATarget":
                            if fields[2] == 'miRNA':
                                addPosType(res, line, "T")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using T type to annotate.")
                                addPosType(res, line, "T")
                                #break

                        elif fields[1] == "CORUM":
                            if fields[2] == 'Complex':
                                addPosType(res, line, "P")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using P type to annotate.")
                                addPosType(res, line, "P")
                                #break

                        elif fields[1] == "Orthologues":
                            if fields[2] == 'S.tuberosum':
                                addPosType(res, line, "N")
                            elif fields[2] in ('A.thaliana'):
                                addPosType(res, line, "N")
                            else:
                                print("IsoAnnotLite can not identify the feature " + str(fields[2]) + " in source " + str(fields[1]) + ", using N type to annotate.")
                                addPosType(res, line, "N")
                                #break

                        else:
                            print("IsoAnnotLite can not identify the source " + str(fields[1]) + ", in line:\n" + line + "\nUsing N type to annotate. You should edit method 'updateGTF' and add your feature type following the same structure than the other features.")
                            addPosType(res, line, "N")
                            #break

                    else:
                        print(fields)
                        print("Error in line (has not 9 fields):\n" + line)
                        break
        
        res.close()

def readGFFandGetData(filenameMod):
    global verbose
    # open annotation file and process all data
    dcTrans = {}
    dcExon = {}
    dcTransFeatures = {}
    dcGenomic = {}
    dcSpliceJunctions = {}
    dcProt = {}
    dcProtFeatures = {}
    dcTranscriptAttributes = {}
    

    dcTransID = {}

    with open(filenameMod, 'r') as f:
        # process all entries - no header line in file
        for line in f:
            if len(line) == 0:
                break
            else:
                if line and line[0] != "#":
                    fields = line.split("\t")
                    if len(fields) == 9:
                        
                        transcript = fields[0]
                        text = fields[8].split(" ")
                        #transcriptID = text[0]
                        #transcriptID = transcriptID[3:-1]

                        if fields[1] == "tappAS":
                            if fields[2] in ["transcript", "gene", "CDS"]:
                                if not dcTrans.get(str(transcript)):
                                    dcTrans.update({str(transcript) : [line]})
                                else:                
                                    dcTrans.update({str(transcript) : dcTrans.get(str(transcript)) + [line]})
                            elif fields[2] in ["exon"]:
                                if not dcExon.get(str(transcript)):
                                    dcExon.update({str(transcript) : [line]})
                                else:                
                                    dcExon.update({str(transcript) : dcExon.get(str(transcript)) + [line]})
                            elif fields[2] in ["genomic"]:
                                if not dcGenomic.get(str(transcript)):
                                    dcGenomic.update({str(transcript) : [line]})
                                else:                
                                    dcGenomic.update({str(transcript) : dcGenomic.get(str(transcript)) + [line]})
                            elif fields[2] in ["splice_junction"]:
                                if not dcSpliceJunctions.get(str(transcript)):
                                    dcSpliceJunctions.update({str(transcript) : [line]})
                                else:                
                                    dcSpliceJunctions.update({str(transcript) : dcSpliceJunctions.get(str(transcript)) + [line]})
                            elif fields[2] in ["protein"]:
                                if not dcProt.get(str(transcript)):
                                    dcProt.update({str(transcript) : [line]})
                                else:                
                                    dcProt.update({str(transcript) : dcProt.get(str(transcript)) + [line]})
                        #Transcript Information
                        elif fields[1] == "TranscriptAttributes":
                            if not dcTranscriptAttributes.get(str(transcript)):
                                dcTranscriptAttributes.update({str(transcript) : [line]})
                            else:                
                                dcTranscriptAttributes.update({str(transcript) : dcTranscriptAttributes.get(str(transcript)) + [line]})
                        #Feature information
                        else:
                            if text[-1].endswith("T\n"):
                                if not dcTransFeatures.get(str(transcript)):
                                    dcTransFeatures.update({str(transcript) : [line]})
                                else:                
                                    dcTransFeatures.update({str(transcript) : dcTransFeatures.get(str(transcript)) + [line]})
                            elif text[-1].endswith("P\n") or text[-1].endswith("G\n") or text[-1].endswith("N\n"):
                                if not dcProtFeatures.get(str(transcript)):
                                    dcProtFeatures.update({str(transcript) : [line]})
                                else:
                                    dcProtFeatures.update({str(transcript) : dcProtFeatures.get(str(transcript)) + [line]})
    
    return dcTrans, dcExon, dcTransFeatures, dcGenomic, dcSpliceJunctions, dcProt, dcProtFeatures, dcTranscriptAttributes

def updateGeneDescriptions(dcTrans, dc_GFF3_raw_annot):
    for SQtrans in dcTrans.keys():
        lines = dcTrans.get(SQtrans)
        cont = 0
        toUpdate = False
        for line in lines:
            fields = line.split("\t")
            if fields[2]=="gene":
                if dc_GFF3_raw_annot.get(SQtrans):
                    linesGFF3 = dc_GFF3_raw_annot.get(SQtrans)
                    for lineGFF3 in linesGFF3:
                        fieldsGFF3 = lineGFF3.split("\t") #line in 2 position
                        if fieldsGFF3[2]=="gene":
                            fields[8] = fieldsGFF3[8]
                            res = '\t'.join(map(str, fields))
                            line = res
                            lines[cont] = line
                            toUpdate = True
            cont = cont + 1
        if toUpdate:
            dcTrans.update({str(SQtrans) : lines})

def generateFinalGFF3(dcTrans, dcExon, dcTransFeatures, dcGenomic, dcSpliceJunctions, dcProt, dcProtFeatures, dcTranscriptAttributes, filename):
    global verbose
    # open new file
    res = open(filename, "w")
    strand = ""
    for SQtrans in dcTrans.keys():
        t = dcTrans.get(SQtrans)
        strand = t[0].split("\t")
        strand = strand[6]
        if t:
            for line in t:
                res.write(line)
            
        tf = dcTransFeatures.get(SQtrans)
        if tf: 
            for line in tf:
                res.write(line)

        g = dcGenomic.get(SQtrans)
        if g:
            for line in g:
                res.write(line)

        e = dcExon.get(SQtrans)
        if e:
            if strand == "+":
                for line in e:
                    res.write(line)
            else:
                for i in range(len(e)-1,-1,-1):
                    res.write(e[i])

        sj = dcSpliceJunctions.get(SQtrans)
        if sj:
            for line in sj:
                res.write(line)

        p = dcProt.get(SQtrans)
        if p:
            for line in p:
                res.write(line)

        pf = dcProtFeatures.get(SQtrans)
        if pf:
            for line in pf:
                res.write(line)

        ta = dcTranscriptAttributes.get(SQtrans)
        if ta:
            for line in ta:
                res.write(line)
    
    res.close()

############
# Parámetros
############

# -GTF (Corrected) de SQANTI3
# -Classification de SQANTI3
# -Junctions de SQANTI3
# -GFF3 de referencia
# output name

def main():
    global USE_GFF3
    global USE_NAME
    global USE_STDOUT
    global ALL_AS_NOVELS
    global INTRONIC
    global STATS
    global SAVE_PROB_TRANSCRIPTS
    global version
    global verbose
    #arguments
    parser = argparse.ArgumentParser(description="IsoAnnotLite " + str(version) + ": Transform SQANTI 3 output files to generate GFF3 to tappAS.")
    parser.add_argument('corrected', help='\t\t*_corrected.gtf file from SQANTI 3 output.') 
    parser.add_argument('classification', help='\t\t*_classification.txt file from SQANTI 3 output.')
    parser.add_argument('junction', help='\t\t*_junctions.txt file from SQANTI 3 output.')
    parser.add_argument('-gff3', help='\t\ttappAS GFF3 file to map its annotation to your SQANTI 3 data (only if you use the same reference genome in SQANTI 3).', required = False)
    parser.add_argument('-o', help='\t\tOutput name for the custom GFF3 file.', required = False)
    parser.add_argument('-stdout', help='\t\Output name where save all the print results (only when -gff3).', required = False)
    parser.add_argument('-novel', help='\t\Annotate transcripts using all gene information instead using only the transcript of reference (just for transcripts with reference).', required = False, action='store_true')
    parser.add_argument('-nointronic', help='\t\Do not annotate intronic features.', required = False, action='store_true')
    parser.add_argument('-statistics', help='\t\Show Feature Level Summary (statistics) [currently not used, by default we show all the statistics results].', required = False, action='store_true')
    parser.add_argument('-saveTranscriptIDs', help='\t\Save problematic transcript IDs in five different files. Transcripts not annotated by positional transference ("file_trans_not_annot_by_PF.txt"). Novel transcripts not annotated by positional transference ("file_novel_not_annot_by_by_PF.txt"). SQ3 reference gene not found in GFF3 annotation ("file_reference_gene_not_annot.txt"). Transcripts not annotated because any features were found in GFF3 annotation ("file_reference_transcript_not_annot.txt"). SQ reference gene field is empty ("file_transcript_wo_gene_ID.txt").', required = False, action='store_true')

    args = parser.parse_args()

    # path and prefix for output files
    args.corrected = os.path.abspath(args.corrected)
    if not os.path.isfile(args.corrected):
        sys.stderr.write("ERROR: '%s' doesn't exist\n" %(args.corrected))
        sys.exit()

    args.classification = os.path.abspath(args.classification)
    if not os.path.isfile(args.classification):
        sys.stderr.write("ERROR: '%s' doesn't exist\n" %(args.classification))
        sys.exit()

    args.junction = os.path.abspath(args.junction)
    if not os.path.isfile(args.junction):
        sys.stderr.write("ERROR: '%s' doesn't exist\n" %(args.junction))
        sys.exit()

    if args.gff3:
        USE_GFF3 = True
        args.gff3 = os.path.abspath(args.gff3)
        if not os.path.isfile(args.gff3):
            sys.stderr.write("ERROR: '%s' doesn't exist\n" %(args.gff3))
            sys.exit()

    if args.o:
        USE_NAME = True
        args.o = "".join(args.o)
        if len(args.o)==0:
            sys.stderr.write("ERROR: -o has not value\n")
            sys.exit()

    if args.stdout:
        USE_STDOUT = True
        args.stdout = "".join(args.stdout)
        if len(args.stdout)==0:
            sys.stderr.write("ERROR: -stdout has not value\n")
            sys.exit()

    if args.novel:
        ALL_AS_NOVELS = True
    else:
        ALL_AS_NOVELS = False

    if args.nointronic:
        INTRONIC = False

    if args.statistics:
        STATS = True
    else: 
        STATS = True #works perfectly

    if args.saveTranscriptIDs:
        SAVE_PROB_TRANSCRIPTS = True
    else:
        SAVE_PROB_TRANSCRIPTS = False

    # Running functionality
    if ALL_AS_NOVELS:
        sys.stdout.write("\n\nRunning IsoAnnot Lite " + str(version) + " (using all gene information) ...\n")
    else:
        sys.stdout.write("\n\nRunning IsoAnnot Lite " + str(version) + "...\n")
    run(args)

def run(args):
    import time
    global USE_GFF3
    global USE_NAME
    global verbose
    
    t1 = time.time()
    #corrected = input("Enter your file name for \"corrected.gtf\" file from SQANTI 3 (with extension): ")
    gtf = args.corrected
    #classification = input("Enter your file name for \"classification.txt\" file from SQANTI 3 (with extension): ")
    classification = args.classification
    #junctions = input("Enter your file name for \"junctions.txt\" file from SQANTI 3 (with extension): ")
    junctions = args.junction
    #GFF3 download from tappAS.org/downloads

    ########################
    # MAPPING SQANTI FILES #
    ########################

    if USE_GFF3:
        gff3 = args.gff3
    
        #File names
        if USE_NAME:
            if ".gff3" in args.o:
                filename = args.o
            else:    
                filename = args.o + ".gff3"
            filenameMod = filename[:-5] + "_mod" + filename[-5:]
        else:
            filename = "tappAS_annot_from_SQANTI3.gff3"
            filenameMod = filename[:-5] + "_mod" + filename[-5:]

        if USE_STDOUT:
            if ".txt" in args.stdout:
                filename_prints = args.stdout
            else:
                filename_prints = args.stdout + ".txt"
        else:
            filename_prints = "output.txt"

        #################
        # START PROCESS #
        #################
        print("Reading reference annotation file and creating data variables...")
        #dc_GFF3 = {trans : [[start,end,line], [start,end,line], ...]}
        #dc_GFF3exonsTrans = {start : [trans, trans, ...]}
        #dc_GFF3transExons = {trans : [[start,end], [start,end]...]}
        #dc_GFF3coding = {trans : [CDSstart, CDSend]}
        #dc_GFF3geneTrans = {gene : [trans1, trans2...]}
        #dc_GFF3_raw_annot = {trans : [line, line...]} #for not annot lines
        #dc_gene_description = {gene : [ID=...]} #description line
        dc_GFF3, dc_GFF3exonsTrans, dc_GFF3transExons, dc_GFF3coding, dc_GFF3strand, dc_GFF3geneTrans, dc_GFF3transGene, dc_GFF3geneNameTrans, dc_GFF3_raw_annot, dc_gene_description = readGFF(gff3) #dc_GFF3exons is sorted

        print("\nReading SQANTI 3 Files and creating an auxiliar GFF...")
        #dc_SQexons = {trans : [[start,end], [start,end]...]}
        #dc_SQcoding = {trans : [CDSstart, CDSend, orf]}
        #dc_SQtransGene = {trans : [gene, category, transAssociated]}
        #dc_SQstrand = {trans : strand}
        dc_SQexons, dc_SQcoding, dc_SQtransGene, dc_SQstrand, dc_geneID2geneName, dc_geneName2geneID = createGTFFromSqanti(gtf, classification, junctions, dc_gene_description, filename) 
        
        #create gene-trans dictionary for SQ file
        dc_SQgeneTrans = createDCgeneTrans(dc_SQtransGene)
        
        print("Transforming CDS local positions to genomic position...")
        #Transformar características a posiciones genómicas
        #dc_GFF3coding = {trans : [exon1, exon2...]}
        dc_SQcoding = transformCDStoGenomic(dc_SQcoding, dc_SQexons, dc_SQstrand)
        dc_GFF3coding = transformCDStoGenomic(dc_GFF3coding, dc_GFF3transExons, dc_GFF3strand)

        print("Transforming feature local positions to genomic position in GFF3...")
        #Transformar características a posiciones genómicas
        dc_GFF3_Genomic = transformTransFeaturesToGenomic(dc_GFF3, dc_GFF3transExons, dc_GFF3coding, dc_GFF3strand)

        print("Generating Transcriptome per each gene...") #dc_GFF3_Genomic
        dc_GFF3Gene_Genomic = getTranscriptomePerGene(dc_SQgeneTrans, dc_GFF3_Genomic) # {gene : [[start,end,line], [start,end,line], ...]}

        print("Mapping transcript features between GFFs...")
        mappingFeatures(dc_SQexons, dc_SQcoding, dc_SQtransGene, dc_SQgeneTrans, dc_SQstrand, dc_GFF3exonsTrans, dc_GFF3transExons, dc_GFF3_Genomic, dc_GFF3Gene_Genomic, dc_GFF3coding, dc_GFF3geneTrans, dc_geneID2geneName, dc_geneName2geneID, dc_GFF3transGene, filename, filename_prints) #edit tappAS_annotation_from_Sqanti file
        #dc_geneID2geneName, dc_geneName2geneID

        print("Adding extra information to GFF3 columns...")
        updateGTF(filename, filenameMod)

        print("Reading GFF3 to sort it correctly...")
        dcTrans, dcExon, dcTransFeatures, dcGenomic, dcSpliceJunctions, dcProt, dcProtFeatures, dcTranscriptAttributes = readGFFandGetData(filenameMod)

        print("Updating Gene Descriptions...")
        updateGeneDescriptions(dcTrans, dc_GFF3_raw_annot)
        #dcTrans = updateGeneDescriptions(dcTrans, dc_GFF3)

        #Remove old files
        os.remove(filename)
        os.remove(filenameMod)

        dcTransFeatures = transformTransFeaturesToLocale(dcTransFeatures, dc_SQexons)

        print("Generating final GFF3...")
        generateFinalGFF3(dcTrans, dcExon, dcTransFeatures, dcGenomic, dcSpliceJunctions, dcProt, dcProtFeatures, dcTranscriptAttributes, filename)

        t2 = time.time()
        time = t2-t1
        print("Time used to generate new GFF3: " + "%.2f" % time + " seconds.\n")

        print("Exportation complete.\nYour GFF3 result is: '" + filename + "'\n")

    #####################
    # JUST SQANTI FILES #
    #####################

    else:
        #File names
        if USE_NAME:
            if ".gff3" in args.o:
                filename = args.o
            else:    
                filename = args.o + ".gff3"
            filenameMod = filename[:-5] + "_mod" + filename[-5:]
        else:
            filename = "tappAS_annot_from_SQANTI3.gff3"
            filenameMod = filename[:-5] + "_mod" + filename[-5:]

        #################
        # START PROCESS #
        #################
        print("\nReading SQANTI 3 Files and creating an auxiliar GFF...")

        #dc_SQexons = {trans : [[start,end], [start,end]...]}
        #dc_SQcoding = {trans : [CDSstart, CDSend, orf]}
        #dc_SQtransGene = {trans : [gene, category, transAssociated]}

        #we do not have dc_gene_descripton - NO-GFF3 Annotation
        dc_gene_description = {}
        dc_SQexons, dc_SQcoding, dc_SQtransGene, dc_SQstrand, dc_geneID2geneName, dc_geneName2geneID = createGTFFromSqanti(gtf, classification, junctions, dc_gene_description, filename) 

        print("Adding extra information to relative columns...")
        updateGTF(filename, filenameMod)

        print("Reading GFF3 to sort it correctly...")
        dcTrans, dcExon, dcTransFeatures, dcGenomic, dcSpliceJunctions, dcProt, dcProtFeatures, dcTranscriptAttributes = readGFFandGetData(filenameMod)

        #Remove old files
        os.remove(filename)
        os.remove(filenameMod)

        print("Generating final GFF3...")
        generateFinalGFF3(dcTrans, dcExon, dcTransFeatures, dcGenomic, dcSpliceJunctions, dcProt, dcProtFeatures, dcTranscriptAttributes, filename)

        t2 = time.time()
        time = t2-t1
        print("Time used to generate new GFF3: " + "%.2f" % time + " seconds.\n")

        print("Exportation complete.\nYour GFF3 result is: '" + filename + "'\n")

if __name__ == "__main__":
    main()
