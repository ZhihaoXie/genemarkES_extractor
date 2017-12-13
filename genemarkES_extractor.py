#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# FileName:  genemarkES_extractor.py
# Author:    Zhihao Xie  \(#^o^)/
# Date:      2017/11/28 9:50
# Version:   v1.0.0
# CopyRight: Copyright ©Zhihao Xie, All rights reserved.

import re, sys, os
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if len(sys.argv) < 4:
    print("usage: python3 %s <genemarkES_out(gff)> <genome_fasta> <out_prefix>" % sys.argv[0])
    sys.exit(1)

gff_file = os.path.abspath(sys.argv[1])
genome = os.path.abspath(sys.argv[2])
output = os.path.abspath(sys.argv[3])

if len(sys.argv) >= 5:
    codon = int(sys.argv[4])
else:
    codon = 1

gene_site = {}
gene_strand = {}
with open(gff_file) as fh:
    for line in fh:
        if re.search(r"^#|^\s*$", line):
            continue
        elif len(line) == 0:
            break
        fields = line.strip().split("\t")
        if fields[2] == "CDS":
            start = fields[3]
            end = fields[4]
            strand = fields[6]
            phase = fields[7]
            if phase == ".":
                phase = 0
            else:
                phase = int(phase)
            gene_id = re.findall(r"gene_id(\s|=)\"?(\w+)\"?;?", fields[-1])[0][1]
            if not gene_id:
                gene_id = re.findall(r"transcript_id(\s|=)\"?(\w+)\"?;?", fields[-1])[0][1]
            gene_site.setdefault(fields[0],{}).setdefault(gene_id,[])
            gene_site[fields[0]][gene_id].append((start,end,phase))
            gene_strand.setdefault(fields[0],{}).setdefault(gene_id, strand)


fasta_hash = {}
for rec in SeqIO.parse(genome, 'fasta'):
    fasta_hash[rec.id] = rec

output_fh = open(output + ".fnn", 'w')
output_fh2 = open(output + ".faa", 'w')
for contig_id in sorted(fasta_hash.keys()):
    if contig_id in gene_site:
        number = 1
        for gene_id in sorted(gene_site[contig_id].keys()):
            gene_str = ""
            protein_str = ""
            for s,e,p in gene_site[contig_id][gene_id]:
                start_site = int(s) - 1
                end_site = int(e)
                phase = p
                cds_str = fasta_hash[contig_id].seq[start_site:end_site]
                gene_str += cds_str
                # protein of cds
                subrecord_cds = SeqRecord(cds_str, id="cds_1", description='')
                if gene_strand[contig_id][gene_id] == "-":
                    rc_seq = subrecord_cds.seq.reverse_complement()
                    subrecord_cds = SeqRecord(rc_seq, id=subrecord_cds.id, description='')
                tmp_cds_str = str(subrecord_cds.seq)[phase:]
                subrecord_cds = SeqRecord(Seq(tmp_cds_str), id=subrecord_cds.id, description='')
                protein_cds = subrecord_cds.seq.translate(table=codon)
                protein_str += str(protein_cds).rstrip("*")
            subrecord = SeqRecord(gene_str.upper(), id=contig_id + "_gene" + str(number), description=str(len(gene_str))+"_nt")
            if gene_strand[contig_id][gene_id] == "-":
                rc_seq = subrecord.seq.reverse_complement()
                subrecord = SeqRecord(rc_seq, id=subrecord.id, description=subrecord.description)
            # 核酸
            SeqIO.write(subrecord, output_fh, 'fasta')
            # 蛋白
            subrecord_prot = SeqRecord(Seq(protein_str, IUPAC.protein), id=contig_id + "_gene" + str(number), description='')
            SeqIO.write(subrecord_prot, output_fh2, 'fasta')
            #protein = subrecord.seq.translate(table=codon)
            #protein_record = SeqRecord(Seq(str(protein).rstrip("*"), IUPAC.protein), id=subrecord.id, description='')
            #SeqIO.write(protein_record, output_fh2, 'fasta')

            number += 1

output_fh2.close()
output_fh.close()
