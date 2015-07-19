# -*- coding: utf-8 -*-
import sys
import re
import csv

'Django libraries and packages'
from django.http import HttpResponse
from django.shortcuts import render, get_object_or_404
from django.conf import settings

'Biopython libraries'
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import Seq

'Graphical Libraries'
import matplotlib.pyplot as plt
import numpy as np
from units import unit

'Trie libraries and packages'
from optparse import OptionParser
from pytrie import SortedStringTrie as trie

'Import external python scripts'
import vcf
import graph_tool as gt
import graph_simple as gs
import graph_common as gc
import graph_differ as gd

'Import defined models'
from .models import Technology, Gene
from .forms import DocumentForm

'TGTCGGCCATATCGAGCCATTGAGCAGCGGCAAGGACAAGACCAGCCTGATCCTGGCCGTGCCCAACCGCGCCGGCGCCGTCTATGACATGCTGGCGCCGATGGCGGCCA'
# Create your views here.
def gene_list(request):
    genes = Gene.objects.all()
    
    return render(request, 'pages/gene_list.html', {'genes': genes})

def gene_detail(request, gene_id):
    
    gene = get_object_or_404(Gene,pk=gene_id) 
    return render(request, 'pages/gene_detail.html', {'gene': gene})

def technology_list(request, gene_id):
    
    gene = get_object_or_404(Gene,pk=gene_id) 
    technologies = gene.technologies.all()
    form = DocumentForm()
    return render(request, 'pages/technology_list.html', {'technologies': technologies, 'form': form})
    
def technology_detail(request, tech_id):
    
    technology = get_object_or_404(Technology,pk=tech_id) 
    
    return render(request, 'pages/technolgy_detail.html', {'technology': technology})

def search(request):
    
    if request.method == 'POST':
        sequence = request.POST.get('sequence')
        if len(sequence) == 0:
            return HttpResponseNotFound('Go back and please input a sequence')
        
        pk = request.POST.get('action', None)
        tech = get_object_or_404(Technology,pk=pk)
        
        'Check if user upload snp-file'
        form = DocumentForm(request.POST, request.FILES) 
        if form.is_valid():
            tech.snp_file = snp_file = request.FILES['docfile']
            tech.save()
        'If no upload nothing happens'
              
        save_folder = gs.draw_simple(tech, sequence)
        tech.barGraph = save_folder
        tech.save()
    
        return render(request, 'pages/technology_detail.html', {'technology': tech})
    
    return HttpResponseNotFound('<h1>Page not found</h1>')

def common_search(request):
    if request.method == 'POST':
        sequence = request.POST.get('sequence')
        technolgies = request.POST.getlist('checks')
        
        'Check if only 2 technologies have been selected'
        if len(technolgies) < 2 or len(technolgies) > 2:
            return HttpResponseNotFound('Go back and select only two technologies please')
        if len(sequence) == 0:
            return HttpResponseNotFound('Go back and please input a sequence')
        
        tech1 = get_object_or_404(Technology,pk=technolgies[0])
        tech2 = get_object_or_404(Technology, pk=technolgies[1])
        
        'Check if user upload snp-file'
        form = DocumentForm(request.POST, request.FILES) 
        if form.is_valid():
            tech1.snp_file = snp_file = request.FILES['docfile']
            tech2.snp_file = snp_file = request.FILES['docfile']
            tech1.save()
            tech2.save()
        'If no upload nothing happens'
        
        save_folder1 = gc.draw_common(tech1, tech2, sequence)
        save_folder2 = gd.draw_differ(tech1,tech2,sequence)
        tech1.commonGraph = save_folder1
        tech1.differGraph = save_folder2
        
        tech2.commonGraph = save_folder1
        tech2.differGraph = save_folder2
        
        tech1.compared_with = tech2.technology
        tech2.compared_with = tech1.technology
        
        tech1.save()      
        tech2.save()
        
        return render(request, 'pages/common_motif.html', {'technology': tech1})
    
    return HttpResponseNotFound('<h1>Page not found</h1>')       
    
#===============================================================================
#  Helper methods retrieve files and create immage
#===============================================================================
