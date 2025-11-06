#!/public/wangm/miniconda3/bin/python
# -*- coding: utf-8 -*-

import pathlib
import time
import json
import shutil
import subprocess
from Bio import SeqIO
from ete3 import NCBITaxa
from Bio.SeqUtils import gc_fraction
from script.function import getblast


workdir = pathlib.Path.cwd()

tmp_dir = pathlib.Path(workdir) / 'tmp'
in_dir = tmp_dir / 'fasta'
gb_dir = tmp_dir / 'gbk'



def rename(runID, infile):
    run_dir = tmp_dir / runID
    run_dir.mkdir(parents=True, exist_ok=True)

    gb_dir.mkdir(parents=True, exist_ok=True)  # ensure gb_dir exists

    i = 1
    id_dict = {}
    newIDfa = run_dir / f'{runID}_newID.fa'

    with open(newIDfa, 'w') as newfa:
        for seq_record in SeqIO.parse(infile, 'fasta'):
            if len(seq_record.id) > 15:
                realID = f'{seq_record.id[:15]}...'
            else:
                realID = seq_record.id
            contigID = f'contig_{i}'
            seq_record.id = contigID
            seqfa = str(seq_record.seq)
            newfa.write(f'>{contigID}\n{seqfa}\n')
            id_dict[contigID] = realID
            i += 1

    return id_dict

def get_time():

	return time.asctime( time.localtime(time.time()) )

def Taxonomy(runID):
    newIDfa = tmp_dir / runID / f'{runID}_newID.fa'
    report = tmp_dir / runID / f'{runID}_kraken.report'
    output = tmp_dir / runID / f'{runID}_kraken.output'
    drawout = tmp_dir / runID / 'kraken.html'

    annocmd = [
        'kraken',
        '--db', 'krakenDB',
        '--report', str(report),
        '--output', str(output),
        str(newIDfa)
    ]

    subprocess.run(annocmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Optional: drawing step (commented in original)
    # drawcmd = [
    #     '/opt/R/3.6.3/bin/Rscript',
    #     './script/sankey.R',
    #     str(report),
    #     str(drawout)
    # ]
    # subprocess.run(drawcmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    spdict = {}
    with open(output, 'r') as taxainfo:
        for line in taxainfo:
            lines = line.strip().split('\t')
            ID, taxid = lines[1], lines[2]
            if taxid == '0':
                spdict[ID] = '-'
            else:
                spname, strainame = get_ranks(taxid)
                spdict[ID] = spname

    return drawout, spdict, report

def get_ranks(taxid):
	ncbi = NCBITaxa()
	lineage = ncbi.get_lineage(taxid)   
	names = ncbi.get_taxid_translator(lineage)
	lineage2ranks = ncbi.get_rank(names)
	ranks2lineage = dict((rank,taxid) for (taxid, rank) in lineage2ranks.items())

	strainid = ''
	if 'species' in ranks2lineage:
		spid = ranks2lineage['species']
		if 'strain' in ranks2lineage:
			strainid = ranks2lineage['strain']
		else:
			strainid = ''
	elif 'genus' in ranks2lineage:
		spid = ranks2lineage['genus']
		strainid = ''
	elif 'phylum' in ranks2lineage:
		spid = ranks2lineage['phylum']
		strainid = ''

	try:
		spname = list(ncbi.get_taxid_translator([spid]).values())[0]
	except:
		spname = '-'

	strainame = ''
	if strainid:
		strainame = list(ncbi.get_taxid_translator([strainid]).values())[0]

	return spname,strainame

def getbase(runID):
    
    newIDfa = tmp_dir / runID / f'{runID}_newID.fa'
    statcmd = ['seqkit', 'stats', '-a', str(newIDfa)]

    # Run command and capture output
    result = subprocess.run(statcmd, capture_output=True, text=True, check=True)
    stats_output = result.stdout.splitlines()

    for line in stats_output:
        lines = line.strip().split()
        if lines[0] != 'file':
            lengt = lines[4]
            count = lines[3]
            n50 = lines[12]

    basefile = tmp_dir / runID / f'{runID}_info.json'
    basedict = {
        'JobID': runID,
        'Submission date': get_time(),
        'Total length': f'{lengt} bp',
        'Contig number': count,
        'Sequence N50': f'{n50} bp'
    }

    with open(basefile, 'w') as gein2:
        json.dump(basedict, gein2, indent=4)

    return basefile

def gc(fasta_file, start, end):
    
    record = SeqIO.read(fasta_file, 'fasta')
    sequence = record.seq[start-1 : end]
    gcs = f'{gc_fraction(sequence) * 100:.2f}'
    
    return gcs

def calculate_gc(fasta_file, start, end, window_size, step_size):

	record = SeqIO.read(fasta_file, "fasta")
	if start == 0:
		start = 1
	sequence = record.seq[start-1:end]

	# windows = []
	gc_contents = []
	pos = []
	j = start/1000 + 0.025
	for i in range(0, len(sequence) - window_size + 1, step_size):
		window = sequence[i:i+window_size]
		gc_content = gc_fraction(window) * 100
		gc_contents.append(gc_content)
		pos.append(round(j, 4))
		j += 0.05

	gcdict = {
		        'xData':pos,
		        'datasets':[{
		            'name':'',
		            'data':gc_contents,
		            'unit':'%',
		            'type':'line',
		            "valueDecimals": 1
		        }]
	    }

	return gcdict

def preanno(runID):
    newIDfa = tmp_dir / runID / f'{runID}_newID.fa'
    anno_fa = tmp_dir / runID / f'{runID}.faa'
    anno_gff = tmp_dir / runID / f'{runID}.gff'

    anno_cmd = [
        'prodigal',
        '-c',
        '-m',
        '-q',
        '-p', 'meta',
        '-f', 'gff',
        '-i', str(newIDfa),
        '-a', str(anno_fa),
        '-o', str(anno_gff)
    ]

    subprocess.run(anno_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def scanf(hmmlist):

	ICEcount = []
	for line in hmmlist:
		if 'MOB' in line:
			ICEcount.append('MOB')
		elif 't4cp' in line or 'tcpA' in line:
			ICEcount.append('t4cp')
		elif 'FA' in line:
			ICEcount.append('T4SS')
		elif line in ['Phage_integrase','UPF0236','Recombinase','rve','TIGR02224','TIGR02249','TIGR02225','PB001819']:
			ICEcount.append('Int')
		else:
			ICEcount.append('T4SS')
	if ICEcount.count('MOB') and ICEcount.count('t4cp') and ICEcount.count('Int') and ICEcount.count('T4SS') >= 5:
		return True
	else:
		return False

def prescan(runID):
    preanno(runID)

    anno_fa = tmp_dir / runID / f'{runID}.faa'
    scanfile = tmp_dir / runID / f'{runID}_prescan'

    scancmd = [
        './tool/hmmscan2',
        '--tblout', str(scanfile),
        './data/ICEscan.hmm',
        str(anno_fa)
    ]

    subprocess.run(scancmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    icedict = {}
    chosen = []

    with open(scanfile, 'r') as outfile:
        for line in outfile:
            if not line.startswith('#'):
                lines = line.strip().split()
                if lines[2] in icedict:
                    continue
                id_parts = lines[2].split('_')
                key = '_'.join(id_parts[0:2])
                if float(lines[4]) < 0.00001:
                    if key in icedict:
                        icedict[key].append(lines[0])
                    else:
                        icedict[key] = [lines[0]]

    for k, v in icedict.items():
        if scanf(v):
            chosen.append(k)

    return chosen

def prokkanno(runID, input_file, threads):
    # annotate the fasta sequence by prokka
    gb_out = gb_dir  # gb_dir is already a pathlib.Path object

    cmd = [
        'prokka',
        str(input_file),
        '--force',
        '--cpus', str(threads),
        '--cdsrnaolap',
        '--prefix', runID,
        '-o', str(gb_out)
    ]

    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
def getgff(runID):
    #open the gff file
    gffile = gb_dir / f'{runID}.gff'
    trnadict = {}
    posdict = {}
    totalnum = 0
    header = ''
    # two dicts for trna pos and position of every CDS.
    with open(gffile,'rt') as gffin:
        for line in gffin.readlines():
            if 'ID=' in line:
                lines = line.strip().split('\t')
                ids = lines[8].split(';')[0].split('=')[1]
                if not header:
                    header = ids.split('_')[0]
                product = lines[8].split('product=')[1]
                pos = [lines[3],lines[4],lines[6],product]
                if lines[2] == 'tRNA' or lines[2] == 'tmRNA':
                    trnadict[ids] = product
                posdict[ids] = pos
                totalnum += 1
    
    return trnadict, posdict, header, totalnum

def getnum(ID):

	return int(ID.split('_')[1].lstrip('0'))

def find_max_distance(numbers):
	max_distance = -1
	max_distance_index = -1

	for i in range(len(numbers) - 1):
		distance = abs(numbers[i] - numbers[i+1])
		if distance > max_distance:
			max_distance = distance
			max_distance_index = i

	if max_distance_index == -1:
		return None

	return numbers[max_distance_index], numbers[max_distance_index + 1]

def pos_tag(pos,posdict,ICE,final,totalnum,dirtag):

	tICE = ICE
	tfinal = final
	for k,v in posdict.items():
		vstart,vend = int(v[0]),int(v[1])
		if int(pos) <= vend:
			if dirtag == 's':
				tICE = getnum(k)
				tfinal = max(1, tICE - 5)
			else:
				if vstart > int(pos):
					tICE = getnum(k) - 1 
				else:
					tICE = getnum(k)
				tfinal = min(totalnum, tICE + 5)
			break
	return tICE, tfinal

def merge_tRNA(runID,ICEdict,DRlist):

	trnadict,posdict,header,totalnum = getgff(runID)
	fICE = getnum(next(iter(ICEdict)))
	eICE = getnum(list(ICEdict.keys())[-1])

	nfICEnum = max(1, fICE - 5)
	neICEnum = min(totalnum, eICE + 5)

	ICEtagnum = [nfICEnum,neICEnum]
	trnalist = []
	for key,value in trnadict.items():
		if nfICEnum <= getnum(key) <= neICEnum:
			ICEtagnum.append(getnum(key))
			trnalist.append(value)

	ICEtagnum.sort()
	finalstart,finalend = find_max_distance(ICEtagnum)

	myDR1 = posdict[zill(header,fICE)][0]
	myDR2 = ''
	myDR3 = ''
	myDR4 = posdict[zill(header,eICE)][1]

	if trnalist:
		if finalstart == nfICEnum:
			eICE = finalend
			finalend = min(totalnum, finalend + 5)
			myDR4 = posdict[zill(header,eICE)][1]					
			for line in DRlist:
				DRs = line.split(',')
				if int(DRs[3]) - int(DRs[0]) > 500000:
					continue
				if int(DRs[3]) - int(DRs[0]) < 5000:
					continue					
				if int(posdict[zill(header,eICE)][0]) < int(DRs[3]) < int(posdict[zill(header,eICE)][1]):
					checktrna = 0
					for key,value in trnadict.items():
						if int(DRs[0]) <= value[0] <= int(DRs[3]) and int(DRs[0]) <= value[1] <= int(DRs[3]):
							checktrna += 1
					if checktrna >= 2:
						break

					fICE,finalstart = pos_tag(DRs[0],posdict,fICE,finalstart,totalnum,'s')
					myDR1 = DRs[0]
					myDR2 = DRs[1]
					myDR3 = DRs[2]
					myDR4 = DRs[3]
					break

		elif finalend == neICEnum:
			fICE = finalstart
			finalstart =  max(1, finalstart - 5)
			myDR1 = posdict[zill(header,fICE)][0]
			for line in DRlist:
				DRs = line.split('|')
				if int(DRs[3]) - int(DRs[0]) > 500000:
					continue	
				if int(DRs[3]) - int(DRs[0]) < 5000:
					continue									
				if int(posdict[zill(header,fICE)][0]) < int(DRs[0]) < int(posdict[zill(header,fICE)][1]):
					checktrna = 0
					for key,value in trnadict.items():
						if int(DRs[0]) <= value[0] <= int(DRs[3]) and int(DRs[0]) <= value[1] <= int(DRs[3]):
							checktrna += 1
					if checktrna >= 2:
						break
					eICE,finalend = pos_tag(DRs[3],posdict,eICE,finalend,totalnum,'e')
					myDR1 = DRs[0]
					myDR2 = DRs[1]
					myDR3 = DRs[2]
					myDR4 = DRs[3]									
					break

	return myDR1,myDR2,myDR3,myDR4,fICE,eICE,finalstart,finalend,posdict,header,trnalist

def get_DR(runID, infile):
    
    DRindex = tmp_dir / runID / f'{runID}_DR'
    # fasta and gbk can both be the input of mkvtree
    mkvtree_cmd = ['mkvtree', '-db', str(infile), '-indexname', str(DRindex), '-dna', '-pl', '-allout']
    vmatch_cmd = ['vmatch', '-l', '15', str(DRindex)]
    subprocess.run(mkvtree_cmd, check = True)
    process = subprocess.run(vmatch_cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out = process.stdout.decode()
    # DRout : 22    0 1228868   D    22    0 1736928   0    8.87e-02     44   100.00
    # length_left, seq_num_i, post_left, strand, length_right, seq_num_j, pos_right, distance, E-value, score, identity
    DRlist = []
    for line in out.strip().split('\n'):
        lines = line.strip().split('\t')
        if not line.startswith('#'):
            DR = [str(int(lines[2]) + 1), str(int(lines[2]) + int(lines[0])), str(int(lines[6]) + 1), 
      str(int(lines[6])+int(lines[4]))]   #left_start, left_end, right_start, right_end
            DRlist.append(','.join(DR))

    return DRlist

def get_ICE(runID, infile):
    
    ICE_dir = tmp_dir / runID / f'{runID}_ICE'
    ICE_res = ICE_dir / 'all_systems.tsv'
    # remove the file folder if it is already exists
    if ICE_dir.exists():
        shutil.rmtree(ICE_dir)
    
    anno_fa = gb_dir / f'{runID}.faa'
    # run macsyfinder
    ICE_cmd = ['macsyfinder', '--db-type', 'ordered_replicon', '--models-dir', 
               './data/macsydata/', '--models', 'ICEscan', 'all', '--replicon-topology', 'linear', 
               '--coverage-profile', '0.3', '--sequence-db', str(anno_fa), '-o', str(ICE_dir)]
    subprocess.run(ICE_cmd, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    
    with open(ICE_res, 'r') as ICEin:
        ICEdict = {}
        infodict = {}
        
        for line in ICEin.readlines():
            if 'Chromosome' in line:
                lines = line.strip().split('\t')
                if 'UserReplicon_IME' not in lines[5]:
                    gbname = lines[1]
                    tags = get_feat(lines[2])
                    mpf = lines[4].split('/')[-1].split('_')[1]
                    
                    if 'Relaxase@' in tags:
                        mob = tags.split('@')[1]
                    else:
                        mob = ''
                    ICEtag = 'ICE'+lines[5].split('_')[-1]
                    
                    ICEdict.setdefault(ICEtag,{})[gbname]=tags
                    if ICEtag not in infodict:
                        infodict[ICEtag] = {'mob': [], 'mpf': []}
                    if mob not in infodict[ICEtag]['mob']:
                        if mob:
                            infodict[ICEtag]['mob'].append(mob)
                    if mpf not in infodict[ICEtag]['mpf']:
                        infodict[ICEtag]['mpf'].append(mpf)	
    
    dictICE = {}
    posdict = {}
    trnalist = []
    header = ''
    DRlist = get_DR(runID,infile)
    for key,value in ICEdict.items():
        myDR1,myDR2,myDR3,myDR4,fICE,eICE,finalstart,finalend,posdict,header,trnalist = merge_tRNA(runID,value,DRlist)
        dictICE[key] = [myDR1,myDR2,myDR3,myDR4,fICE,eICE,finalstart,finalend]
    
    return dictICE,ICEdict,posdict,header,trnalist,infodict

def args(runID):

	return getblast(runID)

def zill(header,num):

	return header+ '_' +str(num).zfill(5)

def get_args(argdict,vfdict,isdict,dfdict,metaldict,popdict,symdict,gene,feature,product):

	feature = [feature]
	product = [product]

	if gene in argdict:
		feature.append('AR')
		product.append(argdict[gene])
	if gene in vfdict:
		feature.append('VF')
		product.append(vfdict[gene])
	if gene in isdict:
		feature.append('IS')
		product.append(isdict[gene])
	if gene in dfdict:
		feature.append('Defense')
		product.append(dfdict[gene])

	if gene in metaldict:
		feature.append('Metal')
		product.append(metaldict[gene])
	if gene in popdict:
		feature.append('Degradation')
		product.append(popdict[gene])
	if gene in symdict:
		feature.append('Symbiosis')
		product.append(symdict[gene])

	feature = '; '.join(list(filter(None, feature)))
	product = '; '.join(list(filter(None, product)))

	return feature,product

def oritseq(runID, regi, infile, start, end):
    
    oritseq = '-'
    fafile = tmp_dir / runID / f'{regi}_fororit.fa'
    
    with open(fafile, 'w') as orif:
        seq = getfa(infile, start, end)
        orif.write('>fororit\n')
        orif.write(seq)

    oriT_Database = pathlib.Path(workdir) / 'data' / 'oriT_db'
    blastn_out = tmp_dir / runID / f'{regi}_oriTout'

    blast_cmd = [
        'blastn',
        '-db', str(oriT_Database),
        '-query', str(fafile),
        '-evalue', '0.01',
        '-word_size', '11',
        '-outfmt', '6 std qlen slen',
        '-num_alignments', '1',
        '-out', str(blastn_out)
    ]

    subprocess.run(blast_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    with open(blastn_out, 'r') as oritout:
        for line in oritout:
            lines = line.strip().split()
            if lines[0]:
                matchl = int(lines[3])
                slen = int(lines[13])
                ident = float(lines[2])
                hvalue = (matchl / slen) * ident
                if hvalue > 0.49:
                    oritseq = getfa(fafile, str(int(lines[6])-1), lines[7])
                    break

    return oritseq

def get_feat(feat):

	featuredict = {
		'Phage_integrase':'Integrase','UPF0236':'Integrase',
        'Recombinase':'Integrase','rve':'Integrase',
        'TIGR02224':'Integrase','TIGR02249':'Integrase',
        'TIGR02225':'Integrase','PB001819':'Integrase'
	}
	if feat in featuredict:
		return 'Integrase@'+feat
	elif 'T4SS_MOB' in feat:
		tag = feat.split('_')[1]
		return 'Relaxase@'+tag
	elif 't4cp' in feat:
		tag = feat.split('_')[1]
		return 'T4CP@' + tag
	elif 'tcpA' in feat:
		tag = feat.split('_')[1]
		return 'T4CP@' + tag
	else:
		return 'T4SS@'+feat.replace('T4SS_','')

def getcolor(feature,product):

	coldict = {'DR':'black','Gene':'#C0C0C0',
		   'Hyp':'#DCDCDC','Integrase':'blue',
		   'Transposase':'yellow','T4SS':'lightpink','T4CP':'orange',
		   'Relaxase':'brown','AR':'red','tRNA':'black',
		   'Flank':'gray','VF':'#ba8448','Defense':'#00B050',
		   'Metal':'#03A89E','Degradation':'#640B0F','Symbiosis':'#FFFFCD'
	}

	namedict = {'Hyp':'Hypothetical protein','Gene':'Other gene',
		    'AR':'Antibiotic resistance gene',
		    'VF':'Virulence factor','Metal':'Metal resistance',
		    'Flank':'Flank region','Defense':'Defense system',
		    'Transposase':'Transposase','Relaxase':'Relaxase',
		    'T4CP':'T4CP','T4SS':'T4SS','Integrase':'Integrase',
		    'Degradation':'Degradation','Symbiosis':'Symbiosis'
	}

	if 'Integrase' in feature:
		feature = 'Integrase'
	elif 'T4SS' in feature:
		feature = 'T4SS'
	elif 'T4CP' in feature:
		feature = 'T4CP'
	elif 'Relaxase' in feature:
		feature = 'Relaxase'
	elif 'IS' in feature:
		feature = 'Transposase'
	elif 'VF' in feature:
		feature = 'VF'
	elif 'AR' in feature:
		feature = 'AR'
	elif 'Defense' in feature:
		feature = 'Defense'
	elif 'Metal' in feature:
		feature = 'Metal'
	elif 'Degradation' in feature:
		feature = 'Degradation'
	elif 'Symbiosis' in feature:
		feature = 'Symbiosis'

	elif feature == 'Flank':
		feature == 'Flank'
	elif feature == '':
		if product == 'hypothetical protein':
			feature = 'Hyp'
		else:
			feature = 'Gene'
	else:
		feature = 'Gene'

	return coldict[feature], namedict[feature]

def gstrand(instra):

	strands = {'+' : 1, '-' : -1}
	return strands[instra]

def getfa(infile,s,e):

	seq_record = SeqIO.read(infile, "fasta")
	sequence = seq_record.seq[int(s):int(e)]
	return str(sequence)

def get_map(sprunID, spdict, id_dict):
    
    final_dir = pathlib.Path(workdir) / 'tmp' / sprunID / 'result'
    js_dir = final_dir / 'js'
    js_dir.mkdir(parents=True, exist_ok=True)

    gcmap = pathlib.Path(workdir) / 'script' / 'js' / 'gcmap.js'
    viewfile = pathlib.Path(workdir) / 'script' / 'js' / 'view.html'

    fasta_file = tmp_dir / sprunID / f'{sprunID}.fa'
    dictICE, ICEdict, posdict, header, trnalist, infodict = get_ICE(sprunID, fasta_file)
    argdict, vfdict, isdict, dfdict, metaldict, popdict, symdict = args(sprunID)

    ICEss = {}
    for key, value in dictICE.items():
        genelist = []
        regi = f'{sprunID}_{key}'
        regijs = f'contig_{sprunID.split("_contig_", 1)[-1]}_{key}'

        genefile = final_dir / f'{regi}_gene.json'
        infofile = final_dir / f'{regi}_info.json'
        gcjson = js_dir / f'{regijs}_gc.js'
        mapfile = js_dir / f'{regijs}.js'
        htmlfile = final_dir / f'{regi}.html'

        [myDR1, myDR2, myDR3, myDR4, fICE, eICE, finalstart, finalend] = value
        
        start = finalstart
        while start < fICE:
            gene = zill(header, start)
            s, e, strand, pro = posdict[gene]
            pos = s + '..' + e + ' [' + strand + '], ' + str(int(e) - int(s) + 1)
        
            feature = 'Flank'
            product = pro
            feature, product = get_args(argdict, vfdict, isdict, dfdict, metaldict, popdict, symdict, gene, feature, product)
            if 'hypothetical protein;' in product:
                product = product.replace('hypothetical protein;', '')
        
            start += 1
            content = {
                'gene': gene,
                'pos': pos,
                'prod': product,
                'featu': feature
            }
            genelist.append(content)
        
        mov = fICE
        while mov <= eICE:
            gene = zill(header, mov)
            s, e, strand, pro = posdict[gene]
            pos = s + '..' + e + ' [' + strand + '], ' + str(int(e) - int(s) + 1)
        
            if gene in ICEdict[key]:
                [feature, pro11] = ICEdict[key][gene].split('@')
            else:
                feature, pro11 = '', ''
        
            if pro11:
                if pro == 'hypothetical protein':
                    product = pro11
                else:
                    product = pro + ', ' + pro11
            else:
                product = pro
        
            feature, product = get_args(argdict, vfdict, isdict, dfdict, metaldict, popdict, symdict, gene, feature, product)
            mov += 1
            content = {
                'gene': gene,
                'pos': pos,
                'prod': product,
                'featu': feature
            }
            genelist.append(content)
        
        while mov <= finalend:
            gene = zill(header, mov)
            s, e, strand, pro = posdict[gene]
            pos = s + '..' + e + ' [' + strand + '], ' + str(int(e) - int(s) + 1)
        
            feature = 'Flank'
            product = pro
            feature, product = get_args(argdict, vfdict, isdict, dfdict, metaldict, popdict, symdict, gene, feature, product)
            if 'hypothetical protein;' in product:
                product = product.replace('hypothetical protein;', '')
        
            mov += 1
            content = {
                'gene': gene,
                'pos': pos,
                'prod': product,
                'featu': feature
            }
            genelist.append(content)
        
        with open(genefile, 'w') as gene_file:
            json.dump(genelist, gene_file, indent=4)
        
        contigID = sprunID.split('_', 1)[1]
        
        sgene = zill(header, fICE)
        egene = zill(header, eICE)
        s1, e1, strand1, pro1 = posdict[sgene]
        s2, e2, strand2, pro2 = posdict[egene]
        if myDR1 == '0':
            myDR1 = '1'
        
        ICEss[regi] = ','.join([myDR1, myDR4, str(fICE), str(eICE)])
        
        host = spdict[contigID]
        gcc = gc(fasta_file, int(myDR1), int(myDR4))
        source = id_dict[contigID]
        
        if myDR2:
            DR1 = getfa(fasta_file, myDR1, myDR2)
            DR2 = getfa(fasta_file, myDR3, myDR4)
            DRw = 'attL:' + myDR1 + '..' + myDR2 + '(' + DR1 + ')  ' + 'attR:' + myDR3 + '..' + myDR4 + '(' + DR2 + ')'
        else:
            DRw = '-'
        
        oritseqs = oritseq(sprunID, regi, fasta_file, myDR1, myDR4)
        # oritdesc = "<br>".join([oritseqs[i:i+63] for i in range(0, len(oritseqs), 63)])
        
        ICEinfo = {
            'Contig source': source,
            'Host Strain': host,
            'GC Content (%)': gcc,
            'Length (bp)': str(int(e2) - int(s1) + 1),
            'oriT seq': oritseqs,
            'DRs': DRw,
            'Relaxase Type': ','.join(infodict[key]['mob']),
            'Mating pair formation systems': ','.join(infodict[key]['mpf']),
            'Close to tRNA': ','.join(trnalist)
        }
        with open(infofile, 'w') as info_file:
            json.dump(ICEinfo, info_file, indent=4)
        
        i = 1
        mapzlist = []
        mapflist = []
        for gene in genelist:
            color, name = getcolor(gene['featu'], gene['prod'])
            start = gene['pos'].split(' ')[0].split('..')[0]
            end = gene['pos'].split(' ')[0].split('..')[1]
            strand = gstrand(gene['pos'].split('[')[1].split(']')[0])
            product = gene['prod']
        
            if product == '':
                product = 'hypothetical protein'
        
            anno = {
                'start': start,
                'end': end,
                'strand': strand,
                'locus_tag': 'M' + str(i),
                'type': 'others',
                'color': color,
                'description': 'Location: ' + gene['pos'].split(' ')[0] + ' (' + gene['pos'].split(' ')[2] + ' bp)<br>Type: ' + name + '<br>Detail: ' + product
            }
            if strand == 1:
                mapzlist.append(anno)
            else:
                mapflist.append(anno)
            i += 1
        
        head = 'var borders = [];\nvar tta_codons = [];\nvar orfs ='
        s = genelist[0]['pos'].split(' ')[0].split('..')[0]
        e = genelist[-1]['pos'].split(' ')[0].split('..')[1]
        
        gcdict = calculate_gc(fasta_file, int(s), int(e), 500, 50)
        with open(gcmap, 'r') as original_file:
            original_content = original_file.read()
        with open(gcjson, 'w') as gein2:
            gein2t = 'var jsonData = ' + str(gcdict) + ';'
            gein2.write(gein2t)
            gein2.write(original_content)
        
        # with open(gcjson,'w') as gein2:
        #     json.dump(gcdict, gein2, indent=4)
        
        maps = str(mapzlist) + ';\nvar orfs2 =' + str(mapflist) + ';\nvar clusterf2 = { start: ' + s + ', end: ' + \
               e + ', idx: 1, orfs: orfs, borders: borders, tta_codons:tta_codons, \
               label: \'\', unordered: true };\nvar clusterr2 = { start: ' + s + ', end: ' + \
               e + ', idx: 2, orfs: orfs2, borders: borders, tta_codons:tta_codons, \
               label: \'\', unordered: true };\nsvgene.drawClusters("' + regijs + '", [clusterf2, clusterr2], 50, 920);'
        with open(mapfile, 'w') as map_file:
            map_file.write(head + maps)
        
        with open(viewfile, 'r') as file:
            file_content = file.read()
        new_content = file_content.replace('XXXX', regijs)
        with open(htmlfile, 'w') as file:
            file.write(new_content)
    
    return ICEss

def copy_files(source_dir, destination_dir):
    
    source_dir = pathlib.Path(source_dir)
    destination_dir = pathlib.Path(destination_dir)

    if source_dir.is_file():
        shutil.copy(source_dir, destination_dir)
    elif source_dir.is_dir():
        destination_dir.mkdir(parents=True, exist_ok=True)
        for item in source_dir.iterdir():
            source_item = item
            destination_item = destination_dir / item.name
            if source_item.is_dir():
                copy_files(source_item, destination_item)
            else:
                shutil.copy2(source_item, destination_item)

def delete_folders_starting_with_keyword(dir, keyword):
    
    for folder in dir.rglob(f'{keyword}*'):
        if folder.is_dir():
            shutil.rmtree(folder)

def getfasta(runID,resultdir,id_dict,key,s,e,stag,etag):
    
    fafile = tmp_dir / runID / f'{runID}.fa'
    faafile = tmp_dir / 'gbk' / f'{runID}.faa'
    
    outfa = resultdir / f'{key}.fa'
    outfaa = resultdir / f'{key}.faa'
    
    seq_record = SeqIO.read(fafile, "fasta")
    with open(outfa, "w") as output_handle1:
        sequence = seq_record.seq[int(s)-1:int(e)]
        ID = '_'.join(seq_record.id.split('_')[-2:])
        seq_record.description = ''
        seq_record.seq = sequence
        seq_record.id = id_dict[ID] + ' ' + s +'-'+e
        SeqIO.write(seq_record, output_handle1, "fasta")
        
    faa_records = SeqIO.parse(faafile, "fasta")
    with open(outfaa, "w") as output_handle2:
        for faa_record in faa_records:
            seq_id = getnum(faa_record.id)
            if int(stag) <= seq_id <= int(etag):
                SeqIO.write(faa_record, output_handle2, "fasta")

def _meta(runID, input_file, threads):
    
    resultdir = workdir / 'result' / runID
    resultdir.mkdir(parents=True, exist_ok=True)
    
    jsback = workdir / 'script' / 'js'
    
    id_dict = rename(runID, input_file)
    chosenfa = prescan(runID)
    drawout,spdict,report = Taxonomy(runID)
#	spdict = {}
    basefile = getbase(runID)
#	copy_files(drawout, resultdir)
    copy_files(basefile, resultdir)
    copy_files(report, resultdir)
    
    ICEsum = tmp_dir / runID / f'{runID}_ICEsum.json'
    newIDfa = tmp_dir / runID / f'{runID}_newID.fa'
    
    i = 1 
    ICEsumlist = []
    for seq_record in SeqIO.parse(newIDfa, "fasta"):
        if seq_record.id in chosenfa:
            sprunID = runID + '_' + seq_record.id
            seqfa = str(seq_record.seq)
            newfolder = tmp_dir / sprunID
            newfolder.mkdir(parents=True, exist_ok=True)
            
            spfa = newfolder / f'{sprunID}.fa'
            with open(spfa, 'w') as outfa:
                outfa.write(f'>{sprunID}\n{seqfa}\n')
            
            final_dir = newfolder / 'result'
            final_dir.mkdir(parents=True, exist_ok=True)
            
            prokkanno(sprunID,spfa)
            ICEss = get_map(sprunID,spdict,id_dict)
            copy_files(final_dir, resultdir)
            
            if ICEss:
                for key,value in ICEss.items():
                    [s,e,stag,etag] = value.split(',')
                    lengt = int(e) - int(s) + 1
                    ICEs = {
						'id' : str(i),
				        'seqid': id_dict[seq_record.id],
				        'species':spdict[seq_record.id],
#						'species':'',
				        'location': s+'..'+e,
				        'length': lengt,
				        'detail': key
				    }
                    ICEsumlist.append(ICEs)
                    getfasta(sprunID,resultdir,id_dict,key,s,e,stag,etag)
                    i += 1 
    
    with open(ICEsum,'w') as ice_file:
        json.dump(ICEsumlist, ice_file, indent=4)
    copy_files(ICEsum, resultdir)
    jsdir = resultdir / 'js'
    copy_files(jsback, jsdir)
    delete_folders_starting_with_keyword(tmp_dir, runID)
