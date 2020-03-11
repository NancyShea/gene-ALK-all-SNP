#！coding:utf-8
#author:XieCuinanCAU
#This script can create file which contain every single nuecleotide\
#change and it's standard mutation ID like in clinvar.2020.3.11.15:14. 


#remember to change 5 places each time 

with open('/Users/Xiecuinan/Desktop/ALK-g38-NP_004295.2.txt') as f:#change1.
    ref = f.read()
global seq
seq = ref.replace('\n','')
global seqList
seqList = list(seq)
global num
num = len(seqList)
global gene 
gene = "ALK"#change2. gene name,like'"GENE"'
global NM 
NM =  "NM_004304.5" #change3.depends on gene,must like'"NM_004304.5(ALK):c.*A>A(p.Met1Met)"'
global locID
locID = 928#change4. geneID,like'928'#fron nm_nnnnnn File,CDS location
global Chr
Chr = "2"#change5.chromesome number
global codon_dict
codon_dict = {
'TTT':'Phe','TCT':'Ser','TAT':'Tyr','TGT':'Cys',
'TTC':'Phe','TCC':'Ser','TAC':'Tyr','TGC':'Cys',
'TTA':'Leu','TCA':'Ser','TAA':'STOP','TGA':'STOP',
'TTG':'Leu','TCG':'Ser','TAG':'STOP','TGG':'Trp',
'CTT':'Leu','CCT':'Pro','CAT':'His','CGT':'Arg',
'CTC':'Leu','CCC':'Pro','CAC':'His','CGC':'Arg',
'CTA':'Leu','CCA':'Pro','CAA':'His','CGA':'Arg',
'CTG':'Leu','CCG':'Pro','CAG':'His','CGG':'Arg',
'ATT':'Ile','ACT':'Thr','AAT':'Asn','AGT':'Ser',
'ATC':'Ile','ACC':'Thr','AAC':'Asn','AGC':'Ser',
'ATA':'Ile','ACA':'Thr','AAA':'Lys','AGA':'Arg',
'ATG':'Met','ACG':'Thr','AAG':'Lys','AGG':'Arg',
'GTT':'Val','GCT':'Ala','GAT':'Asp','GGT':'Gly',
'GTC':'Val','GCC':'Ala','GAC':'Asp','GGC':'Gly',
'GTA':'Val','GCA':'Ala','GAA':'Glu','GGA':'Gly',
'GTG':'Val','GCG':'Ala','GAG':'Glu','GGG':'Gly',
}
#prepare done!

# ALK-g38:Chr:2:
#29193224..29193922,29196770..29196860,
#29197542..29197676,#29207171..29207272,

#29209786..29209878,29213984..29214081,
#29220706..29220835,#29222344..29222408, 
# #29222517..29222607,#29223342..29223528,

#29225461..29225565,#29226922..29227074,
#29227574..29227672,#29228884..29229066,

#29232304..29232448,#29233565..29233696,
#29239680..29239830,#29251105..29251267,
# #29275099..29275227,#29275402..29275496,


#29296888..29297057,29318304..29318404,
#29320751..29320882,#29328350..29328481,

#29383732..29383859,#29531915..29532116,
#29694850..29695014,29717578..29717697,
# #29919993..29920659
def get_number(seq):#GRChr37,copy from g38 and plus a number.#another way is use ncbi remap tool 
	result = []
	g = []
	x = []
	for i in range(29193224,29193922+1,1):
		g1 = str(i)
		g.append(g1)
	for i in range(29196770,29196860+1,1):
		g2 = str(i)
		g.append(g2)
	for i in range(29197542,29197676+1,1):
		g3 = str(i)
		g.append(g3)
	for i in range(29207171,29207272+1,1):
		g4 = str(i)
		g.append(g4)
	for i in range(29209786,29209878+1,1):
		g5 = str(i)
		g.append(g5)
	for i in range(29213984,29214081+1,1):
		g6 = str(i)
		g.append(g6)
	for i in range(29220706,29220835+1,1):
		g7 = str(i)
		g.append(g7)
	for i in range(29222344,29222408+1,1):
		g8 = str(i)
		g.append(g8)
	for i in range(29222517,29222607+1,1):
		g9 = str(i)
		g.append(g9)
	for i in range(29223342,29223528+1,1):
		g10 = str(i)
		g.append(g10)
	for i in range(29225461,29225565+1,1):
		g11 = str(i)
		g.append(g11)
	for i in range(29226922,29227074+1,1):
		g12 = str(i)
		g.append(g12)
	for i in range(29227574,29227672+1,1):
		g13 = str(i)
		g.append(g13)
	for i in range(29228884,29229066+1,1):
		g14 = str(i)
		g.append(g14)
	for i in range(29232304,29232448+1,1):
		g15 = str(i)
		g.append(g15)
	for i in range(29233565,29233696+1,1):
		g16 = str(i)
		g.append(g16)
	for i in range(29239680,29239830+1,1):
		g17 = str(i)
		g.append(g17)
	for i in range(29251105,29251267+1,1):
		g18 = str(i)
		g.append(g18)
	for i in range(29275099,29275227+1,1):
		g19 = str(i)
		g.append(g19)
	for i in range(29275402,29275496+1,1):
		g20 = str(i)
		g.append(g20)
	for i in range(29296888,29297057+1,1):
		g21 = str(i)
		g.append(g21)
	for i in range(29318304,29318404+1,1):
		g22 = str(i)
		g.append(g22)
	for i in range(29320751,29320882+1,1):
		g23 = str(i)
		g.append(g23)
	for i in range(29328350,29328481+1,1):
		g24 = str(i)
		g.append(g24)
	for i in range(29383732,29383859+1,1):
		g25 = str(i)
		g.append(g25)
	for i in range(29531915,29532116+1,1):
		g26 = str(i)
		g.append(g26)
	for i in range(29694850,29695014+1,1):
		g27 = str(i)
		g.append(g27)
	for i in range(29717578,29717697+1,1):
		g28 = str(i)
		g.append(g28)
	for i in range(29919993,29920659+1,1):
		g29 = str(i)
		g.append(g29)
#NM_004304.5(ALK):c.4796C>A (p.Pro1599His)
#GRCh37:Chr2:29416157
#GRCh38:Chr2:29193291
	for i in range(0,len(g)):
		x = int(g[i]) + 222866# -726649
		for j in range(0,4):
			gg37 = str(x)
			result.append(gg37)
	return result
		
def cli_number(seq):#GRChr38,in refseq file download from ncbi-genbank
	result = []
	g = []
	for i in range(29193224,29193922+1,1):
		g1 = str(i)
		g.append(g1)
	for i in range(29196770,29196860+1,1):
		g2 = str(i)
		g.append(g2)
	for i in range(29197542,29197676+1,1):
		g3 = str(i)
		g.append(g3)
	for i in range(29207171,29207272+1,1):
		g4 = str(i)
		g.append(g4)
	for i in range(29209786,29209878+1,1):
		g5 = str(i)
		g.append(g5)
	for i in range(29213984,29214081+1,1):
		g6 = str(i)
		g.append(g6)
	for i in range(29220706,29220835+1,1):
		g7 = str(i)
		g.append(g7)
	for i in range(29222344,29222408+1,1):
		g8 = str(i)
		g.append(g8)
	for i in range(29222517,29222607+1,1):
		g9 = str(i)
		g.append(g9)
	for i in range(29223342,29223528+1,1):
		g10 = str(i)
		g.append(g10)
	for i in range(29225461,29225565+1,1):
		g11 = str(i)
		g.append(g11)
	for i in range(29226922,29227074+1,1):
		g12 = str(i)
		g.append(g12)
	for i in range(29227574,29227672+1,1):
		g13 = str(i)
		g.append(g13)
	for i in range(29228884,29229066+1,1):
		g14 = str(i)
		g.append(g14)
	for i in range(29232304,29232448+1,1):
		g15 = str(i)
		g.append(g15)
	for i in range(29233565,29233696+1,1):
		g16 = str(i)
		g.append(g16)
	for i in range(29239680,29239830+1,1):
		g17 = str(i)
		g.append(g17)
	for i in range(29251105,29251267+1,1):
		g18 = str(i)
		g.append(g18)
	for i in range(29275099,29275227+1,1):
		g19 = str(i)
		g.append(g19)
	for i in range(29275402,29275496+1,1):
		g20 = str(i)
		g.append(g20)
	for i in range(29296888,29297057+1,1):
		g21 = str(i)
		g.append(g21)
	for i in range(29318304,29318404+1,1):
		g22 = str(i)
		g.append(g22)
	for i in range(29320751,29320882+1,1):
		g23 = str(i)
		g.append(g23)
	for i in range(29328350,29328481+1,1):
		g24 = str(i)
		g.append(g24)
	for i in range(29383732,29383859+1,1):
		g25 = str(i)
		g.append(g25)
	for i in range(29531915,29532116+1,1):
		g26 = str(i)
		g.append(g26)
	for i in range(29694850,29695014+1,1):
		g27 = str(i)
		g.append(g27)
	for i in range(29717578,29717697+1,1):
		g28 = str(i)
		g.append(g28)
	for i in range(29919993,29920659+1,1):
		g29 = str(i)
		g.append(g29)
	
	for i in range(0,len(g)):
		x = int(g[i])
		cli= str(x)
		for j in range(0,4):
			result.append(cli)
	return result


def loc_number(seq):
    result = []
    global c
    c = []
    for i in range(0,num):
        for j in range(0,4):
            x = locID + i ############locID from mRNAseq. 
			#Useless "x". Just cuz asked by boss.
			#"x" is nucleotide number in mRNA,928,929...end.
            cn =1+i	#"cn" is AAnumber,1,2,...end
            c.append(str(cn))
            idnm = str(x)
            result.append(idnm)
    return result

def aanumber(seq):
    result = []
    for i in range(0,len(nnn)):
        for j in range(0,12):###########
        	a = i + 1
        	st = str(a)
        	result.append(st)
    return result 

def duplicate(seq):
    result = []
    for i in range (0,num):
        for j in range(0,4):
            result.append(seqList[i])
    return result
def seperate(seq):
    result = []
    global nnn
    nnn = []
    for i in range (0,num,3):
        codon = seq[i:i+3]
        nnn.append(codon)
        #print(nnn)
        for j in range(0,12):
            result.append(seq[i:i+3])
    return result

def translate(seq):
    result = []
    for i in range(0,len(nnn)):
        aa = codon_dict.get (nnn[i],'X')
        for j in range(0,12):
            result.append(aa)
    return  result

def mutation(seq):
    result = []
    global atcg
    atcg = []
    #for i in range(0,len(codon),12):
        #three = codon[i]
    for i in range(0,len(nnn)):#range(0,len(codon),12),sanme with range(0,len(nnn))
        three = nnn[i]
        for j in range(0,len(three)):
            tar = list(three)
            tar[j] = 'A'
            tar1 = ''.join(tar)#join(tar),not join(tar[i])
            
            tar = list(three)#must repeat this words
            tar[j] = 'T'
            tar2 = ''.join(tar)            
            
            tar = list(three)
            tar[j] = 'C'
            tar3 = ''.join(tar)            

            tar = list(three)
            tar[j] = 'G'
            tar4 = ''.join(tar)

            result.append(tar1)
            result.append(tar2)
            result.append(tar3)
            result.append(tar4)                               
    return result

def trans_mut(seq):
    result = []
    for i in range(0,len(new_codon)):
        aa = codon_dict.get (new_codon[i],'X')
        result.append(aa)
    return  result

def last_column(seq):
    
    atcg = "ATCG"*num
    global ATCG
    ATCG = list(atcg)
    result = []
    for i in range(0,len(nt)):
        one = NM# idmn need something find in refseq.txt,used in "NM_XXX"
        two = gene
        three = nt[i]
        four = str(ATCG[i])
        five = AA[i]
        six = str(AAnum[i])
        seven = new_aa[i]
        zuobiao2 = c[i]
        shit = "{}({})".format(one,two)
        shit = shit + ":c.{}".format(zuobiao2)
        shit = shit + "{}>{}(p.{}{}{})".format(three,four,five,six,seven)
        result.append(shit)
    return result

def main():
#ready to use functions!
	
    global idnm
    idnm = loc_number(seq)
    global nt
    nt = duplicate(seq)
    global cod
    cod = seperate(seq)
    global AA
    AA = translate(seq)
    global new_codon
    new_codon = mutation(seq)
    global AAnum
    AAnum = aanumber(seq)
    global new_aa
    
    new_aa = trans_mut(seq)
    shit = last_column(seq)
    g38 = cli_number(seq)
    g38.reverse()
    g37 = get_number(seq)
    g37.reverse()
	
	#print(len(idnm),len(nt),len(cod),len(AAnum),len(shit),len(g38),len(g37))
    #19452 19452 19452 19452 19452 19452 19452
    for i in range(0,len(nt)):
		gg38 = g38[i]
		gg37 = g37[i]
		print idnm[i],nt[i],cod[i],AA[i],ATCG[i],new_codon[i],new_aa[i],AAnum[i],shit[i],"GRCh38:Chr:{}:{}".format(Chr,gg38),"GRCh37:Chr:{}:{}".format(Chr,gg37)
		#print(idnm[i],nt[i],cod[i],AA[i],ATCG[i],new_codon[i],new_aa[i],AAnum[i],shit[i],"GRCh38:Chr:{}:{}".format(Chr,gg38),"GRCh37:Chr:{}:{}".format(Chr,gg37))       
   #print()
             
if __name__ == "__main__":
    main()



    
