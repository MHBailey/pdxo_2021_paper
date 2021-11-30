def captureRestofRow(s,r,c,header): #sheet row and col
    out = []
    for col in range(c+header,s.ncols):
        x = s.cell(r,col)
        value = x.value
        try: value = str(float(value)) #Keep if string or number
        except: pass
        if value != "":
            out.append(value)

    return out

def captureRestofCol(s,r,c,header): #sheet row and col 
    out = []
    for row in range(r+header,s.nrows):
        x = s.cell(row,c)
        value = x.value
        try: value = str(float(value))
        except: pass 
        if value != "":
            out.append(value)
    
    return out 


def captureMiniTable2Anchor(s,r,c, stop):
    out = []
    for col in range(c,s.ncols):
        x = s.cell(r,col)   
        value = x.value
        try: value = str(float(value))
        except: pass 
        if value == stop:
            break
        if value != stop:
            thisCol = captureRestofCol(s,r,col,1)
            if thisCol: #Check to make sure there is data in column
                out.append(captureRestofCol(s,r,col,1))
    
    return out 



def divide_chunks(l, n): #creates a generator 
    # looping till length l 
    for i in range(0, len(l), n):  
        yield l[i:i + n]
 
def sum_chunks(l,n): #creates a generator 
    for i in range(0, len(l), n):
        yield sum(l[i:i+n])

def avg_list(l,w):
    for i in l:
        yield sum(i)/(len(i)*w)

def normalize_chunks(d,l):
    for i in range(0,len(l)):
        for j in l[i]:
            yield(j/d[i])
        

def captureDay0(s,r,c, reps, numconc): #add 1 to starting value to int lists 
    out = []
    for col in range(c,c+reps):
        l = [float(x) for x in captureRestofCol(s,r,col,1)]
        out.append(list(sum_chunks(l,numconc))) #list generator def 
    
    oo = list(zip(*out))
    return list(avg_list(oo,numconc)) 

def captureSolvent(s,r,c, reps, numconc, d0):
    out = []
    for col in range(c,c+reps):
        l = [float(x) for x in captureRestofCol(s,r,col, 1)]
        ll = list(divide_chunks(l,numconc))
        lll = list(normalize_chunks(d0,ll))
        out.append(list(sum_chunks(lll,numconc))) #list generator def 
    
    oo = list(zip(*out))
    return list(avg_list(oo,numconc))

def formatKat(s,r,c,d0,con,solv,sr):
    nconc = len(set(con))

    n = [float(x) for x in captureRestofCol(s,r,c,1)]
    l = list(divide_chunks(n,nconc))
    nn = list(normalize_chunks(d0,l)) #Normalize to Day0
    ll = list(divide_chunks(nn,nconc))
    nnn = list(normalize_chunks(solv,ll)) #Normalize Day0 to solvent
    lll = list(divide_chunks(nnn,nconc))

    raw = n
    dnorm = nn
    snorm = nnn
    #NOTE: added Dec 6. 2019. I need to capture day 0 raw value
    redat0 = np.repeat(d0,nconc)
    resolv = np.repeat(sr,nconc)

    o = list(zip(raw,dnorm,snorm,redat0,resolv))
    oo = list(divide_chunks(o,nconc))
    
    return(oo)



def printALL(A,p,d,ids,c,sn):
    assert(len(A) == len(p)) #means I have the same number of arrays in ALL as plates/drugs replicates
    assert(len(p) == len(d)) #headers match
    
    #NOTE: If more data is needed add column here.
    normalized = ["Raw","Day0","Solvent","RawDay0","Vehicle"]
    dcnt = 0

    for i in range(0,len(p)):
        if dcnt == 4:
            dcnt = 0
        dcnt += 1
        plate = p[i]
        drug = d[i]
        for j in range(0,len(A[i])):
            oid = ids[j][0]
            for k in range(0,len(A[i][j])):
                concentration = c[k]
                for m in range(0,len(A[i][j][k])):
                    norm = normalized[m]
                    value = A[i][j][k][m]
                    variable = "V"+str(dcnt) #For completely long format 
                    tissue = (oid.replace(" ","_")+"_plate_"+str(plate))
                    out = [tissue,drug,str(concentration),norm,sn,variable,str(value)]
                    print("\t".join(out))
    #sys.exit()


def parseKat(f1):
    out = ["Plate","Drug","Vol","Normalized","DrugSet","Variable","Value"]
    print("\t".join(out))

    wb = xlrd.open_workbook(f1)
    for sheet in wb.sheets():
        platenums = []
        targets = []
        concentrations = [] 
        ids = [] 
        reps = 4 
        numconc = 8 #This will throw an error if concentrations are different 
        capturedDay0 = False
        firstDMSO = False
        d0avgs = []
        drugs = []
        solvent = []
        #NOTE: Dec 6 Added for raw value 
        solventR = []

        dmsoCnt = 0
        drugCnt = 0
        excelCol = 0
        sheetname = sheet.name
        ALL = []



        for row in range(0,sheet.nrows):
            for col in range(0,sheet.ncols):
#                print(str(row)+":"+str(col))
                x = sheet.cell(row,col)
                value = x.value
                try: value = str(int(value)) #Keep if string or number
                except: pass
                if "plate" in value.lower():
                    #print(row,col)
                    platenums = captureRestofRow(sheet,row,col,1)
  
                if "target" in value.lower():
                    #print(row,col)
                    targets = captureRestofRow(sheet,row,col,1)
                
                if "identifier" in value.lower():
                    #print(row,col)
                    o = captureMiniTable2Anchor(sheet,row,col,"Day 0") #Dictionary 
                    #print(o)
                    ids = list(zip(*o)) #the * notation is an "UnPACKING" step ;) so cool
                    #print(ids)
                
                if "Day 0" == value and capturedDay0 == False:
                    d0avgs = captureDay0(sheet, row, col, reps, numconc)    
                    capturedDay0 = True

                if "uM" == value:
                    concentrations = captureRestofCol(sheet,row,col,1)
                 
                if "dmso" in value.lower() and firstDMSO == False: #this will finish off the table with the structure in place
                    assert(len(set(concentrations)) == numconc) #This checks the hardcoded data 
                    drugs = [x.lower() for x in  captureRestofRow(sheet,row,col,0)]
                    excelCol = col
                    firstDMSO = True
 
                if "dmso" in value.lower():
                    dmsoCnt += 1
                    if dmsoCnt == 1: 
                        solvent = captureSolvent(sheet, row, col, reps, numconc, d0avgs)
                        solventR = captureDay0(sheet, row, col, reps, numconc)

                    if dmsoCnt == 4:
                        dmsoCnt = 0



                if value.lower() in drugs:
                    drugCnt += 1
                    pos = col-excelCol
                    o = formatKat(sheet,row,col,d0avgs,concentrations,solvent,solventR)
                    ALL.append(o)
                    

    
        printALL(ALL,platenums,drugs,ids,concentrations,sheetname) #NOTE makes sure this is the last step



               
if __name__ == "__main__":
    import sys
    import xlrd
    import datetime
    import numpy as np 
    #datetime.datetime.strptime('24052010', '%d%m%Y').date()
    files = parseKat(sys.argv[1])



