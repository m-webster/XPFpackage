from lxml import html
import requests
import numpy as np
import os.path

dirname = 'C:/Users/mark/Dropbox/PhD/XPFpackage/code examples/code-tables/'

def tableScrape():
    page = requests.get('http://www.codetables.de/TableIII256.php')
    tree = html.fromstring(page.content)
    constr = tree.xpath("//a[@target='bounds']/@href")
    f = open('codeList.txt', 'w')
    f.write('\n'.join(constr))
    f.close
    
def scrapeAll():
    global dirname
    ## get query string
    f = open(dirname + 'codeList.txt', 'r')
    for qstring in f:
        qstring = qstring.replace("\n","")
        pageScrape(qstring)
    f.close  

def pageScrape(qstring):
    global dirname

    ## get n and k to set filename
    qvals = parseQstring(qstring) 
    n,k = qvals['n'],qvals['k']
    print(n,k)
    filename = f'{n}-{k}.txt'
    if not os.path.exists(dirname + filename):
        ## get page content
        mypath = 'http://www.codetables.de/'+qstring 
        page = requests.get(mypath)

        ## write to filename
        f = open(dirname + filename, 'wb')
        f.write(page.content)
        f.close

def pageParse(n,k):
    global dirname
    filename = f'{n}-{k}.txt'
    f = open(dirname + filename, 'r')
    tree = html.fromstring(f.read())
    mystr = tree.xpath('//pre/text()')
    myarr = mystr[0].split("\n\n")
    temp = []
    if len(myarr)> 3:
        mystr = myarr[2]
        for s in ['  ','[',']']:
            mystr = mystr.replace(s,'')
        myarr = mystr.split('\n')
        for i in range(len(myarr)):
            mycols = myarr[i].split('|')
            temp.append([c.split(' ') for c in mycols])
        temp = [[[int(x) for x in c] for c in r] for r in temp]
    return temp

def parseQstring(qstring):
    myarr = qstring.split("?")
    temp = dict()
    if len(myarr) == 0:
        return temp
    myarr = myarr[1].split('&')
    for a in myarr:
        a = a.split("=")
        if len(a) > 0:
            temp[a[0]] = a[1]
    return temp

# n,k=20,6
# sList = pageParse(n,k)
# print(np.array(sList))
# scrapeAll()
