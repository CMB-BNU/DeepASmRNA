#! /usr/bin/env python
# -*- coding: utf-8 -*-
#coding=utf-8
#blast_parser_hsp_2_1000bp_SE" 

Usage="""NAME.py INPUT OUTPUT"""

import sys
import re
from Bio.Blast import NCBIXML#要把Bio放在当前目录
from Bio import SeqIO

#打开两个文件，blast输出结果xml，输出文件
transFile = sys.argv[2]
seqs={}
for seq_record in SeqIO.parse(transFile, "fasta"):
    seqs[seq_record.id]=seq_record.seq._data


def coverage_query_alignment_obj(alignmentObj,e_value_thresh=1000000):#得到目标序列与比对序列的identity

    this_matches = 0#初始化

    for hsp in alignmentObj.hsps:#对所有的hsp进行循环
        if hsp.expect <= e_value_thresh:#对于阈值符合条件的
            this_query = blast_record.query_length#获取序列的长度
            this_matches += hsp.identities#将其中的identity累加
            #print(this_matches)
    if this_matches >0:
        query_percent =  (this_matches * 100. )/ this_query#将identity的总和除以序列长度
        query_percent = round(query_percent,2)
        #print(query_percent)
        subject_percent =  (this_matches * 100. )/ alignment.length#将identity的总和除以序列长度
        subject_percent = round(subject_percent,2)
        #print(subject_percent)

    coverage_percent = max(query_percent, subject_percent)
    #print(query_percent, subject_percent)
    return coverage_percent


def identity_query_alignment_obj(hsp,e_value_thresh=1000000):#得到目标序列与比对序列的identity
    
    if hsp.expect <= e_value_thresh:#将其中的
        this_percent =  hsp.identities * 100/ hsp.align_length#将比对上的长度除以hsp的长度，就是除去gap和mismatch

    return this_percent

#def coverage(hsp):
    
#hspList =[str(qhStart), str(qhEnd), queryName, subjectName, str(shStart), str(shEnd), str(hspLength), str(round(queryIdentityPercent,2)),str(queryLength),str(subjectLength)]  
def adj(step,up,aa,down):
    if type(up) == bytes:
        up = up.decode('utf-8')
    if type(aa) == bytes:
        aa = aa.decode('utf-8')
    if type(down) == bytes:
        down = down.decode('utf-8')
    newup=up=up.upper()
    newaa=aa=aa.upper()
    newdown=down=down.upper()
    if step >-3:
        step=5
    else:
        step=abs(step)
    #if step >=0:
    #    return newup+","+newaa+","+newdown
    gt = [i.start() for i in re.finditer("TG",''.join(reversed(up[-step:])))]#GT所在位置
    ag = [i.start() for i in re.finditer("AG",down[:step])]#AG所在位置
    aa_gt = [i.start() for i in re.finditer("GT",aa[:step])]#GT所在aa的位置
    aa_ag = [i.start() for i in re.finditer("GA",''.join(reversed(aa[-step:])))]#AG所在aa的位置
    #print(''.join(reversed(up[-step:])))
    #print(gt)
    #print(ag)
    #print(aa_gt)
    #print(aa_ag)
    if gt:
        if ag:
            newup = up[:-gt[0]-2]#从头截止到GT
            newdown=down[ag[0]+2:]#从AG到尾
            newaa = up[-gt[0]-2:]+aa+down[:ag[0]+2]#aa加上两端移动的
        else:
            newup = up[:-gt[0]-2]#从头截止到GT
            #newdown=aa[-gt[0]-2:]+down#加上aa的一部分，从头到尾
            #newaa = up[-gt[0]-2:]+aa[:-gt[0]-2]#up端移动的加aa去除down端减少的
            newdown=down#加上aa的一部分，从头到尾
            newaa = up[-gt[0]-2:]+aa#up端移动的加aa去除down端减少的
    else:
        if ag:
            #newup = up+aa[:ag[0]+2]#从头截止到GT
            #newdown=down[ag[0]+2:]#加上aa的一部分，从头到尾
            #newaa = aa[ag[0]+2:]+down[:ag[0]+2]
            newup = up#从头截止到GT
            newdown=down[ag[0]+2:]#加上aa的一部分，从头到尾
            newaa = aa+down[:ag[0]+2]
        else:
            if aa_gt:
                if aa_ag:
                    if len(aa)>4 and aa_ag[0] !=0:
                        newup = up+aa[:aa_gt[0]]#从头截止到aa的GT
                        newdown=aa[-aa_ag[0]:]+down#加上aa的一部分，从头到尾
                        newaa = aa[aa_gt[0]:-aa_ag[0]]#up端移动的加aa去除down端减少的
                else:
                    newup = up+aa[:aa_gt[0]]#从头截止到aa的GT
                    #newdown=down[aa_gt[0]:]#去除移动的距离
                    #newaa = aa[aa_gt[0]:]+down[:aa_gt[0]]#up端移动的加aa去除down端减少的
                    newdown=down#去除移动的距离
                    newaa = aa[aa_gt[0]:]#up端移动的加aa去除down端减少的
            else:
                if aa_ag:
                    if aa_ag[0] !=0:
                        #newup = up[:-aa_ag[0]]#从头截止到aa的GT
                        #newdown=aa[-aa_ag[0]:]+down#去除移动的距离
                        #newaa = up[-aa_ag[0]:]+aa[aa_ag[0]:-aa_ag[0]]#up端移动的加aa去除down端减少的
                        newup = up#从头截止到aa的GT
                        newdown=aa[-aa_ag[0]:]+down#去除移动的距离
                        newaa = aa[:-aa_ag[0]]#up端移动的加aa去除down端减少的
    if gt and aa_ag:
        dong =  [val for val in gt if val in aa_ag]
        if dong:
            newup = up[:-dong[0]-2]#从头截止到GT
            newdown=aa[-dong[0]:]+down#加上aa的一部分，从头到尾
            newaa = up[-dong[0]-2:]+aa[:-dong[0]]#up端移动的加aa去除down端减少的
        else:
            newup = up[:-gt[0]-2]#从头截止到GT
            newdown=aa[-aa_ag[0]:]+down#加上aa的一部分，从头到尾
            newaa = up[-gt[0]-2:]+aa[:-aa_ag[0]]#up端移动的加aa去除down端减少的
            if aa_ag[0] == 0:
                newdown=down#加上aa的一部分，从头到尾
                newaa = up[-gt[0]-2:]+aa#up端移动的加aa去除down端减少的

    if aa_gt and ag:
        dong =  [val for val in ag if val in aa_gt]
        if dong:
            newup = up+aa[:dong[0]]#从头截止到GT
            newdown=down[dong[0]+2:]#加上aa的一部分，从头到尾
            newaa = aa[dong[0]:]+down[:dong[0]+2]#up端移动的加aa去除down端减少的
        else:
            newup = up+aa[:aa_gt[0]]#从头截止到GT
            newdown=down[ag[0]+2:]#加上aa的一部分，从头到尾
            newaa = aa[aa_gt[0]:]+down[:ag[0]+2]#up端移动的加aa去除down端减少的
    if newup[-1:]+newaa[:1] == "GT":
        newup = newup[:-1]
        newaa = "G"+newaa
    if newaa[-1:]+newdown[:1] == "AG":
        newaa = newaa+"G"
        newdown = newdown[1:]
    return newup+","+newaa+","+newdown


def compare(thishsp, lasthsp):#hsp include complenment as follow: Qstart Qend Sstart Send diffslen

    if lasthsp:
        (thisqStart, thisqEnd, thisqName, thissName, thissStart, thissEnd, thisdiffslen, thisidentity, thiscoverage, thisqlen, thisslen) = thishsp
        (lastqStart, lastqEnd, lastqName, lastsName, lastsStart, lastsEnd, lastdiffslen, lastidentity, lastcoverage, lastqlen, lastslen) = lasthsp
        if int(thisdiffslen) > 0 and int(lastdiffslen) > 0:
            if int(thisqEnd) > int(lastqEnd) and int(thissEnd) > int(lastsEnd) :#因为已经排序，所以要排除hsp相互包含的情况，只需判断后者的end比前者小即可
                gapQ = int(thisqStart) - int(lastqEnd)#----- 5 -----
                gapS = int(thissStart) - int(lastsEnd)#-----   100  ------
                mingap = min(gapQ,gapS)
                maxgap = max(gapQ,gapS)
                diffgap = maxgap - mingap
                if mingap <= 1 and mingap >= -10 and diffgap > 1 and diffgap <= 500:#小gap离得近，大gap离得远
                    #return True#改之前
                    #return thisqName+"+"+thissName+"\t"+str(mingap)+"-"+str(diffgap) +"\t"+ thisidentity +"\t" + thiscoverage+"\t"+lastidentity +"\t" + lastcoverage#为了得到可变长度
                    startq = 1
                    starts = 1
                    if int(lastqEnd)-50 > 1:
                        startq = int(lastqEnd)-50
                    if int(lastsEnd)-50 > 1:
                        starts = int(lastsEnd)-50
                    if gapQ > gapS:
                        if gapQ <= 0:
                            newseq = adj(mingap,seqs[thisqName][startq:int(lastqEnd)],seqs[thisqName][int(lastqEnd):int(lastqEnd)+diffgap],seqs[thisqName][int(lastqEnd)+diffgap:int(lastqEnd)+diffgap+50])
                            return thisqName+"+"+thissName+","+newseq
                        else:
                            newseq = adj(mingap,seqs[thisqName][startq:int(lastqEnd)],seqs[thisqName][int(lastqEnd):int(thisqStart)],seqs[thisqName][int(thisqStart):int(thisqStart)+50])#为了序列
                            return thisqName+"+"+thissName+","+newseq#为了序列
                    else:
                        if gapS <= 0:
                            newseq = adj(mingap,seqs[thissName][starts:int(lastsEnd)],seqs[thissName][int(lastsEnd):int(lastsEnd)+diffgap],seqs[thissName][int(lastsEnd)+diffgap:int(lastsEnd)+diffgap+50])
                            return thisqName+"+"+thissName+","+newseq
                        else:
                            newseq = adj(mingap,seqs[thissName][starts:int(lastsEnd)],seqs[thissName][int(lastsEnd):int(thissStart)],seqs[thissName][int(thissStart):int(thissStart)+50])#为了序列
                            return thisqName+"+"+thissName+","+newseq#为了序列
                else:
                    return False #"bufuhetiaojian"
        #if int(thisdiffslen) < 0 and int(lastdiffslen) < 0:
        #   if int(thisqEnd) > int(lastqEnd) and int(thissStart) > int(lastsStart) :#对于负链来说，应该是start比end大：：：：
        #       gapQ = abs(int(thisqStart) - int(lastqEnd))
        #       gapS = abs(int(thissEnd) - int(lastsStart))
        #       if gapQ <= 5 and gapS >= 10 or gapS <= 5 and gapQ >= 10:
        #           return True
        #       else:
        #           return False
        
        #两列序列对比
        #            if gapQ > gapS:
        #                if gapQ <= 0:
        #                    return thisqName+"+"+thissName+","+mingap+"+"+maxgap+"+"+str(diffgap)+"+"+seqs[thisqName][startq:int(lastqEnd)]+","+seqs[thisqName][int(lastqEnd):int(lastqEnd)+diffgap]+","+seqs[thisqName][int(lastqEnd)+diffgap:int(lastqEnd)+diffgap+50]+"\n"+thisqName+"+"+thissName+","+mingap+"+"+maxgap+"+"+str(diffgap)+"+"+seqs[thissName][starts:int(lastsEnd)]+","+seqs[thissName][int(lastsEnd):int(lastsEnd)+50]
        #                else:
        #                    return thisqName+"+"+thissName+","+mingap+"+"+maxgap+"+"+str(diffgap)+"+"+seqs[thisqName][startq:int(lastqEnd)]+","+seqs[thisqName][int(lastqEnd):int(thisqStart)]+","+seqs[thisqName][int(thisqStart):int(thisqStart)+50]+"\n"+thisqName+"+"+thissName+","+mingap+"+"+maxgap+"+"+str(diffgap)+"+"+seqs[thissName][starts:int(lastsEnd)]+","+seqs[thissName][int(lastsEnd):int(lastsEnd)+50]#为了序列
        #            else:
        #                if gapS <= 0:
        #                    return thisqName+"+"+thissName+","+mingap+"+"+maxgap+"+"+str(diffgap)+"+"+seqs[thissName][starts:int(lastsEnd)]+","+seqs[thissName][int(lastsEnd):int(lastsEnd)+diffgap]+","+seqs[thissName][int(lastsEnd)+diffgap:int(lastsEnd)+diffgap+50]+"\n"+thisqName+"+"+thissName+","+mingap+"+"+maxgap+"+"+str(diffgap)+"+"+seqs[thisqName][startq:int(lastqEnd)]+","+seqs[thisqName][int(lastqEnd):int(lastqEnd)+50]
        #                else:
        #                    return thisqName+"+"+thissName+","+mingap+"+"+maxgap+"+"+str(diffgap)+"+"+seqs[thissName][starts:int(lastsEnd)]+","+seqs[thissName][int(lastsEnd):int(thissStart)]+","+seqs[thissName][int(thissStart):int(thissStart)+50]+"\n"+thisqName+"+"+thissName+","+mingap+"+"+maxgap+"+"+str(diffgap)+"+"+seqs[thisqName][startq:int(lastqEnd)]+","+seqs[thisqName][int(lastqEnd):int(lastqEnd)+50]#为了序列
        #正常序列
        #            if gapQ > gapS:
        #                if gapQ <= 0:
        #                    return thisqName+"+"+thissName+","+str(mingap)+"+"+str(maxgap)+"+"+str(diffgap)+"+"+seqs[thisqName][startq:int(lastqEnd)]+","+seqs[thisqName][int(lastqEnd):int(lastqEnd)+diffgap]+","+seqs[thisqName][int(lastqEnd)+diffgap:int(lastqEnd)+diffgap+50]
        #                else:
        #                    return thisqName+"+"+thissName+","+mingap+"+"+maxgap+"+"+str(diffgap)+"+"+seqs[thisqName][startq:int(lastqEnd)]+","+seqs[thisqName][int(lastqEnd):int(thisqStart)]+","+seqs[thisqName][int(thisqStart):int(thisqStart)+50]#为了序列
        #            else:
        #                if gapS <= 0:
        #                    return thisqName+"+"+thissName+","+mingap+"+"+maxgap+"+"+str(diffgap)+"+"+seqs[thissName][starts:int(lastsEnd)]+","+seqs[thissName][int(lastsEnd):int(lastsEnd)+diffgap]+","+seqs[thissName][int(lastsEnd)+diffgap:int(lastsEnd)+diffgap+50]
        #                else:
        #                    return thisqName+"+"+thissName+","+mingap+"+"+maxgap+"+"+str(diffgap)+"+"+seqs[thissName][starts:int(lastsEnd)]+","+seqs[thissName][int(lastsEnd):int(thissStart)]+","+seqs[thissName][int(thissStart):int(thissStart)+50]#为了序列

            
    else:
        return False

def sort_key(s):
    # 排序关键字匹配
    # 匹配开头数字序号
    if s:
        try:
            c = re.findall('^\d+', s)[0]
        except:
            c = -1
        return int(c)

def strsort(alist):
    alist.sort(key=sort_key)
    return alist

#打开两个文件，blast输出结果xml，输出文件
inPutFile=sys.argv[1]
outPutFile = sys.argv[3]
#存到两个文件句柄
result_handle = open(inPutFile,"r")
out_handle = open(outPutFile, "w")
hspsss = open("hsp.txt", "w")
allinfo = open("allinfo.txt", "w")
cover = open("coveragebotname.txt", "w")
#输出文件的文件头
header = ['QueryName', 'SubjectName', 'QhStart', 'QhEnd', 'ShStart', 'ShEnd', 'hsplength', 'QueryIdentity', 'queryLength', 'subjectLength', 'HSP']
out_handle.write("\t".join(header)+'\n')

#读入xml格式的文件，blast结果存在Blast_records中
Blast_records = NCBIXML.parse(result_handle)

for blast_record in Blast_records:
    queryName = blast_record.query#query序列名称
    queryLength = int(blast_record.query_length)#query序列长度
    subjectName = ""#subject序列有好几条，所以先赋空值
    subjectLength = ""
    queryIdentityPercent = 0 #序列的identity
    for alignment, description in zip(blast_record.alignments,blast_record.descriptions):
        hspNum = int(len(alignment.hsps)) #每一个比对，都有很多hsp，High-scoring segment pair。
        subjectName = alignment.title
        subjectNames = subjectName.split(' ')
        subjectName = subjectNames[1]#有些转录本名中间有空格
        subjectLength = int(alignment.length)
        #print(queryName,subjectName,sep=",")
        querycoveragePercent = coverage_query_alignment_obj(alignment)
        #print(querycoveragePercent,queryName,subjectName,sep=",",file=cover)
        breakflag = 0
        if queryName != subjectName and hspNum > 1 and querycoveragePercent > 70:#将不是与自己相比较的去除掉，去除hsp为空的情况
            if queryLength >= 100 and subjectLength >= queryLength:#对比序列长度大于300，100aa对应的长度。设置subject比较长，可以防止重复计数
                allhsp = []#定义将所有hsp都存进来的数组
                lasthsp = []
                thishsp = []
                for hsp in alignment.hsps:
                    queryIdentityPercent = identity_query_alignment_obj(hsp)
                    hspLength = hsp.align_length#hsp的长度
                    qhStart = hsp.query_start
                    qhEnd = hsp.query_end
                    shStart = hsp.sbjct_start
                    shEnd = hsp.sbjct_end
                    if queryIdentityPercent > 90:#4step中关于hsp长度和identity的限制
                        #hspList =[str(qhStart), str(qhEnd), queryName, subjectName, str(shStart), str(shEnd), str(hspLength), str(round(queryIdentityPercent,2)),str(queryLength),str(subjectLength)]  
                        hspList =[str(qhStart), str(qhEnd), queryName, subjectName, str(shStart), str(shEnd), str(hspLength), str(round(queryIdentityPercent,2)), str(round(querycoveragePercent,2)),str(queryLength),str(subjectLength)]   
                        hspstr = '\t'.join(hspList)#把信息合成一个字符串，方便后面排序，制表符分割，方便后续拆解
                        allhsp.append(hspstr)#将所有的hsp的信息字符串加入到allhsp这个数组里面
                    else:
                        breakflag = 1
                allhsp_sort = strsort(allhsp)#按照在序列上的位置，进行排序。从小到大
                print(queryName,subjectName,allhsp_sort,file=hspsss)
                if breakflag:
                    break

                if len(allhsp_sort) > 1:
                    as_num=0
                    for hsp in allhsp_sort:
                        thishsp = hsp.split('\t')
                        figure = compare(thishsp, lasthsp)
                        if figure:
                            as_num=as_num+1
                            thishsp[0],thishsp[1],thishsp[2],thishsp[3] = thishsp[2],thishsp[3],thishsp[0],thishsp[1]
                            lasthsp[0],lasthsp[1],lasthsp[2],lasthsp[3] = lasthsp[2],lasthsp[3],lasthsp[0],lasthsp[1]
                            out_handle.write("\t".join(lasthsp)+'\thsp1\t'+figure+'\n'+"\t".join(thishsp)+'\thsp2\n')
                            thishsp[2],thishsp[3],thishsp[0],thishsp[1] = thishsp[0],thishsp[1],thishsp[2],thishsp[3]
                            lasthsp[2],lasthsp[3],lasthsp[0],lasthsp[1] = lasthsp[0],lasthsp[1],lasthsp[2],lasthsp[3]
                            print(figure,sep=",")
                        lasthsp = thishsp
                    #if as_num :
                        #print(queryName,subjectName,as_num)
                        #outList =[queryName, str(queryLength), subjectName, str(subjectLength), str(qhStart), str(qhEnd), str(shStart), str(shEnd), str(hspNum), str(diffslen), str(round(queryIdentityPercent,2)), str(hspLength)]


result_handle.close()
out_handle.close()
