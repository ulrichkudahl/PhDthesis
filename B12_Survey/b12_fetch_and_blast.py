#! /usr/bin/python -u

def main():
   import sys,os,time
   from Bio import SeqIO
   from Bio.Blast.Applications import NcbiblastpCommandline as blastcmd
   from Bio.Blast import NCBIXML
   from datetime import datetime
   import threading
   import multiprocessing
   
   #define blast folder:
   global blast_loc
   blast_loc = '/applications/ncbi-blast+/ncbi-blast-2.2.30+/bin/'
   
   
   global sqldb
   sqldb = '/home/ujk20/b12_blaster/sqldb/b12_march16.db'
   
   global blastdb
   #blastdb = '/data/public_data/nr_blast/nr'
   blastdb = '/data/public_data/NCBI/NCBI_nr/nr'


   global profile_db
   profile_db = "/data/public_data/NCBI/NCBI_cdd/cdd_3_13/Cdd.pn"

   global blast_threads
   blast_threads = 1
   
   max_cpus = 44
   #multiprocessing.cpu_count() -1
   
   if len(sys.argv) != 3:
      print 'program exited. Propper usage: taxid gene_file'
      
      exit()
   print 'running script'
   
   toplevel = sys.argv[1]
   test_gene_file = sys.argv[2]
   global folder
   folder = '/home/ujk20/Dropbox/Cambridge/programming/b12_blaster/temp/'
   #folder = os.path.dirname(os.path.realpath(__file__)) + '/temp/'
   if not os.path.exists(folder):
      os.makedirs(folder)

   
   names_file = '/home/ujk20/taxdump/new_taxid/names.dmp'
   tax2name = taxfile2dict(names_file)
   toplevel_name = tax2name[toplevel]
   print 'Toplevel ' + toplevel_name + '\t' + str(toplevel)

   startTime = datetime.now()  
   
   tax_file ='/home/ujk20/taxdump/taxdump_220316/nodes.dmp'
   (taxids,ranky,parent2children) = chilly(tax_file,toplevel)
   
   gi_tax_file = '/home/ujk20/taxdump/taxdump_220316/gi_taxid_prot.dmp'
   #gi_taxid_prot.dmp'
   #meso_gi.txt
   #p_denti_gi.dmp
   tax2gi_dict = tax2gilist(taxids,gi_tax_file)
   
   print 'Total time: ' + str(datetime.now() - startTime)

   import subprocess
   reciprocal_blastdb  = folder + 'reciprocal_blast_db'
   subprocess.call('cp '+ test_gene_file + ' ' + reciprocal_blastdb, shell=True)
   subprocess.call(blast_loc+'makeblastdb -dbtype "prot" -in ' + reciprocal_blastdb,shell=True)
   #/applications/ncbi-blast+/ncbi-blast-2.2.30+/bin/blastp
   
   
   while not os.path.isfile(blastdb):
         time.sleep(10)
         print 'sleeping'
   
   print 'Working folder: ' + folder

   counter = 0
   over = 0
   threshold = 1000
   sql_db = sql_reader(taxids)
   #sql_db = {}
   (test_gene_list,gene2path) = gene_routine(test_gene_file,folder)
   
   #### To add an extra gene to all the species already in the system
   #taxids = set(sql_db.keys())
   print 'Tax correcting and file creation started' + 'Total time: ' + str(datetime.now() - startTime) 
   (taxids,tax_gene_dict) = tax_corrector(sql_db,taxids,tax2gi_dict,parent2children,test_gene_list,ranky,threshold,folder)
   print 'Tax correcting and file creation done'  + 'Total time: ' + str(datetime.now() - startTime) 
   
   
   #### To add an extra gene to all the species already in the system
   #taxids = set(sql_db.keys())
   
   
   print 'Number of species ' + str(len(taxids))
   global comb_results
   comb_results = {}
   from multiprocessing import Process, Queue, Manager
   processes = []
   indy=0

   mgr = Manager()
   comb_res2 = mgr.dict()
   
   for taxid in tax_gene_dict.keys():
      #print 'number of active threads ' + str(len(multiprocessing.active_children()))
      while len(multiprocessing.active_children()) > max_cpus:
         time.sleep(1)
      processes.append((taxid,Process(target=blast_analyser,args=(folder,blastdb,taxid,tax_gene_dict[taxid],gene2path,comb_results,reciprocal_blastdb,comb_res2))))
      processes[indy][1].start()
      
      indy+=1
      over += 1
      counter +=1
      fraction = max(len(tax_gene_dict)/100,1)
      if counter % fraction == 0:
         fraction = round(float(counter) *100/ len(tax_gene_dict),2)
         print 'Fraction done ' + str(fraction) + '%' + '\t\t' + 'Total time: ' + str(datetime.now() - startTime)  + ' Number of threads running ' + str(len(multiprocessing.active_children()))
   
      
   for p in processes:
      p[1].join()
   #comb_results = comb_res2
#   print comb_res2
   reskeys = comb_res2.keys()
   for i in reskeys:
      comb_results[i] = comb_res2[i]
#http://eli.thegreenplace.net/2012/01/16/python-parallelizing-cpu-bound-tasks-with-multiprocessing/      
#http://stackoverflow.com/questions/2846653/python-multithreading-for-dummies
#
#def f(name):
#    print 'hello', name
#
#if __name__ == '__main__':
#    p = Process(target=f, args=('bob',))
#    p.start()    
#    p.join()
# http://stackoverflow.com/questions/25627313/multiprocessing-dont-use-all-cpu?lq=1



   pid = os.getpid()
   output = open(folder+'results_'+str(pid)+'.txt','w+')
   output2 = open(folder+'sumerized_results'+str(pid)+'.txt','w+')
   for tx in comb_results:
      if tx in tax2gi_dict:
         number_of_proteins = len(tax2gi_dict[tx])
      else:
         number_of_proteins = 0
      if tx in tax2name:
         namy = tax2name[tx]
      else:
         namy =str(tx)
      
      summay = {'aerobic_synthesis_genes':0,'anaerobic_synthesis_genes':0,'Proteins':number_of_proteins}
      anaerobe_not_in = []
      aerobe_not_in = []
      if isinstance(comb_results[tx],dict):
            for item in comb_results[tx]:
               #print item + '\t' + str(comb_results[tx][item])
               if isinstance(comb_results[tx][item],tuple) and len(comb_results[tx][item]) ==5:
                  tub =comb_results[tx][item]
                  
                  if tub[0] == 'aerobic_synthesis':
                     if tub[3] < 10**-10:
                        summay['aerobic_synthesis_genes'] += 1
                     else:
                        aerobe_not_in.append(item)
                  elif tub[0] == 'anaerobic_synthesis':
                     if tub[3] < 10**-10:
                        summay['anaerobic_synthesis_genes'] += 1
                     else:
                        anaerobe_not_in.append(item)
                  elif tub[0] == 'B12_related':
                     if tub[3] <10**-10:
                        summay[item] = 1
                     else:
                        summay[item] = 0
                  elif tub[0] == 'B12_transport':
                     if tub[3] <10**-10:
                        summay[item] = 1
                        #summay['btub_gi'] = tub[4]
                     else:
                        summay[item] = 0
                  
                  elif tub[0] == 'test':
                     summay[item] = tub[3]
      
      if tx in tax2gi_dict:
         if number_of_proteins > threshold:
            #print 'inside here'
            #print namy + '\t' + str(summay)
            output2.write(namy + '\t' + str(summay)+'\n')
         output.write(namy+'\t'+tx+'\t'+str(len(tax2gi_dict[tx]))+'\n')
      else:
         output.write(namy+'\t'+tx+'\t'+str(0)+'\n')
      for item in comb_results[tx]:
         #print comb_results[tx]
         output.write(str(item) + '\t' +str(comb_results[tx][item])+'\n')
      output.write('########################\n\n')
   output.close()
   output2.close()
   #print comb_results
   sql_writer(comb_results,tax2name,tax2gi_dict)
   
   
   print 'Total time: ' + str(datetime.now() - startTime)
   print 'Species/strians analysed ' + str(over)


def gene_routine(test_gene_file,folder):
   from Bio import SeqIO
   import re,os
   gene_list = []
   gene2path = {}
   for seq_rec in SeqIO.parse(test_gene_file,'fasta'):
      match = re.search(r'gene:(.+)\|pathway:(.+)',seq_rec.id)
      if match:
         gene = match.group(1)
         pathway = match.group(2)
         
         gene2path[gene] = pathway
         
         gene_list.append(gene)
         blastin = folder+str(gene)+'_in.fasta'
         if not os.path.isfile(blastin):
            out = open(blastin,'w')
            out.write('>'+seq_rec.id+'\n'+str(seq_rec.seq)+'\n')
            out.close()
   return (gene_list,gene2path)


def tax_corrector(sql_db,taxids,tax2gi_dict,parent2children,test_gene_list,ranky,threshold,folder):
   
   import os
   taxids2 = []
   
   for id in taxids:
      if ranky[id] == 'species' and id in tax2gi_dict:
         if len(tax2gi_dict[id]) >threshold:
            taxids2.append(id)
         else:
            try:
               children = parent2children[id]
            except:
               children= [0]
            temp_dict = {0:0}
            for child in children:
               if child in tax2gi_dict:
                  temp_dict[child] = len(tax2gi_dict[child])
               else:
                  temp_dict[child] = 0
            key,value = max(temp_dict.iteritems(), key=lambda x:x[1])
            if value > threshold:
               taxids2.append(key)
            else:
               taxids2.append(id)
      elif ranky[id] == 'species':
         taxids2.append(id)
   
   #tax_gene_dict = dict.fromkeys(taxids2,[])
   tax_gene_dict ={}
   for tt in taxids2:
      if tt in tax2gi_dict:
         if len(tax2gi_dict[tt]) > threshold:
            tax_gene_dict[tt] = []
   
   
   for taxid in tax_gene_dict:
      if taxid in sql_db:
         genes_in_db = sql_db[taxid]
      else:
         genes_in_db = []            
      for test_gene in test_gene_list:
         if test_gene not in genes_in_db:
            tax_gene_dict[taxid].append(test_gene)
   
   
   
   
   in_folder = os.listdir(folder)
   files_written = 0
   for t in tax_gene_dict:
      if len(tax_gene_dict) == 0:
         del tax_gene_dict[t]
      gi_list = t+'gi_list.fasta'
      if gi_list not in in_folder:
         files_written +=1
         gi_file = open(folder+gi_list,'w+')
         for gi in tax2gi_dict[t]:
            gi_file.write(str(gi)+'\n')
         gi_file.close()
   print 'GI files written ' + str(files_written)
            
   return(taxids2,tax_gene_dict)



def sql_reader(taxids):
   import sqlite3 as lite
   import os
   
   #dir = '/home/ujk20/b12_blaster'
   dir = os.path.dirname(os.path.realpath(__file__))
   con = lite.connect(sqldb)
   prev = {}
   with con:
      con.row_factory = lite.Row
      cur = con.cursor() 
      cur.execute("SELECT * FROM test")
      rows = cur.fetchall()
      for row in rows:
         tax = str(row['taxid'])
         gene = row['gene']
         if tax in prev:
            prev[tax].append(gene)
         else:
            prev[tax] = [gene]
   return(prev)
   
def sql_writer(res_dic,tax2name,tax2gi_dict):
   import sqlite3 as lite
   import os
   
   dir = os.path.dirname(os.path.realpath(__file__))
   #dir = '/home/ujk20/b12_blaster'
   print  'farfar dir ' +dir
   con = lite.connect(sqldb)
#   headers = ('id int,Type Text, Sci_name TEXT, num_prot int, pathway text, coverage real, identities real, e_val real')
   
   with con:
      cur = con.cursor()
      for tx in res_dic:
         if tx in tax2name:
            sciname = tax2name[tx]
         else:
            sciname = tx
         if tx in tax2gi_dict:
            num_prot = len(tax2gi_dict[tx])
         else:
            num_prot = 0
            tubby = (tx,sciname,num_prot,'na',0,0,float(100))
            #print tubby
            continue
         for gene in res_dic[tx]:
           # print gene + '\t' + str(res_dic[tx][gene])
            pathway = res_dic[tx][gene][0]
            coverage = round(float(res_dic[tx][gene][1]),5)
            identities = round(float(res_dic[tx][gene][2]),5)
            eval = float(res_dic[tx][gene][3])
            gi = res_dic[tx][gene][4]
            reciprocal_gene = str(res_dic[tx][gene][5])
            domains_list = str(res_dic[tx][gene][6])
#            print 'reci far fafr' + reciprocal_gene +'farfar'
            tubby = (tx,sciname,num_prot,gene,pathway,coverage,identities,eval,gi,reciprocal_gene,domains_list)
            #print tubby
            #(taxid int, Sci_name TEXT, num_prot int, gene text ,pathway text, coverage real, identities real, e_val real,GI text)
            
            cur.execute("INSERT INTO test VALUES "+str(tubby))
            
            
            
def blast_analyser(folder,blastdb,taxid,test_gene_list,gene2path,comb_results,reciprocal_blastdb,comb_res2):
   #gene_list =tax_gene_dict[taxid]

   import Bio.Blast.Applications 
   import Bio.Blast.NCBIXML
   import cStringIO
   import re
   import subprocess
   
   results = {}
   for gene in test_gene_list:
      
      resline = ()
      reciprocal_gene = "NA"
      pathway = gene2path[gene]
      results[gene] = (pathway,0,0,100,0000000)
      blastin = folder+str(gene)+'_in.fasta'
      blastout = folder+str(gene)+'_blast_out.fasta'
      gi_list = folder+taxid+'gi_list.fasta'
      blastp_cline = Bio.Blast.Applications.NcbiblastpCommandline(cmd=blast_loc+'blastp',num_threads=blast_threads, max_target_seqs=10, query=blastin, db=blastdb, evalue=10,outfmt=5,gilist=gi_list)   
      stdout, stderr = blastp_cline()      
      blast_xml = cStringIO.StringIO(stdout)
      
#      for l in blast_xml:
 #        print l
      
      blast_record = Bio.Blast.NCBIXML.read(blast_xml)
      query_length = blast_record.query_length
      if len(blast_record.alignments) > 0:
         for align in blast_record.alignments:
            hsp = align.hsps[0]
            coverage = 0
            coverage  = float(hsp.align_length) / max(float(query_length),1)
            #print blast_record.query + '\t' +'eval' + '\t'+ str(hsp.expect) + '\t' + 'alignment_lengt' +   '\t' + str(hsp.align_length) + '\t' + 'input_lenght'+ '\t' + 'coverage' + '\t' +str(coverage) + '\t' + 'Title '+ '\t' + str(align.title)
            match = re.search('gi\|(\w+)\|',align.title)
            if match:
               gi = str(match.group(1))
            else:
               gi = 'NA'
            break

         target_file = folder + "target_seq_"+taxid+".fasta"
         
         get_entry_cmd = blast_loc+"blastdbcmd -db " + blastdb + " -dbtype 'prot'  -entry " + gi + ' > ' +target_file
         output = subprocess.call(get_entry_cmd,shell=True)
         del output
         reciprocal_blast_gene = reciprocal_blast(target_file,reciprocal_blastdb)
         domains_list = domain_finder(target_file,gene)
         resline = (pathway,coverage,float(hsp.identities)/float(hsp.align_length),hsp.expect,gi,reciprocal_blast_gene,domains_list)         
         #print resline
      else:
         resline = (pathway,0,0,100,0000000,"NA",[])      
      results[gene] = resline
   

   comb_res2[taxid] = results
   comb_results[taxid] = results

def domain_finder(target_file,gene):
   import re,Bio.Blast.Applications,Bio.Blast.NCBIXML,cStringIO,itertools
      
   blastp_cline = Bio.Blast.Applications.NcbiblastpCommandline(cmd=blast_loc +'rpsblast',num_threads=blast_threads, query=target_file, db=profile_db, evalue=0.01,outfmt=5)
   stdout, stderr = blastp_cline()
   blast_xml = cStringIO.StringIO(stdout)         
   blast_record = Bio.Blast.NCBIXML.read(blast_xml)
   prev = 0
   domains_list = []
   domain = ""
   eval = 10000
   if len(blast_record.alignments) > 0:
            for align in blast_record.alignments:
               for hsp in align.hsps:
                  if hsp.score <= prev:
                     continue
                  match = re.search('\s(.+?),',align.title)
                  if match:
                     domain = match.group(1)
                     eval = hsp.expect      
               domains_list.append((domain,eval))
   return (domains_list)
   
def  reciprocal_blast(target_file,reciprocal_blastdb):
         import re,Bio.Blast.Applications,Bio.Blast.NCBIXML,sys,cStringIO
         #Reciprocal Matches serach
         
         #Finding the reciprocal result
         blastp_cline = Bio.Blast.Applications.NcbiblastpCommandline(cmd=blast_loc +'blastp',num_threads=1, max_target_seqs=1, query=target_file, db=reciprocal_blastdb, evalue=10,outfmt=5)   
         stdout, stderr = blastp_cline()
         
         blast_xml = cStringIO.StringIO(stdout)
         blast_record = Bio.Blast.NCBIXML.read(blast_xml)
         prev = 0
         
         if len(blast_record.alignments) > 0:
            for align in blast_record.alignments:
               for hsp in align.hsps:
                  if hsp.score <= prev:
                     continue
                  match = re.search('gene:(.+)\|pathway',align.title)
                  if match:
                     gene_reciprocal = match.group(1)
                     return(gene_reciprocal)
                  else:
                     print 'no match - Something is wrong in reciprocal matching -'
         return('NA')
      
def tax2prot(taxid,out_file):
   import os.path
   from Bio import SeqIO
   #if os.path.isfile(out_file):
   #   print 'Already downloaded genome'
   #   return(0)
   

   from Bio import Entrez
   import os
   import sys
   Entrez.email = 'ujk20@cam.ac.uk'
   searchhandle = Entrez.esearch(db='protein',term=str('txid'+taxid),retmode = 'fasta',retmax = 1000000)
   searchresults = Entrez.read(searchhandle)
   protids = searchresults['IdList']
   len_prot = len(protids)
   
   out = open(out_file,'w+')
   for prot in protids:
      out.write(prot+'\n')
   out.close()
   return(len_prot)

def chilly(dump_file,id):   
   dump = open(dump_file,'rU')
   parent2children = {}
   ranky = {}
   for line in dump:
      line_list = line.split('|')
      
      child = str(line_list[0]).strip()
      parent = str(line_list[1]).strip()
      rank = str(line_list[2]).strip()
      
      ranky[child] = rank
      
      if parent in parent2children:
         parent2children[parent].append(child)
      else:
         parent2children[parent] = [child]


   if id in parent2children:
      children = parent2children[id]
   else:
      return([id],ranky,parent2children)
   for chi in children:
      if chi in parent2children:
         children.extend(parent2children[chi])
   children.append(id) 

   end_points = []
   for chi in children:
      if ranky[chi] == 'species':
         end_points.append(chi)
         if chi in parent2children:
            end_points.extend(parent2children[chi])
         
   if ranky[id] == 'species':
      end_points.append(id)
   return(set(end_points),ranky,parent2children)

def tax2gilist(taxids,gi_tax_file):

   tax2gi_dict = {}
   file = open(gi_tax_file,'rU')
   for line in file:
      line_list = line.split()
      tax = line_list[1]
      gi = line_list[0]
      if tax not in taxids:
         continue
      if tax in tax2gi_dict:
         tax2gi_dict[tax].append(gi)
      else:
         tax2gi_dict[tax] = [gi]

   for tax in taxids:
      if tax not in tax2gi_dict:
         tax2gi_dict[tax] = []
   return(tax2gi_dict)


def taxfile2dict(filename):
   import sys,re
   file = open(filename,'rU')
   size_dict = {}
   dict = {}
   for line in file:
      match = re.search('scientific name',line)
      if match:
         ary = line.split('|')
         dict[ary[0].strip()] = ary[1].strip()   
   return(dict) 
    
if __name__ == '__main__':
  main()
