"""@brief
"""
"""@author Sebastien E. Bourban and Noemie Durand
"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
from config import OptionParser,parseConfigFile, parseConfig_DoxygenTELEMAC
from parserFortran import scanSources, parseDoxyWrap, parseArgs
from os import environ, path, rename, sep
from utils import getFileContent, putFileContent
import sys
import re

# _____                             ________________________________
# ____/ Global Regular Expressions /_______________________________/
#
f77comment = re.compile(r'[C!#*]',re.I)
f77continu2 = re.compile(r'(\s{5}\S\s*)(?P<line>.*)',re.I)
f90types = '(CHARACTER|LOGICAL|INTEGER|REAL|COMPLEX|DOUBLE\s*(PRECISION\s*(COMPLEX|)|COMPLEX))\s*?(\**\s*?\d+|\**\(.*?\))?|TYPE\s*\([\w\s,=(*)]*?\)'
typ_args = re.compile(r'\s*?(.*?::)?\s*?(?P<vars>.*?)\s*\Z',re.I)
typ_title = re.compile(r'\s*?(?P<type>(%s))\s*?(?P<after>.*?)\s*\Z'%(f90types),re.I)
brief = re.compile(r'\s*(%s)\s*\S*\s*?(?P<after>.*?)\s*\Z'%('!> @brief'),re.I)

# _____                  ___________________________________________
# ____/ General Toolbox /__________________________________________/
#

def readDOXYGEN(ilines,icount,lname,list,move):

   # ~~ Parsers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   olines = [] 
   doxy_char = '!>  @'; doxy_cont_char = '!>'
   temp_char = '!'    ; cont_char = '!+'
   htmldico = ['<table>','</table>','<tr>','</tr>','<th>','</th>','<td>','</td>','</center>']

   for name,doxy,docs,body,vars in parseDoxyWrap(ilines,icount):
      history = []
      who = list[lname][name]

      #reformats development history
      if doxy['hist'] != []:
         for i in range(len(doxy['hist'])) :
            for format in htmldico :
               doxy['hist'][i] = doxy['hist'][i].replace(format,'')
         i = 0
         while i < len(doxy['hist']) :
            if '<center>' in doxy['hist'][i] :
               release = doxy['hist'][i].strip(doxy_cont_char).strip().strip('<center>').strip()
               if len(release.split('.')) > 1 :
                  release = ' V' + release.split('.')[0] + 'P' + release.split('.')[1]
               date = doxy['hist'][i+1].strip(doxy_cont_char)
               name = doxy['hist'][i+2].strip(doxy_cont_char)
               for j in [0,1,2,3,4,5,6,7,8,9] :   #deletes the telephone number
                  name = name.replace(str(j),'')
               for j in ['....','...','----','---'] :
                  name = name.replace(j,'')
               note = doxy['hist'][i+3].strip(doxy_cont_char)
               i = i+4
               history.append([name,date,release,note])
            else :
               i = i+1

      #writes out docs (with some formatting)                          e.g.
      olines.extend(['!'+20*' '+len(docs[0].strip())*'*'+'\n'])       #!                    *****************
      olines.extend([' '+20*' '+docs[0].strip()+'\n'])                #                     SUBROUTINE CORFON
      olines.extend(['!'+20*' '+len(docs[0].strip())*'*'+'\n'])       #!                    *****************
      olines.extend(['!'])                                            #!
      for i in range(1,len(docs),1) :
         olines.extend([docs[i]])
      olines.extend(['!'])                                            #!
      olines.extend(['!'+71*'*'+'\n'])                                #!***********************************************************************
      line = '! '+lname.upper()
      if history != [] : line = line+'  '+history[0][2]+ \
                                 30*' '+history[0][1]+'\n'
      olines.append(line)                                             #! ARTEMIS   V6P0                                   21/08/2010
      olines.extend(['!'+71*'*'+'\n'+'!\n'])                          #!***********************************************************************

      #adds tags

      # ~~ Function and User Defined Information
      #file = []; hist = []; fcts = []; bugs = []; warn = []; note = []; refs = []; code = []; para = []; resu = []
      for info in ['file','fcts','note','warn','bugs','refs','code'] :
         if doxy[info] != [] :
            tmp = []
            doxy[info][0] = doxy[info][0].replace('<br>','\n'+cont_char)
            doxy[info][0] = doxy[info][0].replace(doxy_char, temp_char)
            tmp.extend(doxy[info][0].split(' ',1))    #formats the tags such that descriptions is 10char. from margin
            line = tmp[0]
            if len(tmp) > 1 :
               if info <> 'refs' :
                  line = line+(10-len(tmp[0]))*' '+tmp[1].strip()
               else :
                  line = line+2*' '+tmp[1].strip()
            olines.append(line)
            for i in range(1,len(doxy[info]),1) :
               doxy[info][i] = doxy[info][i].replace('<br>','\n'+cont_char)
               doxy[info][i] = doxy[info][i].replace(doxy_char, temp_char)
               doxy[info][i] = doxy[info][i].replace(doxy_cont_char, cont_char)
               if doxy[info][i] <> '!endcode' : olines.append(doxy[info][i])
            olines.extend([temp_char+'\n'])
  
      # ~~ Development History
      if history != []:
         #writes them back in reverse order and keeps empty lines when no information is supplied
         #earlier version had if "history[i][1].strip() <> '' :"  which only displayed relevant information
         #it got a bit tricky to doxygenise (recognise a date from a comment)
         #this will hopefully be covered by new script with pop-up window requiring the user to fill-in
         #metadata at the SVN commission stage
         for i in range(len(history)-1,-1,-1) :    #writes them back
            olines.extend(['!history'+2*' '+history[i][0].strip()])
            olines.extend([cont_char+8*' '+history[i][1].strip()])   #if history[i][1].strip() <> '' : 
            olines.extend([cont_char+8*' '+history[i][2].strip()])   #if history[i][2].strip() <> '' : 
            #if history[i][3].strip() <> '' :
            #lines limited to 69char
            history[i][3] = history[i][3].replace('<br>',' '+cont_char+' ') #+cont_char+' ')
            note = history[i][3].strip().split()
            while 1 :
               l = 0; j = 0
               while j <= len(note)-1 and l <= 65 and note[j] != cont_char :
                  l = l + len(note[j]) + 1
                  j = j+1
               if j != len(note) and note[j] != cont_char:   #l>65 broke the while loop
                  j = j-1
               line = cont_char + 3*' '
               for k in range(j) :
                  line = line + note[0] + ' '
                  note.pop(0)
               olines.append(line)
               if note != [] and note[0] == cont_char :      #note[j]=cont_char broke the while loop
                  note.pop(0)
                  olines.append(cont_char)
               if note == [] or (len(note) == 1 and note[0] == '') : break
            olines.extend(['!\n'])
      olines.extend(['!'+71*'~'+'\n'])

      #concatenates vars (with some formatting) and body to re-write fortran code
      #!!! if the variables in the in/out table are not in argument, move description
      #!!! to declarations_[code].f if in "use" (dico move), delete otherwise (local)

      # ~~ Arguments
      for k in vars.keys() :         #eg 'T1,2'
         found = True
         k_list = k.split(',')       #eg ['T1','2']
         for item in k_list :        #eg 'T1'
            if item not in who['args'] :
               found = False
               # ~~ Uses
               if who['vars']['use'] != {}:
                  for u in who['vars']['use']:
                     if item in who['vars']['use'][u]:
                        mfile = list[mdl[u][0]][u]['path']+sep+u.lower()
                        if mfile not in move.keys():
                           move.update({mfile:[]})
                        move[mfile].append([item,temp_char+' '+' '.join(vars[k][1])])
         if not found : del vars[k]

      for v in sorted(vars.keys()):
         name = v + 14*' '
         ino = ['-','-','-']
         if '<' in vars[v][0]: ino[0] = '<'
         if '>' in vars[v][0]: ino[2] = '>'
         line = '!| ' + name[0:15] + '|' + ''.join(ino) + '| ' + '\n!|                |   | '.join(vars[v][1])
         olines.extend([line])
      olines.extend(['!'+71*'~'+'\n'])
      olines.extend(['!\n'])

      # ~~ Main Body of Source Code
      count = 0
      while count < len(body) :
         proc = re.match(f77continu2,body[count])
         if proc: body[count] = '     &'+body[count][6:]   #sets all continuation characters to &
         proc = re.match(f77comment,body[count])
         if proc: body[count] = '!'+body[count][1:]        #sets all comment characters to !
         proc = re.match(brief,body[count])
         if proc: body[count] = '!brief'+proc.group('after') #replaces !> @brief with !brief in modules
         olines.append(body[count])
         count = count+1
                       
   return olines,move

def readTemplate(mfile):

   SrcF = open(mfile.lower()+'.f','r')
   ilines = SrcF.readlines()
   SrcF.close()
   # ~~ Type Declaration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   add = {}
   for i in range(len(ilines)) :
      proc = re.match(typ_title,ilines[i])
      if proc :
         proc = re.match(typ_args,proc.group('after'))
         if proc :
            if proc.group('vars') != None :
               list = parseArgs(proc.group('vars'))
               for j in range(len(move[mfile])) :
                  if move[mfile][j][0] in list :
                     add.update({i:[move[mfile][j][0],move[mfile][j][1]]})

   return ilines,add

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#

__author__="Sebastien Bourban; Noemie Durand"
__date__ ="$19-Jul-2010 08:51:29$"

if __name__ == "__main__":
   debug = False

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Reads config file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   print '\n\nLoading Options and Configurations\n\
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n'
   CFGNAME = ''
   SYSTELCFG = 'systel.cfg'
   if environ.has_key('SYSTELCFG'): SYSTELCFG = environ['SYSTELCFG']
   parser = OptionParser("usage: %prog [options] \nuse -h for more help.")
   parser.add_option("-c", "--configname",
                      type="string",
                      dest="configName",
                      default=CFGNAME,
                      help="specify configuration name, default is randomly found in the configuration file" )
   parser.add_option("-f", "--configfile",
                      type="string",
                      dest="configFile",
                      default=SYSTELCFG,
                      help="specify configuration file, default is systel.cfg" )
   options, args = parser.parse_args()
   if not path.isfile(options.configFile):
      print '\nNot able to get to the configuration file: ' + options.configFile + '\n'
      sys.exit()
   
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# ~~~~ Works for only one configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   cfgname = options.configName
   if options.configName == '':
      cfgname = parseConfigFile(options.configFile).keys()[0]
   if cfgname not in parseConfigFile(options.configFile).keys():
      print '\nNot able to get to find your configuration in the configuration file: ' + options.configFile + '\n'
      sys.exit()

   cfg = parseConfig_DoxygenTELEMAC(cfgname)[cfgname]

   # ~~ Scans all source files to build a relation database ~~
   print '\n\nScanning the source code\n\
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n'
   fic,mdl,sbt,fct,prg,dep,all = scanSources(cfgname,cfg)

   # ~~ Scans all source files to update Doxygen ~~~~~~~~~~~~~~~~
   move = {}
   for mod in fic.keys():
      print '\nWriting out ' + mod.upper() + ' source code following standrad template format\n\
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      for ifile in fic[mod].keys():

         # ~~ Reads the content of the source file ~~~~~~~~~~~
         ilines = getFileContent(ifile)
         # ~~ Updates its Doxygen content ~~~~~~~~~~~~~~~~~~~~
         olines,move = readDOXYGEN(ilines,len(fic[mod][ifile]),mod,all,move)
         # ~~ Renames the original source file as *.hrw
         rename(ifile,path.splitext(ifile)[0]+'.hrw')
         # ~~ Writes the content of the source file ~~~~~~~~~~
         ofile = ifile
         putFileContent(ofile,olines)

   # ~~ Adds the documentation to the relevant modules ~~~~~~~
   for mfile in move.keys():

      # ~~ Reads the content of the template source file ~~~~~
      ilines,add = readTemplate(mfile)
      # ~~ Updates the modules with extra documentation ~~~~~~
      ofile = open(mfile+'.hrw2','w')
      for i in range(len(ilines)) :
         if i in sorted(add) :
            ofile.write(add[i][1]+'\n')
         ofile.write(ilines[i])
      ofile.close()
      rename(mfile+'.hrw',mfile+'.hrw.old')
      rename(mfile+'.hrw2',mfile+'.hrw')

   sys.exit()

