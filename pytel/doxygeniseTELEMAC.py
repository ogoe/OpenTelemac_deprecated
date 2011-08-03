"""@brief
"""
"""@author Sebastien E. Bourban and Noemie Durand
"""
# _____          ___________________________________________________
# ____/ Imports /__________________________________________________/
#
from config import OptionParser,parseConfigFile, parseConfig_DoxygenTELEMAC
from parserFortran import scanSources, parseTemplateWrap, parseVars
from os import environ, path
from utils import getFileContent, putFileContent, createDirectories
import sys
import re

# _____                             ________________________________
# ____/ Global Regular Expressions /_______________________________/
#
f77continu2 = re.compile(r'(\s{5}\S\s*)(?P<line>.*)',re.I)

# _____                  ___________________________________________
# ____/ General Toolbox /__________________________________________/
#

def createDOXYGEN(ilines,lname,icount,list):

   doxy_cont_char = '!>'

   # ~~ Parsers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   olines = []
   for name,doxy,docs,body,vars in parseTemplateWrap(ilines,icount):
      if name=='PLANTE':
         print '\n'.join(ilines[0:3])
         print name, lname
         print list.keys()
         print list['bief'].keys()
      who = list[lname][name]

      #writes some of the doxygen tags documented in code

      # ~~ Function and User Defined Information
      #file = []; hist = []; fcts = []; bugs = []; warn = []; note = []; refs = []; code = []; para = []; resu = []
      for info in ['file','fcts','note','warn','bugs','refs','code'] :
         if doxy[info] != [] :
            olines.extend(['C'+71*'~'+'\n\n'])
            for i in range(len(doxy[info][0])) :
               if doxy[info][0][i] <> doxy_cont_char+'\n' :
                  if i <> len(doxy[info][0])-1 and doxy[info][0][i+1] == doxy_cont_char+'\n' :
                     doxy[info][0][i] = doxy[info][0][i][:len(doxy[info][0][i])-1]+'<br>\n'
                  olines.append(doxy[info][0][i])
            olines.extend(['\n'])
      olines.extend(['C'+71*'~'+'\n\n'])

      # ~~ Final Uses ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if who['uses'] != {}:
         line = '!>  @par Use(s)\n!><br>' + ', '.join(sorted(who['uses'].keys()))
         olines.extend([line])
         olines.extend(['\n'])

      # ~~ Final Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if who['args'] != [] or who['vars']['use'] != {} or who['vars']['cmn'] != [] or who['vars']['dec'] != [] or who['vars']['als'] != []:
         line = '!>  @par Variable(s)\n!>  <br><table>'
         olines.extend([line])

         # ~~ Arguments
         if who['args'] != []:
            line = '!>     <tr><th> Argument(s)\n!>    </th><td> ' + ', '.join(sorted(who['args'])) + '\n!>   </td></tr>'
            olines.extend([line])
         # ~~ Uses
         if who['vars']['use'] != {}:
            line = '!>     <tr><th> Use(s)\n!>    </th><td>'
            for u in who['vars']['use']:
               uv = []
               for v in who['vars']['use'][u]:
                  uv.append('\n!> @link ' + u + '::' + v + ' ' + v + '@endlink')
               line = line + '<hr>\n!> ' + u + ' :<br>' + ', '.join(sorted(uv))
            line = line.replace('<td><hr>\n','<td>\n') + '\n!>   </td></tr>'
            olines.extend([line])

         # ~~ Commons
         if who['vars']['cmn'] != []:
            line = '!>     <tr><th> Common(s)\n!>    </th><td>'
            for u in who['vars']['cmn']:
               line = line + '<hr>\n!> ' + u[0] + ' : ' + ', '.join(u[1])
            line = line.replace('<td><hr>\n','<td>\n') + '\n!>   </td></tr>'
            olines.extend([line])

         # ~~ Declars
         if who['vars']['dec'] != []:
            line = '!>     <tr><th> Internal(s)\n!>    </th><td> ' + ', '.join(sorted(who['vars']['dec'])) + '\n!>   </td></tr>'
            olines.extend([line])

         # ~~ Externals
         #if who['vars']['xtn'] != []:
         #   line = '!>     <tr><th> External(s)\n!>    </th><td> ' + ', '.join(sorted(who['vars']['xtn'])) + '\n!>   </td></tr>'
         #   olines.extend([line])

         # ~~ Aliases
         if who['vars']['als'] != []:
            line = '!>     <tr><th> Alias(es)\n!>    </th><td> ' + ', '.join(sorted(who['vars']['als'])) + '\n!>   </td></tr>'
            olines.extend([line])

         line = '!>     </table>\n'
         olines.extend([line])
         olines.extend(['\n'])

      # ~~ Final Calls ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if who['calls'] != {} or who['functions'] != []:
         line = '!>  @par Call(s)\n!>  <br><table>'
         olines.extend([line])
         if who['calls'] != {}:
            line = '!>     <tr><th> Known(s)\n!>    </th><td> ' + '(), '.join(sorted(who['calls'].keys())) + '()\n!>   </td></tr>'
            olines.extend([line])
         if who['functions'] != []:
            line = '!>     <tr><th> Unknown(s)\n!>    </th><td> ' + ', '.join(sorted(who['functions'])) + '\n!>   </td></tr>'
            olines.extend([line])
         line = '!>     </table>\n'
         olines.extend([line])
         olines.extend(['\n'])

      # ~~ Final Called ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if who['called'] != []:
         line = '!>  @par Called by\n!><br>' + '(), '.join(sorted(who['called'])) + '()\n'
         olines.extend([line])
         olines.extend(['\n'])

      # ~~ Final History ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if doxy['hist'] != []:
         olines.extend(['C'+71*'~'+'\n\n'])
         olines.append('!>  @par Development history')
         olines.append('!>  <br><table>')
         olines.append('!> <tr><th> Release </th><th> Date </th><th> Author </th><th> Notes </th></tr>')
         
         for i in range(len(doxy['hist'])-1,-1,-1) :    #writes them back from most to least recent
            name = ''; date = ''; release = ''; note = ''
            
            name = doxy['hist'][i][0].split('!>  @history')
            if name[0] == '' and len(name) > 1 : name.pop(0)
            name = name[0].strip('\n').strip()   #name can be empty

            date = doxy['hist'][i][1].strip(doxy_cont_char).strip().strip('\n')   #date can be empty
            
            release = doxy['hist'][i][2].strip(doxy_cont_char).strip().strip('\n')   #release can be empty
            if release != '' :
               k = 1
               if release[0].lower() <> 'v' : k = 0
               release = 'v'+release[k]+'p'+release[k+2]   #release is reformatted as v*p*

            for k in range(3,len(doxy['hist'][i]),1) :   #note can be empty
               if doxy['hist'][i][k].strip() != doxy_cont_char :
                  note = note + doxy['hist'][i][k].strip(doxy_cont_char).strip() + ' '
               else :
                  note = note + '<br>'

            olines.append('!>   <tr><td><center> '+release+'     </center>')
            olines.append('!>    </td><td> '+date)
            olines.append('!>    </td><td> '+name)
            olines.append('!>    </td><td> '+note)
            olines.append('!>   </td></tr>')
               
         olines.extend(['!>  </table>\n\n'])

      # ~~ Extended Description
      olines.extend(['C'+71*'~'+'\n\n'])
      if vars.keys() != '' or who['args'] != '':
         for a in who['args']:
            found = False
            for v in vars.keys():
               if a in parseVars(v): found = True
            if not found: vars.update({a:['',['']]})
         olines.extend(['!>  @par Details of primary variable(s)\n!>  <br><table>'])
         olines.extend(['!>\n!>     <tr><th>Name(s)</th><th>(in-out)</th><th>Description</th></tr>'])
         for v in sorted(vars.keys()):
            ino = ['-','-','-']
            if '<' in vars[v][0]: ino[0] = '<'
            if '>' in vars[v][0]: ino[2] = '>'
            line = '!>          <tr><td>' + v + '\n!></td><td>' + ''.join(ino) + '</td><td>' + '\n!>                  '.join(vars[v][1]) + '\n!>    </td></tr>'
            olines.extend([line])
         olines.extend(['!>     </table>'])

      # ~~ Program Name
      olines.extend(['C\nC'+71*'#'+'\nC'])
      for count in range(len(docs)):
         proc = re.match(f77continu2,docs[count])
         if proc: docs[count] = '     &'+docs[count][6:]
      olines.extend(docs)

      olines.extend(['C\n!'+71*'~'])
      # ~~ Variables in Argument
      for v in sorted(vars.keys()):
         name = v + '              '
         ino = ['-','-','-']
         if '<' in vars[v][0]: ino[0] = '<'
         if '>' in vars[v][0]: ino[2] = '>'
         line = '!| ' + name[0:15] + '|' + ''.join(ino) + '| ' + '\n!|                |   | '.join(vars[v][1])
         olines.extend([line])
      olines.extend(['!'+71*'~'+'\n'])
      olines.extend(['!\n'])

      # ~~ Main Body of Source Code
      for count in range(len(body)):
         proc = re.match(f77continu2,body[count])
         if proc: body[count] = '     &'+body[count][6:]   #think has been taken care of by oneoffTELEMAC.py already
      olines.extend(body)

   return olines

# _____             ________________________________________________
# ____/ MAIN CALL  /_______________________________________________/
#

__author__="Sebastien Bourban; Noemie Durand"
__date__ ="$19-Jul-2010 08:51:29$"

if __name__ == "__main__":
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

   debug = False

   # ~~ Scans all source files to build a relation database ~~
   print '\n\nScanning the source code\n\
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n'
   fic,mdl,sbt,fct,prg,dep,all = scanSources(cfgname,cfg,False)

   # ~~ Scans all source files to update Doxygen ~~~~~~~~~~~~~~~~
   for mod in fic.keys():
      print '\nCreating the DOXYGEN Headers for ' + mod + '\n\
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      for ifile in fic[mod].keys():

         # ~~ Read the content of the source file ~~~~~~~~~~~~
         ilines = getFileContent(ifile)
         # ~~ Update its Doxygen content ~~~~~~~~~~~~~~~~~~~~~
         olines = createDOXYGEN(ilines,mod,len(fic[mod][ifile]),all)
         # ~~ Make sure the distination exists ~~~~~~~~~~~~~~~
         ofile = ifile.replace(cfg['TELDIR'],path.join(cfg['TELDIR'],cfgname))
         createDirectories(path.dirname(ofile))
         # ~~ Write the content of the source file ~~~~~~~~~~~
         putFileContent(ofile,olines)

   sys.exit()

