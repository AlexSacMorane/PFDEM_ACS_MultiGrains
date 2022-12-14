# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a report class.
This report is a .txt file with information from the simulation
"""
#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import os
import time

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Report:

#-------------------------------------------------------------------------------

    def __init__(self, Name, Datetime):
        '''
        Defining a report.

            Input :
                itself (a report)
                a name (a string)
                a start time (a datetime class)
            Output :
                Nothing but a file is generated (a file)
        '''
        if type(Name) == str:
            if Name[-4:] != '.txt':
                Name = Name + '.txt'
            self.name = Name
            self.datetimestart = str(Datetime)

            Last_modification = 0
            L_file = os.listdir()
            for file in L_file:
                if file[-3:] == '.py' and file !='main.py':
                    if os.path.getmtime(file) > Last_modification:
                        Last_modification =  os.path.getmtime(file)

            file_to_write = open(self.name,'w')
            file_to_write.write('Last compilation: '+str(time.ctime(Last_modification))+'\n')
            file_to_write.write('Simulation started '+str(Datetime)[:19]+'\n\n')
            file_to_write.close()
        else:
            print('Error')
        self.text_tempo = ''
        self.text_tempo_statut = False

#-------------------------------------------------------------------------------

    def write(self, Text):
        '''
        Write Text (a string) in the report.

            Input :
                itself (a report)
                a text (a string)
            Output :
                Nothing but the report is updated (a file)
        '''

        file_to_write = open(self.name,'a')
        file_to_write.write(Text)
        file_to_write.close()

#-------------------------------------------------------------------------------

    def write_and_print(self, Text_to_write, Text_to_print):
        '''
        Write Text_to_write (a string) in the report and print Text_to_print (a string)

            Input :
                itself (a report)
                    itself (a report)
                two texts (a string)
            Output :
                Nothing but the report is updated (a file) and something is printed (an user interface)
        '''
        file_to_write = open(self.name,'a')
        file_to_write.write(Text_to_write)
        file_to_write.close()
        print(Text_to_print)

#-------------------------------------------------------------------------------

    def tic_tempo(self, Datetime):
        '''
        Save a temporary start time.

        Work with tac_tempo() to compute a time cost for a simulation step.

            Input :
                itself (a report)
                a start time (a datetime class)
            Output :
                Nothing but the report gets an updated attribut (a datetime class)
        '''
        self.datetimestart_tempo = str(Datetime)

#-------------------------------------------------------------------------------

    def tac_tempo(self, Datetime, Step_name):
        '''
        Save a temporary start time.

        Work with tic_tempo() to compute a time cost for a simulation step.

            Input :
                itself (a report)
                a stop time (a datetime class)
            Output :
                Nothing but the report gets an updated attribut (a datetime class) and the report is updated (a file)
        '''

        self.datetimeend_tempo = str(Datetime)

        dyear = int(self.datetimeend_tempo[:4])-int(self.datetimestart_tempo[:4])
        dmonth = int(self.datetimeend_tempo[5:7])-int(self.datetimestart_tempo[5:7])
        dday = int(self.datetimeend_tempo[8:10])-int(self.datetimestart_tempo[8:10])
        dhour = int(self.datetimeend_tempo[11:13])-int(self.datetimestart_tempo[11:13])
        dmin = int(self.datetimeend_tempo[14:16])-int(self.datetimestart_tempo[14:16])
        dsec = int(self.datetimeend_tempo[17:19])-int(self.datetimestart_tempo[17:19])

        dt = dsec + dmin*60 + dhour*60*60 + dday*60*60*24 + dmonth*60*60*24*30.5 + dyear*60*60*24*30.5*365.25 #in sec

        dt_day = dt // (60*60*24) #in dday
        dt = dt % (60*60*24)
        dt_hour = dt // (60*60) #in hour
        dt = dt % (60*60)
        dt_min = dt // 60 # in min
        dt = dt % 60

        file_to_write = open(self.name,'a')
        file_to_write.write('Time spent during "'+Step_name+'" : '+str(dt_day)+' days '+str(dt_hour)+' hours '+str(dt_min)+' min '+str(dt)+' sec\n')
        file_to_write.close()

#-------------------------------------------------------------------------------

    def end(self, Datetime):
     '''
     Work with init() to compute the total time cost of the simulation.

         Input :
             itself (a report)
             a stop time (a datetime class)
         Output :
             Nothing but the file is updated (a file)
     '''
     self.datetimeend = str(Datetime)

     dyear = int(self.datetimeend[:4])-int(self.datetimestart[:4])
     dmonth = int(self.datetimeend[5:7])-int(self.datetimestart[5:7])
     dday = int(self.datetimeend[8:10])-int(self.datetimestart[8:10])
     dhour = int(self.datetimeend[11:13])-int(self.datetimestart[11:13])
     dmin = int(self.datetimeend[14:16])-int(self.datetimestart[14:16])
     dsec = int(self.datetimeend[17:19])-int(self.datetimestart[17:19])

     dt = dsec + dmin*60 + dhour*60*60 + dday*60*60*24 + dmonth*60*60*24*30.5 + dyear*60*60*24*30.5*365.25 #in sec

     dt_day = dt // (60*60*24) #in dday
     dt = dt % (60*60*24)
     dt_hour = dt // (60*60) #in hour
     dt = dt % (60*60)
     dt_min = dt // 60 # in min
     dt = dt % 60

     file_to_write = open(self.name,'a')
     file_to_write.write('\n'+'Simulation ended '+str(Datetime)[:19]+'\n')
     file_to_write.write('Time spent : '+str(dt_day)+' days '+str(dt_hour)+' hours '+str(dt_min)+' min '+str(dt)+' sec\n')
     file_to_write.close()

 #-------------------------------------------------------------------------------
