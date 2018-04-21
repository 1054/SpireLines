#!/usr/bin/env python
# 
# Query the Herschel Science Archive (HSA) by RA Dec
# 



# See HSA documentation:
# http://archives.esac.esa.int/hsa/archive-help/cl/
# see "TAP resources"






# Import python packages

import os, sys, time, numpy, astropy

if sys.version_info < (3, 0):
    import httplib
    import urllib
    from urllib import urlencode
else:
    import http.client as httplib # Python 3
    import urllib.parse # Python 3
    from urllib.parse import urlencode

from xml.dom.minidom import parseString

host = "archives.esac.esa.int"
port = 80
pathinfo = "/hsa/whsa-tap-server/tap/async"



# Read user inputs

input_ra = numpy.nan
input_dec = numpy.nan
input_radius = numpy.nan
input_instrument = []
output_table = ''
output_format = 'csv' # 'votable' 'json'

i = 1
while i < len(sys.argv):
    str_arg = sys.argv[i].lower().replace('--','-').replace('_','-')
    if str_arg == "-out" or str_arg == "-output":
        if i+1 < len(sys.argv):
            i = i+1
            output_table = sys.argv[i]
    elif str_arg == "-instr" or str_arg == "-instrument" or str_arg == "-instrument-id":
        if i+1 < len(sys.argv):
            i = i+1
            input_instrument.append(int(sys.argv[i]))
    elif str_arg == "-spire-photo":
        input_instrument.append((2,4))
        input_instrument.append((2,5))
        input_instrument.append((2,6))
        input_instrument.append((2,7))
    else:
        if numpy.isnan(input_ra):
            input_ra = float(sys.argv[i])
        elif numpy.isnan(input_dec):
            input_dec = float(sys.argv[i])
        elif numpy.isnan(input_radius):
            if float(sys.argv[i]) > 0:
                input_radius = float(sys.argv[i]) / 60.0 # convert input arcmin to degree
    i = i+1

if numpy.isnan(input_ra) or numpy.isnan(input_dec):
    print('Usage: ')
    print('    HSA_tool_TAP_query_by_radec.py RA Dec [Radius_arcmin]')
    print('')
    sys.exit()

if numpy.isnan(input_radius):
    input_radius = 5.0 / 60.0 # 5 arcmin search radius in default

if output_table == '':
    output_table = "HSA_TAP_queried_table.%s"%(output_format)



# Example query selection
# SELECT DISTINCT( bib_code ) ,hsa.publication.bib_code, hsa.publication.title, hsa.publication.author 
# FROM hsa.publication INNER JOIN hsa.v_active_observation ON hsa.publication.observation_id=hsa.v_active_observation.observation_id 
# WHERE ( (hsa.v_active_observation.instrument_oid = '2') ) 
#         AND ( (hsa.v_active_observation.observing_mode_oid = '4') 
#                OR (hsa.v_active_observation.observing_mode_oid = '5') 
#                OR (hsa.v_active_observation.observing_mode_oid = '6') 
#                OR (hsa.v_active_observation.observing_mode_oid = '7') 
#                OR (hsa.v_active_observation.observing_mode_oid = '283') 
#              ) 
#         AND (hsa.v_active_observation.status != 'FAILED' OR hsa.v_active_observation.status IS NULL ) 
#         AND (((hsa.v_active_observation.polygon_fov IS NOT NULL 
#         AND INTERSECTS(CIRCLE('ICRS',308.71805,60.1536778,0.08333333333333333),hsa.v_active_observation.polygon_fov)=1) 
#         OR (hsa.v_active_observation.fov IS NOT NULL AND INTERSECTS(CIRCLE('ICRS',308.71805,60.1536778,0.08333333333333333),hsa.v_active_observation.fov)=1))) 
# ORDER BY hsa.publication.bib_code DESC



# WHAT IS IN THE 'hsa.v_active_observation'
# http://satscm.esac.esa.int/trac/HSA/wiki/HUIcurrent
# hsa.v_active_observation table view
#   observation_oid
#   observation_id
#   ra
#   dec
#   pa
#   od_number
#   target_name
#   start_time
#   end_time
#   duration
#   observer
#   prop_end
#   observing_mode_oid
#   instrument_oid



# BASICS ABOUT ADQL (Astronomical Data Query Language)
# https://www.asterics2020.eu/dokuwiki/lib/exe/fetch.php?media=open:wp4:school3:adql.pdf
#   String concatenation is done using the || operator. Strings also support LIKE that supports
#     patterns. % is “zero or more arbitrary characters”, “exactly one arbitrary character” (like
#     * and ? in shell patterns).
#   





# Prepare query selection

query_selection = "SELECT DISTANCE(POINT('ICRS',ra,dec), POINT('ICRS',%s,%s)) AS dist, * FROM hsa.v_active_observation"%(input_ra,input_dec)

query_condition = "WHERE ( "
query_condition = query_condition + "(1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',%s,%s,%s)))"%(input_ra, input_dec, input_radius)
if len(input_instrument)>0:
    query_condition = query_condition + " AND ( "
    for i in range(len(input_instrument)):
        if i>0: query_condition = query_condition + " OR "
        query_condition = query_condition + "('%d'=instrument_oid AND '%d'=observing_mode_oid)"%(input_instrument[i][0], input_instrument[i][1])
    query_condition = query_condition + " )"
if True:
    query_condition = query_condition + " AND ( "
    query_condition = query_condition + "aor NOT LIKE 'Calibration%%'"
    query_condition = query_condition + " )"
query_condition = query_condition + " )"

query_order = "ORDER BY dist ASC"


# Prepare query param

params = urlencode({\
    "REQUEST": "doQuery", \
    "LANG":    "ADQL", \
    "FORMAT":  output_format, \
    "PHASE":  "RUN", \
    "JOBNAME":  "Any name (optional)", \
    "JOBDESCRIPTION":  "Any description (optional)", \
    "QUERY":  query_selection + " " + query_condition + " " + query_order
    })

headers = {\
    "Content-type": "application/x-www-form-urlencoded", \
    "Accept":       "text/plain" \
    }

connection = httplib.HTTPConnection(host, port)
connection.request("POST",pathinfo,params,headers)

#Status
response = connection.getresponse()
print("Status: " +str(response.status), "Reason: " + str(response.reason))

#Server job location (URL)
location = response.getheader("location")
print("Location: " + location)

#Jobid
jobid = location[location.rfind('/')+1:]
print("Job id: " + jobid)

connection.close()

#-------------------------------------
#Check job status, wait until finished

while True:
    connection = httplib.HTTPConnection(host, port)
    connection.request("GET",pathinfo+"/"+jobid)
    response = connection.getresponse()
    data = response.read()
    #XML response: parse it to obtain the current status
    dom = parseString(data)
    phaseElement = dom.getElementsByTagName('uws:phase')[0]
    phaseValueElement = phaseElement.firstChild
    phase = phaseValueElement.toxml()
    print("Status: " + phase)
    #Check finished
    if phase == 'COMPLETED': break
    #wait and repeat
    time.sleep(0.2)

#print("Data:")
#print(data)

connection.close()

#-------------------------------------
#Get results
connection = httplib.HTTPConnection(host, port)
connection.request("GET",pathinfo+"/"+jobid+"/results/result")
response = connection.getresponse()
data = response.read()
outputFileName = output_table
outputFile = open(outputFileName, "wb")
outputFile.write(data)
outputFile.close()
connection.close()

print("Data saved in: " + outputFileName)






