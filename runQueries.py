# -*- coding: utf-8 -*-
#
# Copyright (C) 2013 Google Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Command-line skeleton application for BigQuery API.
Usage:
  $ python sample.py

You can also get help on all the command-line flags the program understands
by running:

  $ python sample.py --help

"""

import argparse
import httplib2
import os
import sys

from apiclient import discovery
from oauth2client import file
from oauth2client import client
from oauth2client import tools

# Parser for command-line arguments.
parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    parents=[tools.argparser])


# CLIENT_SECRETS is name of a file containing the OAuth 2.0 information for this
# application, including client_id and client_secret. You can see the Client ID
# and Client secret on the APIs page in the Cloud Console:
# <https://cloud.google.com/console#/project/912363036792/apiui>
CLIENT_SECRETS = os.path.join(os.path.dirname(__file__), 'client_secrets.json')

# Set up a Flow object to be used for authentication.
# Add one or more of the following scopes. PLEASE ONLY ADD THE SCOPES YOU
# NEED. For more information on using scopes please see
# <https://developers.google.com/+/best-practices>.
FLOW = client.flow_from_clientsecrets(CLIENT_SECRETS,
  scope=[
      'https://www.googleapis.com/auth/bigquery',
      'https://www.googleapis.com/auth/cloud-platform',
      'https://www.googleapis.com/auth/devstorage.full_control',
      'https://www.googleapis.com/auth/devstorage.read_only',
      'https://www.googleapis.com/auth/devstorage.read_write',
    ],
    message=tools.message_if_missing(CLIENT_SECRETS))


def main(argv):
  # Parse the command-line flags.
  flags = parser.parse_args(argv[1:])
  

  # If the credentials don't exist or are invalid run through the native client
  # flow. The Storage object will ensure that if successful the good
  # credentials will get written back to the file.
  storage = file.Storage('sample.dat')
  credentials = storage.get()
  if credentials is None or credentials.invalid:
    credentials = tools.run_flow(FLOW, storage, flags)

  # Create an httplib2.Http object to handle our HTTP requests and authorize it
  # with our good Credentials.
  http = httplib2.Http()
  http = credentials.authorize(http)

  # Construct the service object for the interacting with the BigQuery API.
  service = discovery.build('bigquery', 'v2', http=http)

  try:
    print "Successfully connected to service."
    
    # Read in file with queries to do:
    filey = fopen('output/anxiety_Feb-04-2014-01-58-33_samples.tab','r')
    keepers = []
    for f in filey:
      keepers.append(f.strip('\n').split('\t'))

    colnames = keepers.pop(0)
    
    print "Running queries for " + str(len(keepers)) + " samples..."

    # First get probe ids for each table? (or comes in result?)
    for k in keepers:
      subid = k[0]
      structure = k[1]
      x = k[11]; y=k[12]; z=k[13]
      # For each of our 12 PACall tables:
      for t in range(0,12): 
        bqCommand = "bq query --destination_table=allen_brain_result." + searchTerm + str(t) + " --append_table \"SELECT * FROM allen_brain_human.pacall" + str(t) + " WHERE specimen_id = " + subid +" AND structure_id = " + structure + " AND mni_x = " + x + " AND mni_y = " + y + " AND mni_z = " + z + "\""
        query_data = {'query':bqCommand}
	query_request = service.jobs()

	# Make a call to the BigQuery API
	query_response = query_request.query(projectId=912363036792,body=query_data).execute()
        print query_response

  except client.AccessTokenRefreshError:
    print ("The credentials have been revoked or expired, please re-run"
      "the application to re-authorize")


# For more information on the BigQuery API you can visit:
#
#   https://developers.google.com/bigquery/docs/overview
#
# For more information on the BigQuery API Python library surface you
# can visit:
#
#   https://developers.google.com/resources/api-libraries/documentation/bigquery/v2/python/latest/
#
# For information on the Python Client Library visit:
#
#   https://developers.google.com/api-client-library/python/start/get_started
if __name__ == '__main__':
  main(sys.argv)
