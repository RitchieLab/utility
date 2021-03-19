#!/usr/bin/env python

import apsw
#import base64
import json
import re
import sys
import urllib
import urllib2


TRIBE_URL = 'https://tribe.greenelab.com'
TRIBE_API_AUTHORIZE = '/oauth2/authorize/'
TRIBE_API_TOKEN = '/oauth2/token/'
TRIBE_API_GENESET = '/api/v1/geneset/'
TRIBE_API_VERSION = '/api/v1/version/'
#TRIBE_CLIENT_ID = 'YutNMQu8gd18GHbDtuNIw3x3kB6tE1x5tUJ0MmRf'
#TRIBE_CLIENT_SECRET = 'j08PwhxishbkY8MQYoKdoHgRqimtKC6q37xKLcAPjd27tCrkIC7kvtwzNGgETRmcOklbUChbmIJ5A1t0pl4KytQeXoQc6oURq8lhE8MvLv5g0A67ZMfuetKEDnSRt368'
TRIBE_CLIENT_ID = '8BKOSPDJG4hc4aGA0Zo5VxVcED6eoaHTjM4ZmzUX'
TRIBE_CLIENT_SECRET = '1OSTjJU5EzzS9DR1kFl02AhRUtwJUH6RXq2iJSeU4qKrT3WL2YMpY91rNgagNOkYGbny7D8BqSvYi4l7Zd5UpS1F1pmfAbHBRGGuL6wspdaFjjX3kfeA9VKJwr1T8hkY'
TRIBE_USER_NAME = 'atf3@ucs.psu.edu'
TRIBE_USER_PASS = 'Ymsl4epyKzbW2C4LqHfwzGx8fb0rWmKzINy1ARDlZvYuFrqgvqkN406pOWu3Zgd4'


def DEBUG_REQUEST(r):
	if False: #TODO
		print r.get_method(), r.get_full_url()
		print "\n".join(": ".join(h) for h in r.header_items())
		print ""
		print r.get_data()
#DEBUG_REQUEST()


def DEBUG_RESPONSE(r):
	if False: #TODO
		print r.getcode(), r.geturl()
		print r.info()
		print ""
		data = r.read()
		print data
		return data
	else:
		return r.read()
#DEBUG_RESPONSE()


if len(sys.argv) < 2:
	print "Usage: %s <loki.db> [source]" % (sys.argv[0],)
	sys.exit(2)


print "Opening LOKI database ..."
conn = apsw.Connection(sys.argv[1])
cursorS = conn.cursor()
cursorU = conn.cursor()
cursorG = conn.cursor()
nsID = dict(cursorS.execute("SELECT namespace,namespace_id FROM `namespace`"))
taxID = max(int(row[0]) for row in cursorS.execute("SELECT value FROM setting WHERE setting='taxonomy'"))
if taxID == 3702:
	species = 'arabidopsis-thaliana'
elif taxID == 559292 or taxID == 4932:
	species = 'saccharomyces-cerevisiae'
elif taxID == 6239:
	species = 'caenorhabditis-elegans'
elif taxID == 7227:
	species = 'drosophila-melanogaster'
elif taxID == 7955:
	species = 'danio-rerio'
elif taxID == 9606:
	species = 'homo-sapiens'
elif taxID == 10090:
	species = 'mus-musculus'
elif taxID == 10116:
	species = 'rattus-norvegicus'
elif taxID == 208964:
	species = 'pseudomonas-aeruginosa'
else:
	print "... ERROR: unsupported taxonomy",taxID
	sys.exit(1)
print "... OK: %s (%d)" % (species,taxID)


print "Requesting TRIBE API token ..."
url = TRIBE_URL + TRIBE_API_TOKEN
headers = {
#	'authorization': 'Basic ' + base64.standard_b64encode(TRIBE_CLIENT_ID + ':' + TRIBE_CLIENT_SECRET)
}
data = urllib.urlencode({
	'client_id'     : TRIBE_CLIENT_ID,
	'client_secret' : TRIBE_CLIENT_SECRET,
	'grant_type'    : 'password',
	'username'      : TRIBE_USER_NAME,
	'password'      : TRIBE_USER_PASS,
})
request = urllib2.Request(url, data, headers)
#DEBUG_REQUEST(request)
try:
	response = urllib2.urlopen(request)
except urllib2.HTTPError as e:
	print e.headers
	raise
#data = DEBUG_RESPONSE(response)
data = response.read()
tribe_token = str(json.loads(data)['access_token'])
response.close()
print "... OK"


print "Uploading genesets ..."
url = TRIBE_URL + TRIBE_API_GENESET
headers = {
	'authorization' : 'OAuth ' + tribe_token,
	'content-type'  : 'application/json',
}
query = {
	'slug'     : '',
	'show_tip' : 'true',
}
geneset = {
	'organism'    : '/api/v1/organism/' + species,
	'public'      : False,
	'xrdb'        : 'Entrez',
	'description' : '',
	'title'       : '',
	'abstract'    : '',
	'annotations' : {},
}
version = {
	'geneset'     : '',
#	'public'      : False,
	'xrdb'        : 'Entrez',
	'description' : 'Via LOKI',
	'annotations' : {},
}
sourceNamespace = {
	'biogrid'  : 'biogrid_id',
	'go'       : 'go_id',
	'kegg'     : 'kegg_id',
	'mint'     : 'mint_id',
	'pfam'     : 'pfam_id',
	'pharmgkb' : 'pharmgkb_id',
	'reactome' : 'reactome_id',
}
sqlS = "SELECT source_id,source FROM `source`"
if len(sys.argv) > 2:
	sqlS += " WHERE source='%s'" % (sys.argv[2],)
for source_id,source_name in cursorS.execute(sqlS):
	if source_name not in sourceNamespace:
		print ("  ERROR: source '%s' has no default namespace" % (source_name,))
		continue
	sqlU = "SELECT group_id,MIN(name) AS slug,label,description FROM `group` JOIN `group_name` USING (group_id) WHERE `group`.source_id=? AND namespace_id=? GROUP BY group_id"
	sqlU += " LIMIT 1" #TODO
	for group_id,group_name,group_label,group_description in cursorU.execute(sqlU, (source_id, nsID[sourceNamespace[source_name]])):
		sqlG = "SELECT bn.name FROM `group_biopolymer` AS gb JOIN biopolymer_name AS bn ON bn.biopolymer_id = gb.biopolymer_id WHERE gb.group_id = ? AND bn.namespace_id = ?"
		genes = dict( (long(entrez_gid),list()) for entrez_gid, in cursorG.execute(sqlG, (group_id,nsID['entrez_gid'])) )
		if not genes:
			continue
		
		slug = re.sub('[^A-Za-z0-9_-]+', '-', (source_name + '-' + group_name + '-' + species)).strip('-').lower()
		print ("  %s ..." % (slug,)),
		
		query['slug'] = slug
		request = urllib2.Request(TRIBE_URL + TRIBE_API_GENESET + '?' + urllib.urlencode(query), None, headers)
		DEBUG_REQUEST(request)
		response = urllib2.urlopen(request)
		data = DEBUG_RESPONSE(response)
		response.close()
		uri = None
		tip = None
		for obj in json.loads(data)['objects']:
			if obj['creator']['username'] == "atf32":
				uri = obj['resource_uri']
				tip = obj['tip']['resource_uri']
		
		if uri and tip:
			print "updating ...",
			url = TRIBE_URL + TRIBE_API_VERSION
			version['geneset'] = uri
			version['parent'] = tip
			version['annotations'] = genes
			data = json.dumps(version)
		else:
			print "creating ...",
			url = TRIBE_URL + TRIBE_API_GENESET
			geneset['slug'] = slug
			geneset['title'] = '(%s) %s' % (source_name,group_label)
			geneset['description'] = group_name
			geneset['abstract'] = group_description or ""
			geneset['annotations'] = genes
			data = json.dumps(geneset)
		request = urllib2.Request(url, data, headers)
		DEBUG_REQUEST(request)
		response = urllib2.urlopen(request)
		data = DEBUG_RESPONSE(response)
		response.close()
		print " OK"
	#foreach rowU
#foreach rowS
print "... OK"
