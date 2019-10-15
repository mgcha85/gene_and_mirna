import urllib.request

print('Beginning file download with urllib2...')

url = 'http://geneontology.org/ontology/go-basic.obo'
urllib.request.urlretrieve(url, '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Software/GO/goatools-master/data/go-basic.obo')
