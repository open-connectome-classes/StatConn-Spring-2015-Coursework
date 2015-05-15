import os
import urllib2
#url = "https://www.dropbox.com/sh/idt3d0gylplyo31/AAB9GMfg30XjCGtrfdku5inFa/wormatlas.graphml?dl=0"
url="https://www.dropbox.com/sh/idt3d0gylplyo31/AABDYXDGy2Q0bSquj8uYKdjka/michaelnorris_wormatlas.graphml?dl=1"

#os.system('wget 
response = urllib2.urlopen(url)
file = response.read()
f = open('data/wormatlas_downloaded.graphml', 'w')
f.write(file)
f.close()
print "downloaded dataset"
