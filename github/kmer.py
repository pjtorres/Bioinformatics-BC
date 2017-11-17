#practice on gut hub, script allows you to get a list of 3 kmers
string='AAAGGG'
n = len('AAAGGG')
k=3
kmers=[]
for i in range(0, n-k+1):
    kmers.append(string[i:i+k])
    print k
print kmers
   
