import sys
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIXML

SHARED_THRES = 2

def add_blast_entry(read, species, score):
    if species not in read:
        read[species] = [0, 0]
    read[species][0] += 1
    if read[species][1] < score:
        read[species][1] = score 

def parse_blast_xml(f):
    blast_recs = NCBIXML.parse(f)
    reads = {}
    try:
        for rec in blast_recs:
            for aln in rec.alignments:
                if rec.query not in reads:
                    reads[rec.query] = {}
                title = ' '.join(aln.title.split(' ')[1:3])
                print title
                for hsp in aln.hsps:
                    add_blast_entry(reads[rec.query], title, float(hsp.score))
    except:
        print "Given XML contain errors, some sequences can be missed"
    return reads

def parse_blast_txt(f):
    reads = {}
    for l in f:
        if l[0] == '#':
            continue
        fields = l.strip().split('\t')
        if len(fields) != 12:
            continue
        if fields[0] not in reads:
            reads[fields[0]] = {}
        add_blast_entry(reads[fields[0]], fields[1], float(fields[-1]))
    return reads

def find_species(reads):
    species = {}
    ids = set()
    for read in reads:
        for subj in reads[read].keys():
            id_fields = subj.split('|')
            ids.add(id_fields[1])
            if id_fields[3] not in species:
                species[id_fields[3]] = []
            species[id_fields[3]].append([read, reads[read][subj]])
    return species, ids

def add_reads(entry, reads):
    for r1 in entry:
        for i, r2 in enumerate(reads):
#            print r1[0] + ' ' + r2[0]
            if r1[0] == r2[0]:
#                print 'eq'
                r1[1] = [max([r1[1][0], r2[1][0]]), max([r1[1][1], r2[1][1]])]
                del reads[i]
                break
    entry += reads
    return entry

def print_species(species):
    for key in sorted(species.keys(), key = lambda k: len(species[k]), reverse=True):
        print key + "\t" + str(len(species[key]))

def define_species(species_id, ids):
    species = {}
    print "Number of different gid: " + str(len(species_id.keys()))
    for i in xrange(0, len(ids), 150):
        Entrez.email = "A.N.Other@example.com"
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", \
                               id=','.join(list(ids)[i:i+150]))
        for subj in SeqIO.parse(handle, 'gb'):
            descr = ' '.join(subj.description.split(' ')[:2])
            if descr not in species:
                species[descr] = []
            species[descr] = add_reads(species[descr], species_id[subj.id])
    return species

def check_and_cluster(species, reads, cluster):
    num_shared = 0
    dist = []
#    print species
    for r2 in reads:
        if next((r for r in cluster[0] if r == r2[0]), None):
            num_shared += 1
        else:
            dist.append(r2[0])
    if (num_shared != 0 and ((len(reads) >= len(cluster[0]) and len(cluster[0]) - num_shared <= SHARED_THRES)
            or (len(cluster[0]) > len(reads) and len(reads) - num_shared <= SHARED_THRES))):
#        print num_shared
        cluster[0] += dist
        cluster[1] += [[species, len(reads)]]
#        print cluster
        return True
    return False

def cluster_species(species):
    clusters = []
    for key in species:
        added = False
        for c in clusters:
            added = check_and_cluster(key, species[key], c)
            if added:
                break
        if not added:
            clusters.append([[r[0] for r in species[key]], [[key, len(species[key])]]])

    return clusters

def print_clusters(clusters):
    for c in clusters:
        print '' 
        for s in sorted(c[1], key = lambda sp: sp[1], reverse = True):
            print s[0] + '\t' + str(s[1])

if __name__ == '__main__':

    f = open(sys.argv[1])
    
    reads = parse_blast_txt(f)
#    print reads
    species, ids = find_species(reads) 
    species = define_species(species, ids) 
#    print_species(species)
    clusters = cluster_species(species)
    print "Clusters:"
    print_clusters(clusters)
