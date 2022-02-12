# Takes filtered chimeric junctions aggregated accross samples input file and
# output site lists
# which can then be filtered to generate integration window bed file

# aggregate chimeric junctions file have the following columns:
#     1. non-shRNA chromosome in the chimeric junction
#     2. coordinate on chromosome in col 1.
#     3. read coverage aggregated across all replicates in which this 
#     non-shRNA+coordinate combination was found

# Many integration site estimates from the chimeric junctions file fall within a 
# few bps of each other. These are collapsed to a single coordinate within 100 bp
# with the highest read depth coordinate as representative

class site:
	def __init__(self,chrom,coordinate,count):
		self.chrom=chrom
		self.coordinate=int(coordinate)
		self.count=int(count)
	
	def get_distance(self,site2):
		if self.chrom!=site2.chrom:
			return 0
		else:
			return abs(self.coordinate-site2.coordinate)
	
	def __str__(self):
		return "Site --> chr: {0}, coordinate: {1}, count: {2}".format(self.chrom,self.coordinate,self.count)

class sitelist:
	def __init__(self):
		self.sites=[] #list of sites
		self.chrom=""
		self.bestcoordinate=0 # coordinate with the most count
		self.totalcount=0
	
	def add_site(self,site):
		self.sites.append(site)
		self.update_best()
	
	def update_best(self):
		if self.chrom=="":
			self.chrom=self.sites[0].chrom
		allcounts=list(map(lambda x:x.count,self.sites))
		max_count=max(allcounts)
		index_max_count=allcounts.index(max_count)
		self.bestcoordinate=self.sites[index_max_count].coordinate
		self.totalcount=sum(allcounts)
		# print(allcounts,max_count,index_max_count,self.bestcoordinate)
	
	def get_distance(self,site):
		if self.chrom!=site.chrom:
			return 0
		else:
			return abs(self.bestcoordinate-site.coordinate)
	
	def __str__(self):
		return_string = "SiteList --> chr: {0}, bestcoordinate: {1}, nelements: {2}, totalcount: {3}".format(self.chrom,self.bestcoordinate,len(self.sites),self.totalcount)
		return_string += "\n"
		for s in self.sites:
			site_return_string = str(s)
			return_string += site_return_string
			return_string += "\n"
		return return_string

# load filtered.chimeric.junctions.aggregate_across_samples
import sys
# filtered chimeric junctions aggregated accross samples input file
data=open(sys.argv[1]).readlines()

list_of_sitelists=[]
nelements=0
for d in data:
	l=d.strip().split()
	chrom=l[0]
	coordinate=l[1]
	count=l[2]
	nelements+=1
	s=site(chrom,coordinate,count)
	if nelements==1:
		newsl=sitelist()
		newsl.add_site(s)
		list_of_sitelists.append(newsl)
		continue
	for sli,sl in enumerate(list_of_sitelists):
		dist=sl.get_distance(s)
		if dist !=0 and dist < 100:
			sl.add_site(s)
			break
	else:
		newsl=sitelist()
		newsl.add_site(s)
		list_of_sitelists.append(newsl)

for sl in list_of_sitelists:
	print(sl)

