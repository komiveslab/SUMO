### Import Commands
library(tools)
library(XML)
library(system)
### /Import Commands


# Create input folder in TPPmaindir
### Set Working Directory
tppmainDir="C:/Inetpub/wwwroot/ISB/data"
tppinput="input"
tppsubDir="2017-04-10"
dir.create(file.path(tppmainDir,tppsubDir))
dir.create(file.path(tppmainDir,tppinput))
setwd(file.path(tppmainDir,tppinput))
getwd()
### /Set Working Directory

#Place files in input directory

### Running post-MSGF+ Workflow
files<-list.files(pattern="\\.mzXML")
for (item in files){
	rawdataname=file_path_sans_ext(item)
	rawDir=file.path(tppsubDir,rawdataname)
	createFolders(rawDir,rawdataname)
	msConvert(rawDir,rawdataname)
	msGF(rawDir,rawdataname,enzyme,instrument,fragmentation,rawdataname,"diGLY")
	file.rename(from=file.path(tppmainDir,tppinput,paste(rawdataname,".raw",sep="")),to=file.path(tppmainDir,rawDir,"1_Spectra",paste(rawdataname,".raw",sep="")))
	file.rename(from=file.path(tppmainDir,tppinput,paste(rawdataname,".mzXML",sep="")),to=file.path(tppmainDir,rawDir,"1_Spectra",paste(rawdataname,".mzXML",sep="")))
	file.rename(from=file.path(tppmainDir,tppinput,paste(rawdataname,".mzid",sep="")),to=file.path(tppmainDir,rawDir,"2_MS-GF+",paste(rawdataname,".mzid",sep="")))
	file.rename(from=file.path(tppmainDir,tppinput,paste(rawdataname,".txt",sep="")),to=file.path(tppmainDir,rawDir,"2_MS-GF+",paste(rawdataname,".txt",sep="")))
	parseMzid(rawDir,rawdataname)
	idConvert(rawDir,rawdataname)
	file.copy(file.path(tppmainDir,rawDir,"1_Spectra",paste(rawdataname,".mzXML",sep="")),file.path(tppmainDir,rawDir,"3_TPP",paste(rawdataname,".mzXML",sep="")))
	pepProphet(rawDir,rawdataname)
	parsePepxml(rawDir,rawdataname)
	ptmProphet(rawDir,rawdataname,rawdataname)
}
for (item in files){
	rawdataname=file_path_sans_ext(item)
	rawDir=file.path(tppsubDir,rawdataname)
	pvalue(rawDir,rawdataname)
}
### /Running post-MSGF+ Workflow



### Running FilterFromAll
pepproph=read.PepProph(input="[xlsoutput]")
filt=filterAll(object=pepproph,name="files_merged_mods.txt",fasta="C:/MS-GF+/database/UniProt.human.20141017.RNFISnr.150contams.fasta")
### /Running FilterFromAll



### Compare Datasets
compare.datasets(uspcontrolsfilt,uspfilt)
### /Compare Datasets



### Compare Objects
compare.objects(sumofilt,uspcontrolsfilt,uspfilt)
### /Compare Objects

### Inset Function Manual
inset=function(protein,site){
	for(object in list(ubobj,ubfilt,sumoobj,nature2016all)){
		if(site %in% object@modsummary[[protein]]){
			print(paste(protein,":",site," in dataset = TRUE",sep=""))
			}
		else{
			print(paste(protein,":",site," in dataset = FALSE",sep=""))
			}
		}
	}
### /Inset Function Manual



### Inset Function apply
insetx=function(x){
	site=x$Site
	protein=x$Protein
	for(object in list(ubobj,ubfilt,sumoobj,nature2016all)){
		if(site %in% object@modsummary[[protein]]){
			print(paste(protein,":",site," in dataset = TRUE",sep=""))
			}
		else{
			print(paste(protein,":",site," in dataset = FALSE",sep=""))
			}
		}
	}
### /Inset Function apply



### Inset Function apply usage
apply(sumodown,1,insetx)
### /Inset Function apply usage


### CombineFilesInFolder

setwd(file.path(tppmainDir,"SENPWalpComb"))
files=list.files()
filesmerge=rbind()
for(item in files){
	print(item)
	filesdelim=read.delim(item,header=T)
	filesmerge=rbind(filesmerge,filesdelim)
}
write.table(filesmerge,"files_merged.xls",row.names=FALSE,sep="\t")
### /CombineFilesInFolder



###COMPARETABLES
##Import
#sumofilt
setwd(file.path(tppmainDir,"2017-04-17"))
sumopepproph=read.PepProph(input=file.path(tppmainDir,"2017-04-17","SUMO_merged12_filt0.9.txt"))
sumofilt=filterAll(object=sumopepproph,name=file.path(tppmainDir,"2017-04-17","SUMO_merged_mods_20170412.txt"),fasta="C:/MS-GF+/database/UniProt.human.20141017.RNFISnr.150contams.fasta",writetsv=FALSE)
#ubfilt
setwd(file.path(tppmainDir,"2017-04-17"))
ubpepproph=read.PepProph(input=file.path(tppmainDir,"2017-04-17","Ub12_merged.txt"))
ubfilt=filterAll(object=ubpepproph,name=file.path(tppmainDir,"2017-04-17","Ub_merged_mods_20170417.txt"),fasta="C:/MS-GF+/database/UniProt.human.20141017.RNFISnr.150contams.fasta",writetsv=FALSE)
#PSP
ubobj=read.psp(file.path("C:/Users/Admin/Desktop/Workflow/OtherDatabases/Ubiquitination_site_dataset"))
sumoobj=read.psp(file.path("C:/Users/Admin/Desktop/Workflow/OtherDatabases/Sumoylation_site_dataset"))
methylobj=read.psp(file.path("C:/Users/Admin/Desktop/Workflow/OtherDatabases/Methylation_site_dataset"))
acetylobj=read.psp(file.path("C:/Users/Admin/Desktop/Workflow/OtherDatabases/Acetylation_site_dataset"))
#Hendricks
natureobj=read.nature(file.path("C:/Users/Admin/Desktop/Workflow/OtherDatabases/nature.txt"))
nature2016=read.nature2016(file.path("C:/Users/Admin/Desktop/Workflow/OtherDatabases/nature2016.txt"))
#SUMO Uniprot
uniprotobj=read.uniprot(file.path("C:/Users/Admin/Desktop/Workflow/OtherDatabases/uniprot-annotation%253A%2528type%253Acrosslnk+SUMO%2529.txt"))
##Run
setwd(file.path(tppmainDir,"2017-04-17"))
compare.tables(sumofilt, natureobj, nature2016all, nature2016pub, nature2016re, nature2016score, ubobj, sumoobj, methylobj, acetylobj, uniprotobj,object2=ubfilt,compare=TRUE,agnostic=FALSE,name="20170417_SUMO_comparetables_re.tsv")
compare.tables(ubfilt, natureobj, nature2016all, nature2016pub, nature2016re, nature2016score, ubobj, sumoobj, methylobj, acetylobj, uniprotobj,compare=FALSE,agnostic=FALSE,name="20170419_Ub_comparetables.tsv")
###/COMPARETABLES