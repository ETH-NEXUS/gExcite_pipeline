import os.path
import sys
import inspect
import copy
import pandas as pd
from snakemake.utils import validate

samples = pd.read_table(config["inputOutput"]["sample_map"]).set_index("sample", drop=False)
validate(samples, "../schema/sample_map.schema.yaml")

# Retrieve hashed samples from the sample map of a given experiment
HashedSamples = samples["sample"].tolist()

# Retrieve filename listing tags for hashed samples.
def getTagFileHashedSamples(wildcards):
    sampleMap = dict(zip(samples['sample'], samples['HashingFile']))
    # Tests if file is existing
    if not os.path.isfile(sampleMap[wildcards.sample]):
        raise ValueError(
        "Sample '%s' does not contain a valid hashing status/file in the sample map!"
        % (wildcards.sample)
        )
    if wildcards.sample not in sampleMap.keys():
        raise ValueError(
            "Sample '%s' not found in the sample map!" % (wildcards.sample)
        )
    return sampleMap[wildcards.sample]

# Retrieve demultiplexed sample names.
def getDemultiplexedSamples():
    HashedSampleNames = samples["sample"].tolist()
    DemultiplexedSamples = []
    for sample in HashedSampleNames:
        tableEntry = samples.loc[samples["sample"]==sample]
        tagFile = tableEntry.loc[sample,"HashingFile"]
        subsampleTable = pd.read_table(tagFile,header=None,index_col=False, sep=",")
        subsamples = subsampleTable.iloc[:,2].tolist()
        DemultiplexedSamples.extend(subsamples)
    return DemultiplexedSamples
        

# Retrieve the number of target cells corresponding to a given sample set (both GEX and ADT)
def getTargetCells(wildcards):
    # Create dictionary with samples and nTargetCells
    sampleMap = dict(zip(samples['sample'], samples['nTargetCells']))
    if wildcards.sample not in sampleMap.keys():
        raise ValueError(
            "Sample '%s' not found in the sample map!" % (wildcards.sample)
        )
    return sampleMap[wildcards.sample]


# Preprend the retrieved target cell with the prefix required for citeseq count.
def getTargetCellsCiteseqCount(wildcards):
    value = getTargetCells(wildcards)
    out = []
    if value == ".":
        out.append("")
    else:
        out.append("-cells " + str(value))
    return out


# Preprend the retrieved target cell with the prefix required for cellranger count.
def getTargetCellsCellranger(wildcards):
    value = getTargetCells(wildcards)
    out = []
    if value == ".":
        out.append("")
    else:
        out.append("--force-cells " + str(value))
    return out


# Retrieve the ADT sequencing run name corresponding to a given sample set
def getSeqRunName(wildcards):
    # Create dictionary with samples and SeqRunName
    sampleMap = dict(zip(samples['sample'], samples['SeqRunName']))
    if wildcards.sample not in sampleMap.keys():
        raise ValueError(
            "Sample '%s' not found in the sample map!" % (wildcards.sample)
        )
    return sampleMap[wildcards.sample]


# Retrieve the ADT feature reference file corresponding to a given sample set
def getFeatRefFile(wildcards):
    # Create dictionary with samples and FeatureRefFile
    sampleMap = dict(zip(samples['sample'], samples['featureReferenceFile']))
    if wildcards.sample not in sampleMap.keys():
        raise ValueError(
            "Sample '%s' not found in the sample map!" % (wildcards.sample)
        )
    return sampleMap[wildcards.sample]

def list_fastqs(base,Rvalue):
    """Returns comma separated list of fastqs."""
    onlyfiles = [base + f for f in listdir(base) if isfile(join(base, f)) and Rvalue in f]
    return onlyfiles
