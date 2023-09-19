import os.path
import sys
import inspect
import copy
import pandas as pd
from snakemake.utils import validate

samples = pd.read_table(config["inputOutput"]["sample_map"]).set_index(
    "sample", drop=False
)
validate(samples, "../schema/sample_map.schema.yaml")

# Retrieve hashed samples from the sample map of a given experiment
HashedSamples = samples["sample"].tolist()


# Retrieve filename listing tags for hashed samples.
def getTagFileHashedSamples(wildcards):
    sampleMap = dict(zip(samples["sample"], samples["HashingFile"]))
    # Tests if file is existing
    if not os.path.isfile(sampleMap[wildcards.sample_set]):
        raise ValueError(
            "Sample '%s' does not contain a valid hashing status/file in the sample map!"
            % (wildcards.sample_set)
        )
    if wildcards.sample_set not in sampleMap.keys():
        raise ValueError(
            "Sample '%s' not found in the sample map!" % (wildcards.sample_set)
        )
    return sampleMap[wildcards.sample_set]


# Retrieve demultiplexed sample names.
def getDemultiplexedSamples(sample):
    HashedSampleNames = samples["sample"].tolist()
    DemultiplexedSamples = []
    for HashedSample in HashedSampleNames:
        if sample == HashedSample:
            tableEntry = samples.loc[samples["sample"] == sample]
            tagFile = tableEntry.loc[sample, "HashingFile"]
            subsampleTable = pd.read_table(
                tagFile, header=None, index_col=False, sep=","
            )
            DemultiplexedSamples = subsampleTable.iloc[:, 2].tolist()
    return DemultiplexedSamples


# Retrieve demultiplexed sample names.
def getCompleteSampleNames():
    HashedSampleNames = samples["sample"].tolist()
    DemultiplexedSamples = []
    for sample in HashedSampleNames:
        tableEntry = samples.loc[samples["sample"] == sample]
        tagFile = tableEntry.loc[sample, "HashingFile"]
        subsampleTable = pd.read_table(tagFile, header=None, index_col=False, sep=",")
        subsamples = subsampleTable.iloc[:, 2].tolist()
        for subsample in subsamples:
            DemultiplexedSamples.append(sample + "." + subsample)
    return DemultiplexedSamples


# Retrieve sample names for a sample without multiplexing
def getSimpleSampleNames():
    SimpleSampleNames = samples["sample"].tolist()
    return SimpleSampleNames


# Retrieve the number of target cells corresponding to a given sample set (both GEX and ADT)
def getTargetCells(wildcards):
    # Create dictionary with samples and nTargetCells
    sampleMap = dict(zip(samples["sample"], samples["nTargetCells"]))
    if wildcards.sample_set not in sampleMap.keys():
        raise ValueError(
            "Sample '%s' not found in the sample map!" % (wildcards.sample_set)
        )
    return sampleMap[wildcards.sample_set]


# Retrieve the number of target cells corresponding to a given (not multiplexed) sample (both GEX and ADT)
def getTargetCells_simple(wildcards):
    # Create dictionary with samples and nTargetCells
    sampleMap = dict(zip(samples["sample"], samples["nTargetCells"]))
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


# Preprend the retrieved target cell with the prefix required for cellranger count.
def getTargetCellsCellranger_simple(wildcards):
    value = getTargetCells_simple(wildcards)
    out = []
    if value == ".":
        out.append("")
    else:
        out.append("--force-cells " + str(value))
    return out


# Retrieve the ADT sequencing run name corresponding to a given sample set
def getSeqRunName(wildcards):
    # Create dictionary with samples and SeqRunName
    sampleMap = dict(zip(samples["sample"], samples["SeqRunName"]))
    if wildcards.sample_set not in sampleMap.keys():
        raise ValueError(
            "Sample '%s' not found in the sample map!" % (wildcards.sample_set)
        )
    return sampleMap[wildcards.sample_set]


# Retrieve the ADT feature reference file corresponding to a given sample set
def getFeatRefFile(wildcards):
    # Create dictionary with samples and FeatureRefFile
    sampleMap = dict(
        zip(samples["sample"], WORKDIR + "/" + samples["featureReferenceFile"])
    )
    if wildcards.sample_set not in sampleMap.keys():
        raise ValueError(
            "Sample '%s' not found in the sample map!" % (wildcards.sample_set)
        )
    return sampleMap[wildcards.sample_set]


# Retrieve the ADT feature reference file corresponding to a given simple (unhashed) sample
def getFeatRefFileSimple(wildcards):
    # Create dictionary with samples and FeatureRefFile
    sampleMap = dict(
        zip(samples["sample"], WORKDIR + "/" + samples["featureReferenceFile"])
    )
    if wildcards.sample not in sampleMap.keys():
        raise ValueError(
            "Sample '%s' not found in the sample map!" % (wildcards.sample)
        )
    return sampleMap[wildcards.sample]


def list_fastqs(base, Rvalue):
    """Returns comma separated list of fastqs."""
    onlyfiles = [
        base + f for f in listdir(base) if isfile(join(base, f)) and Rvalue in f
    ]
    return onlyfiles
