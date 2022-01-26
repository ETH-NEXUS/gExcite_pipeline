import os.path
import sys
import inspect
import copy
import pandas as pd
from snakemake.utils import validate

fail_instantly = False

class Error(object):
    def __init__(self, key, name):
        self.__key = key
        self.__name = name

    def __add__(self, other):
        return self

    def __call__(self, wildcards=None):
        sys.exit(
            """
            ===============================================
            You have not specified '{}' for '{}'
            ===============================================
            """.format(
                self.__key, self.__name
            )
        )

    def __getitem__(self, value):
        return Error(key=self.__key, name=self.__name)


class Config(object):
    def __init__(self, kwargs, name="Config"):
        self.__name = name
        self.__members = {}
        for (key, value) in kwargs.items():
            if isinstance(value, dict):
                self.__members[key] = Config(kwargs=value, name=key)
            else:
                self.__members[key] = value

    def __getitem__(self, key):
        if key in self.__members:
            return self.__members[key]
        else:
            if fail_instantly:
                sys.exit(
                    """
                    ===============================================
                    You have not specified '{}' for '{}'
                    ===============================================
                    """.format(
                        key, self.__name
                    )
                )
            else:
                return Error(key=key, name=self.__name)

samples = pd.read_table(config["inputOutput"]["sample_map"]).set_index("sample", drop=False)
validate(samples, "../schema/sample_map.schema.yaml")

config = Config(config)

# Retrieve hashed samples from the sample map of a given experiment
HashedSamples = samples.loc[samples["HashingStatus"] != "."]
HashedSampleNames = HashedSamples["sample"].tolist()

# Retrieve filename listing tags for hashed samples.
def getTagFileHashedSamples(wildcards):
    HashedSamples = samples.loc[samples["HashingStatus"] != "."]
    sampleMap = dict(zip(HashedSamples['sample'], HashedSamples['HashingStatus']))
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
