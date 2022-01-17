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

NonHashedSamples = samples.loc[samples["HashingStatus"] == "."]
NonHashedSampleNames = NonHashedSamples["sample"].tolist()



# Require all final files, from both hashed and non-hashed samples (if any)
def getInputFiles(wildcards):
    allFiles = []

    # Cellranger output for hashed samples (in preprocessing folder)
    cr_hashed_gex = expand(
        "results/pooled_samples/cellranger_gex/{sample}.features.tsv", sample=HashedSampleNames
    )
    cr_hashed_adt = expand(
        "results/pooled_samples/cellranger_adt/{sample}.features.tsv", sample=HashedSampleNames
    )

#    # Cellranger output for non-hashed samples (also in preprocessing folder)
#    cr_nonhashed_gex = expand(
#        "results/cellranger_gex/{sample}.features.tsv", sample=getNonHashedSampleNames()
#    )
#    cr_nonhashed_adt = expand(
#        "results/cellranger_adt/{sample}.features.tsv", sample=getNonHashedSampleNames()
#    )

    # Hashing framework: symlinks to cellranger output (in pooled sample folders)
    root_hashed_gex = expand(
        "results/pooled_samples/cellranger_gex/{sample}.features.tsv",
        sample=HashedSampleNames,
    )
    root_hashed_adt = expand(
        "results/pooled_samples/cellranger_adt/{sample}.features.tsv",
        sample=HashedSampleNames,
    )
    # Citeseq output ( in pooled analysis folder)
    citeseq = expand(
        "results/pooled_samples/citeseq_count/{sample}.run_report.yaml",
        sample=HashedSampleNames,
    )
    # Hashing output ( in pooled analysis folder)
    hashing = expand(
        "results/pooled_samples/hashing_analysis/{sample}.complete_hashing.txt",
        sample=HashedSampleNames,
    )

    # Add all strings to list of required files. Empty will be skipped
    # Cellranger results
    for f in cr_hashed_gex:
        allFiles.append(f)
    for f in cr_hashed_adt:
        allFiles.append(f)
#    for f in cr_nonhashed_gex:
#        allFiles.append(f)
#    for f in cr_nonhashed_adt:
#        allFiles.append(f)
    # Symlinks to analysis folder
    for f in root_hashed_gex:
        allFiles.append(f)
    for f in root_hashed_adt:
        allFiles.append(f)
    # Citeseq
    for f in citeseq:
        allFiles.append(f)
    # Hashing
    for f in hashing:
        allFiles.append(f)
    return allFiles


# Retrieve filename listing tags for non-hashed samples.
def getTagFileHashedSamples(wildcards):
    if not "SAMPLEMAPPING" in globals():
        return ["NOMAPPINGFILE"]
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ["NOMAPPINGFILE"]
    sampleMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            if line.startswith("sample"):
                continue
            if line.strip() != "":
                lineSplit = line.strip().split()
                sample = lineSplit[0].strip()
                status = lineSplit[1].strip()
                if status not in ["NH", "."]:
                    if not os.path.isfile(status):
                        raise ValueError(
                            "Sample '%s' does not contain a valid hashing status/file in the sample map!"
                            % (sample)
                        )
                if sample in sampleMap.keys():
                    raise ValueError(
                        "Sample '%s' is not unique in the sample map!"
                        % (wildcards.sample)
                    )
                sampleMap[sample] = status
    if wildcards.sample not in sampleMap.keys():
        raise ValueError(
            "Sample '%s' not found in the sample map!" % (wildcards.sample)
        )
    return sampleMap[wildcards.sample]



# Retrieve the number of target cells corresponding to a given sample set (both GEX and ADT)
def getTargetCells(wildcards):
    if not "SAMPLEMAPPING" in globals():
        return ["NOMAPPINGFILE"]
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ["NOMAPPINGFILE"]
    sampleMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            if line.startswith("sample"):
                continue
            if line.strip() != "":
                lineSplit = line.strip().split()
                sample = lineSplit[0].strip()
                nCells = lineSplit[3].strip()
                if sample in sampleMap.keys():
                    raise ValueError(
                        "Sample '%s' is not unique in the sample map!"
                        % (wildcards.sample)
                    )
                sampleMap[sample] = "."
                if nCells != ".":
                    sampleMap[sample] = int(nCells)

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
    if not "SAMPLEMAPPING" in globals():
        return ["NOMAPPINGFILE"]
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ["NOMAPPINGFILE"]
    sampleMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            if line.startswith("sample"):
                continue
            if line.strip() != "":
                lineSplit = line.strip().split()
                sample = lineSplit[0].strip()
                seqRun = lineSplit[2].strip()
                if sample in sampleMap.keys():
                    raise ValueError(
                        "Sample '%s' is not unique in the sample map!"
                        % (wildcards.sample)
                    )
                if seqRun == ".":
                    raise ValueError(
                        "Sample '%s' does not contain a sequencing run name in the sample map!"
                        % (wildcards.sample)
                    )
                sampleMap[sample] = seqRun

    if wildcards.sample not in sampleMap.keys():
        raise ValueError(
            "Sample '%s' not found in the sample map!" % (wildcards.sample)
        )
    return sampleMap[wildcards.sample]


# Retrieve the ADT feature reference file corresponding to a given sample set
def getFeatRefFile(wildcards):
    if not "SAMPLEMAPPING" in globals():
        return ["NOMAPPINGFILE"]
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ["NOMAPPINGFILE"]
    sampleMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            if line.startswith("sample"):
                continue
            if line.strip() != "":
                lineSplit = line.strip().split()
                sample = lineSplit[0].strip()
                featFile = lineSplit[4].strip()
                if sample in sampleMap.keys():
                    raise ValueError(
                        "Sample '%s' is not unique in the sample map!"
                        % (wildcards.sample)
                    )
                if featFile == ".":
                    raise ValueError(
                        "Sample '%s' does not contain a feature reference file in the sample map!"
                        % (wildcards.sample)
                    )
                sampleMap[sample] = featFile

    if wildcards.sample not in sampleMap.keys():
        raise ValueError(
            "Sample '%s' not found in the sample map!" % (wildcards.sample)
        )
    return sampleMap[wildcards.sample]

def list_fastqs(base,Rvalue):
    """Returns comma separated list of fastqs."""
    onlyfiles = [base + f for f in listdir(base) if isfile(join(base, f)) and Rvalue in f]
    return onlyfiles
