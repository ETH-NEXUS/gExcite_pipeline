import os.path
import sys
import inspect
import copy

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


config = Config(config)


def getSampleNames():
    output = (
        []
    )  # [samplename.replace(FASTQDIR,'').replace('/','')for samplename in glob.glob(FASTQDIR + '*/')]
    if output == []:
        if not "SAMPLEMAPPING" in globals():
            return ["NOMAPPINGFILE"]
        try:
            open(SAMPLEMAPPING, "r")
        except IOError:
            return ["NOMAPPINGFILE"]
        sampleMap = dict()
        with open(SAMPLEMAPPING, "r") as f:
            for line in f:
                if line.strip() != "":
                    lineSplit = line.strip().split()
                    sample = lineSplit[1]
                    if not (sample in output):
                        output.append(sample)
    return output


# Retrieve hashed samples from the sample map of a given experiment
def getHashedSampleNames():
    output = (
        []
    )  # [samplename.replace(FASTQDIR,'').replace('/','')for samplename in glob.glob(FASTQDIR + '*/')]
    if output == []:
        if not "SAMPLEMAPPING" in globals():
            return ["NOMAPPINGFILE"]
        try:
            open(SAMPLEMAPPING, "r")
        except IOError:
            return ["NOMAPPINGFILE"]
        sampleMap = dict()
        with open(SAMPLEMAPPING, "r") as f:
            for line in f:
                if line.strip() != "":
                    lineSplit = line.strip().split()
                    sample = lineSplit[1].strip()
                    status = lineSplit[2].strip()
                    if status not in ["NH", "."]:
                        if not os.path.isfile(status):
                            raise ValueError(
                                "Sample '%s' does not contain a valid hashing status/file in the sample map!"
                                % (sample)
                            )
                    # Only include sample in output string if it is hashed ie. "H"
                    if os.path.isfile(status):
                        if not (sample in output):
                            output.append(sample)
    return output


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
            if line.strip() != "":
                lineSplit = line.strip().split()
                sample = lineSplit[1].strip()
                status = lineSplit[2].strip()
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


# Retrieve non-hashed samples from the sample map of a given experiment
def getNonHashedSampleNames():
    output = (
        []
    )  # [samplename.replace(FASTQDIR,'').replace('/','')for samplename in glob.glob(FASTQDIR + '*/')]
    if output == []:
        if not "SAMPLEMAPPING" in globals():
            return ["NOMAPPINGFILE"]
        try:
            open(SAMPLEMAPPING, "r")
        except IOError:
            return ["NOMAPPINGFILE"]
        sampleMap = dict()
        with open(SAMPLEMAPPING, "r") as f:
            for line in f:
                if line.strip() != "":
                    lineSplit = line.strip().split()
                    sample = lineSplit[1].strip()
                    status = lineSplit[2].strip()
                    if status not in ["NH", "H", "."]:
                        if not os.path.isfile(status):
                            raise ValueError(
                                "Sample '%s' does not contain a valid hashing status in the sample map!"
                                % (sample)
                            )
                    # Only include sample in output string if it is non-hashed ie. "NH"
                    if status in ["NH", "."]:
                        if not (sample in output):
                            output.append(sample)
    return output


def getExperimentNames():
    output = (
        []
    )  # [samplename.replace(FASTQDIR,'').replace('/','')for samplename in glob.glob(FASTQDIR + '*/')]
    if output == []:
        if not "SAMPLEMAPPING" in globals():
            return ["NOMAPPINGFILE"]
        try:
            open(SAMPLEMAPPING, "r")
        except IOError:
            return ["NOMAPPINGFILE"]
        sampleMap = dict()
        with open(SAMPLEMAPPING, "r") as f:
            for line in f:
                if line.strip() != "":
                    lineSplit = line.strip().split()
                    exp = lineSplit[0]
                    if not (exp in output):
                        output.append(exp)
    return output


def getSampleNamesFromExperimentNames(wildcards):
    if not "SAMPLEMAPPING" in globals():
        return ["NOMAPPINGFILE"]
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ["NOMAPPINGFILE"]
    expMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            if line.strip() != "":
                lineSplit = line.strip().split()
                exp = lineSplit[0]
                sample = lineSplit[1]
                sampleType = lineSplit[2]
                tpoint = lineSplit[3]
                if exp not in expMap.keys():
                    expMap[exp] = []
                expMap[exp].append(sample)
    return expMap[wildcards.experiment]


def checkFilesAgainstSampleNames(files, sampleNames):
    finalFiles = []
    for f in files:
        for name in sampleNames:
            if name + "/" == f[0 : len(name + "/")]:
                finalFiles.append(f)

    return finalFiles


def getSingleFastqFiles(SAMPLENAMES):
    files = [
        file.replace(FASTQDIR, "").replace(".fastq.gz", "")
        for file in glob.glob(FASTQDIR + "*/SINGLEEND/*.fastq.gz")
    ]
    if files == []:
        files = [
            file.replace(FASTQDIR, "").replace(".fastq", "")
            for file in glob.glob(FASTQDIR + "*/SINGLEEND/*.fastq")
        ]

    return checkFilesAgainstSampleNames(files, SAMPLENAMES)


# return [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq.gz')]
# return [file.replace(FASTQDIR, '').replace('.fastq','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq')]


def getPairedFastqFiles(SAMPLENAMES):
    files = [
        file.replace(FASTQDIR, "").replace(".fastq.gz", "")
        for file in glob.glob(FASTQDIR + "*/PAIREDEND/*R[12].fastq.gz")
    ]
    if files == []:
        files = [
            file.replace(FASTQDIR, "").replace(".fastq", "")
            for file in glob.glob(FASTQDIR + "*/PAIREDEND/*R[12].fastq")
        ]

    return checkFilesAgainstSampleNames(files, SAMPLENAMES)


# return [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*R[12].fastq.gz')]
# return [file.replace(FASTQDIR, '').replace('.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*R[12].fastq')]


def getPairedFastqFilesWithoutR(SAMPLENAMES):
    files = [
        file.replace(FASTQDIR, "").replace("_R1.fastq.gz", "")
        for file in glob.glob(FASTQDIR + "*/PAIREDEND/*_R1.fastq.gz")
    ]
    if files == []:
        files = [
            file.replace(FASTQDIR, "").replace("_R1.fastq", "")
            for file in glob.glob(FASTQDIR + "*/PAIREDEND/*_R1.fastq")
        ]

    return checkFilesAgainstSampleNames(files, SAMPLENAMES)


# return [file.replace(FASTQDIR, '').replace('_R1.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq.gz')]
# return [file.replace(FASTQDIR, '').replace('_R1.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq')]

## Deprecated
# def getNormalTumorFiles():
#     if not 'SAMPLEMAPPING' in globals():
#         return ['NOMAPPINGFILE']
#     try:
#         open(SAMPLEMAPPING, "r")
#     except IOError:
#         return ['NOMAPPINGFILE']
#     output = []
#     sampleMap = dict()
#     with open(SAMPLEMAPPING, "r") as f:
#         for line in f:
#             if line.strip() != "":
#                 lineSplit = line.strip().split()
#                 exp = lineSplit[0]
#                 sample = lineSplit[1]
#                 sampleType = lineSplit[2]
#                 tpoint = lineSplit[3]
#                 if exp not in sampleMap.keys():
#                     sampleMap[exp] = dict()
#                 if tpoint not in sampleMap[exp].keys():
#                     sampleMap[exp][tpoint] = dict()
#                 if sampleType not in sampleMap[exp][tpoint].keys():
#                     sampleMap[exp][tpoint][sampleType] = []
#                 sampleMap[exp][tpoint][sampleType].append(sample)
#     for expKey, expValue in sampleMap.items():
#         for tpointKey, tpointValue in expValue.items():
#             if 'T' in tpointValue and 'N' in tpointValue:
#                 for sampleTumor in tpointValue['T']:
#                     for sampleNormal in tpointValue['N']:
#                         output.append(sampleTumor + '_vs_' + sampleNormal)
#     return output

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
            if line.strip() != "":
                lineSplit = line.strip().split()
                sample = lineSplit[1].strip()
                nCells = lineSplit[4].strip()
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
            if line.strip() != "":
                lineSplit = line.strip().split()
                sample = lineSplit[1].strip()
                seqRun = lineSplit[3].strip()
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
            if line.strip() != "":
                lineSplit = line.strip().split()
                sample = lineSplit[1].strip()
                featFile = lineSplit[5].strip()
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
