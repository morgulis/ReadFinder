# ReadFinder
Efficient short nucleotide read filtering application.

## Prerequisites

ReadFinder needs Boost C++ libraries (version 1.62 or later) installed in order 
to compile successfully. In addition, in order to enable direct access to 
NCBI SRA archive, NCBI NGS and VDB libraries must be installed. 
The sections below describe how to install these prerequisites.

In the following we assume that all the commands are performed under a
directory pointed to by environment variable `ROOT`, so the first command
to do is to change into that directory.

```
> cd $ROOT
```

### Quick build

Download the build script from github.

```
> wget https://github.com/morgulis/ReadFinder/blob/master/src/scripts/build.sh?raw=true -O build.sh
```

Build ReadFinder (note: the order of options is important).

```
> bash ./build.sh [--no-boost] [--no-ngs] $PWD/build
```

`--no-boost` can be specified if there is a system wide installation of boost version at least 1.63.
`--no-ngs` can be specified to disable support for direct access to SRA.

The last argument must be an absolute path to the build directory.

If the build is successful, then `readfinder` executable will appear in
`build/readfinder/build/readfinder/` sub-directory of `$ROOT`.

### Boost

If Boost development packages are not already installed on the system, or if
Boost version is too old, the following describes how to install it locally.

Assuming we are in `$ROOT` directory, create subdirectory for Boost distribution.

```
> mkdir -p boost/download boost/src boost/install
```

Download the latest stable version (1.72.0 at the time of this writing)

```
> pushd boost/download
> wget https://dl.bintray.com/boostorg/release/1.72.0/source/boost_1_72_0.tar.bz2
> popd
```

Now `$ROOT/boost/download` contains file `boost_1_72_0.tar.bz2`.

Next, unpack boost into `$ROOT/boost/src` directory.

```
> pushd boost/src
> tar xvjf ../download/boost_1_72_0.tar.bz2
```

This puts boost source distribution under `$ROOT/boost/src/boost_1_72_0`.

Next build boost and install libraries, and headers under
`$ROOT/boost/install`.

```
> cd boost_1_72_0
> ./bootstrap.sh --prefix=$ROOT/boost/install
> ./b2 install
> popd
```

`$ROOT/boost/install/` now contains two subdirectories: `include/`, with
boost headers, and `lib/` with libraries.

### NGS

To install NGS software, first create directories for sources and binaries.

```
> cd $ROOT
> mkdir -p NGS/ngs-src NGS/ngs NGS/vdb-src NGS/vdb NGS/build
```

#### Install NGS software.

Configure and pull NGS SDK git repository into `$ROOT/NGS/ngs-src`.

```
> pushd NGS/ngs-src
> git init
> git remote add -f origin https://github.com/ncbi/ngs.git
> git config core.sparseCheckout true
> echo "ngs-sdk" >> .git/info/sparse-checkout
> git pull origin master
```

At this point `$ROOT/NGS/ngs-src` will contain subdirectory `ngs-sdk` with
NGS SDK sources.

Next step is to build NGS SDK and install it in `$ROOT/NGS/ngs`.

```
> cd ngs-sdk
> ./configure --prefix=$ROOT/NGS/ngs --build-prefix=$ROOT/NGS/build
> make
> make install
> popd
```
The library and headers are now installed under `$ROOT/NGS/ngs/`.

#### Install VDB software

Clone ncbi-vdb repository from github.

```
> pushd $ROOT/NGS/vdb-src
> git clone https://github.com/ncbi/ncbi-vdb.git
```
At this point `$ROOT/NGS/vdb-src` will contain subdirectory `ncbi-vdb` with
NGS SDK sources.

Next step is to build VDB and install it in `$ROOT/NGS/vdb`.

```
> cd ncbi-vdb
> ./configure --prefix=$ROOT/NGS/vdb --build-prefix=$ROOT/NGS/build --with-ngs-sdk-prefix=$ROOT/NGS/ngs
> make
> make install
> popd
```

The library and headers are now installed under `$ROOT/NGS/vdb/`.

## ReadFinder installation

Clone readfinder git repository.

```
> mkdir -p $ROOT/readfinder/build
> pushd $ROOT/readfinder
> git clone https://github.com/morgulis/ReadFinder
> cd build
```

ReadFinder sources are now installed under `$ROOT/readfinder/ReadFinder`.

If boost is installed in `$ROOT/boost/install`, then prepend the cmake commands below 
with `Boost_DIR=$ROOT/boost/install` and also execute:

```
export LD_LIBRARY_PATH=$ROOT/boost/install/lib:$LD_LIBRARY_PATH
```

To configure build without NGS:

```
> cmake -DCMAKE_BUILD_TYPE=Release ../ReadFinder/src
```

To configure build with NGS support:

```
> NGS_ROOT=$ROOT/NGS/ngs VDB_ROOT=$ROOT/NGS/vdb cmake -DENABLE_NGS=static -DCMAKE_BUILD_TYPE=Release ../ReadFinder/src
```

After build is configured, run `make`.

```
> make
> popd
```

ReadFinder executable should be in `$ROOT/readfinder/build/readfinder/` directory.

## Usage

Readfinder searches for exact matches of 21-mer words between reads and reference
sequences and selects reads that are covered by such matches within a narrow band
of diagonals up to a user specified threshold. The output is fasta formatted file
of selected reads.

### Readfinder invocation and command line options

```
readfinder [global options] [common options] [action specific options]
```

#### Global options

-----------------------------------
```--version [-v]```

Print application version information to the standard output and exit.
All other options are ignored.
  
-----------------------------------
```--help [-h]``` 

Print application usage help message to the standard output and exit.
All other options, except `--version [-v]` are ignored.
  
-----------------------------------  
```--action [-A] <action_name>```

Default value: `search`.

Action to perform. Possible values are `mkdb` and `search`.
`mkdb` action prepares a reference database. `search` action scans the
reads for matches to the reference database.
  
#### Common options

The following options work for both `mkdb` and `search` actions.

-----------------------------------
```--log-file <file_name>```

Redirect program log to the specified file. By default log goes to the
standard output.

-----------------------------------
```--trace-level <level>```

Default value: `warning`.

Sets the severity threshold for log messages. Possible values are:
`quiet`, `info`, `warning`, `error`.

-----------------------------------
```--threads [-t] <int>```

Default value: `0`

Max number of worker threads to use. The value of `0` means selects
max hardware treads supported by the machine.

-----------------------------------
```--quiet```

Disable progress reporting.

#### Options specific to `mkdb` action

-----------------------------------
```--input [-i] <file_name>```

Name of the Fasta file containing reference sequences. If the name
ends in `.gz` the file is assummed compressed by gzip.

-----------------------------------
```--output [-o] <name>```

Base name for database files. ReadFinder database consists of several
files. This parameter is the common prefix for the database file names.

-----------------------------------
```--mkidx```

Add word index to the database. By default the index is created on the
fly during the search. However this can reduce performance for large
references. In this case this flags allows to pre-compute the index
and store it along with other database files. During the search, the
index is automatically loaded from disk, if present.

-----------------------------------
```--mkws```

Add word bitset to the database. This option must be specified if
`--prescreen` option is used for search. This allows for large
speed and memory optimization, when most of the reads are not expected
to match the reference.

#### Options specific to `search` action

-----------------------------------
```--db [-d] <name>```

Base name for the database (previousely created by `mkdb` action).

-----------------------------------
```--input [-i] <name>```

Input specification. `name` is either a file name or SRA accession,
depending on the value of `--in-fmt` option. If SRA accession is
specified, then paired run will induce paired search; non-paired run
will induce non-paired search. If a file is specified, then the 
search is non-paired. Use `--first-mate` and `--second-mate`
options to specify input files for paired search.

-----------------------------------
```--first-mate [-1] <file_name>```

Name of the file containing first mates of the input reads.

-----------------------------------
```--second-mate [-2] <file_name>```

Name of the file containing second mates of the input reads.

-----------------------------------
```--in-fmt [-F] <format_name>```

Default value: `fasta` or `zfasta`, if input file name ends in `.gz`

The following values are supported:
  - `fasta` --- uncompressed fasta file(s)
  - `zfasta` --- gzip compressed fasta file(s)
  - `fastq` --- uncompressed fastq file(s)
  - `zfastq` --- gzip compressed fastq file(s)
  - `sra` --- SRA accession.

In the case of paired input, mixed formats are not supported.

-----------------------------------
```--memory [-M] <positive int>```

Memory limit in megabytes. This is hint. The actual memory usage can
be slightly higher. Note that there must be at least enough space for
the in-memory part of the database, otherwise the program will exit
with error. By default memory limit is set to 90% of the available
machine memmory.

-----------------------------------
```--keep-db-loaded```

If this flag is specified, the whole database index is kept pre-loaded
in memory. This can help performance in setups with slow disk or
network access.

-----------------------------------
```--batch-size [-b] <positive_int>```

Default value: `5000000`

Process input in batches of this size.

-----------------------------------
```--first-read <non-negative_int>```

Default value: `0`

Index of the first read to process (0-based).

-----------------------------------
```--last-read <non-negative_int>```

Default value: number of reads - 1

Index of the last read to process (0-based).

-----------------------------------
```--output [-o] <file_name>```

Default: standard output

Name of the output Fasta file containing matching reads.

-----------------------------------
```--continuous-bases [-C] <non-negative_int>```

Default: `21`

Minimum length of coninuous exact match needed to select the
sequence as matching. A sequence is selected if it satisfies
`--continuous-bases`, `--covered-bases` constraint,
and `--coverage-ratio` constraint.

-----------------------------------
```--covered-bases [-B] <non-negative_int>```

Default: `0`

Minimum number of bases covered by matches within a diagonal
band (see `--max-diag-delta` option) needed to select the
sequence as matching. A sequence is selected if it satisfies
`--continuous-bases`, `--covered-bases` constraint, and
`--coverage-ratio` constraint.

-----------------------------------
```--coverage-ratio [-c] <float>```

Default value: `0.5`

The value must be between `0.0` and `1.0`. The minimum
proportion of bases covered by matches within a diagonal
band (see `--max-diag-delta` option) needed to select the
sequence as matching. A sequence is selected if it satisfies
`--continuous-bases`, `--covered-bases` constraint, and
`--coverage-ratio` constraint.

-----------------------------------
```--max-diag-delta [-D] <non-negative_int>```

Default value: `4`

Matches used for `--covered-bases` and `--coverage-ratio` tests
must on diagonals that are at most this far from each other.

-----------------------------------
```--paired```

This flag should be used to interpret a single fasta or fastq
file as containing paired reads. The input file should contain
even number of sequences. Consecutive pairs of sequences constitute
a paired read. The ids of the read mates must be the same or
differ only in the suffix: ".1" for the first mate, ".2" for the
second mate. This flag has no effect for sra input type or when
two input files are used.

-----------------------------------
```--per-mate-marks [-m]```

By default, if a mate in a paired read is matched, both mate
sequences appear in the output. If this flag is specified,
only the mate that matches will appear in the output.

-----------------------------------
```--prescreen```

If this flag is specified, a special check is performed to
quickly filter out reads that can't satisfy matching criteria.
This optimization should only be used when it is known that
very small ratio of reads will match. The database must be
created with `--mkws` option in order to use this feature.


## Examples

Assume reference sequence is in file `ref.fa`. To create the 
database out of it:

```readfinder -A mkdb -i ref.fa -o ref```

This will generate the following files:
- `ref.dat`
- `ref.ids`
- `ref.off`

To create additional reference index (especially useful for large references):

```readfinder -A mkdb -i ref.fa -o ref --mkidx```

This will generate the following additional files:
- `ref.idm`
- `ref.idx`

To create word bitset needed for read pre-screening:

```readfinder -A mkdb -i ref.fa -o ref --mkidx --mkws```

This adds another file `ref.ws`.

To match SRA accession SRR6399765 against `ref` database in batches of 100M reads,
pre-screening the reads, matching at least 40 bp and 60% of the mate sequence,
selecting only matching mates, and sending results to SRR6399765.out.fa:

```readfinder -A search -d ref -F sra -i SRR6399765 -m --prescreen -b 100000000 -B 40 -c 0.6 -o SRR6399765.out.fa```

To do the same when reads are in compressed fastq files SRR6399765_1.fq.gz and SRR6399765_2.fq.gz:

```readfinder -A search -d ref -F fastq -1 SRR6399765_1.fq.gz -2 SRR6399765_2.fq.gz -m --prescreen -b 100000000 -B 40 -c 0.6 -o SRR6399765.out.fa```


