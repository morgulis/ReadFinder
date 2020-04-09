# ReadFinder
Efficient short nucleotide read filtering application.

## Prerequisites

ReadFinder needs Boost C++ libraries (version 1.62 or later) installed in order 
to compile successfully. In addition, in order to enable direct access to 
NCBI SRA archive, NCBI NGS and VDB libraries must be installed. 
The sections below describe how to install these prerequisites.

In the following we assume that all the commands are performed under a
directory pointed to by environment variable `ROOT`. So the first command
to do is to change into that directory.

```
> cd $ROOT
```

### Boost

If Boost development packages are not already installed on the system, or if
Boost version is too old, the following describes how to install it locally.

Assuming we are in `$ROOT` directory, create subdirectory for boost distribution.

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

If boost is installed in `$ROOT/boost/install`, then prepend the following command 
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
> NGS_ROOT=$ROOT/NGS/ngs VDB_ROOT=$ROOT/NGS/vdb cmake -DENABLE_NGS=Y -DCMAKE_BUILD_TYPE=Release ../ReadFinder/src
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

```--version [-v]```

Print application version information to the standard output and exit.
All other options are ignored.
  
```--help [-h]``` 

Print application usage help message to the standard output and exit.
All other options, except `--version [-v]` are ignored.
  
```--action [-A] <action_name>```

Default value: `search`.

Action to perform. Possible values are `mkdb` and `search`.
  
#### Common options

#### Options specific to `mkdb` action

#### Options specific to `search` action

## Examples
