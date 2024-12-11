#!/bin/bash

# This script will attempt to build ReadFinder software and its prerequisites.

USAGE="
USAGE build.sh [--no-boost] [--no-ngs] install-root-dir

    --no-boost          -- BOOST libraries will not be built,
                           in this case BOOST headers and libraries
                           should be installed in standard system
                           locations

    --no-ngs            -- NGS libraries will not be built,
                           in this case NCBI SRA access will be disabled

    install-root-dir    -- all build/installation work will be done under
                           this directory; the directory will be created;
                           it's an error if the directory already exists

    NOTES:

        1. install-root-dir must be an absolute path;
        2. to build BOOST wget must be installed and accessible via PATH;
        3. to build NGS and ReadFinder GNU make and git must be installed
            and accessible via PATH.
        4. to build ReadFinder cmake must be installed and accecssible via PATH
"

################################################################################
#
# initial setup
#

BOOST_DIST="https://archives.boost.io/release/1.85.0/source/boost_1_85_0.tar.bz2"

if [[ $# < 1 ]] ; then echo "$USAGE" ; exit 1 ; fi

boost=1
ngs=1

if [[ "$1" == "--no-boost" ]] ; then
    boost=0 ;
    shift ;
    if [[ $# < 1 ]] ; then echo "$USAGE" ; exit 1 ; fi
fi

if [[ "$1" == "--no-ngs" ]] ; then
    ngs=0 ;
    shift ;
    if [[ $# < 1 ]] ; then echo "$USAGE" ; exit 1 ; fi
fi

echo "build boost: $boost"
echo "build ngs: $ngs"

if [[ -e $1 ]] ; then
    echo "ERROR: $1 exists" ;
    echo "$USAGE" ;
    exit 1 ;
fi

root=$1

case $root in
    /*) echo "building in $root" ;;
    *)  echo "ERROR: last argument must be an absolute path" ; echo "$USAGE" ; exit 1 ;;
esac

mkdir $root || { echo "can't create $root" ; exit 1 ; }
pushd $root > /dev/null

################################################################################
#
# build boost, if requested
#
if [[ $boost == 1 ]] ; then
    echo "building boost" ;
    BOOSTDIR="boost" ;
    mkdir -p "$BOOSTDIR/download" "$BOOSTDIR/src" "$BOOSTDIR/install" || \
        { echo "can't create $BOOSTDIR or its subdirectory" ; exit 1 ; } ;

    # download boost
    #
    pushd $BOOSTDIR/download > /dev/null ;
    echo $BOOST_DIST ;
    wget $BOOST_DIST || { echo "can't download boost libraries" ; exit 1 ; } ;
    echo "downloaded $BOOST_DIST" ;
    popd > /dev/null ;

    # unpack boost
    #
    pushd $BOOSTDIR/src > /dev/null ;
    tar xvjf ../download/boost_*.tar.bz2 || \
        { echo "unpack of boost failed" ; exit 1 ; } ;
    echo "boost unpacked" ;

    # build boost
    #
    pushd boost_* > /dev/null ;
    ./bootstrap.sh --with-libraries=program_options,iostreams \
                   --prefix=$root/$BOOSTDIR/install || \
        { echo "error booststraping boost" ; exit 1 ; } ;
    ./b2 install ;
    echo "built $BOOST_DIST" ;
    popd > /dev/null ;
    popd > /dev/null ;
fi

################################################################################
#
# Prepare to build NGS and VDB
#
if [[ $ngs == 1 ]] ; then
    NGSDIR="NGS" ;
    mkdir -p $NGSDIR/ngs-src $NGSDIR/ngs \
             $NGSDIR/vdb-src $NGSDIR/vdb \
             $NGSDIR/sratk-src $NGSDIR/build || \
        { echo "can't create $NGSDIR or its subdirectory" ; exit 1; } ;

################################################################################
#
# Build NGS libraries
#
    pushd $NGSDIR/ngs-src > /dev/null ;
    git init || { echo "git init failed" ; exit 1 ; } ;
    git remote add -f origin https://github.com/ncbi/ngs.git || \
        { echo "failed to add remote https://github.com/ncbi/ngs.git" ; exit 1 ; } ;
    git config core.sparseCheckout true || \
        { echo "git config failed" ; exit 1 ; } ;
    echo "ngs-sdk" >> .git/info/sparse-checkout ;
    git pull origin master || { echo "failed to pull from remote" ; exit 1 ; } ;
    cd ngs-sdk ;
    ./configure --prefix=$root/$NGSDIR/ngs --build-prefix=$root/$NGSDIR/build || \
        { echo "NGS configure failed" ; exit 1 ; } ;
    make || { echo "NGS make failed" ; exit 1 ; } ;
    make install || { echo "NGS install failed" ; exit 1 ; } ;
    popd > /dev/null ;

################################################################################
#
# Build VDB libraries
#
    pushd $NGSDIR/vdb-src > /dev/null ;
    git clone https://github.com/ncbi/ncbi-vdb.git || \
        { echo "failed to clone ncbi-vdb" ; exit 1 ; } ;
    cd ncbi-vdb ;
    ./configure --prefix=$root/$NGSDIR/vdb --build-prefix=$root/$NGSDIR/build || \
        { echo "ncbi-vdb configure failed" ; exit 1 ; } ;
    make || { echo "ncbi-vdb make failed" ; exit 1 ; } ;
    make install || { echo "ncbi-vdb install failed" ; exit 1 ; } ;
    popd > /dev/null ;

################################################################################
#
# Build SRA tools libraries
#
    pushd $NGSDIR/sratk-src > /dev/null ;
    git clone https://github.com/ncbi/sra-tools.git || \
        { echo "failed to clone sra-tools" ; exit 1 ; } ;
    cd sra-tools ;
    ./configure --prefix=$root/$NGSDIR/vdb \
                --with-ncbi-vdb-prefix=$root/$NGSDIR/vdb \
                --build-prefix=$root/$NGSDIR/build || \
        { echo "sra-tools configure failed" ; exit 1 ; } ;
    make || { echo "sra-tools make failed" ; exit 1 ; } ;
    make install || { echo "sra-tools install failed" ; exit 1 ; } ;
    popd > /dev/null ;

fi

################################################################################
#
# Build ReadFinder
#
mkdir -p $root/readfinder/build || \
    { echo "can't create $root/readfinder/build" ; exit 1 ; }
pushd $root/readfinder > /dev/null
git clone https://github.com/morgulis/ReadFinder || \
    { echo "can't clone from https://github.com/morgulis/ReadFinder" ; exit 1 ; }
cd build
CMAKE_OPTS=""
if [[ $boost == 1 ]] ; then
    export Boost_DIR=$root/$BOOSTDIR/install ;
    export LD_LIBRARY_PATH=$Boost_DIR/lib:$LD_LIBRARY_PATH ;
fi
if [[ $ngs == 1 ]]; then
    export NGS_ROOT=$root/$NGSDIR/ngs ;
    export VDB_ROOT=$root/$NGSDIR/vdb ;
    CMAKE_OPTS=" -DENABLE_NGS=static " ;
fi
cmake $CMAKE_OPTS -DCMAKE_BUILD_TYPE=Release ../ReadFinder/src || \
    { echo "cmake failed for readfinder" ; exit 1; }
make || { echo "readfinder make failed" ; exit 1 ; }
popd > /dev/null

popd > /dev/null

