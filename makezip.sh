#!/bin/bash

pwd=`pwd`
thisdir=`basename $pwd`

cd ..

find $thisdir -not -name "makezip.sh" -not -wholename "*/.git*"|zip ${thisdir}/alm${1}.zip -@

