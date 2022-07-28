#!/bin/sh
export CMSSW_TOOLS=$CMSSW_RELEASE_BASE/config/toolbox/$SCRAM_ARCH/tools/selected
echo $CMSSW_TOOLS
export TKONLINESW_BASE=$(cat $CMSSW_TOOLS/tkonlinesw.xml | grep name=\"TKONLINESW_BASE\" | grep default | sed 's/\" default=\"/=/g' | sed 's/<environment name=\"//g' | sed 's/\"\/>//g' | sed 's/TKONLINESW_BASE=//g' | tr -d " ")
echo $TKONLINESW_BASE
export XERCES_C_BASE=$(cat $CMSSW_TOOLS/xerces-c.xml | grep name=\"XERCES_C_BASE\" | grep default | sed 's/\" default=\"/=/g' | sed 's/<environment name=\"//g' | sed 's/\"\/>//g' | sed 's/XERCES_C_BASE=//g' | tr -d " ")
echo $XERCES_C_BASE
export ORACLE_BASE=$(cat $CMSSW_TOOLS/oracle.xml | grep name=\"ORACLE_BASE\" | grep default | sed 's/\" default=\"/=/g' | sed 's/<environment name=\"//g' | sed 's/\"\/>//g' | sed 's/ORACLE_BASE=//g' | tr -d " ")
echo $ORACLE_BASE
export PYTHON3_BASE=$(cat $CMSSW_TOOLS/python3.xml | grep name=\"PYTHON3_BASE\" | grep default | sed 's/\" default=\"/=/g' | sed 's/<environment name=\"//g' | sed 's/\"\/>//g' | sed 's/PYTHON3_BASE=//g' | tr -d " ")
echo $PYTHON3_BASE

swig -c++ -python -classic AccessDb.i

g++ -D DATABASE -D CMS_TK_64BITS -L $XERCES_C_BASE/lib -L $TKONLINESW_BASE/lib -I $PYTHON3_BASE/include/python3.8 -L $ORACLE_BASE/lib -I $ORACLE_BASE/include/ -I $TKONLINESW_BASE/include/ -I $XERCES_C_BASE/include/ -lxerces-c -locci -lclntsh  -lFed9UUtils -lICUtils -lFed9UDeviceFactory -llibDeviceDescriptions -std=c++0x -c AccessDb.cc AccessDb_wrap.cxx -fPIC

g++ -shared AccessDb.o AccessDb_wrap.o  -o _AccessDb.so -L $XERCES_C_BASE/lib -L $TKONLINESW_BASE/lib -L $ORACLE_BASE/lib -I $TKONLINESW_BASE/include/  -I $ORACLE_BASE/include/ -I $XERCES_C_BASE/include/ -lFed9UUtils -lICUtils -lFed9UDeviceFactory -lDeviceDescriptions -lxerces-c -locci -lclntsh
 
