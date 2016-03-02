#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-MacOSX
CND_DLIB_EXT=dylib
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/BedFactory.o \
	${OBJECTDIR}/src/FastaFactory.o \
	${OBJECTDIR}/src/KmersFactory.o \
	${OBJECTDIR}/src/SNPFactory.o \
	${OBJECTDIR}/src/SVMPredict.o \
	${OBJECTDIR}/src/berror.o \
	${OBJECTDIR}/src/bmemory.o \
	${OBJECTDIR}/src/bstring.o \
	${OBJECTDIR}/src/chebyshev.o \
	${OBJECTDIR}/src/features.o \
	${OBJECTDIR}/src/gamma.o \
	${OBJECTDIR}/src/lgamma.o \
	${OBJECTDIR}/src/lgammacor.o \
	${OBJECTDIR}/src/phyper.o \
	${OBJECTDIR}/src/stirlerr.o \
	${OBJECTDIR}/src/svm.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/bedFactoryTest \
	${TESTDIR}/TestFiles/fastaFactoryTest \
	${TESTDIR}/TestFiles/kmerFactoryTest \
	${TESTDIR}/TestFiles/peakTest \
	${TESTDIR}/TestFiles/phyperTest

# Test Object Files
TESTOBJECTFILES= \
	${TESTDIR}/tests/BedFactoryTest.o \
	${TESTDIR}/tests/FastaFactoryTest.o \
	${TESTDIR}/tests/KmerFactoryTest.o \
	${TESTDIR}/tests/PeakTest.o \
	${TESTDIR}/tests/PhyperTest.o

# C Compiler Flags
CFLAGS=-g -Wall

# CC Compiler Flags
CCFLAGS=-g -Wall
CXXFLAGS=-g -Wall

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ../../lib/libkm.a

../../lib/libkm.a: ${OBJECTFILES}
	${MKDIR} -p ../../lib
	${RM} ../../lib/libkm.a
	${AR} -rv ../../lib/libkm.a ${OBJECTFILES} 
	$(RANLIB) ../../lib/libkm.a

${OBJECTDIR}/src/BedFactory.o: src/BedFactory.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/BedFactory.o src/BedFactory.cpp

${OBJECTDIR}/src/FastaFactory.o: src/FastaFactory.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/FastaFactory.o src/FastaFactory.cpp

${OBJECTDIR}/src/KmersFactory.o: src/KmersFactory.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/KmersFactory.o src/KmersFactory.cpp

${OBJECTDIR}/src/SNPFactory.o: src/SNPFactory.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/SNPFactory.o src/SNPFactory.cpp

${OBJECTDIR}/src/SVMPredict.o: src/SVMPredict.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/SVMPredict.o src/SVMPredict.cpp

${OBJECTDIR}/src/berror.o: src/berror.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/berror.o src/berror.c

${OBJECTDIR}/src/bmemory.o: src/bmemory.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bmemory.o src/bmemory.c

${OBJECTDIR}/src/bstring.o: src/bstring.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bstring.o src/bstring.c

${OBJECTDIR}/src/chebyshev.o: src/chebyshev.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chebyshev.o src/chebyshev.c

${OBJECTDIR}/src/features.o: src/features.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/features.o src/features.c

${OBJECTDIR}/src/gamma.o: src/gamma.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/gamma.o src/gamma.c

${OBJECTDIR}/src/lgamma.o: src/lgamma.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/lgamma.o src/lgamma.c

${OBJECTDIR}/src/lgammacor.o: src/lgammacor.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/lgammacor.o src/lgammacor.c

${OBJECTDIR}/src/phyper.o: src/phyper.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/phyper.o src/phyper.c

${OBJECTDIR}/src/stirlerr.o: src/stirlerr.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/stirlerr.o src/stirlerr.c

${OBJECTDIR}/src/svm.o: src/svm.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/svm.o src/svm.cpp

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-tests-subprojects .build-conf ${TESTFILES}
.build-tests-subprojects:

${TESTDIR}/TestFiles/bedFactoryTest: ${TESTDIR}/tests/BedFactoryTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/bedFactoryTest -Wl,-S $^ ${LDLIBSOPTIONS} 

${TESTDIR}/TestFiles/fastaFactoryTest: ${TESTDIR}/tests/FastaFactoryTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/fastaFactoryTest $^ ${LDLIBSOPTIONS} 

${TESTDIR}/TestFiles/kmerFactoryTest: ${TESTDIR}/tests/KmerFactoryTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/kmerFactoryTest $^ ${LDLIBSOPTIONS} 

${TESTDIR}/TestFiles/peakTest: ${TESTDIR}/tests/PeakTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/peakTest $^ ${LDLIBSOPTIONS} 

${TESTDIR}/TestFiles/phyperTest: ${TESTDIR}/tests/PhyperTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/phyperTest $^ ${LDLIBSOPTIONS} 


${TESTDIR}/tests/BedFactoryTest.o: tests/BedFactoryTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/BedFactoryTest.o tests/BedFactoryTest.cpp


${TESTDIR}/tests/FastaFactoryTest.o: tests/FastaFactoryTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -Iincludes -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/FastaFactoryTest.o tests/FastaFactoryTest.cpp


${TESTDIR}/tests/KmerFactoryTest.o: tests/KmerFactoryTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -Iincludes -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/KmerFactoryTest.o tests/KmerFactoryTest.cpp


${TESTDIR}/tests/PeakTest.o: tests/PeakTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -Iincludes -I. -std=c++14 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/PeakTest.o tests/PeakTest.cpp


${TESTDIR}/tests/PhyperTest.o: tests/PhyperTest.c 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -Iincludes -I. -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/PhyperTest.o tests/PhyperTest.c


${OBJECTDIR}/src/BedFactory_nomain.o: ${OBJECTDIR}/src/BedFactory.o src/BedFactory.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/BedFactory.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Iincludes -std=c++14 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/BedFactory_nomain.o src/BedFactory.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/BedFactory.o ${OBJECTDIR}/src/BedFactory_nomain.o;\
	fi

${OBJECTDIR}/src/FastaFactory_nomain.o: ${OBJECTDIR}/src/FastaFactory.o src/FastaFactory.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/FastaFactory.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Iincludes -std=c++14 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/FastaFactory_nomain.o src/FastaFactory.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/FastaFactory.o ${OBJECTDIR}/src/FastaFactory_nomain.o;\
	fi

${OBJECTDIR}/src/KmersFactory_nomain.o: ${OBJECTDIR}/src/KmersFactory.o src/KmersFactory.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/KmersFactory.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Iincludes -std=c++14 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/KmersFactory_nomain.o src/KmersFactory.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/KmersFactory.o ${OBJECTDIR}/src/KmersFactory_nomain.o;\
	fi

${OBJECTDIR}/src/SNPFactory_nomain.o: ${OBJECTDIR}/src/SNPFactory.o src/SNPFactory.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/SNPFactory.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Iincludes -std=c++14 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/SNPFactory_nomain.o src/SNPFactory.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/SNPFactory.o ${OBJECTDIR}/src/SNPFactory_nomain.o;\
	fi

${OBJECTDIR}/src/SVMPredict_nomain.o: ${OBJECTDIR}/src/SVMPredict.o src/SVMPredict.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/SVMPredict.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Iincludes -std=c++14 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/SVMPredict_nomain.o src/SVMPredict.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/SVMPredict.o ${OBJECTDIR}/src/SVMPredict_nomain.o;\
	fi

${OBJECTDIR}/src/berror_nomain.o: ${OBJECTDIR}/src/berror.o src/berror.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/berror.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iincludes -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/berror_nomain.o src/berror.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/berror.o ${OBJECTDIR}/src/berror_nomain.o;\
	fi

${OBJECTDIR}/src/bmemory_nomain.o: ${OBJECTDIR}/src/bmemory.o src/bmemory.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/bmemory.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iincludes -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bmemory_nomain.o src/bmemory.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/bmemory.o ${OBJECTDIR}/src/bmemory_nomain.o;\
	fi

${OBJECTDIR}/src/bstring_nomain.o: ${OBJECTDIR}/src/bstring.o src/bstring.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/bstring.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iincludes -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bstring_nomain.o src/bstring.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/bstring.o ${OBJECTDIR}/src/bstring_nomain.o;\
	fi

${OBJECTDIR}/src/chebyshev_nomain.o: ${OBJECTDIR}/src/chebyshev.o src/chebyshev.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chebyshev.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iincludes -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chebyshev_nomain.o src/chebyshev.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chebyshev.o ${OBJECTDIR}/src/chebyshev_nomain.o;\
	fi

${OBJECTDIR}/src/features_nomain.o: ${OBJECTDIR}/src/features.o src/features.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/features.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iincludes -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/features_nomain.o src/features.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/features.o ${OBJECTDIR}/src/features_nomain.o;\
	fi

${OBJECTDIR}/src/gamma_nomain.o: ${OBJECTDIR}/src/gamma.o src/gamma.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/gamma.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iincludes -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/gamma_nomain.o src/gamma.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/gamma.o ${OBJECTDIR}/src/gamma_nomain.o;\
	fi

${OBJECTDIR}/src/lgamma_nomain.o: ${OBJECTDIR}/src/lgamma.o src/lgamma.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/lgamma.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iincludes -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/lgamma_nomain.o src/lgamma.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/lgamma.o ${OBJECTDIR}/src/lgamma_nomain.o;\
	fi

${OBJECTDIR}/src/lgammacor_nomain.o: ${OBJECTDIR}/src/lgammacor.o src/lgammacor.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/lgammacor.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iincludes -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/lgammacor_nomain.o src/lgammacor.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/lgammacor.o ${OBJECTDIR}/src/lgammacor_nomain.o;\
	fi

${OBJECTDIR}/src/phyper_nomain.o: ${OBJECTDIR}/src/phyper.o src/phyper.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/phyper.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iincludes -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/phyper_nomain.o src/phyper.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/phyper.o ${OBJECTDIR}/src/phyper_nomain.o;\
	fi

${OBJECTDIR}/src/stirlerr_nomain.o: ${OBJECTDIR}/src/stirlerr.o src/stirlerr.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/stirlerr.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iincludes -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/stirlerr_nomain.o src/stirlerr.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/stirlerr.o ${OBJECTDIR}/src/stirlerr_nomain.o;\
	fi

${OBJECTDIR}/src/svm_nomain.o: ${OBJECTDIR}/src/svm.o src/svm.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/svm.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -O2 -Iincludes -std=c++14 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/svm_nomain.o src/svm.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/src/svm.o ${OBJECTDIR}/src/svm_nomain.o;\
	fi

# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${TESTDIR}/TestFiles/bedFactoryTest || true; \
	    ${TESTDIR}/TestFiles/fastaFactoryTest || true; \
	    ${TESTDIR}/TestFiles/kmerFactoryTest || true; \
	    ${TESTDIR}/TestFiles/peakTest || true; \
	    ${TESTDIR}/TestFiles/phyperTest || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ../../lib/libkm.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
